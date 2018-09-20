# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from copy import deepcopy

# http://www.numpy.org/
import numpy as np

# atomman imports
from ..compatibility import range
try:
    from .cythonized import nlist_cython
except:
    cython_imported = False
else:
    cython_imported = True
from .dvect import dvect

def nlist(system, cutoff, cmult=1, code=None, initialsize=20, deltasize=10):
    """
    Calculates a neighbor list for all atoms in a System taking periodic
    boundaries into account.
    
    Parameters
    ----------
    system : atomman.System 
        The system to calculate the neighbor list for.
    cutoff : float
        Radial cutoff distance for identifying neighbors.
    cmult : int, optional
        Parameter associated with the binning routine. Default value is most
        likely the fastest.
    code: str, optional
        Option for specifying which underlying code function to use:
        - 'cython' uses the version of the function built in cython (faster).
        - 'python' uses the purely python version.
        Default is 'cython' if the code can be imported, otherwise 'python'.
    initialsize : int, optional
        The number of neighbor positions to initially assign to each atom.
        Default value is 20.
    deltasize : int, optional
        Specifies the number of extra neighbor positions to allow each atom
        when the number of neighbors exceeds the underlying array size.
        Default value is 10.
        
    Returns
    -------
    numpy.ndarray
        Two dimensional array where the first dimension corresponds to each
        atom.  In each row (second dimension), the first term is the atom's
        coordination number, and the following terms from 1 to coordination+1
        are the list of neighboring atoms (by index).
    """
    
    if code is None:
        if cython_imported is True:
            return nlist_cython(system, cutoff, cmult=cmult, code=code,
                                initialsize=initialsize, deltasize=deltasize)
        elif cython_imported is False:
            return nlist_python(system, cutoff, cmult=cmult, code=code,
                                initialsize=initialsize, deltasize=deltasize)
    
    elif code == 'cython':
        if cython_imported is True:
            return nlist_cython(system, cutoff, cmult=cmult, code=code,
                                initialsize=initialsize, deltasize=deltasize)
        else:
            raise ValueError('cython version of nlist not loaded')
    
    elif code == 'python':
        return nlist_python(system, cutoff, cmult=cmult, code=code,
                            initialsize=initialsize, deltasize=deltasize)
    
    else:
        raise ValueError("Invalid code style: only 'cython' and 'python' allowed")

def nlist_python(system, cutoff, cmult=1, code=None, initialsize=20,
                 deltasize=10):
    """
    Calculates a neighbor list for all atoms in a System taking periodic
    boundaries into account.
    
    Parameters
    ----------
    system : atomman.System 
        The system to calculate the neighbor list for.
    cutoff : float
        Radial cutoff distance for identifying neighbors.
    cmult : int, optional
        Parameter associated with the binning routine. Default value is most
        likely the fastest.
    code: str, optional
        Option for specifying which underlying code function of dvect to use:
        - 'cython' uses the version of the function built in cython (faster).
        - 'python' uses the purely python version.
        Default is 'cython' if the code can be imported, otherwise 'python'.
    initialsize : int, optional
        The number of neighbor positions to initially assign to each atom.
        Default value is 20.
    deltasize : int, optional
        Specifies the number of extra neighbor positions to allow each atom
        when the number of neighbors exceeds the underlying array size.
        Default value is 10.
    """
    
    # Extract parameters
    natoms = system.natoms
    vects = system.box.vects
    origin = system.box.origin
    pos = system.atoms.view['pos']
    pbc = system.pbc
    
    # Build neighbors array
    neighbors = np.zeros((natoms, initialsize + 1), dtype='int64')
    
    # Determine orthogonal super-box that fully encompasses the system
    corners = origin + np.array([[0, 0, 0], 
                                 vects[0], 
                                 vects[1], 
                                 vects[2],
                                 vects[0] + vects[1],
                                 vects[0] + vects[2],
                                 vects[1] + vects[2],
                                 vects[0] + vects[1] + vects[2]])
    supermin = corners.min(axis=0) - 1.01 * cutoff
    supermax = corners.max(axis=0) + 1.01 * cutoff
    
    # Construct bins
    binsize = cutoff / cmult
    xbins = np.arange(supermin[0], supermax[0] + binsize, binsize)
    ybins = np.arange(supermin[1], supermax[1] + binsize, binsize)
    zbins = np.arange(supermin[2], supermax[2] + binsize, binsize)

    # Build xyz box index for each atom
    xindex = np.digitize(pos[:, 0], xbins) - 1
    yindex = np.digitize(pos[:, 1], ybins) - 1
    zindex = np.digitize(pos[:, 2], zbins) - 1
    xyzindex = np.hstack((xindex[:, np.newaxis], yindex[:, np.newaxis], zindex[:, np.newaxis]))
    
    # Relate atom's id to xyz_index
    atomindex = np.arange(natoms, dtype='int64')
    
    # Identify all bins with real atoms
    realbins = unique_rows(xyzindex)
    
    # Create iterators based on pbc
    check = [(0,), (0,), (0,)]
    for i in range(3):
        if pbc[i]:
            check[i] = (-1, 0, 1) 
    
    # Construct list of ghost atoms in the super-box
    ghostpos = []
    for x in check[0]:
        for y in check[1]:
            for z in check[2]:
                if x == 0 and y == 0 and z == 0:
                    pass
                else:
                    newpos = x * vects[0] + y * vects[1] + z * vects[2] + pos
                    bool_check = np.empty_like(newpos, dtype=bool)
                    for i in range(3):
                        bool_check[:, i] = np.all([newpos[:, i] > supermin[i], newpos[:, i] < supermax[i]], axis=0)
                    
                    newindex = np.where(np.all(bool_check, axis=1))[0]
                    try:
                        ghostpos = np.vstack((ghostpos, newpos[newindex]))
                        ghostindex = np.hstack((ghostindex, newindex))
                    except:
                        ghostpos = newpos[newindex]
                        ghostindex = newindex
    
    # Append xyzindex and atomindex lists with ghost atoms
    if len(ghostpos) > 0:
        xindex = np.digitize(ghostpos[:, 0], xbins) - 1
        yindex = np.digitize(ghostpos[:, 1], ybins) - 1
        zindex = np.digitize(ghostpos[:, 2], zbins) - 1
        xyzghostindex = np.hstack((xindex[:, np.newaxis], 
                                   yindex[:, np.newaxis],
                                   zindex[:, np.newaxis]))
        xyzindex = np.vstack((xyzindex, xyzghostindex))
        atomindex = np.hstack((atomindex, ghostindex))
    
    # Assign atoms and ghost atoms to xyz bins
    xyzbins = np.zeros((len(xbins), len(ybins), len(zbins), 41), dtype='int64')
    for i in range(len(atomindex)):
        x, y, z = xyzindex[i]
        xyzbins[x, y, z, 0] += 1
        try:
            xyzbins[x, y, z, xyzbins[x, y, z, 0]] = atomindex[i]
        except:
            oldsize = len(xyzbins[0, 0, 0])
            newbins = np.zeros((len(xbins), len(ybins), len(zbins), oldsize+10),
                                dtype='int64')
            newbins[:, :, :, :oldsize] = xyzbins[:, :, :, :oldsize]
            xyzbins = newbins
            xyzbins[x, y, z, xyzbins[x, y, z, 0]] = atomindex[i]
    
    # Iterate over all bins with real atoms
    for realbin in realbins:
        x, y, z = realbin
        
        # shortlist = all atoms in current bin
        shortlist = xyzbins[x, y, z, 1:xyzbins[x, y, z, 0] + 1]
        longlist = deepcopy(shortlist)
        
        # Add all atoms in half of the nearby bins to longlist
        for dx, dy, dz in box_iter(cmult):
            dc = xyzbins[x + dx, y + dy, z + dz, 0]
            longlist = np.hstack((longlist, xyzbins[x + dx, y + dy, z + dz, 1:dc + 1]))
        
        # Compare all atoms in shortlist to mediumlist (which is longlist starting with u+1).
        for u in range(len(shortlist)):
            upos = pos[shortlist[u]]
            uindex = shortlist[u]
            mediumlist = longlist[u + 1:]
            vpos = pos[mediumlist]
            try:
                d = np.linalg.norm(system.dvect(upos, vpos, code=code), axis=1)
            except:
                d = np.linalg.norm([system.dvect(upos, vpos, code=code)], axis=1)
            vlist = mediumlist[np.where(d < cutoff)]
            for vindex in vlist:
                neighbors = __append_neighbor(neighbors, uindex, vindex, deltasize)
    
    
    # Sort each atom's neighbors
    for i in range(len(neighbors)):
        neighbors[i, 1 : neighbors[i, 0]+1] = np.sort(neighbors[i, 1 : neighbors[i, 0]+1])
    
    return neighbors

def __append_neighbor(n, a, b, deltasize=10):
    """Adds atom ids a and b to each other's list of neighbors."""
    if b not in n[a, 1:n[a, 0]+1] and a != b:
        n[a, 0] += 1
        n[b, 0] += 1
        try:
            n[a, n[a, 0]] = b
            n[b, n[b, 0]] = a
        except:
            old_size = len(n[0])
            newlist = np.zeros((len(n), old_size + deltasize), dtype='int64')
            newlist[:, :old_size] = n[:, :old_size]
            n = newlist
            n[a, n[a, 0]] = b
            n[b, n[b, 0]] = a
    return n

def unique_rows(a):
    """Takes two-dimensional array a and returns only the unique rows."""
    return np.unique(a.view(np.dtype((np.void, a.dtype.itemsize*a.shape[1])))).view(a.dtype).reshape(-1, a.shape[1])

def box_iter(num):
    """
    Iterates over unique x,y,z combinations where each x,y,z can equal +-num.
    Excludes 0,0,0 and equal but opposite combinations (e.g. x,y,z used, -x,-y,-z excluded).
    """
    z = -num
    end = False
    while z <= 0:
        y = -num
        while y <= num:
            x = -num
            while x <= num:
                if x == 0 and y == 0 and z == 0:
                    end = True
                    break
                yield x, y, z
                x+=1
            if end is True:
                break
            y+=1
        if end is True:
            break
        z+=1