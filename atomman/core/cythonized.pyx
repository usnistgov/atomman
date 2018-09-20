# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from copy import deepcopy

# http://cython.org/
import cython
cimport cython

# http://www.numpy.org/
import numpy as np
cimport numpy as np

def dvect_cython(pos_0, pos_1, box, pbc):
    """
    Computes the shortest distance between pos_0 and pos_1 using box 
    dimensions and accounting for periodic boundaries.
    
    Parameters
    ----------
    pos_0 : numpy.ndarray, list, or tuple
        Absolute Cartesian vector position(s) to use as reference point(s).
    pos_1 : numpy.ndarray, list, or tuple
        Absolute Cartesian vector position(s) to find relative to pos_0.
    box : atomman.Box
        Defines the system/box dimensions
    pbc : list, tuple, or numpy.ndarray of bool.
        Three Boolean values indicating which of the three box vectors are
        periodic (True means periodic).
    """
    # Convert pos_0 to numpy array with proper dimensions
    pos_0 = np.asarray(pos_0, dtype='float64')
    if pos_0.ndim == 0:
        raise TypeError('Invalid pos_0')
    if pos_0.ndim == 1:
        pos_0 = pos_0[np.newaxis, :]
    
    # Convert pos_1 to numpy array with proper dimensions
    pos_1 = np.asarray(pos_1, dtype='float64')
    if pos_1.ndim == 0:
        raise TypeError('Invalid pos_1')
    if pos_1.ndim == 1:
        pos_1 = pos_1[np.newaxis, :]

    # Broadcast to compatible lengths
    if len(pos_0) == 1:
        pos_0 = np.broadcast_to(pos_0, pos_1.shape)
    elif len(pos_1) == 1:
        pos_1 = np.broadcast_to(pos_1, pos_0.shape)
    elif len(pos_0) != len(pos_1):
        raise ValueError('Incompatible pos lengths')
        
    # Extract box vectors
    bvects = box.vects
        
    # Call the cython function
    return cdvect(pos_0, pos_1, bvects, pbc)

@cython.boundscheck(False)
cdef np.ndarray[np.float64_t, ndim=2] cdvect(
                                             np.ndarray[np.float64_t, ndim=2] pos_0, 
                                             np.ndarray[np.float64_t, ndim=2] pos_1, 
                                             np.ndarray[np.float64_t, ndim=2] bvects, 
                                             pbc):
    """
    Computes the shortest distance between pos_0 and pos_1 using box 
    dimensions and accounting for periodic boundaries.
    
    Parameters
    ----------
    pos_0 : numpy.ndarray
        Absolute Cartesian vector position(s) to use as reference point(s).
    pos_1 : numpy.ndarray
        Absolute Cartesian vector position(s) to find relative to pos_0.
    bvects : numpy.ndarray
        3x3 array defining the system/box dimensions.
    pbc : list, tuple, or numpy.ndarray of bool.
        Three Boolean values indicating which of the three box vectors are
        periodic (True means periodic).
    """
    
    # Define parameters
    cdef Py_ssize_t ni = pos_0.shape[0]
    cdef Py_ssize_t nj = 3
    
    cdef Py_ssize_t i, j, x, y, z, xl, xu, yl, yu, zl, zu
        
    cdef np.ndarray[np.float64_t, ndim=2] d = np.empty_like(pos_0)
    cdef np.ndarray[np.float64_t] test = np.empty(3)
    cdef np.float64_t test2, d2
    
    # Create iterators based on pbc
    if pbc[0]:
        xl, xh = -1, 2
    else:
        xl, xh = 0, 1
    if pbc[1]:
        yl, yh = -1, 2
    else:
        yl, yh = 0, 1
    if pbc[2]:
        zl, zh = -1, 2
    else:
        zl, zh = 0, 1
    
    # Loop over all pos
    for i in range(ni):
        # Compute pos_1 - pos_0 
        for j in range(nj):
            d[i,j] = pos_1[i,j] - pos_0[i,j] 
        
        # Loop over all periodic boundary conditions
        for x in range(xl, xh):
            for y in range(yl, yh):
                for z in range(zl, zh):
                    
                    # Compute pos_1 - pos_0 + boundary image shifts
                    for j in range(nj):
                        test[j] = (pos_1[i,j] - pos_0[i,j] 
                                   + x * bvects[0,j] 
                                   + y * bvects[1,j] 
                                   + z * bvects[2,j])
                    
                    # Replace d if new vector is smaller
                    test2 = test[0]*test[0]+test[1]*test[1]+test[2]*test[2]
                    d2 = d[i,0]*d[i,0]+d[i,1]*d[i,1]+d[i,2]*d[i,2]
                    if test2 < d2:
                        for j in range(nj):
                            d[i,j] = test[j]
    
    return d

def nlist_cython(system, float cutoff, int cmult=1, code=None, initialsize=20,
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
        Option for specifying which underlying code function of dvect to use
        (ignored here).
    initialsize : int, optional
        The number of neighbor positions to initially assign to each atom.
        Default value is 20.
    deltasize : int, optional
        Specifies the number of extra neighbor positions to allow each atom
        when the number of neighbors exceeds the underlying array size.
        Default value is 10.
    """
    
    pos = system.atoms.pos
    vects = system.box.vects
    origin = system.box.origin
    pbc = system.pbc
    return cnlist(pos, vects, origin, pbc, cutoff, cmult=cmult,
                  initialsize=initialsize, deltasize=deltasize)

cdef cnlist(np.ndarray[np.float64_t, ndim=2] pos,
            np.ndarray[np.float64_t, ndim=2] vects,
            np.ndarray[np.float64_t, ndim=1] origin,
            pbc, float cutoff, int cmult=1,
            size_t initialsize=20,
            size_t deltasize=10):
    
    # Extract and define parameters
    cdef size_t maxneighbors = initialsize
    cdef size_t maxatomsperbin = 40
    cdef bint end, new
    cdef size_t natoms = len(pos)
    cdef float cutoff2 = cutoff**2
    cdef float binsize
    
    cdef Py_ssize_t x, y, z, xl, xh, yl, yh, zl, zh, c, dc, dx, dy, dz
    
    cdef size_t numxbins, numybins, numzbins, uindex, vindex, i, u
    
    cdef np.ndarray[np.float64_t, ndim=1] supermin, supermax, xbins, ybins, zbins, upos, dvect2
    cdef np.ndarray[np.float64_t, ndim=2] corners, vpos, dvect
    
    cdef np.ndarray[np.int64_t, ndim=1] xindex, yindex, zindex, shortlist, mediumlist, longlist
    cdef np.ndarray[np.int64_t, ndim=2] xyzindex, xyzghostindex, newneighbors
    cdef np.ndarray[np.int64_t, ndim=4] xyzbins
    
    cdef np.ndarray[np.float64_t, ndim=2] ghostpos = np.empty((0, 3))
    cdef np.ndarray[np.int64_t, ndim=1] ghostindex = np.empty(0, dtype='int64')
    cdef np.ndarray[np.int64_t, ndim=2] neighbors = np.zeros((natoms, maxneighbors+1), dtype='int64')

    # Determine orthogonal superbox that fully encompasses the system
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
    numxbins = len(xbins)
    numybins = len(ybins)
    numzbins = len(zbins)
    
    # Build xyz box index for each atom
    xindex = np.digitize(pos[:, 0], xbins) - 1
    yindex = np.digitize(pos[:, 1], ybins) - 1
    zindex = np.digitize(pos[:, 2], zbins) - 1
    xyzindex = np.hstack((xindex[:, np.newaxis], yindex[:, np.newaxis], zindex[:, np.newaxis]))
    
    # Relate atom's id to xyz_index
    atomindex = np.arange(natoms, dtype='int64')
    
    # Identify all bins with real atoms
    realbins = unique_rows2(xyzindex)
    
    # Create iterators based on pbc
    if pbc[0]:
        xl, xh = -1, 2
    else:
        xl, xh = 0, 1
    if pbc[1]:
        yl, yh = -1, 2
    else:
        yl, yh = 0, 1
    if pbc[2]:
        zl, zh = -1, 2
    else:
        zl, zh = 0, 1
    
    # Construct list of ghost atoms in the superbox
    for x in range(xl, xh):
        for y in range(yl, yh):
            for z in range(zl, zh):
                if x == 0 and y == 0 and z == 0:
                    pass
                else:
                    newpos = x * vects[0] + y * vects[1] + z * vects[2] + pos
                    bool_check = np.empty_like(newpos, dtype=bool)
                    for i in range(3):
                        bool_check[:, i] = np.all([newpos[:, i] > supermin[i], newpos[:, i] < supermax[i]], axis=0)
                    
                    newindex = np.where(np.all(bool_check, axis=1))[0]
                    ghostpos = np.vstack((ghostpos, newpos[newindex]))
                    ghostindex = np.hstack((ghostindex, newindex))
    
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
    xyzbins = np.zeros((numxbins, numybins, numzbins, maxatomsperbin+1), dtype='int64')
    for i in range(len(atomindex)):
        x, y, z = xyzindex[i]
        c = xyzbins[x, y, z, 0] + 1
        
        # Increase size if needed
        if c == maxatomsperbin:
            newbins = np.zeros((numxbins, numybins, numzbins, maxatomsperbin+11), dtype='int64')
            newbins[:, :, :, :maxatomsperbin+1] = xyzbins[:, :, :, :maxatomsperbin+1]
            xyzbins = newbins
            maxatomsperbin += 10
        
        xyzbins[x, y, z, 0] = c
        xyzbins[x, y, z, c] = atomindex[i]
    
    # Iterate over all bins with real atoms
    for i in range(len(realbins)):
        x, y, z = realbins[i]
        c = xyzbins[x, y, z, 0]
        
        # shortlist = all atoms in current bin
        shortlist = xyzbins[x, y, z, 1:c+1]
        longlist = deepcopy(shortlist)
        
        # Add all atoms in half of the nearby bins to longlist
        end = False
        for dz in range(-cmult, cmult+1):
            for dy in range(-cmult, cmult+1):
                for dx in range(-cmult, cmult+1):
                    if dx == 0 and dy == 0 and dz == 0:
                        end = True
                        break
                    dc = xyzbins[x + dx, y + dy, z + dz, 0]
                    longlist = np.hstack((longlist, xyzbins[x + dx, y + dy, z + dz, 1:dc+1]))
                if end:
                    break
            if end:
                break
        
        # Compare all atoms in shortlist to longlist.
        for u in range(len(shortlist)):
            upos = pos[shortlist[u]]
            uindex = shortlist[u]
            mediumlist = longlist[u + 1:]
            vpos = pos[mediumlist]
            
            # Compare distances to cutoff
            dvect = cdvect(np.broadcast_to(upos, (len(vpos), 3)), vpos, vects, pbc)
            dvect2 = dvect[:, 0]**2 + dvect[:, 1]**2 + dvect[:, 2]**2
            for vindex in mediumlist[dvect2 < cutoff2]:
            
                # Test if new neighbor
                if uindex != vindex:
                    new = True
                    for i in range(neighbors[uindex, 0]):
                        if neighbors[uindex, i+1] == vindex:
                            new = False
                            break
                    if new:
                        # Increase coordination
                        neighbors[uindex, 0] += 1
                        neighbors[vindex, 0] += 1
                        
                        # Extend system size if needed
                        if neighbors[uindex, 0] > maxneighbors or neighbors[vindex, 0] > maxneighbors:
                            newneighbors = np.zeros((natoms, maxneighbors + deltasize + 1), dtype='int64')
                            newneighbors[:, :maxneighbors+1] = neighbors[:, :maxneighbors+1]
                            neighbors = newneighbors
                            maxneighbors += deltasize
                        
                        # Assign neighbors
                        neighbors[uindex, neighbors[uindex, 0]] = vindex
                        neighbors[vindex, neighbors[vindex, 0]] = uindex
    
    # Sort each atom's neighbors
    for i in range(len(neighbors)):
        neighbors[i, 1 : neighbors[i, 0]+1] = np.sort(neighbors[i, 1 : neighbors[i, 0]+1])
    
    return neighbors

cdef unique_rows2(np.ndarray[np.int64_t, ndim=2] a):  
    """Takes two-dimensional array a and returns only the unique rows."""
    return np.unique(a.view(np.dtype((np.void, a.dtype.itemsize*a.shape[1])))).view(a.dtype).reshape(-1, a.shape[1])  
