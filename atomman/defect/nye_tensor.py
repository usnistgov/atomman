# coding: utf-8
# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
import warnings

# http://www.numpy.org/
import numpy as np

# atomman imports
from .. import NeighborList
from ..compatibility import range
from ..tools import axes_check

def nye_tensor(system, p_vectors, theta_max = 27, axes=None, neighbors=None,
               cutoff=None):
    """
    Computes strain properties and Nye tensor for a defect containing system.
    
    Parameters
    ----------
    system : atomman.System
        The atomic system to compute the per-atom strain properties and Nye
        tensor for.
    p_vectors : array-like object
        List(s) of radial distance vectors between each atom and its nearest
        neighbors in a perfect crystal setting.  If one list of p_vectors is
        given, then it is applied to all atoms.
    theta_max : float, optional
        The maximum theta angle in degrees to use when searching for matches
        between p vectors and q vectors.  Optimum values are dependent on the
        crystal structure. Default value is 27, which is the original value
        used for fcc crystals.
    axes : array-like object, optional
        3x3 array of right-handed orthogonal axes.  If given, will be used to
        transform the p_vectors before computing the Nye tensor.
    neighbors : atomman.NeighborList, optional
        The neighbor list associated with system to use.  Either neighbors
        or cutoff must be given, or system must have a neighbors attribute.
    cutoff : float
        Cutoff distance for computing a neighbor list for system.  Either
        neighbors or cutoff must be given, or system have a neighbors
        attribute.
    
    Returns
    -------
    dict
        Contains the per-atom properties 'strain', 'strain_invariant_1',
        'strain_invariant_2', 'angular_velocity', and 'Nye_tensor'.
    """
    # Neighbor list setup
    if neighbors is not None:
        assert cutoff is None, 'neighbors and cutoff cannot both be given'
    elif cutoff is not None:
        neighbors = NeighborList(system=system, cutoff=cutoff)
    elif hasattr(system, 'neighbors'):
        neighbors = system.neighbors
    else:
        raise ValueError('neighbors or cutoff is required')
    
    # Broadcast a single p_vectors list to all atoms
    if len(p_vectors) == 1:
        p_vectors = np.broadcast_to(p_vectors, (system.natoms, len(p_vectors[0]), 3))
    elif len(p_vectors) != system.natoms:
        p_vectors = np.broadcast_to(p_vectors, (system.natoms, len(p_vectors), 3))
    
    # Transform p_vectors if axes is given
    if axes is not None:
        p_vectors = np.inner(p_vectors, axes_check(axes))
    
    # Get cos of theta_max
    cos_theta_max = np.cos(theta_max * np.pi / 180)
    
    # Define epsilon array
    eps = np.array([[[ 0, 0, 0],[ 0, 0, 1],[ 0,-1, 0]],
                    [[ 0, 0,-1],[ 0, 0, 0],[ 1, 0, 0]],
                    [[ 0, 1, 0],[-1, 0, 0],[ 0, 0, 0]]])
    
    # Identify largest number of nearest neighbors
    nmax = neighbors.coord.max()
    
    # Initialize variables
    strain = np.empty((system.natoms, 3, 3))
    inv1 = np.empty(system.natoms)
    inv2 = np.empty(system.natoms)
    inv3 = np.empty(system.natoms)
    ang_vel = np.empty(system.natoms)
    nye = np.empty((system.natoms, 3, 3))
    P = np.zeros((nmax, 3))
    Q = np.zeros((nmax, 3))
    G = np.empty((system.natoms, 3, 3))
    gradG = np.empty((3, 3, 3))
    
    # Calculate correspondence tensor, G, and strain data for each atom
    for i in range(system.natoms):
        p = np.asarray(p_vectors[i])
        if p.ndim == 1:
            p = np.array([p])
        p_mags = np.linalg.norm(p, axis=1)
        r1 = p_mags.min()
        
        # Calculate radial neighbor vectors, q
        q = system.dvect(i, neighbors[i])
        if q.ndim == 1:
            q = np.array([q])
        q_mags = np.linalg.norm(q, axis=1)
        
        # Calculate cos_thetas between all p's and q's.
        cos_thetas = (np.dot(p, q.T) /q_mags ).T / p_mags
        
        # Identify best index matches
        index_pairing = cos_thetas.argmax(1)
        
        # Exclude values where theta is greater than theta_max
        index_pairing[cos_thetas.max(1) < cos_theta_max] = -1
        
        # Search for duplicate index_pairings
        u, u_count = np.unique(index_pairing, return_counts=True)
        
        # Check if the particular p has already been assigned to another q
        for n in range(len(q)):
            if index_pairing[n] >=0:
                for k in range(n):
                    if index_pairing[n] == index_pairing[k]:
                        nrad = abs(r1 - q_mags[n])
                        krad = abs(r1 - q_mags[k])
                        # Remove the p-q pair that is farther from r1
                        if nrad < krad:
                            index_pairing[k]=-1
                        else:
                            index_pairing[n]=-1
        
        # Construct reduced P, Q matrices from p-q pairs
        c = 0
        for n in range(len(q)):
            if index_pairing[n] >= 0:
                Q[c] = q[n]
                P[c] = p[index_pairing[n]]
                c+=1
                
        # Compute lattice correspondence tensor, G, from P and Q
        if c == 0:
            G[i] = np.identity(3)
            warnings.warn('An atom lacks pair sets. Check neighbor list')
        else:
            G[i], resid, rank, s = np.linalg.lstsq(Q[:c], P[:c])
        
        # Compute strain properties from G
        strain[i] = ((np.identity(3) - G[i]) + (np.identity(3) - G[i]).T) / 2.
        inv1[i] = strain[i,0,0] + strain[i,1,1] + strain[i,2,2]
        inv2[i] = (strain[i,0,0] * strain[i,1,1] 
                 + strain[i,0,0] * strain[i,2,2] 
                 + strain[i,1,1] * strain[i,2,2] 
                - strain[i,0,1]**2 - strain[i,0,2]**2 - strain[i,1,2]**2)
        inv3[i] = np.linalg.det(strain[i])
        rot = ((np.identity(3) - G[i]) - (np.identity(3) - G[i]).T) / 2.
        ang_vel[i] = (rot[0,1]**2 + rot[0,2]**2 + rot[1,2]**2)**0.5
    
    # Construct the gradient tensor of G, gradG for each atom
    for i in range(system.natoms):
        Q = system.dvect(i, neighbors[i])
        if Q.ndim == 1:
            Q = np.array([Q])
        dG = G[neighbors[i]] - G[i]
        for x in range(3):
            gradG[x,:] = np.linalg.lstsq(Q, dG[:,x,:])[0].T
        
        # Use gradG to calculate the nye tensor
        nye[i] = -np.einsum('ijm,ikm->jk', eps, gradG)
    
    return {'strain':strain, 'strain_invariant_1':inv1,
            'strain_invariant_2':inv2, 'strain_invariant_3':inv3,
            'angular_velocity':ang_vel, 'Nye_tensor':nye}