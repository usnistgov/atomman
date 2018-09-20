# coding: utf-8
# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

# atomman imports
from .. import NeighborList
from ..compatibility import range

def slip_vector(system_0, system_1, neighbors=None, cutoff=None):
    """
    Compute the slip vectors for all atoms.  Note that this differs from the
    original formulation in that it is not normalized by number of slipped
    neighbors.
    
        s_i = -Σ_j d_ij(t) - d_ij(0)
    
    where j is neighbor atoms of atom i, and d_ij() is vector distance between
    atoms i and j at time t.
    
    Parameters
    ----------
    system_0 : atomman.system
        The base/reference system to use.
    system_1 : atomman.system
        The defect/current system to use.
    neighbors : atomman.NeighborList, optional
        The neighbor list associated with system_0 to use.  Either neighbors
        or cutoff must be given, or system_0 must have a neighbors attribute.
    cutoff : float, optional
        Cutoff distance for computing a neighbor list for system_0.  Either
        neighbors or cutoff must be given, or system_0 have a neighbors
        attribute.
    
    Returns
    -------
    numpy.ndarray
        The computed slip vectors.
    """
    # Check that the two systems have the same number of atoms
    if system_0.natoms != system_1.natoms:
        raise ValueError('systems have different number of atoms')
    
    # Neighbor list setup
    if neighbors is not None:
        assert cutoff is None, 'neighbors and cutoff cannot both be given'
    elif cutoff is not None:
        neighbors = NeighborList(system=system_0, cutoff=cutoff)
    elif hasattr(system_0, 'neighbors'):
        neighbors = system_0.neighbors
    else:
        raise ValueError('neighbors or cutoff is required')
    
    # Calculate the slip vector
    slip = np.zeros((system_0.natoms, 3))
    for i in range(system_0.natoms):
        slip[i] = -np.sum(system_1.dvect(i, neighbors[i])
                        - system_0.dvect(i, neighbors[i]), axis=0)
    
    return slip