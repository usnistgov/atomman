# coding: utf-8
# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

# atomman imports
from .. import NeighborList
from ..compatibility import range

def nye_tensor_p(system, neighbors=None, cutoff=None):
    """
    Generates a list of p vectors for each atom to be used by the nye_tensor()
    function.  The p vectors correspond to the radial distance vectors between
    an atom and each of its neighbors within a cutoff distance in a "perfect
    crystal" reference state.
    
    Parameters
    ----------
    system : atomman.system
        The base/reference system to use.  This should be a defect-free
        perfect crystal system with atom ids directly corresponding to atoms
        in any system that you want to analyze with the Nye tensor.
    neighbors : atomman.NeighborList, optional
        The neighbor list associated with system to use.  Either neighbors
        or cutoff must be given, or system must have a neighbors attribute.
    cutoff : float
        Cutoff distance for computing a neighbor list for system.  Either
        neighbors or cutoff must be given, or system have a neighbors
        attribute.
    
    Returns
    -------
    numpy.ndarray
        The list of p distance vectors for each atom in system.
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
    
    p = []
    for i in range(system.natoms):
        p.append(system.dvect(i, neighbors[i]))
    
    return np.asarray(p)