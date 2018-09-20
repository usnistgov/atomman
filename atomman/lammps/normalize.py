# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from copy import deepcopy

# http://www.numpy.org/
import numpy as np

def normalize(system, return_transform=False):
    """
    The normalize function takes any arbitrary system and transforms it to
    be compatible with LAMMPS.  In particular, LAMMPS systems must have:
    1. Right-handed box vectors.
    2. avect = [lx, 0.0, 0.0]
    3. bvect = [xy, ly,  0.0]
    4. cvect = [xz, yz,  lz]
    5. All atoms initially inside the box dimensions.
    Note: large box tilt factors are not adjusted with this function.
    As such, the LAMMPS command 'box tilt large' may be needed.
    
    Parameters
    ----------
    system : atomman.System
        The system to normalize.
    return_transform : bool, optional
        Indicates if the transformation matrix used during the normalization
        is to be returned as well.  Default value is False
        
    Returns
    -------
    newsystem : atomman.System
        A new system that has been normalized.
    transform : np.ndarray
        The transformation matrix associated with the normalization.  Returned
        if return_transform is True.
    """
    
    # Create a copy of the system 
    system = deepcopy(system)
    
    # Swap cvector direction if box is left-handed.
    if np.dot(np.cross(system.box.avect, system.box.bvect), system.box.cvect) < 0:
        system.box_set(avect=system.box.avect,
                       bvect=system.box.bvect,
                       cvect=-system.box.cvect,
                       origin=system.box.origin + system.box.cvect)
    
    # vects represents Cartesian vectors before transformation
    vects = deepcopy(system.box.vects)
    
    # Rebuild box using a, b, c, alpha, beta, gamma
    system.box_set(a=system.box.a, b=system.box.b, c=system.box.c,
                   alpha=system.box.alpha, beta=system.box.beta, gamma=system.box.gamma,
                   scale=True)
    
    # Wrap to ensure all atoms are within the box
    system.wrap()
    
    # Update the transformation matrix
    transformation = np.linalg.lstsq(vects, system.box.vects)[0].T
    test1 = np.linalg.norm(transformation, axis=1)
    
    assert np.allclose(test1, np.array([1., 1., 1.])), '%f %f %f' % tuple(test1)
    assert np.isclose(transformation[0].dot(transformation[1]), 0.0)
    assert np.isclose(transformation[0].dot(transformation[2]), 0.0)
    assert np.isclose(transformation[1].dot(transformation[2]), 0.0)
    
    if return_transform:
        return system, transformation
    else:
        return system