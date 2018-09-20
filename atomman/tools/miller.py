# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

from .crystalsystem import ishexagonal

__all__ = ['vector3to4', 'vector4to3', 'vectortocartesian']

def vector3to4(indices):
    """
    Converts 3-term Miller [uvw] indices to 4-term hexagonal [uvtw]
    Miller-Bravias indices. Note that vectors will be normalized to
    smallest integer representations.
    
    Parameters
    ----------
    indices : np.ndarray of int
        (..., 3) array of Miller crystallographic indices.
   
    Returns
    -------
    np.ndarray of int
        (..., 4) array of Miller-Bravais crystallographic indices.
    
    Raises
    ------
    AssertionError
        If indices dimensions are not (..., 3).
    """
    # Verify compatible indices
    indices = np.asarray(indices, dtype='int64')
    assert indices.shape[-1] == 3, 'Invalid index dimensions'
    
    # Transform 
    newindices = np.empty(indices.shape[:-1] + (4,), dtype='int64')
    newindices[..., 0] = 2 * indices[..., 0] - indices[..., 1]
    newindices[..., 1] = 2 * indices[..., 1] - indices[..., 0]
    newindices[..., 2] = -(newindices[..., 0] + newindices[..., 1])
    newindices[..., 3] = 3 * indices[..., 2]
    
    # Remove common factors
    for i in range(np.abs(newindices).max(), 1, -1):
        test = np.remainder(newindices, i).sum(axis=-1) == 0
        try:
            newindices[test] = newindices[test] / i
        except:
            # Exception for older numpy versions
            if test:
                newindices[:] = newindices / i
    return newindices
    
def vector4to3(indices):
    """
    Converts 4-term hexagonal Miller-Bravias [uvtw] indices to 3-term 
    Miller [uvw] indices. Note that vectors will be normalized to 
    smallest integer representations.
    
    Parameters
    ----------
    indices : np.ndarray of int
        (..., 4) array of Miller-Bravais crystallographic indices.
   
    Returns
    -------
    np.ndarray of int
        (..., 3) array of Miller crystallographic indices.
   
    Raises
    ------    
    AssertionError
        If indices dimensions are not (..., 4), or if h+k+i != 0.
    """
    # Verify compatible indices
    indices = np.asarray(indices, dtype='int64')
    assert indices.shape[-1] == 4, 'Invalid index dimensions'
    assert np.allclose(indices[...,:3].sum(axis=-1), 0.0), 'Invalid indices: h+k+i != 0'
    
    # Transform
    newindices = np.empty(indices.shape[:-1] + (3,), dtype='int64')
    newindices[..., 0] = 2 * indices[..., 0] + indices[..., 1]
    newindices[..., 1] = 2 * indices[..., 1] + indices[..., 0]
    newindices[..., 2] = indices[..., 3]
    
    # Remove common factors
    for i in range(np.abs(newindices).max(), 1, -1):
        test = np.remainder(newindices, i).sum(axis=-1)==0
        try:
            newindices[test] = newindices[test] / i
        except:
            # Exception for older numpy versions
            if test:
                newindices[:] = newindices / i
    return newindices

def vectortocartesian(indices, box):
    """
    Converts crystal indices to Cartesian vectors relative
    to a given lattice box. 
    
    Parameters
    ----------
    indices : np.ndarray of int
        (..., 3) array of [uvw] Miller crystallographic indices or 
        (..., 4) array of [uvtw] Miller-Bravais crystallographic indices.
    box : atomman.Box
        Box that defines the lattice cell vectors to use. 
   
    Returns
    -------
    np.ndarray of float
        (..., 3) array of Cartesian vectors.
        
    Raises
    ------
    AssertionError
        If indices dimensions are not (..., 3) or (..., 4), or if
        hexagonal indices given with non-hexagonal box.
    """
    indices = np.array(indices, dtype='int64')
    
    # Convert 4-term Miller-Bravais to standard 3-term indices
    if indices.shape[-1] == 4:
        assert ishexagonal(box), 'Hexagonal indices given with non-hexagonal box'
        
        indices = vector4to3(indices)
    
    assert indices.shape[-1] == 3, 'Invalid index dimensions'
        
    return (np.outer(indices[...,0], box.vects[0]) 
          + np.outer(indices[...,1], box.vects[1]) 
          + np.outer(indices[...,2], box.vects[2]))