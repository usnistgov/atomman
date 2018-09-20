# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

def axes_check(axes, tol=1e-8):
    """
    Checks that given axes are orthogonal and right-handed.
    
    Parameters
    ----------
    axes : list of list or np.ndarray
        A 3x3 array of axes vectors.
    tol : float, optional
        Tolerance to use in checking if axes are orthogonal and right-handed.
        Default value is 1e-8.
    
    Returns
    -------
    np.ndarray
        A 3x3 array of the corresponding unit vectors of axes.
    
    Raises
    ------
    ValueError
        If axes are not orthogonal or right-handed.
    """
    # Convert to numpy.ndarray if needed
    axes = np.asarray(axes)
    assert axes.shape == (3,3), 'Invalid axes shape'
    
    # Compute unit axes
    uaxes = (axes.T / np.linalg.norm(axes, axis=1)).T
    
    # Check if orthogonal
    if not np.allclose(np.dot(uaxes, uaxes.T), np.identity(3), atol=tol):
        raise ValueError('axes are not orthogonal')
    
    # Check if right-handed
    if not np.allclose(np.cross(uaxes[0], uaxes[1]), uaxes[2], atol=tol):
        raise ValueError('axes are not right-handed')
    
    return uaxes