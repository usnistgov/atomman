# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

# atomman imports
try:
    from .cythonized import dvect_cython
except:
    cython_imported = False
else:
    cython_imported = True
    
from ..compatibility import range
    
def dvect(pos_0, pos_1, box, pbc, code=None):
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
    code: str, optional
        Option for specifying which underlying code function to use:
        - 'cython' uses the version of the function built in cython (faster).
        - 'python' uses the purely python version.
        Default is 'cython' if the code can be imported, otherwise 'python'.
    """
    
    if code is None:
        if cython_imported is True:
            return dvect_cython(pos_0, pos_1, box, pbc)
        elif cython_imported is False:
            return dvect_python(pos_0, pos_1, box, pbc)
    
    elif code == 'cython':
        if cython_imported is True:
            return dvect_cython(pos_0, pos_1, box, pbc)
        else:
            raise ValueError('cython version of dvect not loaded')
    
    elif code == 'python':
        return dvect_python(pos_0, pos_1, box, pbc)
    
    else:
        raise ValueError("Invalid code style: only 'cython' and 'python' allowed")
    

def dvect_python(pos_0, pos_1, box, pbc):
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
    
    # Convert positions to np.arrays
    pos_0 = np.asarray(pos_0)
    pos_1 = np.asarray(pos_1)
    if pos_0.ndim == 0:
        raise TypeError('Invalid pos_0')
    if pos_1.ndim == 0:
        raise TypeError('Invalid pos_1')
    
    # Get box values
    avect = box.avect
    bvect = box.bvect
    cvect = box.cvect 
    
    # Compute the non-periodic distances between pos_0 and pos_1
    delta = pos_1 - pos_0
    if delta.ndim == 1:
        delta = delta[np.newaxis]
    
    # Create iterators based on pbc
    check = [(0,), (0,), (0,)]
    for i in range(3):
        if pbc[i]:
            check[i] = (-1, 0, 1)
    
    # Add all combinations of system vectors to delta to identify shortest d vector(s)
    d = delta.copy()
    for x in check[0]:
        for y in check[1]:
            for z in check[2]:
                test = delta + (x * avect + y * bvect + z * cvect)
                d = np.where(d.T[0]**2+d.T[1]**2+d.T[2]**2 <= test.T[0]**2+test.T[1]**2+test.T[2]**2, d.T, test.T).T
    
    if len(d) == 1:
        return d[0]
    else:
        return d