# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

__all__ = ['iscubic', 'ishexagonal', 'istetragonal', 'isrhombohedral',
           'isorthorhombic', 'ismonoclinic', 'istriclinic', 'identifyfamily']

def iscubic(box, rtol=1e-05, atol=1e-08):
    """
    Tests if a box is consistent with a standard cubic cell:
    a = b = c
    alpha = beta = gamma = 90
    
    Parameters
    ----------
    box : atomman.Box
        The box object to test.
    rtol : float, optional
        Relative tolerance for testing box parameters. Default value is 1e-5.
    atol : float, optional
        Absolute tolerance for testing box parameters. Default value is 1e-8.
        
    Returns
    -------
    bool
        True if box is a standard cubic cell, False otherwise.
    """
    return (np.isclose(box.a, box.b, atol=atol, rtol=rtol)
            and np.isclose(box.a, box.c, atol=atol, rtol=rtol)
            and np.isclose(box.alpha, 90.0, atol=atol, rtol=rtol)
            and np.isclose(box.beta, 90.0, atol=atol, rtol=rtol)
            and np.isclose(box.gamma, 90.0, atol=atol, rtol=rtol))

def ishexagonal(box, rtol=1e-05, atol=1e-08):
    """
    Tests if a box is consistent with a standard hexagonal cell:
    a = b != c
    alpha = beta = 90
    gamma = 120
    
    Parameters
    ----------
    box : atomman.Box
        The box object to test.
    rtol : float, optional
        Relative tolerance for testing box parameters. Default value is 1e-5.
    atol : float, optional
        Absolute tolerance for testing box parameters. Default value is 1e-8.
        
    Returns
    -------
    bool
        True if box is a standard hexagonal cell, False otherwise.
    """
    return (np.isclose(box.a, box.b, atol=atol, rtol=rtol)
            and np.isclose(box.alpha, 90.0, atol=atol, rtol=rtol)
            and np.isclose(box.beta, 90.0, atol=atol, rtol=rtol)
            and np.isclose(box.gamma, 120.0, atol=atol, rtol=rtol))

def istetragonal(box, rtol=1e-05, atol=1e-08):
    """
    Tests if a box is consistent with a standard tetragonal cell:
    a = b != c
    alpha = beta = gamma = 90
    
    Parameters
    ----------
    box : atomman.Box
        The box object to test.
    rtol : float, optional
        Relative tolerance for testing box parameters. Default value is 1e-5.
    atol : float, optional
        Absolute tolerance for testing box parameters. Default value is 1e-8.
        
    Returns
    -------
    bool
        True if box is a standard tetragonal cell, False otherwise.
    """
    return (np.isclose(box.a, box.b, atol=atol, rtol=rtol)
            and not np.isclose(box.a, box.c, atol=atol, rtol=rtol)
            and np.isclose(box.alpha, 90.0, atol=atol, rtol=rtol)
            and np.isclose(box.beta, 90.0, atol=atol, rtol=rtol)
            and np.isclose(box.gamma, 90.0, atol=atol, rtol=rtol))

def isrhombohedral(box, rtol=1e-05, atol=1e-08):
    """
    Tests if a box is consistent with a standard rhombohedral cell:
    a = b = c
    alpha = beta = gamma != 90
    
    Parameters
    ----------
    box : atomman.Box
        The box object to test.
    rtol : float, optional
        Relative tolerance for testing box parameters. Default value is 1e-5.
    atol : float, optional
        Absolute tolerance for testing box parameters. Default value is 1e-8.
        
    Returns
    -------
    bool
        True if box is a standard rhombohedral cell, False otherwise.
    """
    return (np.isclose(box.a, box.b, atol=atol, rtol=rtol)
            and np.isclose(box.a, box.c, atol=atol, rtol=rtol)
            and np.isclose(box.alpha, box.beta, atol=atol, rtol=rtol)
            and np.isclose(box.alpha, box.gamma, atol=atol, rtol=rtol)
            and not np.isclose(box.alpha, 90.0, atol=atol, rtol=rtol))

def isorthorhombic(box, rtol=1e-05, atol=1e-08):
    """
    Tests if a box is consistent with a standard orthorhombic cell:
    a != b != c
    alpha = beta = gamma = 90
    
    Parameters
    ----------
    box : atomman.Box
        The box object to test.
    rtol : float, optional
        Relative tolerance for testing box parameters. Default value is 1e-5.
    atol : float, optional
        Absolute tolerance for testing box parameters. Default value is 1e-8.
        
    Returns
    -------
    bool
        True if box is a standard orthorhombic cell, False otherwise.
    """
    return (not np.isclose(box.a, box.b, atol=atol, rtol=rtol)
            and not np.isclose(box.a, box.c, atol=atol, rtol=rtol)
            and np.isclose(box.alpha, 90.0, atol=atol, rtol=rtol)
            and np.isclose(box.beta, 90.0, atol=atol, rtol=rtol)
            and np.isclose(box.gamma, 90.0, atol=atol, rtol=rtol))

def ismonoclinic(box, rtol=1e-05, atol=1e-08):
    """
    Tests if a box is consistent with a standard monoclinic cell:
    a != b != c
    alpha = gamma = 90
    beta != 90
    
    Parameters
    ----------
    box : atomman.Box
        The box object to test.
    rtol : float, optional
        Relative tolerance for testing box parameters. Default value is 1e-5.
    atol : float, optional
        Absolute tolerance for testing box parameters. Default value is 1e-8.
        
    Returns
    -------
    bool
        True if box is a standard monoclinic cell, False otherwise.
    """
    return (not np.isclose(box.a, box.b, atol=atol, rtol=rtol)
            and not np.isclose(box.a, box.c, atol=atol, rtol=rtol)
            and np.isclose(box.alpha, 90.0, atol=atol, rtol=rtol)
            and not np.isclose(box.beta, 90.0, atol=atol, rtol=rtol)
            and np.isclose(box.gamma, 90.0, atol=atol, rtol=rtol))

def istriclinic(box, rtol=1e-05, atol=1e-08):
    """
    Tests if a box is consistent with a standard triclinic cell:
    a != b != c
    alpha != 90
    beta != 90
    gamma != 90
    
    Parameters
    ----------
    box : atomman.Box
        The box object to test.
    rtol : float, optional
        Relative tolerance for testing box parameters. Default value is 1e-5.
    atol : float, optional
        Absolute tolerance for testing box parameters. Default value is 1e-8.
        
    Returns
    -------
    bool
        True if box is a standard triclinic cell, False otherwise.
    """
    return (not np.isclose(box.a, box.b, atol=atol, rtol=rtol)
            and not np.isclose(box.a, box.c, atol=atol, rtol=rtol)
            and not np.isclose(box.alpha, box.beta, atol=atol, rtol=rtol)
            and not np.isclose(box.alpha, box.gamma, atol=atol, rtol=rtol))

def identifyfamily(box, rtol=1e-05, atol=1e-08):
    """
    Tests if a box is consistent with a standard representation
    of a crystal system cell.
    
    Parameters
    ----------
    box : atomman.Box
        The box object to test.
    rtol : float, optional
        Relative tolerance for testing box parameters. Default value is 1e-5.
    atol : float, optional
        Absolute tolerance for testing box parameters. Default value is 1e-8.
        
    Returns
    -------
    str or None
        'cubic', 'hexagonal', 'tetragonal', 'rhombohedral', 'orthorhombic',
        'monoclinic' or 'triclinic' if it matches any. None if no matches.
        
    Raises
    ------
    ValueError
        If box is not consistent with a standard cell.
    """
    if iscubic(box, rtol=1e-05, atol=1e-08):
        return 'cubic'
    elif ishexagonal(box, rtol=1e-05, atol=1e-08):
        return 'hexagonal'
    elif istetragonal(box, rtol=1e-05, atol=1e-08):
        return 'tetragonal'
    elif isrhombohedral(box, rtol=1e-05, atol=1e-08):
        return 'rhombohedral'
    elif isorthorhombic(box, rtol=1e-05, atol=1e-08):
        return 'orthorhombic'
    elif ismonoclinic(box, rtol=1e-05, atol=1e-08):
        return 'monoclinic'
    elif istriclinic(box, rtol=1e-05, atol=1e-08):
        return 'triclinic'
    else:
        None