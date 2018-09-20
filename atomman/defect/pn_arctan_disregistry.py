# coding: utf-8
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

# atomman imports
from ..compatibility import int

def pn_arctan_disregistry(xmax=None, xstep=None, xnum=None,
                          burgers=np.array([1.0, 0.0, 0.0]),
                          halfwidth=1, normalize=True, shift=True):
    """
    Computes the classic Peierls-Nabarro arctan disregistry for an array of
    points x.
        
        δ(x) = b / π * arctan(x / ξ) + b / 2
    
    Parameters
    ----------
    xmax : float or None, optional
        Maximum value of x to use.  Minimum value is taken as -xmax.  At least
        2 of xmax, xstep, and xnum must be not None.  Default value is None.
    xstep : float or None, optional
        Step size to use between each x value.  At least 2 of xmax, xstep, and
        xnum must be not None.  Default value is None.
    xnum : int or None, optional
        Number of x values to use.  At least 2 of xmax, xstep, and xnum must
        be not None.  Default value is None.
    burgers : numpy.ndarray, optional
        The Burgers vector for the dislocation. Default value is [1, 0, 0].
    halfwidth : float, optional
        The dislocation halfwidth to use. Default value is 1.
    normalize : bool, optional
        If True (default), the disregistry values will be scaled such that the
        two endpoints differ by exactly one Burgers vector.
    shift : bool, optional
        If True (default), the disregistry will range [0, 0, 0] to burgers.
        If False, the disregistry will range from -burgers to burgers.
        
    Returns
    -------
    x : numpy.ndarray
        The x-coordinates for the disregistry values.
    disregistry : numpy.ndarray
        The disregistry vector at each x-coordinate.
    """
    # Generate missing x parameters
    if xmax is None:
        if xstep is None or xnum is None:
            raise ValueError('At least two parameters must be given')
        xmax = xstep * (xnum - 1) / 2
    elif xstep is None:
        if xnum is None:
            raise ValueError('At least two parameters must be given')
        xstep = (2 * xmax) / (xnum - 1)
    elif xnum is None:
        xnum = ((2 * xmax) / xstep) + 1
    
    # Round xnum to int if needed
    if isinstance(xnum, float):
        if np.isclose(xnum, int(xnum + 0.5)):
            xnum = int(xnum + 0.5)
        else:
            raise ValueError('Invalid parameters: xnum or ((2 * xmax) / xstep) not an integer.')
    
    # Generate x and validate
    x, dx = np.linspace(-xmax, xmax, xnum, retstep=True)
    if not np.isclose(dx, xstep):
        raise ValueError('Incompatible parameters: xmax = xstep * (xnum - 1) / 2')
    
    burgers = np.asarray(burgers)
    
    # δ(x) = b / π * arctan(x / ξ) + b / 2
    disregistry = np.outer(np.arctan(x / halfwidth), burgers / np.pi) + burgers / 2
    
    if normalize is True:
        disregistry = disregistry - disregistry[0]
        disregistry = disregistry * np.linalg.norm(burgers) / np.linalg.norm(disregistry[-1])
    
    if shift is False:
        disregistry -= burgers / 2
    
    return x, disregistry