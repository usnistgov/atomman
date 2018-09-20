# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

from .. import displacement

def disregistry(basesystem, dislsystem, spreadvector=[1.0, 0.0, 0.0],
                planenormal=[0.0, 1.0, 0.0], planepos=[0.0, 0.0, 0.0]):
    """
    Computes the disregistry profile for a dislocation system.  Assumes that
    the dislocation line is along the z-axis and the slip plane is the y=0
    plane.  These assumptions match how atomman.defect.Stroh generates
    dislocations.
    
    Parameters
    ----------
    basesystem : atomman.System
        A perfect reference system with atoms directly corresponding to atoms
        in dislsystem.
    dislsystem : atomman.System
        A dislocation-containing system.
    spreadvector : array-like object, optional
        Unit vector defining the direction associated with the coord-coordinate
        spreading along the slip plane.  Default value is [1, 0, 0]
        (Cartesian coord-coordinates).
    planenormal : array-like object, optional
        Unit vector defining the normal to the slip plane.  Must be
        perpendicular to spreadvector.  Default value is [0, 1, 0] (Cartesian
        y-axis).
    planepos : array-like object, optional
        A position on the slip plane so that the plane can be fully defined.
        The slip plane position should fall between two planes of atoms.
        Default value is [0,0,0].
    
    Returns
    -------
    coord : numpy.ndarray
        The (N,) array of unique coord-coordinates (atomic columns) neighboring
        the slip plane.
    disregistry : numpy.ndarray
        A (N, 3) array of the dislocation's disregistry at each coord.
    """
    # Handle spreadvector, planenormal and planepos
    spreadvector = np.asarray(spreadvector, dtype=float)
    planenormal = np.asarray(planenormal, dtype=float)
    planepos = np.asarray(planepos, dtype=float)
    if not np.isclose(spreadvector.dot(planepos), 0.0, atol=1e-8, rtol=0.0):
        raise ValueError('spreadvector and planepos must be perpendicular')
    
    # Extract pos from basesystem and compute atomic displacements
    basepos = basesystem.atoms.pos
    disp = displacement(basesystem, dislsystem)
    
    # Transform basepos and planepos to calculation coordinates
    allx = np.dot(basepos, spreadvector)
    ally = np.dot(basepos, planenormal)
    midy = np.dot(planepos, planenormal)
    
    # Identify atoms just above and below slip plane
    uniquey = np.unique(ally)
    abovey = uniquey[uniquey > midy].min()
    belowy = uniquey[uniquey < midy].max()
    if np.isclose(abovey, belowy):
        raise ValueError('planepos must fall between atomic planes')
    #uniquey = np.unique(basepos[:, 1])
    #abovey = uniquey[uniquey > 0].min()
    #belowy = uniquey[uniquey < 0].max()
    
    # Get coordinates of atoms just above and below slip plane
    abovex = allx[np.isclose(ally, abovey)]
    belowx = allx[np.isclose(ally, belowy)]
    #abovex = basepos[np.isclose(basepos[:, 1], abovey)][:, 0]
    #belowx = basepos[np.isclose(basepos[:, 1], belowy)][:, 0]
    
    # Identify unique coordinates
    uabovex = np.unique(abovex)
    ubelowx = np.unique(belowx)
    coord = np.union1d(uabovex, ubelowx)
    
    # Get displacements of atoms just above and below slip plane
    abovedisp = disp[np.isclose(ally, abovey)]
    belowdisp = disp[np.isclose(ally, belowy)]
    #abovedisp = disp[np.isclose(basepos[:, 1], abovey)]
    #belowdisp = disp[np.isclose(basepos[:, 1], belowy)]
    
    # Average displacements for the same coordinates
    abovedispmean = np.empty((len(uabovex), 3))
    for i, ix in enumerate(uabovex):
        abovedispmean[i] = abovedisp[np.isclose(abovex, ix)].mean(axis=0)
    belowdispmean = np.empty((len(ubelowx), 3))
    for i, ix in enumerate(ubelowx):
        belowdispmean[i] = belowdisp[np.isclose(belowx, ix)].mean(axis=0)
    
    # Linearly interpolate displacement values to all coord
    abovedispinterp = np.vstack([np.interp(coord, uabovex, abovedispmean[:,0]),
                                 np.interp(coord, uabovex, abovedispmean[:,1]),
                                 np.interp(coord, uabovex, abovedispmean[:,2])]).T
    belowdispinterp = np.vstack([np.interp(coord, ubelowx, belowdispmean[:,0]),
                                 np.interp(coord, ubelowx, belowdispmean[:,1]),
                                 np.interp(coord, ubelowx, belowdispmean[:,2])]).T
    
    # Calculate disregistry
    disregistry = abovedispinterp - belowdispinterp
    
    return coord, disregistry