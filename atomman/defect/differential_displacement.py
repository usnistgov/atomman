# coding: utf-8
# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

# https://matplotlib.org/
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm

# atomman imports
from ..tools import axes_check
from .. import NeighborList

def differential_displacement(system_0, system_1, burgers, plotxaxis='x',
                              plotyaxis='y', xlim=None, ylim=None, zlim=None,
                              neighbors=None, cutoff=None,
                              component='standard', axes=None,
                              plot_scale=1, save_file=None, show=True):
    """
    Generates a differential displacement plot for characterizing dislocations.
    
    Parameters
    ----------
    system_0 : atomman.system
        The base/reference system to use.
    system_1 : atomman.system
        The defect/current system to use.
    burgers : array-like object
        The dislocation's Burgers vector.
    plotxaxis : str or array-like object, optional
        Indicates the Cartesian direction associated with the system's atomic
        coordinates to align with the plotting x-axis.  Values are either 3D
        unit vectors, or strings 'x', 'y', or 'z' for the Cartesian axes
        directions.  plotxaxis and plotyaxis must be orthogonal.  Default value
        is 'x' = [1, 0, 0].
    plotyaxis : str or array-like object, optional
        Indicates the Cartesian direction associated with the system's atomic
        coordinates to align with the plotting y-axis.  Values are either 3D
        unit vectors, or strings 'x', 'y', or 'z' for the Cartesian axes
        directions.  plotxaxis and plotyaxis must be orthogonal.  Default value
        is 'y' = [0, 1, 0].
    xlim : tuple, optional
        The minimum and maximum coordinates along the plotting x-axis to
        include in the fit.  Values are taken in the specified length_unit.
        If not given, then the limits are set based on min and max atomic
        coordinates along the plotting axis.
    ylim : tuple, optional
        The minimum and maximum coordinates along the plotting y-axis to
        include in the fit.  Values are taken in the specified length_unit.
        If not given, then the limits are set based on min and max atomic
        coordinates along the plotting axis.
    zlim : tuple, optional
        The minimum and maximum coordinates normal to the plotting axes
        (i.e. plotxaxis X plotyaxis) to include in the fit.  Values are taken
        in the specified length_unit.  The optimum zlim should encompass only
        a single periodic slice.  If not given, then the limits are set
        based on min and max atomic coordinates along the axis.
    neighbors : atomman.NeighborList, optional
        The neighbor list associated with system_0 to use.  Either neighbors
        or cutoff must be given, or system_0 must have a neighbors attribute.
    cutoff : float, optional
        Cutoff distance for computing a neighbor list for system_0.  Either
        neighbors or cutoff must be given, or system_0 have a neighbors
        attribute.
    component : str, optional
        Indicates the style of the calculation to use.  'standard' (default)
        plots the differential displacements between atoms in the Burgers
        direction.  'xy' plots the differential displacements within the xy
        plotting plane.  This is useful for screw dislocations with localized
        non-screw components.
    axes : arraly-like object, optional
        3x3 transformation array to apply to the Burgers vector to make it
        correspond to the system's orientation.
    plot_scale : float, optional
        Scalar for multiplying the magnitude of the differential displacement
        arrows.  Default value is 1 (no scaling).
    save_file : str, optional
        If given then the plot will be saved to a file with this name.
    show : bool, optional
        Flag for showing the figure. Default is True.
    
    """
    # Transform burgers using axes
    if axes is not None:
        T = axes_check(axes)
        burgers = T.dot(burgers)
    
    
    # Neighbor list setup
    if neighbors is not None:
        assert cutoff is None, 'neighbors and cutoff cannot both be given'
    elif cutoff is not None:
        neighbors = NeighborList(system=system_0, cutoff=cutoff)
    elif hasattr(system_0, 'neighbors'):
        neighbors = system_0.neighbors
    else:
        raise ValueError('neighbors or cutoff is required')
    
    # Give numeric values for str plot axis terms
    if plotxaxis == 'x':
        plotxaxis = [1.0,0.0,0.0]
    elif plotxaxis == 'y':
        plotxaxis = [0.0,1.0,0.0]
    elif plotxaxis == 'z':
        plotxaxis = [0.0,0.0,1.0]
    
    if plotyaxis == 'x':
        plotyaxis = [1.0,0.0,0.0]
    elif plotyaxis == 'y':
        plotyaxis = [0.0,1.0,0.0]
    elif plotyaxis == 'z':
        plotyaxis = [0.0,0.0,1.0]
    
    # Build transformation matrix, T, from plot axes.
    plotxaxis = np.asarray(plotxaxis, dtype='float64')
    plotyaxis = np.asarray(plotyaxis, dtype='float64')
    T = axes_check([plotxaxis, plotyaxis, np.cross(plotxaxis, plotyaxis)])
    
    # Extract positions and transform using T
    pos_0 = np.inner(system_0.atoms.pos, T)
    
    # Transform burgers using plot axes and separate magnitude and direction
    burgers = T.dot(burgers)
    burgers_mag = np.linalg.norm(burgers)
    burgers_uvect = burgers / burgers_mag
    
    # Set plot limits
    if xlim is None:
        xlim = (pos_0[:, 0].min(), pos_0[:, 0].max())
    if ylim is None:
        ylim = (pos_0[:, 1].min(), pos_0[:, 1].max())
    if zlim is None:
        zlim = (pos_0[:, 2].min(), pos_0[:, 2].max())
    
    # Identify only atoms in xlim, ylim, zlim
    in_bounds = ((pos_0[:,0] > xlim[0] - 5.) &
                 (pos_0[:,0] < xlim[1] + 5.) &
                 (pos_0[:,1] > ylim[0] - 5.) &
                 (pos_0[:,1] < ylim[1] + 5.) &
                 (pos_0[:,2] > zlim[0]) &
                 (pos_0[:,2] < zlim[1]))
    
    # Initial plot setup and parameters
    fig, ax1, = plt.subplots(1, 1, squeeze=True, figsize=(7,7), dpi=72)
    ax1.axis([xlim[0], xlim[1], ylim[0], ylim[1]])
    atom_circle_radius = burgers_mag / 10
    arrow_width_scale = 1. / 200.
    
    # Loop over all atoms i in plot range
    for i in np.arange(system_0.natoms)[in_bounds]:
        
        # Plot a circle for atom i
        color = cm.hsv((pos_0[i, 2] - zlim[0]) / (zlim[1] - zlim[0]))
        ax1.add_patch(mpatches.Circle(pos_0[i, :2], atom_circle_radius, fc=color, ec='k'))
        
        # Compute distance vectors between atom i and its neighbors for both systems
        dvectors_0 = np.inner(system_0.dvect(int(i), neighbors[i]), T)
        dvectors_1 = np.inner(system_1.dvect(int(i), neighbors[i]), T)
    
        # Compute differential displacement vectors
        dd_vectors = dvectors_1 - dvectors_0
        
        # Compute center point positions for the vectors
        arrow_centers = pos_0[i] + dvectors_0 / 2
        
        # Plot standard differential displacement component
        if component == 'standard':
            # Compute unit distance vectors
            uvectors_0 = dvectors_0 / np.linalg.norm(dvectors_0, axis=1)[:,np.newaxis]
            
            # Compute component of the dd_vector parallel to the burgers vector
            dd_components = dd_vectors.dot(burgers_uvect)
            dd_components[dd_components > burgers_mag / 2] -= burgers_mag
            dd_components[dd_components < -burgers_mag / 2] += burgers_mag
            
            # Scale arrow lengths and vectors
            arrow_lengths = uvectors_0 * dd_components[:,np.newaxis] * plot_scale
            arrow_widths = arrow_width_scale * dd_components * plot_scale
        
            # Plot the arrows
            for center, length, width in zip(arrow_centers, arrow_lengths,
                                             arrow_widths):
                if width > 1e-7:
                    ax1.quiver(center[0], center[1], length[0], length[1],
                               pivot='middle', angles='xy', scale_units='xy',
                               scale=1, width=width)
        
        # Plot xy differential displacement component
        elif component == 'xy':
            # Scale arrow lengths and vectors
            arrow_lengths = dd_vectors[:, :2] * plot_scale
            arrow_lengths[dd_vectors[:, 2] > 0] *= -1
            arrow_widths = arrow_width_scale * (arrow_lengths[:,0]**2
                                              + arrow_lengths[:,1]**2)**0.5
            
            # Plot the arrows
            for center, length, width in zip(arrow_centers, arrow_lengths, arrow_widths):
                if width > 1e-7:
                    ax1.quiver(center[0], center[1], length[0], length[1], width=width,
                               pivot='middle', angles='xy', scale_units='xy', scale=1)
    
    if save_file is not None:
        plt.savefig(save_file, dpi=800)
    if show == False:
        plt.close(fig)
    plt.show()