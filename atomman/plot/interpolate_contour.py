# coding: utf-8
# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

# https://www.scipy.org/
from scipy.interpolate import griddata

# https://matplotlib.org/
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib import cm

# atomman imports
import atomman.unitconvert as uc
from ..tools import axes_check, vect_angle
from ..compatibility import inttype

def interpolate_contour(system, name, property=None, index=None, magnitude=False,
                        plotxaxis='x', plotyaxis='y', xlim=None, ylim=None,
                        zlim=None, xbins=200, ybins=200, dots=True, czero=True,
                        save=False, show=True, length_unit='angstrom',
                        property_unit=None, cmap='jet'):
    """
    Creates a contour plot of a system's per-atom properties by interpolating
    properties between atoms.
    
    Parameters
    ----------
    system : atomman.System
        The system with the per-atom property that you want to plot.
    name : str
        The name of the per-atom property that you want to plot.
    property : array-like object, optional
        Values for the per-atom property to plot.  If not given, values will
        be taken as the "name" property of system.
    index : int or tuple, optional
        Specifies which component of a multidimensional property to plot.  Not
        needed if the property is scalar.
    magnitude : bool, optional
        If True, plots the per-atom magnitude of a vector property.  Cannot be
        combined with index.  Default value is False.
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
        in the specified length_unit.  If not given, then the limits are set
        based on min and max atomic coordinates along the axis.    
    xbins : int, optional
        Specifies the number of interpolation bins to use along the plotting 
        x-axis.  Default value is 200.
    ybins : int, optional
        Specifies the number of interpolation bins to use along the plotting
        y-axis.  Default value is 200.
    dots : bool, optional
        If True, then the positions of the atoms are shown as circles.
        Default value is True.
    czero : bool, optional
        If True, the range of property values will be centered around zero,
        i.e. cmax = -cmin.  If False, cmax and cmin will be independently
        selected using the property values.  Default value is True.
    save : bool, optional
        If True, the generated plot will be saved to "name.png".  Default
        value is False.
    show : bool, optional
        If True, matplotlib.pyplot.show() is called.  Default value is True.
    length_unit : str, optional
        The units of length to use for the plotting x- and y-axes.  Default
        value is 'angstrom'.
    property_unit : str or None, optional
        The units to use for the property value being plotted.  Default value
        is None, in which no unit conversion is applied.
    cmap : str, optional
        The name of the matplotlib colormap to use.  Default value is 'jet'.
    
    Returns
    -------
    intsum : float
        The area integrated sum of the property over the plotted region.
    avsum : float
        The average property value taken across all plotting bins.
    """
    
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
    pos = uc.get_in_units(np.inner(system.atoms.pos, T), length_unit)
    
    # Set plot limits
    if xlim is None:
        xlim = (pos[:, 0].min(), pos[:, 0].max())
    if ylim is None:
        ylim = (pos[:, 1].min(), pos[:, 1].max())
    if zlim is None:
        zlim = (pos[:, 2].min(), pos[:, 2].max())
    
    # Extract property values
    if property is None:
        property = system.atoms.view[name]
    
    # Handle index
    if index is not None:
        assert magnitude is False, 'index and magnitude cannot be combined'
        if isinstance(index, inttype):
            index = [index]
        else:
            index = list(index)
        for i in index:
            name += '[' + str(i+1) + ']'
        property = property[[Ellipsis] + index]
    
    # Handle magnitude
    elif magnitude is True:
        property = np.linalg.norm(property, axis=1)
        name += '_mag'
    
    assert property.shape == (system.natoms, ), 'property to plot must be a per-atom scalar'
    
    # Identify only atoms in xlim, ylim, zlim
    in_bounds = ((pos[:,0] > xlim[0] - 5.) &
                 (pos[:,0] < xlim[1] + 5.) &
                 (pos[:,1] > ylim[0] - 5.) &
                 (pos[:,1] < ylim[1] + 5.) &
                 (pos[:,2] > zlim[0]) &
                 (pos[:,2] < zlim[1]))
    
    # Extract plotting coordinates and values
    x = pos[in_bounds, 0]
    y = pos[in_bounds, 1]
    v = uc.get_in_units(property[in_bounds], property_unit)
    
    # Generate interpolation grid
    grid, xedges, yedges = grid_interpolate_2d(x, y, v, xbins=xbins, ybins=ybins, range=[xlim, ylim])
    
    # Compute intsum and avsum values
    intsum = np.sum(grid)
    avsum = intsum / (xbins-1) / (ybins-1)
    intsum = intsum * (xlim[1] - xlim[0]) / (xbins-1) * (ylim[1] - ylim[0]) / (ybins-1)
    
    # Generate a pretty figure of the grid
    fig = prettygrid(grid, xedges, yedges, cmap=cmap, propname=name, czero=czero)
    
    # Add dots
    if dots:
        adddots(x, y, xedges, yedges)
    
    # Save, show and close figure
    if save is True:
        plt.savefig(name + '.png', dpi=800)
    if show is True:
        plt.show()
    plt.close(fig)
    
    return intsum, avsum
    
def grid_interpolate_2d(x, y, v, xbins=50, ybins=50, range=None):
    """
    Generates 2D grid of property values by interpolating between measured
    values.
    
    Parameters
    ----------
    x : array-like object
        The list of x-coordinate values.
    y : array-like object
        The list of y-coordinate values.
    v : array-like object
        The list of property values associated with each (x, y) set.
    xbins : int, optional
        The number of bins to use in interpolating between the x coordinates.
        Default value is 50.
    ybins : int, optional
        The number of bins to use in interpolating between the y coordinates.
        Default value is 50.
    range : list, tuple or array-like object
        2x2 list of the [[xmin, xmax], [ymin, ymax]] coordinates for the bins.
        If not given, will use [[x.min(), x.max()], [y.min(), y.max()]].
    
    Returns
    -------
    grid : numpy.ndarray
        The interpolated values in a grid map array.
    xedges : list
        The min, max values of x associated with the grid.
    yedges : list
        The min, max values of y associated with the grid.
    """
    # Handle range and bins options
    if range is None:
        range=[[x.min() ,x.max()], [y.min(), y.max()]]
    
    # Generate 1D interpolation points
    xi = np.linspace(range[0][0], range[0][1], num=xbins)
    yi = np.linspace(range[1][0], range[1][1], num=ybins)
    
    # Generate 2D grid points
    x0, y0 = np.meshgrid(xi, yi)
    
    # Interpolate values to grid points
    grid = griddata((x, y), v, (x0, y0))
    
    return grid, range[0], range[1]

def prettygrid(grid, xedges, yedges, cmap='jet', propname='', czero=True,
               scale=1):
    """
    Generates pretty-looking 2D image maps for grid data using matplotlib.
    
    Parameters
    ----------
    grid : array-like object
        The grid of property values to plot
    xedges : list
        The min, max values of x associated with the grid.
    yedges : list
        The min, max values of y associated with the grid.
    cmap : str, optional
        The name of the matplotlib colormap to use.  Default value is 'jet'.
    propname : str, optional
        Name of the property being printed.  This is used as the plot title.
    czero : bool, optional
        If True, the range of property values will be centered around zero,
        i.e. cmax = -cmin.  If False, cmax and cmin will be independently
        selected using the property values.  Default value is True.
    scale : int, optional
        Scaling factor to apply to min/max property values.  Default value is
        1 (no scaling).
    
    Returns
    -------
    matplotlib.Figure
        The generated figure object.
    """
    
    # Get cmap from str name
    if isinstance(cmap, str):
        cmap = cm.get_cmap(cmap)
    
    # Set vmin = -vmax if czero is True
    if czero is True:
        vmax = abs(grid).max()
        if vmax != 0.0:
            vrounder = np.floor(np.log10(vmax))
            vmax = np.around(2 * vmax / 10.**vrounder) * 10.**vrounder / 2.
        else:
            vmax = 1e-15
        vmin = -vmax
    
    # Set vmin, vmax if czero is False
    else:
        vmax = grid.max()
        vmin = grid.min()

        if abs(grid).max() != 0.0:
            vrounder = np.floor(np.log10(grid.max() - grid.min()))

            vmax = np.around(2 * vmax / 10.**vrounder) * 10.**vrounder / 2.
            vmin = np.around(2 * vmin / 10.**vrounder) * 10.**vrounder / 2.

            if vmax == vmin:
                if vmax > 0:
                    vmin = 0
                else:
                    vmax = 0
        
        else:
            vmax = 1e-15
            vmin = -1e-15
    
    # Scale values if needed
    vmin *= scale
    vmax *= scale
    
    # Generate plot
    fig = plt.figure(figsize=(7.7, 7), dpi=72)
    plt.imshow(grid, extent = [xedges[0], xedges[-1], yedges[-1], yedges[0]], 
               cmap=cmap, norm=plt.Normalize(vmax=vmax, vmin=vmin))
    
    # Add pretty colorbar
    vticks = np.linspace(vmin, vmax, 11, endpoint=True)
    cbar = plt.colorbar(fraction=0.05, pad = 0.05, ticks=vticks)
    cbar.ax.tick_params(labelsize=15)
   
    # Make axis values and title pretty
    plt.xlim(xedges[0], xedges[-1])
    plt.xticks(size=15)
    plt.ylim(yedges[0], yedges[-1])
    plt.yticks(size=15)
    plt.title(propname, size=30)
    
    return fig
    
def adddots(x, y, xedges, yedges):
    """
    Overlays circles onto an active 2D plot to show actual atomic positions.
    
    Parameters
    ----------
    x : array-like object
        List of x-coordinates for the atoms/dots.
    y : array-like object
        List of y-coordinates for the atoms/dots.
    xedges : list
        The [xmin, xmax] values associated with the plot.
    yedges : list
        The [ymin, ymax] values associated with the plot.
    """
    # Get current plotting axis
    ax = plt.gca()
    
    # Set linewidth based on system dimensions
    syswidth = max([abs(xedges[-1]-xedges[0]), abs(yedges[-1]-yedges[0])])
    linewidth = 60. / syswidth
    
    # Add circles at each (x,y) position
    for xi, yi in zip(x, y):
        point = mpatches.Circle((xi, yi), 0.3, ls='solid', lw=linewidth)
        ax.add_artist(point)
        point.set_facecolor('none')
        point.set_edgecolor('k')