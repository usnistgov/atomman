# Standard Python imports
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from collections import OrderedDict
import os

# http://pandas.pydata.org/
import pandas as pd

# http://www.numpy.org/
import numpy as np

# https://www.scipy.org/
from scipy.interpolate import griddata, Rbf

# http://matplotlib.org/
import matplotlib.pyplot as plt

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# atomman imports
from .. import Box
from .. import unitconvert as uc
from ..compatibility import stringtype

class GammaSurface(object):
    """
    Class for representing gamma surfaces, i.e., generalized stacking faults.
    """
    
    def __init__(self, model=None, box=None, a1vect=None, a2vect=None,
                 a1=None, a2=None, E_gsf=None, delta=None):
        """
        Class initializer.
        
        Parameters
        ----------
        model : str, file-like object, DataModelDict
            XML/JSON data model containing the stacking fault information
        box : atomman.Box
            Defines box dimensions allowing conversion between crystal and
            Cartesian vectors.
        """
        if model is not None:
            try:
                assert box is None
                assert a1vect is None
                assert a2vect is None
                assert a1 is None
                assert a2 is None
                assert E_gsf is None
                assert delta is None
            except:
                raise TypeError('model cannot be given with any other parameter')
            else:
                self.model(model=model)
        
        else:
            self.a1vect = a1vect
            self.a2vect = a2vect
            self.box = box
            self.setdata(a1, a2, E_gsf, delta=delta)
    
    @property
    def data(self):
        """pandas.DataFrame : The raw data."""
        return self.__data
    
    @property
    def a1vect(self):
        """numpy.ndarray : Crystal vector for the a1 direction."""
        return self.__a1vect
    
    @a1vect.setter
    def a1vect(self, value):
        if isinstance(value, stringtype):
            value = value.split()
        value = np.asarray(value, dtype=float)
        if value.shape != (3,):
            raise ValueError('a1vect must be a 3D vector')
        self.__a1vect = value
    
    @property
    def a2vect(self):
        """numpy.ndarray : Crystal vector for the a1 direction."""
        return self.__a2vect
    
    @a2vect.setter
    def a2vect(self, value):
        if isinstance(value, stringtype):
            value = value.split()
        value = np.asarray(value, dtype=float)
        if value.shape != (3,):
            raise ValueError('a1vect must be a 3D vector')
        self.__a2vect = value
    
    @property
    def box(self):
        """
        atomman.Box : Defines box dimensions for conversion between crystal
                      and Cartesian vectors.
        """
        return self.__box
    
    @box.setter
    def box(self, value):
        if not isinstance(value, Box):
            raise TypeError('box must be an atomman.Box')
        self.__box = value
    
    def setdata(self, a1, a2, E_gsf, delta=None):
        """
        Sets stacking fault data.
        
        Parameters
        ----------
        a1vect : array-like object
            Crystal vector for the a1 direction.
        a2vect : array-like object
            Crystal vector for the a1 direction.
        a1 : array-like object
            List of fractional coordinates along a1vect corresponding to the
            E_gsf (and delta) values.
        a2 : array-like object
            List of fractional coordinates along a2vect corresponding to the
            E_gsf (and delta) values.
        E_gsf : array-like object
            List of generalized stacking fault energies for the positions
            associated with the corresponding (a1, a2) fractional coordinates.
        delta : array-like object, optional
            List of change in displacements normal to the fault plane for the
            positions associated with the corresponding (a1, a2) fractional
            coordinates.
        """
        data = OrderedDict()
        data['a1'] = a1
        data['a2'] = a2
        data['E_gsf'] = E_gsf
        if delta is not None:
            data['delta'] = delta
        self.__data = pd.DataFrame(data)
        self.fit()
    
    def fit(self):
        """
        Defines the interpolation functions from the raw data.
        """
        self.__E_gsf_fit = Rbf(self.data.a1, self.data.a2, self.data.E_gsf)
        if 'delta' in self.data:
            self.__delta_fit = Rbf(self.data.a1, self.data.a2, self.data.delta)
    
    def model(self, model=None, length_unit='angstrom',
              energyperarea_unit='eV/angstrom^2'):
        """
        Return or set DataModelDict representation of the gamma surface.
        
        Parameters
        ----------
        model : str, file-like object or DataModelDict, optional
            XML/JSON content to extract gamma surface energy from. If not
            given, model content will be generated.
        length_unit : str, optional
            Units to report delta displacement values in when a new model is
            generated. Default value is 'angstrom'.
        energyperarea_unit : str, optional
            Units to report fault energy values in when a new model is
            generated.  Default value is 'mJ/m^2'.
        
        Returns
        -------
        DataModelDict
            A dictionary containing the stacking fault data of the
            GammaSurface object.  Returned if model is not given.
        """
        # Set values if model given
        if model is not None:
            model = DM(model).find('stacking-fault-map')
            
            # Read in box, a1vect and a2vect
            self.box = Box(avect = model['box']['avect'],
                           bvect = model['box']['bvect'],
                           cvect = model['box']['cvect'])
            
            self.a1vect = model['shift-vector-1']
            self.a2vect = model['shift-vector-2']
            
            # Read in stacking fault data
            gsf = model.find('stacking-fault-relation')
            data = OrderedDict()
            data['a1'] = gsf['shift-vector-1-fraction']
            data['a2'] = gsf['shift-vector-2-fraction']
            data['E_gsf'] = uc.value_unit(gsf['energy'])
            try:
                data['delta'] = uc.value_unit(gsf['plane-separation'])
            except:
                pass
            self.__data = pd.DataFrame(data)
            self.fit()
        
        # Generate model
        else:
            model = DM()
            model['stacking-fault-map'] = sfm = DM()
            sfm['box'] = DM()
            sfm['box']['avect'] = list(self.box.avect)
            sfm['box']['bvect'] = list(self.box.bvect)
            sfm['box']['cvect'] = list(self.box.cvect)
            sfm['shift-vector-1'] = list(self.a1vect)
            sfm['shift-vector-2'] = list(self.a2vect)
            sfm['stacking-fault-relation'] = sfr = DM()
            sfr['shift-vector-1-fraction'] = list(self.__data['a1'])
            sfr['shift-vector-2-fraction'] = list(self.__data['a2'])
            sfr['energy'] = uc.model(self.__data['E_gsf'], energyperarea_unit)
            if 'delta' in self.data:
                sfr['plane-separation'] = uc.model(self.__data['delta'], length_unit)
            
            return model
    
    def a12_to_pos(self, a1, a2):
        """
        Conversion function from normalized a1, a2 coordinates to Cartesian
        positions.
        
        Parameters
        ----------
        a1 : float(s)
            Fractional distance(s) along a1 vector.
        a2 : float(s)
            Fractional distance(s) along a2 vector.
        
        Returns
        -------
        np.array
            3D Cartesian position vector(s).
        """
        # Extract data
        a1vect = self.a1vect
        a2vect = self.a2vect
        boxvects = self.box.vects
        
        # Convert a1vect and a2vect from crystal to Cartesian coordinates
        a1vect = np.dot(a1vect, boxvects)
        a2vect = np.dot(a2vect, boxvects)
        
        return np.outer(a1, a1vect) + np.outer(a2, a2vect)
        
    def pos_to_xy(self, pos, xvect=None):
        """
        Conversion function from Cartesian positions to plotting x, y
        coordinates.
        
        Parameters
        ----------
        pos: np.array
            3D Cartesian position vector(s).
        xvect : np.array, optional
            Cartesian vector corresponding to the plotting x-axis. If None (default), this is
            taken as a1vect.
        
        Returns
        -------
        x : float(s)
            Plotting x coordinate(s).
        y : float(s)
            Plotting y coordinate(s).
        """
        # Extract data
        a1vect = self.a1vect
        a2vect = self.a2vect
        
        if xvect is None:
            xvect = a1vect
        zvect = np.cross(a1vect, a2vect)
        if not np.isclose(np.dot(xvect, zvect), 0.0):
            raise ValueError('xvect must be in plane defined by a1vect and a2vect')
        yvect = np.cross(zvect, xvect)
        
        # Build transformation tensor
        transform = np.array([xvect, yvect, zvect])
        transform = (transform.T / np.linalg.norm(transform, axis=1)).T
        
        # Transform coords
        pos = transform.dot(pos.T).T
        
        # Save x, y to gammasurface.data
        return pos[...,0], pos[...,1]
    
    def a12_to_xy(self, a1, a2, xvect=None):
        """
        Conversion function from normalized a1, a2 coordinates to plotting x, y
        coordinates.
        
        Parameters
        ----------
        a1 : float(s)
            Fractional distance(s) along a1 vector.
        a2 : float(s)
            Fractional distance(s) along a2 vector.
        xvect : np.array, optional
            Cartesian vector corresponding to the plotting x-axis. If None (default), this is
            taken as a1vect.
        
        Returns
        -------
        x : float(s)
            Plotting x coordinate(s).
        y : float(s)
            Plotting y coordinate(s).
        """
        pos = self.a12_to_pos(a1, a2)
        return self.pos_to_xy(pos, xvect=xvect)
    
    def pos_to_a12(self, pos):
        """
        Conversion function from Cartesian positions to normalized a1, a2
        coordinates.
        
        Parameters
        ----------
        pos : np.array
            3D Cartesian position vector(s).
        
        Returns
        -------
        a1 : float(s)
            Fractional distance(s) along a1 vector.
        a2 : float(s)
            Fractional distance(s) along a2 vector.
        """
        # Extract data
        a1vect = self.a1vect
        a2vect = self.a2vect
        boxvects = self.box.vects
        
        # Convert a1vect and a2vect from crystal to Cartesian coordinates
        a1vect = np.dot(a1vect, boxvects)
        a2vect = np.dot(a2vect, boxvects)
        a3vect = np.cross(a1vect, a2vect)
        coeffs = np.array([np.array([a1vect, a2vect, a3vect]).T])
        
        a123 = np.linalg.solve(coeffs, pos)
        
        assert np.allclose(a123[...,2], 0.0, atol=1e-6), np.abs(a123[...,2]).max()
        return a123[...,0], a123[...,1]
    
    def xy_to_pos(self, x, y, xvect=None):
        """
        Conversion function from plotting x, y coordinates to Cartesian
        positions.
        
        Parameters
        ----------
        x : float(s)
            Plotting x coordinate(s).
        y : float(s)
            Plotting y coordinate(s).
        xvect : np.array, optional
            Cartesian vector corresponding to the plotting x-axis. If None
            (default), this is taken as a1vect.
        
        Returns
        -------
        pos: np.array
            3D Cartesian position vector(s).
        """
        # Extract data
        a1vect = self.a1vect
        a2vect = self.a2vect
        
        if xvect is None:
            xvect = a1vect
        zvect = np.cross(a1vect, a2vect)
        if not np.isclose(np.dot(xvect, zvect), 0.0):
            raise ValueError('xvect must be in plane defined by a1vect and a2vect')
        yvect = np.cross(zvect, xvect)
        
        pos = np.outer(x, [1,0,0]) + np.outer(y, [0,1,0])
        
        # Build transformation tensor
        transform = np.array([xvect, yvect, zvect])
        transform = (transform.T / np.linalg.norm(transform, axis=1)).T
        transform = np.linalg.inv(transform)
        
        # Transform coords
        pos = transform.dot(pos.T).T
        
        return pos
    
    def xy_to_a12(self, x, y, xvect=None):
        """
        Conversion function from plotting x, y coordinates to normalized a1, a2
        coordinates.
        
        Parameters
        ----------
        x : float(s)
            Plotting x coordinate(s).
        y : float(s)
            Plotting y coordinate(s).
        xvect : np.array, optional
            Cartesian vector corresponding to the plotting x-axis. If None
            (default), this is taken as a1vect.
        
        Returns
        -------
        a1 : float(s)
            Fractional distance(s) along a1 vector.
        a2 : float(s)
            Fractional distance(s) along a2 vector.
        """
        pos = self.xy_to_pos(x, y, xvect=xvect)
        return self.pos_to_a12(pos)
    
    def E_gsf(self, **kwargs):
        """
        Returns values for generalized stacking fault energy interpolated from
        the raw data.  Values can be obtained relative to a1, a2 normalized
        coordinates, x, y plotting coordinates, or pos Cartesian coordinates.
        
        Parameters
        ----------
        a1 : float(s), optional
            Fractional distance(s) along a1 vector.
        a2 : float(s), optional
            Fractional distance(s) along a2 vector.
        pos : np.array, optional
            3D Cartesian position vector(s).
        x : float(s), optional
            Plotting x coordinate(s).
        y : float(s), optional
            Plotting y coordinate(s).
        xvect : np.array, optional
            Cartesian vector corresponding to the plotting x-axis. If None
            (default), this is taken as a1vect.
        """
        # Convert x, y to a1, a2
        if 'x' in kwargs:
            x = kwargs.pop('x')
            y = kwargs.pop('y')
            xvect = kwargs.pop('xvect', None)
            assert len(kwargs) == 0, 'Unknown/incompatible arguments given'
            a1, a2 = self.xy_to_a12(x, y, xvect=xvect)
        
        # Convert pos to a1, a2
        elif 'pos' in kwargs:
            pos = kwargs.pop('pos')
            assert len(kwargs) == 0, 'Unknown/incompatible arguments given'
            a1, a2 = self.pos_to_a12(pos)
        
        # Get a1, a2 from kwargs
        else:
            a1 = np.array(kwargs.pop('a1'))
            a2 = np.array(kwargs.pop('a2'))
            assert len(kwargs) == 0, 'Unknown/incompatible arguments given'
        
        # Wrap all values within 0.0 < a1, a2 < 1.0
        while np.any(a1 > 1.0): 
            a1[a1 > 1.0] -= 1.0
        while np.any(a1 < 0.0): 
            a1[a1 < 0.0] += 1.0
        while np.any(a2 > 1.0): 
            a2[a2 > 1.0] -= 1.0
        while np.any(a2 < 0.0): 
            a2[a2 < 0.0] += 1.0
        
        return self.__E_gsf_fit(a1, a2)
    
    def delta(self, **kwargs):
        """
        Returns values for generalized stacking fault energy interpolated from
        the raw data.  Values can be obtained relative to a1, a2 normalized
        coordinates, x, y plotting coordinates, or pos Cartesian coordinates.
        
        Parameters
        ----------
        a1 : float(s), optional
                Fractional distance(s) along a1 vector.
        a2 : float(s), optional
            Fractional distance(s) along a2 vector.
        pos : np.array, optional
            3D Cartesian position vector(s).
        x : float(s), optional
            Plotting x coordinate(s).
        y : float(s), optional
            Plotting y coordinate(s).
        xvect : np.array, optional
            Cartesian vector corresponding to the plotting x-axis. If None
            (default), this is taken as a1vect.
        """
        if 'delta' not in self.data:
            raise ValueError('delta data not provided')
        
        # Convert x, y to a1, a2
        if 'x' in kwargs:
            x = kwargs.pop('x')
            y = kwargs.pop('y')
            xvect = kwargs.pop('xvect', None)
            assert len(kwargs) == 0, 'Unknown/incompatible arguments given'
            a1, a2 = self.xy_to_a12(x, y, xvect=xvect)
        
        # Convert pos to a1, a2
        elif 'pos' in kwargs:
            pos = kwargs.pop('pos')
            assert len(kwargs) == 0, 'Unknown/incompatible arguments given'
            a1, a2 = self.pos_to_a12(pos)
        
        # Get a1, a2 from kwargs
        else:
            a1 = np.array(kwargs.pop('a1'))
            a2 = np.array(kwargs.pop('a2'))
            assert len(kwargs) == 0, 'Unknown/incompatible arguments given'
        
        # Wrap all values within 0.0 < a1, a2 < 1.0
        while np.any(a1 > 1.0):
            a1[a1 > 1.0] -= 1.0
        while np.any(a1 < 0.0):
            a1[a1 < 0.0] += 1.0
        while np.any(a2 > 1.0):
            a2[a2 > 1.0] -= 1.0
        while np.any(a2 < 0.0):
            a2[a2 < 0.0] += 1.0
        
        return self.__delta_fit(a1, a2)
    
    def E_gsf_surface_plot(self, normalize=False, smooth=True, xvect=None,
                         length_unit='angstrom', energyperarea_unit='eV/angstrom^2',
                         numx=100, numy=100, xsize=None, ysize=None,
                         cmap=None):
        """
        Creates a 2D surface plot from the stacking fault energy values.
        
        Parameters
        ----------
        normalize : bool, optional
            Flag indicating if axes are Cartesian (False, default) or
            normalized by a1, a2 vectors (True).
        smooth : bool, optional
            If True (default), then plot shows smooth interpolated values.
            If False, plot shows nearest raw data values.
        xvect : numpy.array, optional
            Crystal vector to align with the plotting x-axis for 
            non-normalized plots.  If not given, this is taken as the object's
            a1vect.
        length_unit : str, optional
            The unit of length to display non-normalized axes values in.
            Default value is 'angstrom'.
        energyperarea_unit : str, optional
            The unit of energy per area to display the stacking fault energies
            in. Default value is 'eV/angstrom^2'.
        numx : int, optional
            The number of plotting points to use along the x-axis.  Default
            value is 100.
        numy : int, optional
            The number of plotting points to use along the y-axis.  Default
            value is 100.
        xsize : float or None, optional
            Base size in inches for the plot's x dimensions.  If None
            (default), the value is scaled relative to ysize based on the axes
            vectors.  If both size parameters are None, the larger of the two
            is set to 10.
        ysize : float or None, optional
            Base size in inches for the plot's y dimensions.  If None
            (default), the value is scaled relative to xsize based on the axes
            vectors.  If both size parameters are None, the larger of the two
            is set to 10.
        cmap : str
            matplotlib colormap name to use.
            
        Returns
        -------
        matplotlib.figure
        """
        # Generate grids of a1, a2 values from numx, numy
        x_grid, y_grid = np.meshgrid(np.linspace(0, 1, numx),
                                     np.linspace(0, 1, numy))
        
        # Generate grid of values either with or without interpolation
        if smooth is True:
            C = self.E_gsf(a1=x_grid, a2=y_grid)
        else:
            C = griddata(self.data[['a1', 'a2']].values, self.data.E_gsf.values,
                         (x_grid, y_grid), method='nearest')
        
        # Convert units of C using energyperarea_unit
        C = uc.get_in_units(C, energyperarea_unit)
        
        # Set parameters for normalized plots
        if normalize is True:
            yscale = 1
            #yscale = np.linalg.norm(self.a2vect) / np.linalg.norm(self.a1vect)
            xlabel = '$a_1$ = ' + str(self.a1vect)
            ylabel = '$a_2$ = ' + str(self.a2vect)
        
        # Set parameters for absolute plots
        else:
            shape = x_grid.shape
            x_grid, y_grid = self.a12_to_xy(x_grid.flatten(), y_grid.flatten(),
                                            xvect=xvect)
            x_grid.shape = shape
            y_grid.shape = shape
            x_grid = uc.get_in_units(x_grid, length_unit)
            y_grid = uc.get_in_units(y_grid, length_unit)
            yscale = (y_grid.max()-y_grid.min()) / (x_grid.max() - x_grid.min())
            xlabel = 'x (' + length_unit + ')'
            ylabel = 'y (' + length_unit + ')'
        
        # Set default xsize or ysize if needed
        if xsize is None and ysize is None:
            if yscale < 1:
                xsize = 10
            else:
                ysize = 10
        
        # Set figsize
        # xscale accounts for colorbar
        # yscale makes x,y values correctly spaced
        xscale = 1.175
        if xsize is None:
            figsize = (xscale * ysize / yscale, ysize)
        elif ysize is None:
            figsize = (xscale * xsize, xsize * yscale)
        else:
            figsize = (xscale * xsize, ysize)
        
        # Generate plot
        fig = plt.figure(figsize=figsize)
        plt.pcolormesh(x_grid, y_grid, C, cmap=cmap)
        plt.xlabel(xlabel, fontsize='x-large')
        plt.ylabel(ylabel, fontsize='x-large')
        cbar = plt.colorbar(aspect=40, fraction=0.1)
        cbar.ax.set_ylabel('$E_{gsf}$ (' + energyperarea_unit + ')',
                           fontsize='x-large')
        
        return fig
    
    def E_gsf_line_plot(self, vect=None, num=100, xsize=10, ysize=6,
                      length_unit='angstrom', energyperarea_unit='eV/angstrom^2'):
        """
        Generates a line plot for the interpolated generalized stacking fault
        energy along a specified crystallographic vector in the (a1, a2) plane.
        
        Parameters
        ----------
        vect : numpy.array, optional
            Crystallographic vector to plot the gsf along.  Must be in the
            plane defined by the GammaSurface object's a1vect and a2vect
            vectors.
        num : int, optional
            The number of points to evaluate the generalized stacking fault
            energy for.  Default value is 100.
        xsize : float
            The size in inches to make the figure's x dimensions.  Default
            value is 10.
        ysize : float
            The size in inches to make the figure's y dimensions.  Default
            value is 6.
        length_unit : str, optional
            The unit of length to display the x-axis coordinates in.
            Default value is 'angstrom'.
        energyperarea_unit : str, optional
            The unit of energy per area to display the stacking fault energies
            in. Default value is 'eV/angstrom^2'.
        
        Returns
        -------
        matplotlib.figure
        """
        # Extract data
        a1vect = self.a1vect
        a2vect = self.a2vect
        boxvects = self.box.vects
        
        # Assign default vect if needed
        if vect is None:
            vect = a1vect
        vect = np.asarray(vect)
        
        # Convert vect, a1vect and a2vect from crystal to Cartesian coordinates
        svect = np.dot(vect, boxvects)
        a1vect = np.dot(a1vect, boxvects)
        a2vect = np.dot(a2vect, boxvects)
        
        # Assert vect is in (a1,a2) plane: vect * (a1vect x a2vect) == [0,0,0]
        assert np.allclose(np.dot(svect, np.cross(a1vect, a2vect)), 
                           np.zeros(3))
        
        # Generate intermediate pos points
        pos = np.empty((num, 3))
        pos[:, 0] = np.linspace(0, svect[0], num)
        pos[:, 1] = np.linspace(0, svect[1], num)
        pos[:, 2] = np.linspace(0, svect[2], num)
        
        # Evaluate interpolated energy and distance along x
        E = uc.get_in_units(self.E_gsf(pos=pos), energyperarea_unit)
        x = uc.get_in_units(np.linalg.norm(pos, axis=1), length_unit)
        
        # Create plot
        fig = plt.figure(figsize=(xsize, ysize))
        plt.plot(x, E, '-k')
        plt.xlabel('$x$ along ' + str(vect) + ' (' + length_unit + ')',
                   fontsize='x-large')
        plt.ylabel('$E_{gsf}$ (' + energyperarea_unit + ')', fontsize='x-large')
        plt.xlim(0, x.max())
        plt.ylim(E.min(), None)
        
        return fig
    
    def delta_surface_plot(self, normalize=False, smooth=True, xvect=None,
                         length_unit='angstrom',
                         numx=100, numy=100, xsize=None, ysize=None,
                         cmap=None):
        """
        Creates a 2D surface plot from the delta planar displacement values.
        
        Parameters
        ----------
        normalize : bool, optional
            Flag indicating if axes are Cartesian (False, default) or
            normalized by a1, a2 vectors (True).
        smooth : bool, optional
            If True (default), then plot shows smooth interpolated values.
            If False, plot shows nearest raw data values.
        xvect : numpy.array, optional
            Crystal vector to align with the plotting x-axis for 
            non-normalized plots.  If not given, this is taken as the object's
            a1vect.
        length_unit : str, optional
            The unit of length to display non-normalized axes values and delta
            displacements in.  Default value is 'angstrom'.
        numx : int, optional
            The number of plotting points to use along the x-axis.  Default
            value is 100.
        numy : int, optional
            The number of plotting points to use along the y-axis.  Default
            value is 100.
        xsize : float or None, optional
            Base size in inches for the plot's x dimensions.  If None
            (default), the value is scaled relative to ysize based on the axes
            vectors.  If both size parameters are None, the larger of the two
            is set to 10.
        ysize : float or None, optional
            Base size in inches for the plot's y dimensions.  If None
            (default), the value is scaled relative to xsize based on the axes
            vectors.  If both size parameters are None, the larger of the two
            is set to 10.
        cmap : str
            matplotlib colormap name to use.
        
        Returns
        -------
        matplotlib.figure
        """
        if 'delta' not in self.data:
            raise ValueError('delta data not provided')
        
        # Generate grids of a1, a2 values from numx, numy
        x_grid, y_grid = np.meshgrid(np.linspace(0, 1, numx),
                                     np.linspace(0, 1, numy))
        
        # Generate grid of values either with or without interpolation
        if smooth is True:
            C = self.delta(a1=x_grid, a2=y_grid)
        else:
            C = griddata(self.data[['a1', 'a2']].values, self.data.delta.values,
                         (x_grid, y_grid), method='nearest')
        
        # Convert units of C using length_unit
        C = uc.get_in_units(C, length_unit)
        
        # Set parameters for normalized plots
        if normalize is True:
            yscale = 1
            #yscale = np.linalg.norm(self.a2vect) / np.linalg.norm(self.a1vect)
            xlabel = '$a_1$ = ' + str(self.a1vect)
            ylabel = '$a_2$ = ' + str(self.a2vect)
        
        # Set parameters for absolute plots
        else:
            shape = x_grid.shape
            x_grid, y_grid = self.a12_to_xy(x_grid.flatten(), y_grid.flatten(),
                                            xvect=xvect)
            x_grid.shape = shape
            y_grid.shape = shape
            x_grid = uc.get_in_units(x_grid, length_unit)
            y_grid = uc.get_in_units(y_grid, length_unit)
            yscale = (y_grid.max()-y_grid.min()) / (x_grid.max() - x_grid.min())
            xlabel = 'x (' + length_unit + ')'
            ylabel = 'y (' + length_unit + ')'
        
        # Set default xsize or ysize if needed
        if xsize is None and ysize is None:
            if yscale < 1:
                xsize = 10
            else:
                ysize = 10
        
        # Set figsize
        # xscale accounts for colorbar
        # yscale makes x,y values correctly spaced
        xscale = 1.175
        if xsize is None:
            figsize = (xscale * ysize / yscale, ysize)
        elif ysize is None:
            figsize = (xscale * xsize, xsize * yscale)
        else:
            figsize = (xscale * xsize, ysize)
        
        # Generate plot
        fig = plt.figure(figsize=figsize)
        plt.pcolormesh(x_grid, y_grid, C, cmap=cmap)
        plt.xlabel(xlabel, fontsize='x-large')
        plt.ylabel(ylabel, fontsize='x-large')
        cbar = plt.colorbar(aspect=40, fraction=0.1)
        cbar.ax.set_ylabel('$\delta_{gsf}$ (' + length_unit + ')',
                           fontsize='x-large')
        
        return fig
    
    def delta_line_plot(self, vect=None, num=100, xsize=10, ysize=6,
                      length_unit='angstrom'):
        """
        Generates a line plot for the interpolated delta planar shift values
        along a specified crystallographic vector in the (a1, a2) plane.
        
        Parameters
        ----------
        vect : numpy.array, optional
            Crystallographic vector to plot the delta shift along.  Must be in
            the plane defined by the GammaSurface object's a1vect and a2vect
            vectors.
        num : int, optional
            The number of points to evaluate the delta planar shift for.
            Default value is 100.
        xsize : float
            The size in inches to make the figure's x dimensions.  Default
            value is 10.
        ysize : float
            The size in inches to make the figure's x dimensions.  Default
            value is 6.
        length_unit : str, optional
            The unit of length to display non-normalized axes values and delta
            displacements in.  Default value is 'angstrom'.
        
        Returns
        -------
        matplotlib.figure
        """
        if 'delta' not in self.data:
            raise ValueError('delta data not provided')
        
        # Extract data
        a1vect = self.a1vect
        a2vect = self.a2vect
        boxvects = self.box.vects
        
        # Assign default vect if needed
        if vect is None:
            vect = a1vect
        vect = np.asarray(vect)
        
        # Convert vect, a1vect and a2vect from crystal to Cartesian coordinates
        svect = np.dot(vect, boxvects)
        a1vect = np.dot(a1vect, boxvects)
        a2vect = np.dot(a2vect, boxvects)
        
        # Assert vect is in (a1,a2) plane: vect * (a1vect x a2vect) == [0,0,0]
        assert np.allclose(np.dot(svect, np.cross(a1vect, a2vect)), 
                           np.zeros(3))
        
        # Generate intermediate pos points
        pos = np.empty((num, 3))
        pos[:, 0] = np.linspace(0, svect[0], num)
        pos[:, 1] = np.linspace(0, svect[1], num)
        pos[:, 2] = np.linspace(0, svect[2], num)
        
        # Evaluate interpolated delta shifts and distance along x
        E = uc.get_in_units(self.delta(pos=pos), length_unit)
        x = uc.get_in_units(np.linalg.norm(pos, axis=1), length_unit)
        
        # Create plot
        fig = plt.figure(figsize=(xsize, ysize))
        plt.plot(x, E, '-k')
        plt.xlabel('$x$ along ' + str(vect) + ' (' + length_unit + ')',
                   fontsize='x-large')
        plt.ylabel('$\delta_{gsf}$ (' + length_unit + ')', fontsize='x-large')
        plt.xlim(0, x.max())
        plt.ylim(E.min(), None)
        
        return fig