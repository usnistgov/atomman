# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from collections import OrderedDict
from copy import deepcopy

# http://www.numpy.org/
import numpy as np

# https://pandas.pydata.org/
import pandas as pd

# atomman imports
from . import Atoms, Box, dvect, NeighborList
from ..lammps import normalize as lmp_normalize
from ..compatibility import iteritems, range, inttype, stringtype
from ..tools import indexstr, miller
from .. import dump

class System(object):
    """
    A representation of an atomic system.  This combines Atoms with Box and
    adds methods and attributes involving both.
    """
    
    def __init__(self, atoms=Atoms(), box=Box(), pbc=(True, True, True),
                 scale=False, symbols=()):
        """
        Initialize a System by joining an am.Atoms and am.Box instance.
        
        Parameters
        ----------
        atoms : atomman.Atoms, optional
            The underlying Atoms object to build system around.
        box : atomman.Box, optional
            The underlying box object to build system around.
        pbc : tuple or list of bool, optional
            Indicates which of the dimensions related to the three box vectors
            are periodic.  Default value is (True, True, True).
        scale : bool, optional
            If True, atoms.pos will be scaled relative to the box.  Default
            value is False.
        symbols : tuple, optional
            A list of the element symbols for each atom atype.  If len(symbols)
            is less than natypes, then missing values will be set to None.
            Default value is (,), i.e. all values set to None.
        """
        
        # Check data types
        if not isinstance(atoms, Atoms):
            raise TypeError('Invalid atoms type')
        if not isinstance(box, Box):
            raise TypeError('Invalid box type')
        if not isinstance(scale, bool):
            raise TypeError('Invalid scale type')
        
        # Set properties
        self.__atoms = atoms
        self.__box = box
        self.pbc = pbc
        self.__transformation = np.identity(3)
        
        if isinstance(symbols, stringtype):
            symbols = (symbols,)
        assert len(symbols) <= self.natypes
        self.__symbols = tuple(symbols)
        
        # Scale pos if needed
        if scale is True:
            self.atoms_prop('pos', value=atoms.pos, scale=True)
    
    def __str__(self):
        """str : The string representation of a system."""
        return '\n'.join([str(self.box),
                          'natoms = ' + str(self.natoms),
                          'natypes = ' + str(self.natypes),
                          'symbols = ' + str(self.symbols),
                          'pbc = ' + str(self.pbc),
                          str(self.atoms)])
    
    def __len__(self):
        return self.natoms
    
    @property
    def atoms(self):
        """atomman.Atoms : underlying Atoms object."""
        return self.__atoms
    
    @property
    def natoms(self):
        """int : The number of atoms in the Atoms class."""
        return self.__atoms.natoms
    
    @property
    def atypes(self):
        """tuple : List of unique int atom types."""
        return self.__atoms.atypes
    
    @property
    def natypes(self):
        """int : The number of atom types in the Atoms class."""
        return self.__atoms.natypes
    
    @property
    def box(self):
        """atommman.Box : underlying Box object."""
        return self.__box
    
    @property
    def pbc(self):
        """list of bool : The periodic boundary condition settings."""
        return self.__pbc
        
    @pbc.setter
    def pbc(self, value):
        pbc = np.asarray(value, dtype=bool)
        assert pbc.shape == (3,), 'invalid pbc entry' 
        self.__pbc = pbc
    
    @property
    def symbols(self):
        """tuple : The element model symbols associated with each atype."""
        
        # Check if there are symbols for all atypes
        if len(self.__symbols) != self.natypes:
            symbols = [None for x in range(self.natypes)]
            for i in range(len(self.__symbols)):
                symbols[i] = self.__symbols
            self.__symbols = tuple(symbols)
        
        return self.__symbols
    
    @symbols.setter
    def symbols(self, value):
        if isinstance(value, stringtype):
            value = (value,)
        assert len(value) == self.natypes, 'length of symbols does not match natypes'
        self.__symbols = tuple(value)
    
    def atoms_prop(self, key=None, index=None, value=None, a_id=None, scale=False):
        """
        Extends Atoms.prop() by adding a scale argument.
        
        Parameters
        ----------
        key : str, optional
            Per-atom property name.
        index : int, list, slice, optional
            Index of atoms.
        value : any, optional
            Property values to assign.
        a_id : int, optional
            Integer atom index.  Left in for backwards compatibility.
        scale : bool, optional
            Flag indicating if values should be scaled/unscaled. If True:
            - Values being retrieved are scaled from absolute Cartesian to box relative vectors.
            - Values being set are unscaled from box relative to absolute Cartesian vectors.
            Default value is False.
        """
        
        # Check that scale is bool
        if not isinstance(scale, bool):
            raise TypeError('Invalid scale type')
        
        # Call atoms.prop() for scale is False
        if scale is False:
            if value is None:
                return self.atoms.prop(key=key, index=index, a_id=a_id)
            else:
                self.atoms.prop(key=key, index=index, value=value, a_id=a_id)
        
        # Handle property scaling if scale is True
        else:
            # Handle a_id
            if a_id is not None:
                if index is not None:
                    raise ValueError('a_id and index cannot both be given')
                index = a_id
            
            # Get values if value is None
            if value is None:
                # If no key, return copy of atoms (all or a slice) with scaled pos
                if key is None:
                    if index is None:
                        newatoms = deepcopy(self.atoms)
                    else:
                        newatoms = deepcopy(self.atoms[index])
                    newatoms.pos = self.scale(newatoms.pos)
                    return newatoms
                
                # If key is given, return scaled property values
                else:
                    if index is None:
                        value = self.atoms.view[key]
                    else:
                        value = self.atoms.view[key][index]
                    return self.scale(value)
            
            # Set values if value is given
            else:
                # If no key, unscale pos of Atoms value and set to atoms (or a slice)
                if key is None:
                    if not isinstance(value, Atoms):
                        raise TypeError('If key is None, value must be instance of atomman.Atoms')
                    
                    value.pos = self.unscale(value.pos)
                    if index is None:
                        self.atoms[:] = value
                    else:
                        self.atoms[index] = value
                
                # If key is given, unscale value and set to property
                else:
                    value = self.unscale(value)
                    if index is None:
                        self.atoms.view[key] = value
                    else:
                        self.atoms.view[key][index] = value
    
    def atoms_df(self, scale=False):
        """
        Extends Atoms.df() by adding a scale argument.
        
        Parameters
        ----------
        scale : bool or list, optional
            Indicates if/which per-atom properties are to be scaled to box
            relative values.  If False (default), no properties will be scaled.
            If True, pos will be scaled.  Multiple properties can be scaled by
            providing a list of the property names to scale.
            
        Returns
        -------
        pandas.DataFrame
            Tabulated DataFrame version of the atomic data.
        """
        # Handle scale values
        if scale is True:
            scale = ['pos']
        elif scale is False:
            scale = []
        elif not isinstance(scale, list):
            scale = [scale]
        
        # Initialize new dictionary of values
        values = OrderedDict()
        for key in self.atoms.view.keys():
            
            value = self.atoms.view[key]
            if key in scale:
                value = self.scale(value)
            
            # Flatten multidimensional arrays
            for index, istr in indexstr(self.atoms.view[key].shape[1:]):
                newkey = key + istr
                
                # Copy values over
                if index == ():
                    values[newkey] = value
                else:
                    values[newkey] = value[(Ellipsis, ) + index]
        
        # Return DataFrame
        return pd.DataFrame(values)
    
    def box_set(self, **kwargs):
        """
        Extends box.set() with a scale argument.
        
        Parameters
        ----------
        scale : bool
            If True, the scaled (box-relative) positions remain unchanged.
            Default value is False (absolute positions unchanged).
        other kwargs : any
            Any other kwargs allowed by box.set().
        """
        
        # Pop scale
        scale = kwargs.pop('scale', False)
        if not isinstance(scale, bool):
            raise TypeError('Invalid scale type')
        
        # Hold scaled positions constant
        if scale is True:
            spos = self.atoms_prop('pos', scale=True)
            self.box.set(**kwargs)
            self.atoms_prop('pos', value=spos, scale=True)
        
        # Call box.set without scaling
        else:
            self.box.set(**kwargs)
    
    def scale(self, value):
        """
        Scales 3D vectors from absolute Cartesian coordinates to relative box
        coordinates.
        
        Parameters
        ----------
        value : numpy.ndarray
            Values to scale.
        """
        
        # Retrieve parameters
        value = np.asarray(value, dtype=float)
        vects = self.box.vects
        inverse = np.linalg.inv(vects)
        origin = self.box.origin
        
        # Convert
        return (value - origin).dot(inverse)
          
    def unscale(self, value):
        """
        Unscales 3D vectors from relative box coordinates to absolute
        Cartesian coordinates.
        
        Parameters
        ----------
        value : numpy.ndarray
            Values to unscale.
        """
        
        # Retrieve parameters
        value = np.asarray(value, dtype=float)
        vects = self.box.vects
        origin = self.box.origin
        
        # Convert
        return value.dot(vects) + origin
        
    def wrap(self):
        """
        Wrap atoms around periodic boundaries and extend non-periodic
        boundaries such that all atoms are within the box.
        """
        
        # mins and maxs are box dimensions relative to box vectors, i.e 0 to 1
        mins = np.array([0.0, 0.0, 0.0])
        maxs = np.array([1.0, 1.0, 1.0])
        
        # Retrieve scaled pos
        spos = self.atoms_prop('pos', scale=True)
        
        # Loop over three pbc directions
        for i in range(3):
            
            # Wrap atoms across periodic boundaries
            if self.pbc[i]:
                spos[:, i] -= np.floor(spos[:, i])
            
            # Shift min and max to encompass atoms across non-periodic bounds
            else:
                min = spos[:, i].min()
                max = spos[:, i].max()
                if min <= mins[i]: mins[i] = min - 0.001
                if max >= maxs[i]: maxs[i] = max + 0.001
        
        # Unscale spos and save to pos
        self.atoms_prop('pos', value=spos, scale=True)
        
        # Modify box vectors and origin by new min and max
        origin = self.box.origin + mins.dot(self.box.vects) 
        avect = self.box.avect * (maxs[0] - mins[0])
        bvect = self.box.bvect * (maxs[1] - mins[1])
        cvect = self.box.cvect * (maxs[2] - mins[2])
        self.box_set(avect=avect, bvect=bvect, cvect=cvect, origin=origin)
    
    def dvect(self, pos_0, pos_1, code=None):
        """
        Computes the shortest distance between pos_0 and pos_1 using box 
        dimensions and accounting for periodic boundaries.
        
        Parameters
        ----------
        pos_0 : numpy.ndarray or index
            Absolute Cartesian vector position(s) to use as reference point(s).
            If the value can be used as an index, then self.atoms.pos[pos_0]
            is taken.
        pos_1 : numpy.ndarray or index
            Absolute Cartesian vector position(s) to find relative to pos_0.
            If the value can be used as an index, then self.atoms.pos[pos_1]
            is taken.
        code: str, optional
            Option for specifying which underlying code function to use:
            - 'cython' uses the version of the function built in cython (faster).
            - 'python' uses the purely python version.
            Default is 'cython' if the code can be imported, otherwise 'python'.
        """
        # Test if pos_0 and pos_1 can be used as numpy array indices
        try:
            pos_0 = self.atoms.pos[pos_0]
        except:
            pos_0 = np.asarray(pos_0)
        try:
            pos_1 = self.atoms.pos[pos_1]
        except:
            pos_1 = np.asarray(pos_1)
        
        # Call dvect using self's box and pbc
        vects = dvect(pos_0, pos_1, self.box, self.pbc, code=code)
        if len(vects) == 1:
            return vects[0]
        else:
            return vects
    
    def neighborlist(self, **kwargs):
        """
        Builds a neighbor list for the system.  The resulting NeighborList
        object is saved to the object as attribute 'neighbors'.
        
        Parameters
        ----------
        cutoff : float, optional
            Radial cutoff distance for identifying neighbors.  Must be given if
            model is not given.
        model : str or file-like object, optional
            Gives the file path or content to load.  If given, no other
            parameters are allowed.
        cmult : int, optional
            Parameter associated with the binning routine.  Default value is most
            likely the fastest.
        code: str, optional
            Option for specifying which underlying code function of nlist to use:
            - 'cython' uses the version of the function built in cython (faster).
            - 'python' uses the purely python version.
            Default is 'cython' if the code can be imported, otherwise 'python'.
        initialsize : int, optional
            The number of neighbor positions to initially assign to each atom.
            Default value is 20.
        deltasize : int, optional
            Specifies the number of extra neighbor positions to allow each atom
            when the number of neighbors exceeds the underlying array size.
            Default value is 10.
            
        Returns
        -------
        atomman.NeighborList
            The compiled list of neighbors.
        """
        if 'system' in kwargs:
            raise KeyError("Parameter 'system' not allowed")
        else:
            kwargs['system'] = self
        return NeighborList(**kwargs)
    
    def supersize(self, a_size, b_size, c_size):
        """
        Creates a larger system from a given system by replicating it along the
        system's box vectors.
        
        The multiplier values \\*_size are taken to be integer tuples (m, n) where
        m <= 0 and n >= 0.  The system multiplication works such that if n = -m,
        then the seed system's origin will be at the center of the new system.
        If only one integer is given, then it is assigned to m or n depending on
        its sign, and the other value is taken to be 0.
        
        Parameters
        ----------
        a_size : int or tuple of int
            Single int or two integers specifying replication along the avect
            direction.
        b_size -- int or tuple of int
            Single int or two integers specifying replication along the bvect
            direction.
        c_size -- int or tuple of int
            Single int or two integers specifying replication along the cvect
            direction.
        
        Returns
        -------
        atomman.System
            A new system created by replicating the given seed system according to
            the \\*_size parameters.
          
        """
        # Extract parameters
        sizes = [a_size, b_size, c_size]
        mults = np.array([0, 0, 0], dtype=int)
        vects = self.box.vects
        origin = self.box.origin
        spos = self.atoms_prop('pos', scale=True)
        
        # Check the *_size values
        for i in range(3):
            
            # Change single int to tuple of two int
            if isinstance(sizes[i], inttype):
                if sizes[i] > 0:
                    sizes[i] = (0, sizes[i])
                elif sizes[i] < 0:
                    sizes[i] = (sizes[i], 0)
            
            elif isinstance(sizes[i], tuple):
                if True:
                    assert len(sizes[i]) == 2, str(len(sizes[i]))
                    assert isinstance(sizes[i][0], inttype), str(sizes[i][0])
                    assert sizes[i][0] <= 0, str(sizes[i][0])
                    assert isinstance(sizes[i][1], inttype), str(sizes[i][1])
                    assert sizes[i][1] >= 0, str(sizes[i][1])
                else:
                    raise TypeError('Invalid system multipliers')
            else:
                raise TypeError('Invalid system multipliers')
            
            # Calculate full multipliers 
            mults[i] = sizes[i][1] - sizes[i][0]
            if mults[i] == 0:
                raise ValueError('Cannot multiply system dimension by zero')
                
            # Scale box and first set of positions accordingly
            spos[:,i] /= mults[i]
            origin += vects[i] * sizes[i][0]
            vects[i] *= mults[i]
            
        # Initialize new Box and Atoms
        box = Box(vects=vects, origin=origin)
        natoms = self.natoms * mults[0] * mults[1] * mults[2]
        atoms = Atoms(natoms=natoms)
        
        # Copy over all property values (except pos) using numpy broadcasting
        for key in self.atoms_prop():
            if key == 'pos':
                continue
        
            # Get old array
            old = self.atoms.view[key]
            
            # Create new array and broadcast old to it
            new = np.empty((mults[0] * mults[1] * mults[2],) + old.shape, dtype = old.dtype)
            new[:] = old
            
            # Reshape new and save to atoms
            new_shape = new.shape
            new_shape = (new_shape[0] * new_shape[1],) + new_shape[2:]
            atoms.view[key] = np.array(new.reshape(new_shape))
        
        # Expand spos using broadcasting
        new_spos = np.empty((mults[0] * mults[1] * mults[2],) + spos.shape)
        new_spos[:] = spos
        
        # Reshape spos
        new_shape = new_spos.shape
        new_shape = (new_shape[0]*new_shape[1],) + new_shape[2:]
        new_spos = new_spos.reshape(new_shape)
        
        # Use broadcasting to create arrays to add to spos
        test = np.empty(mults[0] * self.natoms)
        test.shape = (self.natoms, mults[0])
        test[:] = np.arange(mults[0])
        x = test.T.flatten()

        test = np.empty(mults[1] * len(x))
        test.shape = (len(x), mults[1])
        test[:] = np.arange(mults[1])
        y = test.T.flatten()
        test.shape = (mults[1], len(x))
        test[:] = x
        x = test.flatten()

        test = np.empty(mults[2] * len(x))
        test.shape = (len(x), mults[2])
        test[:] = np.arange(mults[2])
        z = test.T.flatten()
        test.shape = (mults[2], len(x))
        test[:] = x
        x = test.flatten()
        test[:] = y
        y = test.flatten()
        
        # xyz is displacement values to add to spos
        xyz = (np.hstack((x[:, np.newaxis], y[:, np.newaxis], z[:, np.newaxis])) 
             * np.array([1 / mults[0], 1 / mults[1], 1 / mults[2]]))
        
        # Save pos values, return new System
        atoms.view['pos'] = new_spos + xyz
        
        return System(box=box, atoms=atoms, scale=True, symbols=self.symbols)
    
    def rotate(self, uvws, tol=1e-5, return_transform=False):
        """
        Transforms a System representing a periodic crystal cell from a standard
        orientation to a specified orientation. Note: if hexagonal indices are
        given, the vectors will be reduced to the smallest uvw integer
        representation.
        
        Parameters
        ----------
        uvws : numpy.ndarray of int
            A (3, 3) array of the Miller crystal vectors or a (3, 4) array of
            Miller-Bravais hexagonal crystal vectors to use in transforming the
            system.
        tol : float, optional
            Tolerance parameter used in rounding atomic positions near the
            boundaries to the boundary values.  In box-relative coordinates, any
            atomic positions within tol of 0 or 1 will be rounded to 0 or 1,
            respectively.  Default value is 1e-5.
        return_transform : bool, optional
            Indicates if the transformation matrix associated with the
            rotation is returned.  Default value is False.
        
        Returns
        -------
        atomman.System
            A new fully periodic system rotated and transformed according to the
            uvws crystal vectors.
        transform : np.ndarray
            The transformation matrix associated with the rotation.
            Returned if return_transform is True.
        """
        
        # Check parameters
        try:
            uvws = np.asarray(uvws, dtype='int64')
           
            if uvws.shape == (3, 4):
                uvws = miller.vector4to3(uvws)
            assert uvws.shape == (3, 3)
        except:
            raise ValueError('Invalid uvws crystal indices')
        
        # Get natoms and volume of system
        natoms = self.natoms
        volume = np.abs(self.box.avect.dot(np.cross(self.box.bvect,
                                                    self.box.cvect)))
        
        # Convert uvws to Cartesian units and compute new volume and natoms
        newvects = miller.vectortocartesian(uvws, box=self.box)
        newvolume = np.abs(newvects[0].dot(np.cross(newvects[1], newvects[2])))
        newnatoms = int(round(newvolume / volume) * natoms)
        
        # Check new values
        if newnatoms == 0:
            raise ValueError('New box has no atoms/volume: vectors are parallel or planar')
        
        # Identify box corners of new system wrt uvws
        corners = np.empty((8,3), dtype='int64')
        corners[0] = np.zeros(3)
        corners[1] = uvws[0]
        corners[2] = uvws[1]
        corners[3] = uvws[2]
        corners[4] = uvws[0] + uvws[1]
        corners[5] = uvws[0] + uvws[2]
        corners[6] = uvws[1] + uvws[2]
        corners[7] = uvws[0] + uvws[1] + uvws[2]
        
        # Create a supercell of system that contains all box corners
        a_mults = (corners[:,0].min()-1, corners[:,0].max()+1)
        b_mults = (corners[:,1].min()-1, corners[:,1].max()+1)
        c_mults = (corners[:,2].min()-1, corners[:,2].max()+1)
        system2 = self.supersize(a_mults, b_mults, c_mults)
        
        # Change system.box.vects to newvects
        system2.box_set(vects=newvects, scale=False)
        
        # Round atom positions near box boundaries to the boundaries
        spos = system2.atoms_prop('pos', scale=True)
        spos[np.isclose(spos, 0.0, atol=tol)] = 0.0
        spos[np.isclose(spos, 1.0, atol=tol)] = 1.0
        
        # Identify all atoms whose positions are 0 <= x < 1
        aindex = np.where(((spos[:, 0] >= 0.0) & (spos[:, 0] < 1.0)
                         & (spos[:, 1] >= 0.0) & (spos[:, 1] < 1.0)
                         & (spos[:, 2] >= 0.0) & (spos[:, 2] < 1.0)))
        if len(aindex[0]) != newnatoms:
            raise ValueError('Filtering failed: ' + str(newnatoms) +
                             'atoms expected, ' + str(len(aindex[0])) + ' found')
        
        # Make newsystem by cutting out all atoms in system2 outside boundaries
        newsystem = System(atoms=system2.atoms[aindex], box=system2.box, symbols=self.symbols)
        
        # Return normalized system
        return newsystem.normalize(return_transform=return_transform)
    
    def normalize(self, style='lammps', return_transform=False):
        """
        Normalizes a system's box vectors and atom positions to be compatible
        with simulation codes.
        
        Parameters
        ----------
        style : str, optional
            Indicates the normalization style to use.  Default (and only
            current option) is 'lammps'.
        return_transform : bool, optional
            Indicates if the transformation matrix associated with the
            normalization is returned.  Default value is False.
        Returns
        -------
        newsystem : atomman.System
            A new system that has been normalized.
        transform : np.ndarray
            The transformation matrix associated with the normalization.
            Returned if return_transform is True.
        """
        if style == 'lammps':
            return lmp_normalize(self, return_transform=return_transform)
        else:
            raise ValueError("Unknown style (only 'lammps' is currently supported)")
    
    def dump(self, style, **kwargs):
        """
        Convert a System to another format.
    
        Parameters
        ----------
        style : str
            Indicates the format of the content to dump the atomman.System as.
        kwargs
            Any extra keyword arguments to pass to the underlying dump methods.
            
        Returns
        -------
        str, object or tuple
            Any content returned by the underlying dump methods.
        """
        return dump(style, self, **kwargs)