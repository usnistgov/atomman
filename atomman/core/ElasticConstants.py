# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from copy import deepcopy

# http://www.numpy.org/
import numpy as np

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# atomman imports
from ..tools import axes_check
import atomman.unitconvert as uc
from ..compatibility import iteritems, range

class ElasticConstants(object):
    """Class for storing and converting elastic constant values"""
    
    def __init__(self, **kwargs):
        """
        Initilizes an ElasticConstants instance from one of the parameter
        options.
        
        Parameters
        ----------
        Cij : numpy.ndarray
            (6, 6) array of Voigt representation of elastic stiffness.
        Sij : numpy.ndarray
            (6, 6) array of Voigt representation of elastic compliance.
        Cij9 : numpy.ndarray
            (9, 9) array representation of elastic stiffness.
        Cijkl : numpy.ndarray
            (3, 3, 3, 3) array representation of elastic stiffness.
        Sijkl : numpy.ndarray
            (3, 3, 3, 3) array representation of elastic compliance.
        model : DataModelDict, string, or file-like object 
            Data model containing elastic constants.
        C11, C12, ... C66 : float
            Individual components of Cij for a standardized representation:
            isotropic: C11, C12 (see below for more options)
            cubic: C11, C12, C44
            hexagonal: C11, C12, C13, C33, C44, C66 (2*C66=C11-C12)
            tetragonal: C11, C12, C13, C16, C33, C44, C66 (C16 optional)
            rhombohedral: C11, C12, C13, C14, C15, C33, C44, C66 (2*C66=C11-C12, C15 optional)
            orthorhombic: C11, C12, C13, C22, C23, C33, C44, C55, C66
            monoclinic: C11, C12, C13, C15, C22, C23, C25, C33, C35, C44, C46, C55, C66
            triclinic: all Cij where i <= j
        M, lambda, mu, E, nu, K, C11, C12, C44 : float
            Individual isotropic elastic moduli (exactly two must be given).
        """
        # Initialize for no arguments
        if len(kwargs) == 0:
            self.__c_ij = np.zeros((6,6), dtype='float64')
        
        # Initialize for matrix arguments
        elif 'Cij' in kwargs:
            assert len(kwargs) == 1, 'Cij cannot be specified with other keyword arguments'
            self.Cij = kwargs['Cij']
        elif 'Sij' in kwargs:
            assert len(kwargs) == 1, 'Sij cannot be specified with other keyword arguments'
            self.Sij = kwargs['Sij']
        elif 'Cij9' in kwargs:
            assert len(kwargs) == 1, 'Cij9 cannot be specified with other keyword arguments'
            self.Cij9 = kwargs['Cij9']
        elif 'Cijkl' in kwargs:
            assert len(kwargs) == 1, 'Cijkl cannot be specified with other keyword arguments'
            self.Cijkl = kwargs['Cijkl']
        elif 'Sijkl' in kwargs:
            assert len(kwargs) == 1, 'Sijkl cannot be specified with other keyword arguments'
            self.Sijkl = kwargs['Sijkl']
        
        # Initialize using data model
        elif 'model' in kwargs:
            self.model(**kwargs)
        
        # Initialize for individually specified parameters
        elif len(kwargs) == 2:
            self.isotropic(**kwargs)
        elif len(kwargs) == 3:
            self.cubic(**kwargs)
        elif len(kwargs) == 5:
            self.hexagonal(**kwargs)
        elif len(kwargs) == 6 or len(kwargs) == 7:
            if 'C14' in kwargs:
                self.rhombohedral(**kwargs)
            else:
                self.tetragonal(**kwargs)
        elif len(kwargs) == 8:
            self.rhombohedral(**kwargs)
        elif len(kwargs) == 9:
            self.orthorhombic(**kwargs)
        elif len(kwargs) == 13:
            self.monoclinic(**kwargs)
        elif len(kwargs) == 21:
            self.triclinic(**kwargs)
        else:
            raise TypeError('Invalid argument keywords')
    
    def __str__(self):
        """Calling string returns str(self.Cij)."""
        return str(self.Cij)
    
    @property
    def Cij(self):
        """The stiffness constants in Voigt 6x6 format"""
        return deepcopy(self.__c_ij)
    
    @Cij.setter
    def Cij(self, value):
        value = np.asarray(value, dtype='float64')
        assert value.shape == (6,6),  'Cij must be 6x6'
        
        # Zero out near-zero terms
        assert value.max() > 0.0, 'Cij values not valid'
        value[np.isclose(value/value.max(), 0.0, atol=1e-9)] = 0.0
        
        # Check symmetry
        for i in range(6):
            for j in range(i):
                assert np.isclose(value[i,j], value[j,i], atol=1e-9), '6x6 matrix not symmetric!' 
        self.__c_ij = value
    
    @property
    def Sij(self):
        """The compliance constants in Voigt 6x6 format"""
        return np.linalg.inv(self.Cij)
    
    @Sij.setter
    def Sij(self, value):
        value = np.asarray(value, dtype='float64')
        assert value.shape == (6,6),  'Sij must be 6x6'
        self.Cij = np.linalg.inv(value)
    
    @property
    def Cij9(self):
        """The stiffness constants in 9x9 format"""
        c = self.Cij
        return np.array([[c[0,0],c[0,1],c[0,2],c[0,3],c[0,4],c[0,5],c[0,3],c[0,4],c[0,5]],
                         [c[1,0],c[1,1],c[1,2],c[1,3],c[1,4],c[1,5],c[1,3],c[1,4],c[1,5]],
                         [c[2,0],c[2,1],c[2,2],c[2,3],c[2,4],c[2,5],c[2,3],c[2,4],c[2,5]],
                         [c[3,0],c[3,1],c[3,2],c[3,3],c[3,4],c[3,5],c[3,3],c[3,4],c[3,5]],
                         [c[4,0],c[4,1],c[4,2],c[4,3],c[4,4],c[4,5],c[4,3],c[4,4],c[4,5]],
                         [c[5,0],c[5,1],c[5,2],c[5,3],c[5,4],c[5,5],c[5,3],c[5,4],c[5,5]],
                         [c[3,0],c[3,1],c[3,2],c[3,3],c[3,4],c[3,5],c[3,3],c[3,4],c[3,5]],
                         [c[4,0],c[4,1],c[4,2],c[4,3],c[4,4],c[4,5],c[4,3],c[4,4],c[4,5]],
                         [c[5,0],c[5,1],c[5,2],c[5,3],c[5,4],c[5,5],c[5,3],c[5,4],c[5,5]]])
    
    @Cij9.setter
    def Cij9(self, value):
        value = np.asarray(value, dtype='float64')
        assert value.shape == (9,9), 'Cij9 must be 9x9'
        
        # Check symmetry
        for i in range(6, 9):
            for j in range(9):
                assert value[i,j] == value[i-3, j]
                assert value[j,i] == value[j, i-3]
        self.Cij = value[:6, :6]
    
    @property
    def Cijkl(self):
        """The stiffness constants in 3x3x3x3 format"""
        c = self.Cij
        return np.array([[[[c[0,0],c[0,5],c[0,4]], [c[0,5],c[0,1],c[0,3]], [c[0,4],c[0,3],c[0,2]]],
                          [[c[5,0],c[5,5],c[5,4]], [c[5,5],c[5,1],c[5,3]], [c[5,4],c[5,3],c[5,2]]],
                          [[c[4,0],c[4,5],c[4,4]], [c[4,5],c[4,1],c[4,3]], [c[4,4],c[4,3],c[4,2]]]], 
                        
                         [[[c[5,0],c[5,5],c[5,4]], [c[5,5],c[5,1],c[5,3]], [c[5,4],c[5,3],c[5,2]]],
                          [[c[1,0],c[1,5],c[1,4]], [c[1,5],c[1,1],c[1,3]], [c[1,4],c[1,3],c[1,2]]],
                          [[c[3,0],c[3,5],c[3,4]], [c[3,5],c[3,1],c[3,3]], [c[3,4],c[3,3],c[3,2]]]],
                         
                         [[[c[4,0],c[4,5],c[4,4]], [c[4,5],c[4,1],c[4,3]], [c[4,4],c[4,3],c[4,2]]],
                          [[c[3,0],c[3,5],c[3,4]], [c[3,5],c[3,1],c[3,3]], [c[3,4],c[3,3],c[3,2]]],
                          [[c[2,0],c[2,5],c[2,4]], [c[2,5],c[2,1],c[2,3]], [c[2,4],c[2,3],c[2,2]]]]])
    
    @Cijkl.setter
    def Cijkl(self, value):
        c = np.asarray(value, dtype='float64')
        assert c.shape == (3,3,3,3),  'Cijkl must be 3x3x3x3'
        assert c.max() > 0.0, 'Cij values not valid'
        # Check symmetry
        indexes = np.array([[0,0], [1,1], [2,2], [1,2], [0,2], [0,1]], dtype=int)
        for ij in range(6):
            for kl in range(ij, 6):
                i, j, k, l = indexes[ij,0], indexes[ij,1], indexes[kl,0], indexes[kl,1]
                assert np.isclose(c[i,j,k,l], c[j,i,k,l])
                assert np.isclose(c[i,j,k,l], c[j,i,l,k])
                assert np.isclose(c[i,j,k,l], c[k,l,j,i])
                assert np.isclose(c[i,j,k,l], c[l,k,j,i])
                assert np.isclose(c[i,j,k,l], c[i,j,l,k])
                assert np.isclose(c[i,j,k,l], c[k,l,i,j])
                assert np.isclose(c[i,j,k,l], c[l,k,i,j])
        
        self.Cij = np.array([[c[0,0,0,0], c[0,0,1,1], c[0,0,2,2], c[0,0,1,2], c[0,0,0,2], c[0,0,0,1]],
                             [c[1,1,0,0], c[1,1,1,1], c[1,1,2,2], c[1,1,1,2], c[1,1,0,2], c[1,1,0,1]],
                             [c[2,2,0,0], c[2,2,1,1], c[2,2,2,2], c[2,2,1,2], c[2,2,0,2], c[2,2,0,1]],
                             [c[1,2,0,0], c[1,2,1,1], c[1,2,2,2], c[1,2,1,2], c[1,2,0,2], c[1,2,0,1]],
                             [c[0,2,0,0], c[0,2,1,1], c[0,2,2,2], c[0,2,1,2], c[0,2,0,2], c[0,2,0,1]],
                             [c[0,1,0,0], c[0,1,1,1], c[0,1,2,2], c[0,1,1,2], c[0,1,0,2], c[0,1,0,1]]])
    
    @property
    def Sijkl(self):
        """The compliance constants in 3x3x3x3 format"""
        s = self.Sij
        s[3:,:] = s[3:,:]/2.
        s[:,3:] = s[:,3:]/2.
        return np.array([[[[s[0,0],s[0,5],s[0,4]], [s[0,5],s[0,1],s[0,3]], [s[0,4],s[0,3],s[0,2]]],
                          [[s[5,0],s[5,5],s[5,4]], [s[5,5],s[5,1],s[5,3]], [s[5,4],s[5,3],s[5,2]]],
                          [[s[4,0],s[4,5],s[4,4]], [s[4,5],s[4,1],s[4,3]], [s[4,4],s[4,3],s[4,2]]]],
                         
                         [[[s[5,0],s[5,5],s[5,4]], [s[5,5],s[5,1],s[5,3]], [s[5,4],s[5,3],s[5,2]]],
                          [[s[1,0],s[1,5],s[1,4]], [s[1,5],s[1,1],s[1,3]], [s[1,4],s[1,3],s[1,2]]],
                          [[s[3,0],s[3,5],s[3,4]], [s[3,5],s[3,1],s[3,3]], [s[3,4],s[3,3],s[3,2]]]],
                         
                         [[[s[4,0],s[4,5],s[4,4]], [s[4,5],s[4,1],s[4,3]], [s[4,4],s[4,3],s[4,2]]],
                          [[s[3,0],s[3,5],s[3,4]], [s[3,5],s[3,1],s[3,3]], [s[3,4],s[3,3],s[3,2]]],
                          [[s[2,0],s[2,5],s[2,4]], [s[2,5],s[2,1],s[2,3]], [s[2,4],s[2,3],s[2,2]]]]])
    
    @Sijkl.setter
    def Sijkl(self, value):
        s = np.asarray(value, dtype='float64')
        assert s.shape == (3,3,3,3),  'Sijkl must be 3x3x3x3'
        
        # Check symmetry
        indexes = np.array([[0,0], [1,1], [2,2], [1,2], [0,2], [0,1]], dtype=int)
        for ij in range(6):
            for kl in range(ij, 6):
                i, j, k, l = indexes[ij,0], indexes[ij,1], indexes[kl,0], indexes[kl,1]
                assert np.isclose(s[i,j,k,l], s[j,i,k,l])
                assert np.isclose(s[i,j,k,l], s[j,i,l,k])
                assert np.isclose(s[i,j,k,l], s[k,l,j,i])
                assert np.isclose(s[i,j,k,l], s[l,k,j,i])
                assert np.isclose(s[i,j,k,l], s[i,j,l,k])
                assert np.isclose(s[i,j,k,l], s[k,l,i,j])
                assert np.isclose(s[i,j,k,l], s[l,k,i,j]) 
        
        self.Sij = np.array([[   s[0,0,0,0],    s[0,0,1,1],    s[0,0,2,2], 2.*s[0,0,1,2], 2.*s[0,0,0,2], 2.*s[0,0,0,1]],
                             [   s[1,1,0,0],    s[1,1,1,1],    s[1,1,2,2], 2.*s[1,1,1,2], 2.*s[1,1,0,2], 2.*s[1,1,0,1]],
                             [   s[2,2,0,0],    s[2,2,1,1],    s[2,2,2,2], 2.*s[2,2,1,2], 2.*s[2,2,0,2], 2.*s[2,2,0,1]],
                             [2.*s[1,2,0,0], 2.*s[1,2,1,1], 2.*s[1,2,2,2], 4.*s[1,2,1,2], 4.*s[1,2,0,2], 4.*s[1,2,0,1]],
                             [2.*s[0,2,0,0], 2.*s[0,2,1,1], 2.*s[0,2,2,2], 4.*s[0,2,1,2], 4.*s[0,2,0,2], 4.*s[0,2,0,1]],
                             [2.*s[0,1,0,0], 2.*s[0,1,1,1], 2.*s[0,1,2,2], 4.*s[0,1,1,2], 4.*s[0,1,0,2], 4.*s[0,1,0,1]]])
    
    def transform(self, axes, tol=1e-8):
        """
        Transforms the elastic constant matrix based on the supplied axes.
        
        Parameters
        ----------
        axes : numpy.ndarray
            (3, 3) array giving three right-handed orthogonal vectors to use
            for transforming.
        tol : float, optional
            Relative tolerance to use in identifying near-zero terms.
            
        Returns
        -------
        ElasticConstants
            A new ElasticConstants object that has been transformed.
        """
        axes = np.asarray(axes, dtype='float64')
        T = axes_check(axes)
        
        Q = np.einsum('km,ln->mnkl', T, T)
        C = np.einsum('ghij,ghmn,mnkl->ijkl', Q, self.Cijkl, Q)
        C[abs(C / C.max()) < tol] = 0.0
        
        return ElasticConstants(Cijkl=C)
    
    def isotropic(self, **kwargs):
        """
        Set values with two independent isotropic moduli.
        
        Parameters
        ----------
        C11 : float, optional
            C11 component of Cij.
        C12 : float, optional
            C12 component of Cij.
        C44 : float, optional
            C44 component of Cij.
        M : float, optional
            P-wave modulus(Equivalent to C11).
        lambda : float, optional
            Lame's first parameter (Equivalent to C12).
        mu : float, optional
            Shear modulus (Equivalent to C44).
        E : float, optional
            Young's modulus
        nu : float, optional
            Poisson's ratio
        K : float, optional
            Bulk modulus
        """
        
        try:
            # Handle equivalent terms
            if 'M' in kwargs:
                kwargs['C11'] = kwargs.pop('M')
            if 'lambda' in kwargs:
                kwargs['C12'] = kwargs.pop('lambda')
            if 'mu' in kwargs:
                kwargs['C44'] = kwargs.pop('mu')
            
            # Check len of kwargs
            assert len(kwargs) == 2
            
            # Pop and convert terms
            if 'C11' in kwargs:
                c11 = kwargs.pop('C11')
                
                if 'C12' in kwargs:
                    c12 = kwargs.pop('C12')
                    c44 = (c11 - c12) / 2
                
                else:
                    if 'C44' in kwargs:
                        c44 = kwargs.pop('C44')
                    
                    elif 'E' in kwargs: 
                        E = kwargs.pop('E')
                        S = (E**2 + 9 * c11**2 - 10 * E * c11)**0.5
                        c44 = (3 * c11 + E - S) / 8
                    
                    elif 'nu' in kwargs: 
                        nu = kwargs.pop('nu')
                        c44 = c11 * (1 - 2 * nu) / (2 * (1 - nu))
                    
                    elif 'K' in kwargs: 
                        K = kwargs.pop('K')
                        c44 = 3 * (c11 - K) / 4
                    
                    c12 = c11 - 2 * c44
            
            else:
                if 'C12' in kwargs:
                    c12 = kwargs.pop('C12')
                    
                    if 'C44' in kwargs:
                        c44 = kwargs.pop('C44')
                    
                    elif 'E' in kwargs:
                        E = kwargs.pop('E')
                        R = (E**2 + 9 * c12**2 + 2 * E * c12)**0.5
                        c44 = (E - 3 * c12 + R) / 4
                    
                    elif 'nu' in kwargs:
                        nu = kwargs.pop('nu')
                        c44 = c12 * (1 - 2 * nu) / (2 * nu)
                    
                    elif 'K' in kwargs:
                        K = kwargs.pop('K')
                        c44 = 3 * (K - c12) / 2
                
                elif 'C44' in kwargs:
                    c44 = kwargs.pop('C44')
                    
                    if 'E' in kwargs:
                        E = kwargs.pop('E')
                        c12 = c44 * (E - 2 * c44) / (3 * c44 - E)
                    
                    elif 'nu' in kwargs:
                        nu = kwargs.pop('nu')
                        c12 = 2 * c44 * nu / (1 - 2 * nu)
                    
                    elif 'K' in kwargs:
                        K = kwargs.pop('K')
                        c12 = K - 2 * c44 / 3
                
                elif 'E' in kwargs: 
                    E = kwargs.pop('E')
                    
                    if 'nu' in kwargs: 
                        nu = kwargs.pop('nu')
                        c12 = E * nu / ((1 + nu) * (1 - 2 * nu))
                        c44 = E / (2 * (1 + nu))
                    
                    elif 'K' in kwargs: 
                        K = kwargs.pop('K')
                        c12 = 3 * K * (3 * K - E) / (9 * K - E)
                        c44 = 3 * K * E / (9 * K - E)
                
                elif 'nu' in kwargs:
                    nu = kwargs.pop('nu')
                    
                    if 'K' in kwargs: 
                        K = kwargs.pop('K')
                        c12 = 3 * K * nu / (1 + nu)
                        c44 = 3 * K * (1 - 2 * nu) / (2 * (1 + nu))
                
                c11 = c12 + 2 * c44
        
        except:
            raise TypeError('isotropic style takes two unique keyword arguments of (C11=M, C12=lambda, C44=mu, E, nu, K)')
        
        # Build Cij array
        self.Cij = np.array([[c11, c12, c12, 0.0, 0.0, 0.0],
                             [c12, c11, c12, 0.0, 0.0, 0.0],
                             [c12, c12, c11, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, c44, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, c44, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, c44]]) 
    
    def cubic(self, **kwargs):
        """
        Set values with three independent cubic moduli.
        
        Parameters
        ----------
        C11 : float
            C11 component of Cij.
        C12 : float
            C12 component of Cij.
        C44 : float
            C44 component of Cij.
        """
        
        try:
            # Check len of kwargs
            assert len(kwargs) == 3
            
            # Pop required independent terms
            c11 = kwargs.pop('C11')
            c12 = kwargs.pop('C12')
            c44 = kwargs.pop('C44')
        except:
            raise TypeError('cubic style takes keyword arguments C11, C12, and C66')
        
        # Build Cij array
        self.Cij = np.array([[c11, c12, c12, 0.0, 0.0, 0.0],
                             [c12, c11, c12, 0.0, 0.0, 0.0],
                             [c12, c12, c11, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, c44, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, c44, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, c44]])
    
    def hexagonal(self, **kwargs):
        """
        Set values with five independent hexagonal moduli.
        (2 * C66 = C11 - C12)
        
        Parameters
        ----------
        C11 : float, optional
            C11 component of Cij.
        C12 : float, optional
            C12 component of Cij.
        C13 : float
            C13 component of Cij.
        C33 : float
            C33 component of Cij.
        C44 : float
            C44 component of Cij.
        C66 : float, optional
            C66 component of Cij.
        """
        
        try:
            # Check len of kwargs
            assert len(kwargs) >= 5 and len(kwargs) <=6
            
            # Pop required independent terms
            c33 = kwargs.pop('C33')
            c13 = kwargs.pop('C13')
            c44 = kwargs.pop('C44')
            
            # Pop required dependent terms
            if 'C11' in kwargs and 'C12' in kwargs:
                c11 = kwargs.pop('C11')
                c12 = kwargs.pop('C12')
                c66 = (c11 - c12) / 2
                
                # Check if redundant C66 is given
                if 'C66' in kwargs:
                    assert np.isclose(c66, kwargs['C66'])
                    c66 = kwargs.pop('C66')
            
            elif 'C11' in kwargs and 'C66' in kwargs:
                c11 = kwargs.pop('C11')
                c66 = kwargs.pop('C66')
                c12 = c11 - 2 * c66
            
            elif 'C12' in kwargs and 'C66' in kwargs:
                c12 = kwargs.pop('C12')
                c66 = kwargs.pop('C66')
                c11 = 2 * c66 + c12
            else:
                assert False
        except:
            raise TypeError('hexagonal style takes keyword arguments C33, C13, C44, and two of (2*C66=C11-C12)')
        
        # Build Cij array
        self.Cij = np.array([[c11, c12, c13, 0.0, 0.0, 0.0],
                             [c12, c11, c13, 0.0, 0.0, 0.0],
                             [c13, c13, c33, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, c44, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, c44, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, c66]])
    
    def rhombohedral(self, **kwargs):
        """
        Set values with six or seven independent rhombohedral moduli.
        (2 * C66 = C11 - C12)
        
        Parameters
        ----------
        C11 : float, optional
            C11 component of Cij.
        C12 : float, optional
            C12 component of Cij.
        C13 : float
            C13 component of Cij.
        C14 : float
            C33 component of Cij.
        C15 : float, optional
            C44 component of Cij.
        C33 : float
            C66 component of Cij.
        C44 : float
            C66 component of Cij.
        C66 : float, optional
            C66 component of Cij.
        """
        
        try:
            # Check len of kwargs
            assert len(kwargs) >= 6 and len(kwargs) <= 8
            
            # Pop required independent terms
            c33 = kwargs.pop('C33')
            c13 = kwargs.pop('C13')
            c14 = kwargs.pop('C14')
            c44 = kwargs.pop('C44')
            
            # Pop required dependent terms
            if 'C11' in kwargs and 'C12' in kwargs:
                c11 = kwargs.pop('C11')
                c12 = kwargs.pop('C12')
                c66 = (c11 - c12) / 2
                
                # Check if redundant C66 is given
                if 'C66' in kwargs:
                    assert np.isclose(c66, kwargs['C66'])
                    c66 = kwargs.pop('C66')
            
            elif 'C11' in kwargs and 'C66' in kwargs:
                c11 = kwargs.pop('C11')
                c66 = kwargs.pop('C66')
                c12 = c11 - 2 * c66
            
            elif 'C12' in kwargs and 'C66' in kwargs:
                c12 = kwargs.pop('C12')
                c66 = kwargs.pop('C66')
                c11 = 2 * c66 + c12
            else:
                assert False
            
            # Check for optional term
            if len(kwargs) == 0:
                c15 = 0.0
            else:
                c15 = kwargs.pop('C15')
                assert len(kwargs) == 0
        except:
            raise TypeError('rhombohedral style takes keyword arguments C33, C13, C14, C44, at least two of (2*C66=C11-C12) and optional C15')
        
        # Build Cij array
        self.Cij = np.array([[c11, c12, c13, c14, c15, 0.0],
                             [c12, c11, c13,-c14,-c15, 0.0],
                             [c13, c13, c33, 0.0, 0.0, 0.0],
                             [c14,-c14, 0.0, c44, 0.0,-c15],
                             [c15,-c15, 0.0, 0.0, c44, c14],
                             [0.0, 0.0, 0.0,-c15, c14, c66]])
    
    def tetragonal(self, **kwargs):
        """
        Set values with six or seven independent tetragonal moduli.
        
        Parameters
        ----------
        C11 : float
            C11 component of Cij.
        C12 : float
            C12 component of Cij.
        C13 : float
            C13 component of Cij.
        C16 : float, optional
            C16 component of Cij.
        C33 : float
            C33 component of Cij.
        C44 : float
            C44 component of Cij.
        C66 : float
            C66 component of Cij.
        """
        
        try:
            # Check len of kwargs
            assert len(kwargs) == 6 or len(kwargs) == 7
            
            # Pop required independent terms
            c11 = kwargs.pop('C11')
            c33 = kwargs.pop('C33')
            c12 = kwargs.pop('C12')
            c13 = kwargs.pop('C13')
            c44 = kwargs.pop('C44')
            c66 = kwargs.pop('C66')
            
            # Check for optional term
            if len(kwargs) == 0:
                c16 = 0.0
            else:
                c16 = kwargs.pop('C16')
                assert len(kwargs) == 0
        except:
            raise TypeError('tetragonal style takes keyword arguments C11, C33, C12, C13, C44, C66, and optional C16')
        
        # Build Cij array
        self.Cij = np.array([[c11, c12, c13, 0.0, 0.0, c16],
                             [c12, c11, c13, 0.0, 0.0,-c16],
                             [c13, c13, c33, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, c44, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, c44, 0.0],
                             [c16,-c16, 0.0, 0.0, 0.0, c66]])
    
    def orthorhombic(self, **kwargs):
        """
        Set values with nine independent orthorhombic moduli.
        
        Parameters
        ----------
        C11 : float
            C11 component of Cij.
        C12 : float
            C12 component of Cij.
        C13 : float
            C13 component of Cij.
        C22 : float
            C22 component of Cij.
        C23 : float
            C23 component of Cij.
        C33 : float
            C33 component of Cij.
        C44 : float
            C44 component of Cij.
        C55 : float
            C55 component of Cij.
        C66 : float
            C66 component of Cij.
        """
        
        try:
            # Check len of kwargs
            assert len(kwargs) == 9
            
            # Set required independent terms
            c11 = kwargs['C11']
            c22 = kwargs['C22']
            c33 = kwargs['C33']
            c12 = kwargs['C12']
            c13 = kwargs['C13']
            c23 = kwargs['C23']
            c44 = kwargs['C44']
            c55 = kwargs['C55']
            c66 = kwargs['C66']
        except:
            raise TypeError('orthorhombic style takes keyword arguments C11, C22, C33, C12, C13, C23, C44, C55, C66')
        
        # Build Cij array
        self.Cij = np.array([[c11, c12, c13, 0.0, 0.0, 0.0],
                             [c12, c22, c23, 0.0, 0.0, 0.0],
                             [c13, c23, c33, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, c44, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, c55, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, c66]])
    
    def monoclinic(self, **kwargs):
        """
        Set values with thirteen independent monoclinic moduli.
        
        Parameters
        ----------
        C11 : float
            C11 component of Cij.
        C12 : float
            C12 component of Cij.
        C13 : float
            C13 component of Cij.
        C15 : float
            C15 component of Cij.
        C22 : float
            C22 component of Cij.
        C23 : float
            C23 component of Cij.
        C25 : float
            C25 component of Cij.
        C33 : float
            C33 component of Cij.
        C35 : float
            C35 component of Cij.
        C44 : float
            C44 component of Cij.
        C46 : float
            C46 component of Cij.
        C55 : float
            C55 component of Cij.
        C66 : float
            C66 component of Cij.
        """
        
        try:
            # Check len of kwargs
            assert len(kwargs) == 13
            
            # Set required independent terms
            c11 = kwargs['C11']
            c12 = kwargs['C12']
            c13 = kwargs['C13']
            c15 = kwargs['C15']
            c22 = kwargs['C22']
            c23 = kwargs['C23']
            c25 = kwargs['C25']
            c33 = kwargs['C33']
            c35 = kwargs['C35']
            c44 = kwargs['C44']
            c46 = kwargs['C46']
            c55 = kwargs['C55']
            c66 = kwargs['C66']
        except:
            raise TypeError('monoclinic style takes keyword arguments C11, C12, C13, C15, C22, C23, C25, C33, C35, C44, C46, C55, C66')
        
        # Build Cij array
        self.Cij = np.array([[c11, c12, c13, 0.0, c15, 0.0],
                             [c12, c22, c23, 0.0, c25, 0.0],
                             [c13, c23, c33, 0.0, c35, 0.0],
                             [0.0, 0.0, 0.0, c44, 0.0, c46],
                             [c15, c25, c35, 0.0, c55, 0.0],
                             [0.0, 0.0, 0.0, c46, 0.0, c66]])
    
    def triclinic(self, **kwargs):
        """
        Set values with twenty one independent triclinic moduli
        
        Parameters
        ----------
        C11 : float
            C11 component of Cij.
        C12 : float
            C12 component of Cij.
        C13 : float
            C13 component of Cij.
        C14 : float
            C14 component of Cij.
        C15 : float
            C15 component of Cij.
        C16 : float
            C16 component of Cij.
        C22 : float
            C22 component of Cij.
        C23 : float
            C23 component of Cij.
        C24 : float
            C24 component of Cij.
        C25 : float
            C25 component of Cij.
        C26 : float
            C26 component of Cij.
        C33 : float
            C33 component of Cij.
        C34 : float
            C34 component of Cij.
        C35 : float
            C35 component of Cij.
        C36 : float
            C36 component of Cij.
        C44 : float
            C44 component of Cij.
        C45 : float
            C45 component of Cij.
        C46 : float
            C46 component of Cij.
        C55 : float
            C55 component of Cij.
        C56 : float
            C56 component of Cij.
        C66 : float
            C66 component of Cij.
        """
        
        try:
            # Check len of kwargs
            assert len(kwargs) == 21
            
            # Set required independent terms
            c11 = kwargs['C11']
            c12 = kwargs['C12']
            c13 = kwargs['C13']
            c14 = kwargs['C14']
            c15 = kwargs['C15']
            c16 = kwargs['C16']
            c22 = kwargs['C22']
            c23 = kwargs['C23']
            c24 = kwargs['C24']
            c25 = kwargs['C25']
            c26 = kwargs['C26']
            c33 = kwargs['C33']
            c34 = kwargs['C34']
            c35 = kwargs['C35']
            c36 = kwargs['C36']
            c44 = kwargs['C44']
            c45 = kwargs['C45']
            c46 = kwargs['C46']
            c55 = kwargs['C55']
            c56 = kwargs['C56']
            c66 = kwargs['C66']
        except:
            raise TypeError('triclinic style takes keyword arguments of all Cij where i <= j')
        
        # Build Cij array
        self.Cij = np.array([[c11, c12, c13, c14, c15, c16],
                             [c12, c22, c23, c24, c25, c26],
                             [c13, c23, c33, c34, c35, c36],
                             [c14, c24, c34, c44, c45, c46],
                             [c15, c25, c35, c45, c55, c56],
                             [c16, c26, c36, c46, c56, c66]])
    
    def model(self, model=None, unit=None, crystal_system='triclinic'):
        """
        Return or set DataModelDict representation of the elastic constants.
        
        Parameters
        ----------
        model : DataModelDict, string, or file-like object, optional
            Data model containing exactly one 'elastic-constants' branch to
            read.
        unit : str, optional
            Units or pressure to save values in when building a data model.
            Default value is None (no conversion).
        crystal_system : str, optional
            Indicates the crystal system representation to use when building a
            data model. Default value is 'triclinic' (save all values in Cij).
        
        Returns
        -------
        DataModelDict
            If model is not given as a parameter.
        """
        
        # Set values if model given
        if model is not None:
            
            # Find elastic-constants element
            model = DM(model).find('elastic-constants')
            
            # Read in values
            c_dict = {}
            for C in model['C']:
                key = 'C' + C['ij'][0] + C['ij'][2]
                c_dict[key] = uc.value_unit(C['stiffness'])
            self.Cij = ElasticConstants(**c_dict).Cij 
        
        # Return DataModelDict if model not given
        else:
            
            # Start creating data model
            model = DM()
            model['elastic-constants'] = DM()
            model['elastic-constants']['C'] = C = []
            
            # Get values in the specified unit
            c = self.Cij
            c_dict = DM()
            
            if crystal_system == 'cubic':
                c_dict['1 1'] = (c[0,0] + c[1,1] + c[2,2]) / 3
                c_dict['1 2'] = (c[0,1] + c[0,2] + c[1,2]) / 3
                c_dict['4 4'] = (c[3,3] + c[4,4] + c[5,5]) / 3
            
            elif crystal_system == 'hexagonal':
                c_dict['1 1'] = (c[0,0] + c[1,1]) / 2
                c_dict['3 3'] = c[2,2]
                c_dict['1 2'] = (c[0,1] + (c[0,0] - 2*c[5,5])) / 2
                c_dict['1 3'] = (c[0,2] + c[1,2]) / 2
                c_dict['4 4'] = (c[3,3] + c[4,4]) / 2
            
            elif crystal_system == 'tetragonal':
                c_dict['1 1'] = (c[0,0] + c[1,1]) / 2
                c_dict['3 3'] = c[2,2]
                c_dict['1 2'] = c[0,1]
                c_dict['1 3'] = (c[0,2] + c[1,2]) / 2
                c_dict['1 6'] = (c[0,5] - c[1,5]) / 2
                c_dict['4 4'] = (c[3,3] + c[4,4]) / 2
                c_dict['6 6'] = c[5,5]
            
            elif crystal_system == 'rhombohedral':
                c_dict['1 1'] = (c[0,0] + c[1,1]) / 2
                c_dict['3 3'] = c[2,2]
                c_dict['1 2'] = (c[0,1] + (c[0,0] - 2*c[5,5])) / 2
                c_dict['1 3'] = (c[0,2] + c[1,2]) / 2
                c_dict['1 4'] = (c[0,3] - c[1,3]) / 2
                c_dict['1 5'] = (c[0,4] - c[1,4] - c[3,5]) / 3
                c_dict['4 4'] = (c[3,3] + c[4,4]) / 2
            
            elif crystal_system == 'orthorhombic':
                c_dict['1 1'] = c[0,0]
                c_dict['2 2'] = c[1,1]
                c_dict['3 3'] = c[2,2]
                c_dict['1 2'] = c[0,1]
                c_dict['1 3'] = c[0,2]
                c_dict['2 3'] = c[1,2]
                c_dict['4 4'] = c[3,3]
                c_dict['5 5'] = c[4,4]
                c_dict['6 6'] = c[5,5]
            
            elif crystal_system == 'triclinic':
                c_dict['1 1'] = c[0,0]
                c_dict['1 2'] = c[0,1]
                c_dict['1 3'] = c[0,2]
                c_dict['1 4'] = c[0,3]
                c_dict['1 5'] = c[0,4]
                c_dict['1 6'] = c[0,5]
                c_dict['2 2'] = c[1,1]
                c_dict['2 3'] = c[1,2]
                c_dict['2 4'] = c[1,3]
                c_dict['2 5'] = c[1,4]
                c_dict['2 6'] = c[1,5]
                c_dict['3 3'] = c[2,2]
                c_dict['3 4'] = c[2,3]
                c_dict['3 5'] = c[2,4]
                c_dict['3 6'] = c[2,5]
                c_dict['4 4'] = c[3,3]
                c_dict['4 5'] = c[3,4]
                c_dict['4 6'] = c[3,5]
                c_dict['5 5'] = c[4,4]
                c_dict['5 6'] = c[4,5]
                c_dict['6 6'] = c[5,5]
            
            else:
                raise ValueError('Invalid crystal_system: ' + crystal_system)
            
            for ij, value in iteritems(c_dict):
                cij = DM()
                cij['stiffness'] = uc.model(value, unit)
                cij['ij'] = ij
                C.append(cij)
            return model
    
    def bulk(self, style='Hill'):
        """
        Returns a bulk modulus estimate.
        
        Parameters
        ----------
        style : str
            Indicates which style of estimate to use.  Default value is 'Hill'.
            - 'Hill' -- Hill estimate (average of Voigt and Reuss).
            - 'Voigt' -- Voigt estimate. Uses Cij.
            - 'Reuss' -- Reuss estimate. Uses Sij.
        """
        if style == 'Hill':
            return (self.bulk('Voigt') + self.bulk('Reuss')) / 2
        
        elif style == 'Voigt':
            c = self.Cij
            return ( (c[0,0] + c[1,1] + c[2,2]) + 2*(c[0,1] + c[1,2] + c[0,2]) ) / 9
        
        elif style == 'Reuss':
            s = self.Sij
            return 1 / ( (s[0,0] + s[1,1] + s[2,2]) + 2*(s[0,1] + s[1,2] + s[0,2]) )
        
        else:
            raise ValueError('Unknown estimate style')
    
    def shear(self, style='Hill'):
        """
        Returns a shear modulus estimate.
        
        Parameters
        ----------
        style : str
            Indicates which style of estimate to use.  Default value is 'Hill'.
            - 'Hill' -- Hill estimate (average of Voigt and Reuss).
            - 'Voigt' -- Voigt estimate. Uses Cij.
            - 'Reuss' -- Reuss estimate. Uses Sij.
        """
        if style == 'Hill':
            return (self.shear('Voigt') + self.shear('Reuss')) / 2
        
        elif style == 'Voigt':
            c = self.Cij
            return ( (c[0,0] + c[1,1] + c[2,2]) - (c[0,1] + c[1,2] + c[0,2]) + 3*(c[3,3] + c[4,4] + c[5,5]) ) / 15
        
        elif style == 'Reuss':
            s = self.Sij
            return 15 / ( 4*(s[0,0] + s[1,1] + s[2,2]) - 4*(s[0,1] + s[1,2] + s[0,2]) + 3*(s[3,3] + s[4,4] + s[5,5]) )
        
        else:
            raise ValueError('Unknown estimate style')