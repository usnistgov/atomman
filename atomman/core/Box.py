# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from copy import deepcopy

# http://www.numpy.org/
import numpy as np

# atomman imports
from ..tools import vect_angle

class Box(object):
    """
    A representation of a triclinic (parallelepiped) box.
    """
    
    def __init__(self, **kwargs):
        """
        Initilizes a Box.  If parameters besides origin are given they
        must completely define the box.  Allowed parameter sets are:
        
        - no parameters -> box is set to square unit box with origin = [0,0,0].
        
        - origin. -> Only origin is changed (same as setting origin directly).
        
        - vects, (and origin).
        
        - avect, bvect, cvect, (and origin).
        
        - a, b, c, (alpha, beta, gamma, and origin).
        
        - lx, ly, lz, (xy, xz, yz, and origin).
        
        - xlo, xhi, ylo, yhi, zlo, zhi, (xy, xz, and yz).
        
        See the description of class methods and attributes for more details
        on the allowed parameters.
        
        """
        self.__vects = np.eye(3, dtype='float64')
        self.__origin = np.zeros(3, dtype='float64')
        
        if len(kwargs) > 0:
            self.set(**kwargs)
    
    @property
    def vects(self):
        """numpy.ndarray : Array containing all three box vectors.  Can be set directly."""
        return deepcopy(self.__vects)
    
    @vects.setter
    def vects(self, value):
        self.__vects[:] = value
        
        #Zero out near zero terms
        self.__vects[np.isclose(self.__vects/abs(self.__vects).max(), 0.0, atol=1e-9)] = 0.0
    
    @property
    def origin(self):
        """numpy.ndarray : Box origin position where vects are added to define the box.  Can be set directly."""
        return deepcopy(self.__origin)
    
    @origin.setter
    def origin(self, value):
        self.__origin[:] = value
    
    @property
    def avect(self):
        """numpy.ndarray : Vector associated with the a box dimension."""
        return self.vects[0]
    
    @property
    def bvect(self):
        """numpy.ndarray : Vector associated with the b box dimension."""
        return self.vects[1]
    
    @property
    def cvect(self):
        """numpy.ndarray : Vector associated with the c box dimension."""
        return self.vects[2]
    
    @property
    def a(self):
        """float : The a lattice parameter (magnitude of avect)."""
        return (self.__vects[0,0]**2 + self.__vects[0,1]**2 + self.__vects[0,2]**2)**0.5
    
    @property
    def b(self):
        """float : The b lattice parameter (magnitude of avect)."""
        return (self.__vects[1,0]**2 + self.__vects[1,1]**2 + self.__vects[1,2]**2)**0.5
    
    @property
    def c(self):
        """float : The c lattice parameter (magnitude of avect)."""
        return (self.__vects[2,0]**2 + self.__vects[2,1]**2 + self.__vects[2,2]**2)**0.5 
    
    @property
    def alpha(self):
        """float : The alpha lattice angle in degrees (angle between bvect and cvect)."""
        return vect_angle(self.__vects[1], self.__vects[2])
    
    @property
    def beta(self):
        """float : The beta lattice angle in degrees (angle between avect and cvect)."""
        return vect_angle(self.__vects[0], self.__vects[2])
    
    @property
    def gamma(self):
        """float : The gamma lattice angle in degrees (angle between avect and bvect)."""
        return vect_angle(self.__vects[0], self.__vects[1])
    
    @property
    def lx(self):
        """float : LAMMPS lx box length (avect[0] for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__vects[0,0]
    
    @property
    def ly(self):
        """float : LAMMPS ly box length (bvect[1] for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__vects[1,1]
    
    @property
    def lz(self):
        """float : LAMMPS lz box length (cvect[2] for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__vects[2,2]
    
    @property
    def xy(self):
        """float : LAMMPS xy box tilt factor (bvect[0] for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__vects[1,0]
    
    @property
    def xz(self):
        """float : LAMMPS xz box tilt factor (cvect[0] for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__vects[2,0]
    
    @property
    def yz(self):
        """float : LAMMPS yz box tilt factor (cvect[1] for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__vects[2,1]
    
    @property
    def xlo(self):
        """float : LAMMPS xlo box lo term (origin[0] for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__origin[0]
    
    @property
    def ylo(self):
        """float : LAMMPS ylo box lo term (origin[1] for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__origin[1]
    
    @property
    def zlo(self):
        """float : LAMMPS zlo box lo term (origin[2] for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__origin[2]
    
    @property
    def xhi(self):
        """float : LAMMPS xhi box hi term (origin[0] + lx for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__origin[0] + self.__vects[0,0]
    
    @property
    def yhi(self):
        """float : LAMMPS yhi box hi term (origin[1] + ly for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__origin[1] + self.__vects[1,1]
    
    @property
    def zhi(self):
        """float : LAMMPS zhi box hi term (origin[2] + lz for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__origin[2] + self.__vects[2,2]
    
    def __str__(self):
        """
        The string representation of the box.  Lists the three vectors and origin.
        """
        return '\n'.join(['avect =  [%6.3f, %6.3f, %6.3f]' % (self.__vects[0,0], self.__vects[0,1], self.__vects[0,2]),
                          'bvect =  [%6.3f, %6.3f, %6.3f]' % (self.__vects[1,0], self.__vects[1,1], self.__vects[1,2]),
                          'cvect =  [%6.3f, %6.3f, %6.3f]' % (self.__vects[2,0], self.__vects[2,1], self.__vects[2,2]),
                          'origin = [%6.3f, %6.3f, %6.3f]' % (self.__origin[0],  self.__origin[1],  self.__origin[2])])
    
    def set(self, **kwargs):
        """
        Sets a Box's dimensions.  If parameters besides origin are given they
        must completely define the box.  Allowed parameter sets are:
        
        - no parameters -> box is set to square unit box with origin = [0,0,0].
        
        - origin. -> Only origin is changed (same as setting origin directly).
        
        - vects, (and origin).
        
        - avect, bvect, cvect, (and origin).
        
        - a, b, c, (alpha, beta, gamma, and origin).
        
        - lx, ly, lz, (xy, xz, yz, and origin).
        
        - xlo, xhi, ylo, yhi, zlo, zhi, (xy, xz, and yz).
        
        See the description of class methods and attributes for more details
        on the allowed parameters.
        """
        
        # Set default values if no kwargs given
        if len(kwargs) == 0:
            self.vects = np.eye(3)
            self.origin = np.zeros(3)
        
        # Set directly if vects given
        elif 'vects' in kwargs:
            vects = kwargs.pop('vects')
            origin = kwargs.pop('origin', [0.0, 0.0, 0.0])
            assert len(kwargs) == 0, 'Invalid arguments'
            self.vects = vects
            self.origin = origin
        
        # Call set_vectors if vect inputs given
        elif 'avect' in kwargs:
            self.set_vectors(**kwargs)
        
        # Call set_lengths if length inputs given
        elif 'lx' in kwargs:
            self.set_lengths(**kwargs)
        
        # Call set_hi_los if hi/lo inputs given
        elif 'xlo' in kwargs:
            self.set_hi_los(**kwargs)
        
        # Call set_abc if vector magnitudes are given
        elif 'a' in kwargs:
            self.set_abc(**kwargs)
        
        # Set only origin if given alone
        elif 'origin' in kwargs:
            origin = kwargs.pop('origin')
            assert len(kwargs) == 0, 'Invalid arguments'
            self.origin = origin
        
        else:
            raise TypeError('Invalid arguments')
    
    def set_vectors(self, avect, bvect, cvect, origin=[0.0, 0.0, 0.0]):
        """
        Set the box using the three box vectors.
        
        Parameters
        ----------
        avect : numpy.ndarray
            The 3D vector for the a box dimension.
        bvect : numpy.ndarray
            The 3D vector for the b box dimension.
        cvect : numpy.ndarray
            The 3D vector for the c box dimension.
        origin : numpy.ndarray, optional
            The 3D vector for the box origin position.  Default value is
            (0,0,0).
        """
        # Combine avect, bvect and cvect into vects and set directly
        self.vects = [avect, bvect, cvect]
        self.origin = origin
        
    def set_abc(self, a, b, c, alpha=90.0, beta=90.0, gamma=90.0, origin=[0.0, 0.0, 0.0]):
        """
        Set the box using crystal cell lattice parameters and angles.
        
        Parameters
        ----------
        a : float
            The a lattice parameter.
        b : float
            The b lattice parameter.
        c : float
            The c lattice parameter.
        alpha : float, optional
            The alpha lattice angle in degrees (angle between b and c vectors).
            Default value is 90.0.
        beta : float, optional
            The beta lattice angle in degrees (angle between a and c vectors).
            Default value is 90.0.
        gamma : float, optional
            The gamma lattice angle in degrees (angle between a and b vectors).
            Default value is 90.0.
        origin : numpy.ndarray, optional
            The 3D vector for the box origin position.  Default value is
            (0,0,0).
        """
         
        # Convert to lx, ly, lz, xy, xz, yz
        lx = a
        xy = b * np.cos(gamma * np.pi / 180)
        xz = c * np.cos(beta * np.pi / 180)
        ly = (b**2 - xy**2)**0.5
        yz = (b * c * np.cos(alpha * np.pi / 180) - xy * xz) / ly
        lz = (c**2 - xz**2 - yz**2)**0.5

        # Call set_lengths
        self.set_lengths(lx=lx, ly=ly, lz=lz, xy=xy, xz=xz, yz=yz, origin=origin)
        
    def set_lengths(self, lx, ly, lz, xy=0.0, xz=0.0, yz=0.0, origin=[0.0, 0.0, 0.0]):
        """
        Set the box using LAMMPS box lengths and tilt factors.
        
        Parameters
        ----------
        lx : float
            The LAMMPS box length in the x direction.
        ly : float
            The LAMMPS box length in the y direction.
        lz : float
            The LAMMPS box length in the z direction.
        xy : float, optional
            The LAMMPS box tilt factor in the xy direction.  Default value is 
            0.0.
        xz : float, optional
            The LAMMPS box tilt factor in the xz direction.  Default value is 
            0.0.
        yz : float, optional
            The LAMMPS box tilt factor in the yz direction.  Default value is 
            0.0.
        origin : numpy.ndarray, optional
            The 3D vector for the box origin position.  Default value is
            (0,0,0).
        """
        
        assert lx > 0 and ly > 0 and lz > 0, 'box lengths must be positive'
        
        # Construct vects array
        self.vects = [[lx, 0.0, 0.0],
                      [xy, ly,  0.0],
                      [xz, yz,  lz]]
        self.origin = origin
    
    def set_hi_los(self, xlo, xhi, ylo, yhi, zlo, zhi, xy=0.0, xz=0.0, yz=0.0):
        """
        Set the box using LAMMPS box hi's, lo's and tilt factors.
        
        Parameters
        ----------
        xlo : float
            The LAMMPS xlo box lo term.
        xhi : float
            The LAMMPS xhi box hi term.
        ylo : float
            The LAMMPS ylo box lo term.
        yhi : float
            The LAMMPS yhi box hi term.
        zlo : float
            The LAMMPS zlo box lo term.
        zhi : float
            The LAMMPS zhi box hi term.
        xy : float, optional
            The LAMMPS box tilt factor in the xy direction.  Default value is 
            0.0.
        xz : float, optional
            The LAMMPS box tilt factor in the xz direction.  Default value is 
            0.0.
        yz : float, optional
            The LAMMPS box tilt factor in the yz direction.  Default value is 
            0.0.
        """
        
        # Convert to hi and lo term to lengths and origin
        lx = xhi - xlo
        ly = yhi - ylo
        lz = zhi - zlo
        origin = [xlo, ylo, zlo]
        
        # Call set_lengths
        self.set_lengths(lx=lx, ly=ly, lz=lz, xy=xy, xz=xz, yz=yz, origin=origin)
    
    def is_lammps_norm(self):
        """
        Tests if box is compatible with LAMMPS.
        Note: large box tilt factors not checked.  The LAMMPS command
        'box tilt large' may be needed to run LAMMPS.
        """
        return (self.__vects[0,1] == 0.0
            and self.__vects[0,2] == 0.0
            and self.__vects[1,2] == 0.0
            and self.__vects[0,0] > 0.0
            and self.__vects[1,1] > 0.0
            and self.__vects[2,2] > 0.0)