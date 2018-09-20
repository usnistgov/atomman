# coding: utf-8
# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from copy import deepcopy

# http://www.numpy.org/
import numpy as np

# atomman imports
from ..tools import axes_check

class Stroh(object):
    """
    Class for solving the Eshelby anisotropic solution for a straight
    dislocation or crack using the Stroh method.
    """
    
    def __init__(self, C, burgers, axes=None, m=[1,0,0], n=[0,1,0], tol=1e-8):
        """
        Initialize an instance of the Stroh class and computes the Stroh
        solution.
        
        Parameters
        ----------
        C : atomman.ElasticConstants
            The medium's elastic constants.
        burgers : array-like object
            The dislocation's Cartesian Burgers vector.
        axes : array-like object, optional
            3x3 set of rotational axes for the system. If given, C and burgers
            will be transformed using axes.
        m : array-like object, optional
            The m unit vector for the solution.  m, n, and u (dislocation
            line) should be right-hand orthogonal.  Default value is [1,0,0]
            (x-axis).
        n : array-like object, optional
            The n unit vector for the solution.  m, n, and u (dislocation
            line) should be right-hand orthogonal.  Default value is [0,1,0]
            (y-axis).
        tol : float
            Tolerance parameter used to round off near-zero values.  Default
            value is 1e-8.
        """
        
        self.solve(C, burgers, axes=axes, m=m, n=n, tol=tol)

    def solve(self, C, burgers, axes=None, m=[1,0,0], n=[0,1,0], tol=1e-8):
        """
        Computes the Stroh solution.
        
        Parameters
        ----------
        C : atomman.ElasticConstants
            The medium's elastic constants.
        burgers : array-like object
            The dislocation's Cartesian Burgers vector.
        axes : array-like object, optional
            3x3 set of rotational axes for the system. If given, C and burgers
            will be transformed using axes.
        m : array-like object, optional
            The m unit vector for the solution.  m, n, and u (dislocation line)
            should be right-hand orthogonal.  Default value is [1, 0, 0]
            (x-axis).
        n : array-like object, optional
            The n unit vector for the solution.  m, n, and u (dislocation line)
            should be right-hand orthogonal.  Default value is [0, 1, 0]
            (y-axis).
        tol : float
            Tolerance parameter used to round off near-zero values.  Default
            value is 1e-8.
        """
        # Convert burgers, m, n to numpy arrays if needed
        burgers = np.asarray(burgers, dtype='float64')
        m = np.asarray(m, dtype='float64')
        n = np.asarray(n, dtype='float64')
        
        # Transform burgers and C
        if axes is not None:
            T = axes_check(axes)
            burgers = T.dot(burgers)
            C = C.transform(axes)
        
        # Pull out full 3x3x3x3 elastic constants matrix
        Cijkl = C.Cijkl
        
        # Matrices of Cijkl constants used to construct N
        mm = np.einsum('i,ijkl,l', m, Cijkl, m)
        mn = np.einsum('i,ijkl,l', m, Cijkl, n)
        nm = np.einsum('i,ijkl,l', n, Cijkl, m)
        nn = np.einsum('i,ijkl,l', n, Cijkl, n)
        
        # The four 3x3 matrices that represent the quadrants of N
        NB = -np.linalg.inv(nn)
        NA = NB.dot(nm)
        NC = mn.dot(NA) + mm
        ND = mn.dot(NB)
        
        # N is the 6x6 array, where the eigenvalues are the roots p
        # and the eigenvectors give A and L
        N =  np.array(np.vstack((np.hstack((NA, NB)), np.hstack((NC, ND)))))
        
        # Calculate the eigenvectors and eigenvalues
        eig = np.linalg.eig(N)
        p = eig[0]
        eigvec = np.transpose(eig[1])
        
        # Separate the eigenvectors into A and L
        A = np.array([eigvec[0,:3], eigvec[1,:3], eigvec[2,:3], eigvec[3,:3], eigvec[4,:3], eigvec[5,:3]])
        L = np.array([eigvec[0,3:], eigvec[1,3:], eigvec[2,3:], eigvec[3,3:], eigvec[4,3:], eigvec[5,3:]])
        
        # Calculate normalization factor k
        k = 1. / (2. * np.einsum('si,si->s', A, L))
        
        # Calculation verification checks
        try:
            assert np.allclose(np.einsum('s,si,sj->ij', k,A,L), np.identity(3, dtype='complex128'), atol=tol)
            assert np.allclose(np.einsum('s,si,sj->ij', k,A,A), np.zeros((3,3), dtype='complex128'), atol=tol)
            assert np.allclose(np.einsum('s,si,sj->ij', k,L,L), np.zeros((3,3), dtype='complex128'), atol=tol)
            assert np.allclose(np.einsum('s,t,si,ti->st', k**.5,k**.5,A,L) + np.einsum('s,t,ti,si->st', k**.5,k**.5,A,L),
                                         np.identity(6, dtype='complex128'), atol = tol)
        except:
            raise ValueError('Stroh checks failed!')
        
        # Assign property values
        self.__burgers = burgers
        self.__Cijkl = Cijkl
        self.__tol = tol
        self.__p = p
        self.__A = A
        self.__L = L
        self.__k = k
        self.__m = m
        self.__n = n
    
    @property
    def p(self):
        """numpy.ndarray : p eigenvalues"""
        return deepcopy(self.__p)
    
    @property
    def A(self):
        """numpy.ndarray : A eigenvectors"""
        return deepcopy(self.__A)
    
    @property
    def L(self):
        """numpy.ndarray : L eigenvectors"""
        return deepcopy(self.__L)
    
    @property
    def k(self):
        """numpy.ndarray : k normalization factors"""
        return deepcopy(self.__k)
    
    @property
    def K_coeff(self):
        """float : The energy coefficient"""
        
        # K = b_i K_ij b_j / (b_k b_k)
        return (self.__burgers.dot(self.K_tensor.dot(self.__burgers)) 
                / self.__burgers.dot(self.__burgers))
    
    @property
    def K_tensor(self):
        """numpy.ndarray : The anisotropic energy coefficient tensor"""
        
        # ii is imaginary unit
        ii = np.array([1.j])
        
        # updn is alternating +-
        updn = np.array([1, -1, 1, -1, 1, -1])
        
        # Compute K_tensor
        K = ii * np.einsum('s,s,si,sj->ij', updn, self.k, self.L, self.L)
        
        # Round away near-zero terms
        K = np.real_if_close(K, tol=self.__tol)
        K[np.isclose(K / K.max(), 0.0, atol=self.__tol)] = 0.0
        
        return K
    
    @property
    def preln(self):
        """float : The pre-ln strain energy factor"""
        
        # a = b_i K_ij b_j / (4 π)
        return self.__burgers.dot(self.K_tensor.dot(self.__burgers)) / (4 * np.pi)
    
    @property
    def m(self):
        return self.__m
        
    @property
    def n(self):
        return self.__n
    
    def displacement(self, pos):
        """
        Compute the position-dependent anisotropic displacement.
        
            u_i = 1 / (2 π i) (Σ_a +- k_a A_ai (L_aj*burgers_j) ln(η_a))
        
        Parameters
        ----------
        pos : array-like object
            3D vector position(s).
        
        Returns
        -------
        numpy.ndarray
            The computed 3D vector displacements at all given points.
        """
        # ii is imaginary unit
        ii = np.array([1.j])
        
        # updn is alternating +-
        updn = np.array([1, -1, 1, -1, 1, -1])
        
        # Convert pos to eta coordinate
        eta = self.eta(pos)
        
        # Compute the displacements
        disp = 1 / (2 * np.pi * ii) * np.einsum('a,a,ai,a,na->ni',
                                                updn, self.k, self.A,
                                                self.L.dot(self.__burgers),
                                                np.log(eta))
        
        # Round away near-zero terms
        disp = np.real_if_close(disp, tol=self.__tol)
        
        # Reduce single-value solutions
        if np.asarray(pos).ndim == 1:
            return disp[0]
        else:
            return disp
    
    def stress(self, pos):
        """
        Compute the position-dependent anisotropic stresses.
        
            σ_ij = 1 / (2 π i) (Σ_a +- k_a C_ijkl mpn_al A_ak (L_am*burgers_m) / η_a)
        
        Parameters
        ----------
        pos : array-like object
            3D vector position(s).
        
        Returns
        -------
        numpy.ndarray
            The computed 3x3 stress states at all given points.
        """
        # ii is imaginary unit
        ii = np.array([1.j])
        
        # updn is alternating +-
        updn = np.array([1, -1, 1, -1, 1, -1])
        
        # Convert pos to eta coordinate
        eta = self.eta(pos)
        
        # Compute mpn factor
        mpn = self.__m + np.outer(self.p, self.__n)
        
        # Compute the stresses
        stress = 1 / (2 * np.pi * ii) * np.einsum('a,a,ijkl,al,ak,a,na->nij',
                                                  updn, self.k, self.__Cijkl,
                                                  mpn, self.A,
                                                  self.L.dot(self.__burgers),
                                                  1/eta)
        
        # Round away near-zero terms
        stress = np.real_if_close(stress, tol=self.__tol)
        
        if np.asarray(pos).ndim == 1:
            return stress[0]
        else:
            return stress
    
    def eta(self, pos):
        """
        Compute the eta coordinates based on positions, p, m and n.  Used by
        displacement() and stress().
        
            η_a = x_i m_i + p_a x_j n_j
        
        Parameters
        ----------
        pos : array-like object
            3D vector position(s).
        
        Returns
        -------
        numpy.ndarray
            The computed eta factor at all given points.
        """
        
        x = np.dot(pos, self.__m)
        y = np.dot(pos, self.__n)
        
        return (x + np.outer(self.p, y)).T