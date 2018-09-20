# coding: utf-8
# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
import os
from copy import deepcopy
import warnings

# http://www.numpy.org/
import numpy as np

# https://www.scipy.org/
from scipy.optimize import minimize

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# atomman imports
import atomman.unitconvert as uc
from ..tools import axes_check
from .GammaSurface import GammaSurface

class SDVPN(object):
    """
    Class representation of the semidiscrete variational Peierls-Nabarro
    dislocation model.
    """
    
    def __init__(self, *args, **kwargs):
        """
        Class initializer. Calls either load() or solve() based on the given
        parameters.
        
        Parameters
        ----------
        *args : list 
            Any unnamed arguments passed on to set().
        *kwargs : dict 
            Any named arguments to pass on.  If 'model' is in kwargs, then
            model() is called.  Else, *args and **kwargs are saved as class
            properties, and then solve() is called.
        """
        
        # Load solution from existing model
        if 'model' in kwargs:
            gamma = kwargs.pop('gamma', None)
            if len(kwargs) == 1 and len(args) == 0:
                self.load(kwargs['model'], gamma=gamma)
            else:
                raise ValueError('model keyword cannot be given with calculation parameters')
        
        # Solve new problem
        else:
            # Set default parameters
            self.__set(*args, **kwargs)
            
            # Call solve
            self.solve()
    
    def __set(self, x, disregistry, gamma, axes, K_tensor, **kwargs):
        """
        Sets default values for class parameters.
        
        Parameters
        ----------
        x : numpy.ndarray
            An array of shape (N) giving the x coordinates corresponding to
            the disregistry solution.
        disregistry : numpy.ndarray
            A (N,3) array giving the initial disregistry vector guess at each
            x coordinate.
        gamma : atomman.defect.GammaSurface
            The gamma surface (stacking fault map) to use for computing the
            misfit energy.
        axes : numpy.ndarray
            A (3,3) array defining the crystal directions of the dislocation
            system. Used to orient the disregistry vector to the gamma surface
            in computing the misfit energy.
        K_tensor : numpy.ndarray
            A (3,3) array giving the anisotropic elastic energy coefficients
            corresponding to the dislocation system's orientation.  Can be
            computed with atomman.defect.Stroh.
        tau : numpy.ndarray, optional
            A (3,3) array giving the stress tensor to apply to the system
            using the stress energy term.  Only the xy, yy, and yz components
            are used.  Default value is all zeros.
        alpha : list of float, optional
            The alpha coefficient(s) used by the nonlocal energy term.  Default
            value is [0.0].
        beta : numpy.ndarray, optional
            The (3,3) array of beta coefficient(s) used by the surface energy
            term.  Default value is all zeros.
        cutofflongrange : float, optional
            The cutoff distance to use for computing the long-range energy.
            Default value is 1000 angstroms.
        burgers : numpy.ndarray, optional
            The (3,) array of the dislocation's Burgers vector relative to the
            dislocation system.  Used only by the long-range energy.  Default
            value is all zeros (long-range energy will be excluded).
        fullstress : bool, optional
            Flag indicating which stress energy algorithm to use.  Default
            value is True.
        cdiffelastic : bool, optional
            Flag indicating if the dislocation density for the elastic energy
            component is computed with central difference (True) or simply
            neighboring values (False).  Default value is False.
        cdiffsurface : bool, optional
            Flag indicating if the dislocation density for the surface energy
            component is computed with central difference (True) or simply
            neighboring values (False).  Default value is True.
        cdiffstress : bool, optional
            Flag indicating if the dislocation density for the stress energy
            component is computed with central difference (True) or simply
            neighboring values (False).  Only matters if fullstress is True.
            Default value is False.
        min_method : str, optional
            The scipy.optimize.minimize method to use.  Default value is
            'Powell'.
        min_options : dict, optional
            Any options to pass on to scipy.optimize.minimize. Default value
            is {}.
        """
        
        # Set mandatory parameters
        self.__x = x
        self.__disregistry = disregistry
        self.__gamma = gamma
        self.__axes = axes
        self.__K_tensor = K_tensor
        
        # Set optional keyword parameters
        self.__tau = kwargs.get('tau', np.zeros((3,3)))
        self.__alpha = kwargs.get('alpha', [0.0])
        if not isinstance(self.__alpha, list):
            self.__alpha = [self.__alpha]
        self.__beta = kwargs.get('beta', np.zeros((3,3)))
        self.__cutofflongrange = kwargs.get('cutofflongrange',
                                            uc.set_in_units(1000, 'angstrom'))
        self.__burgers = kwargs.get('burgers', np.zeros(3))
        self.__fullstress = kwargs.get('fullstress', True)
        self.__cdiffelastic = kwargs.get('cdiffelastic', False)
        self.__cdiffsurface = kwargs.get('cdiffsurface', True)
        self.__cdiffstress = kwargs.get('cdiffstress', False)
        self.__min_method = kwargs.get('min_method', 'Powell')
        self.__min_options = kwargs.get('min_options', {})
    
    def solve(self, **kwargs):
        """
        Solves the semidiscrete variational Peierls-Nabarro dislocation
        disregistry through energy minimization using the set class
        properties.  All parameters are optional keyword arguments that
        replace the class properties.
        
        Parameters
        ----------
        x : numpy.ndarray, optional
            An array of shape (N) giving the x coordinates corresponding to
            the disregistry solution.
        disregistry : numpy.ndarray, optional
            A (N,3) array giving the initial disregistry vector guess at each
            x coordinate.
        gamma : atomman.defect.GammaSurface, optional
            The gamma surface (stacking fault map) to use for computing the
            misfit energy.
        axes : numpy.ndarray, optional
            A (3,3) array defining the crystal directions of the dislocation
            system. Used to orient the disregistry vector to the gamma surface
            in computing the misfit energy.
        K_tensor : numpy.ndarray, optional
            A (3,3) array giving the anisotropic elastic energy coefficients
            corresponding to the dislocation system's orientation.  Can be
            computed with atomman.defect.Stroh.
        tau : numpy.ndarray, optional
            A (3,3) array giving the stress tensor to apply to the system
            using the stress energy term.  Only the xy, yy, and yz components
            are used.
        alpha : list of float, optional
            The alpha coefficient(s) used by the nonlocal energy term.
        beta : numpy.ndarray, optional
            The (3,3) array of beta coefficient(s) used by the surface energy
            term.
        cutofflongrange : float, optional
            The cutoff distance to use for computing the long-range energy.
        burgers : numpy.ndarray, optional
            The (3,) array of the dislocation's Burgers vector relative to the
            dislocation system.  Used only by the long-range energy.
        fullstress : bool, optional
            Flag indicating which stress energy algorithm to use.
        cdiffelastic : bool, optional
            Flag indicating if the dislocation density for the elastic energy
            component is computed with central difference (True) or simply
            neighboring values (False).
        cdiffsurface : bool, optional
            Flag indicating if the dislocation density for the surface energy
            component is computed with central difference (True) or simply
            neighboring values (False).
        cdiffstress : bool, optional
            Flag indicating if the dislocation density for the stress energy
            component is computed with central difference (True) or simply
            neighboring values (False).  Only matters if fullstress is True.
        min_method : str, optional
            The scipy.optimize.minimize method to use.
        min_options : dict, optional
            Any options to pass on to scipy.optimize.minimize.
        """
        
        # Define subfunctions
        def decompose(d):
            """Breaks disregistry into components for minimization"""
            d13 = np.concatenate([d[1:-1, 0], d[1:-1, 2]])
            first = d[0]
            last = d[-1]
            return d13, first, last

        def recompose(d13, first, last):
            """Reassembles the disregistry components"""
            half = int(len(d13)/2)
            d = np.zeros((half+2, 3))
            d[0] = first
            d[-1] = last
            d[1:-1, 0] = d13[:half]
            d[1:-1, 2] = d13[half:]
            return d

        def min_func(d13, first, last):
            """Function for minimizing"""
            disregistry = recompose(d13, first, last)
            return self.total_energy(disregistry=disregistry)
        
        # Handle parameters
        self.__x = kwargs.pop('x', self.__x)
        self.__disregistry = kwargs.pop('disregistry', self.__disregistry)
        self.__gamma = kwargs.pop('gamma', self.__gamma)
        self.__axes = kwargs.pop('axes', self.__axes)
        self.__K_tensor = kwargs.pop('K_tensor', self.__K_tensor)
        self.__tau = kwargs.pop('tau', self.__tau)
        self.__alpha = kwargs.pop('alpha', self.__alpha)
        if not isinstance(self.__alpha, list):
            self.__alpha = [self.__alpha]
        self.__beta = kwargs.pop('beta', self.__beta)
        self.__cutofflongrange = kwargs.pop('cutofflongrange', 
                                            self.__cutofflongrange)
        self.__burgers = kwargs.pop('burgers', self.__burgers)
        self.__fullstress = kwargs.pop('fullstress', self.__fullstress)
        self.__cdiffelastic = kwargs.pop('cdiffelastic', self.__cdiffelastic)
        self.__cdiffsurface = kwargs.pop('cdiffsurface', self.__cdiffsurface)
        self.__cdiffstress = kwargs.pop('cdiffstress', self.__cdiffstress)
        self.__min_method = kwargs.pop('min_method', self.__min_method)
        self.__min_options = kwargs.pop('min_options', self.__min_options)
        
        # Solve disregistry
        d13, first, last = decompose(self.disregistry)
        res = minimize(min_func, d13, args=(first, last),
                       method=self.min_method, options=self.min_options)
        self.__disregistry = recompose(res.x, first, last)
        
        self.__res = res
    
    @property
    def x(self):
        """numpy.ndarray : The x coordinates."""
        return deepcopy(self.__x)
    
    @property
    def disregistry(self):
        """numpy.ndarray : The disregistry vector for each x coordinate."""
        return deepcopy(self.__disregistry)
    
    @property
    def gamma(self):
        """atomman.defect.GammaSurface : The stacking fault map."""
        return deepcopy(self.__gamma)
    
    @property
    def axes(self):
        """numpy.ndarray : Axes for orienting the disregistry to gamma."""
        return deepcopy(self.__axes)
    
    @property
    def K_tensor(self):
        """numpy.ndarray : Dislocation energy coefficient tensor for system."""
        return deepcopy(self.__K_tensor)
    
    @property
    def tau(self):
        """numpy.ndarray : The applied stress state."""
        return deepcopy(self.__tau)
    
    @property
    def alpha(self):
        """list of float : Coefficients for nonlocal energy correction."""
        return deepcopy(self.__alpha)
    
    @property
    def beta(self):
        """numpy.ndarray : Coefficients for gradient energy correction."""
        return deepcopy(self.__beta)
    
    @property
    def cutofflongrange(self):
        """float : Cutoff distance for long-range elastic energy."""
        return self.__cutofflongrange
    
    @property
    def burgers(self):
        """numpy.ndarray : Full Burgers vector for long-range elastic energy."""
        return deepcopy(self.__burgers)
    
    @property
    def fullstress(self):
        """bool : Flag indicating which stress algorithm was used."""
        return self.__fullstress
    
    @property
    def cdiffelastic(self):
        """bool : Flag indicating if elastic energy used central difference for computing the dislocation density."""
        return self.__cdiffelastic
    
    @property
    def cdiffsurface(self):
        """bool : Flag indicating if surface energy used central difference for computing the dislocation density."""
        return self.__cdiffsurface
    
    @property
    def cdiffstress(self):
        """bool : Flag indicating if stress energy used central difference for computing the dislocation density."""
        return self.__cdiffstress
    
    @property
    def min_method(self):
        """str : scipy.optimize.minimize method used."""
        return self.__min_method
    
    @property
    def min_options(self):
        """dict : scipy.optimize.minimize options used."""
        return self.__min_options
        
    @property
    def res(self):
        """OptimizeResult : scipy.optimize.minimize result."""
        return self.__res
    
    def disldensity(self, x=None, disregistry=None, cdiff=False):
        """
        Computes the dislocation density as the numerical derivative of
        disregistry with respect to x.  Uses either neighboring values
        
            ρ[i] = (δ[i] - δ[i-1]) / (x[i] - x[i-1])
        
        or central difference
        
            ρ[i] = (δ[i+1] - δ[i-1]) / (x[i+1] - x[i-1])
        
        Parameters
        ----------
        cdiff : bool
            Flag indicating how to compute the derivative.  A value of False
            (default) compares nearest values while a value of True compares
            next-nearest (i.e., central difference).
        
        Returns
        -------
        newx : numpy.array
            The x positions corresponding to the dislocation density values.
        rho : numpy.array
            The computed dislocation density.
        """
        # Default values are class properties
        if x is None:
            x = self.__x
        if disregistry is None:
            disregistry = self.__disregistry
        
        # Extract values
        d = disregistry
        
        if cdiff is False:
            # ρ[i] = (δ[i] - δ[i-1]) / (x[i] - x[i-1])
            rho = ((d[1:] - d[:-1]).T / (x[1:] - x[:-1])).T
            
            # newx is all x except the first
            newx = x[1:]
        
        elif cdiff is True:
            # ρ[i] = (δ[i+1] - δ[i-1]) / (x[i+1] - x[i-1])
            rho = ((d[2:] - d[:-2]).T / (x[2:] - x[:-2])).T
            
            # newx is all x except the first and last
            newx = x[1:-1]
        else:
            raise TypeError('cdiff must be bool')
            
        return (newx, rho)
    
    def misfit_energy(self, x=None, disregistry=None):
        """
        Computes the misfit energy for the disregistry using the stored gamma
        surface
        
            E_misfit = Σ γ(δ)Δx
        
        Parameters
        ----------
        x : numpy.ndarray, optional
            x-coordinates.  Default value is the stored x-coordinates.
        disregistry : numpy.ndarray, optional
            (N, 3) shaped array of disregistry vectors at each x-coordinate.
            Default value is the stored disregistry values.
            
        Returns
        -------
        float
            The misfit energy for the dislocation.
        """
        # Default values are class properties
        if x is None:
            x = self.__x
        if disregistry is None:
            disregistry = self.__disregistry
        
        # Extract values
        d = disregistry
        dx = x[1] - x[0]
        axes = self.__axes
        gamma = self.__gamma
        
        # Strip out y-component of disregistry and transform for gamma
        disreg = np.vstack([d[:,0], np.zeros(len(d)), d[:,2]]).T
        T = axes_check(axes)
        pos = np.inner(disreg, T.T)
        
        # Σ γ(δ)Δx
        return dx * gamma.E_gsf(pos=pos).sum()
    
    def elastic_energy(self, x=None, disregistry=None):
        """
        Computes the short-range configuration-dependent elastic energy term
        for the dislocation based on the dislocation density and K_tensor.
        
            E_elastic = 1/(4π) Σ_i Σ_j χ(i,j,Δx) K_lm ρ_l[i] ρ_m[j]
            
            χ(i,j,Δx) = (3/2) Δx² + ψ(i-1,j-1,Δx) + ψ(i,j,Δx) - ψ(i,j-1,Δx) - ψ(j,i-1,Δx)
            
            ψ(i,j,Δx) = (1/2) (i-j)² Δx² ln(\|i-j\|Δx)
        
        Parameters
        ----------
        x : numpy.ndarray, optional
            x-coordinates.  Default value is the stored x-coordinates.
        disregistry : numpy.ndarray, optional
            (N, 3) shaped array of disregistry vectors at each x-coordinate.
            Default value is the stored disregistry values.
        
        Returns
        -------
        float
            The elastic energy for the dislocation.
        """
        # Define subfunctions
        def chi(i, j, dx):
            """
            Computes the chi subfunction:
                χ(i,j,Δx) = (3/2) Δx² + ψ(i-1,j-1,Δx) + ψ(i,j,Δx)
                                      - ψ(i,j-1,Δx) - ψ(j,i-1,Δx)
            """
            return 3./2. * dx**2 + (psi(i-1, j-1, dx) + psi(i, j, dx)
                                    - psi(i, j-1, dx) - psi(j, i-1, dx))

        def psi(i, j, dx):
            """
            Computes the psi subfunction:
                ψ(i,j,Δx) = (1/2) (i-j)² Δx² ln(|i-j|Δx)
            """
            # Suppress NaN runtime warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                
                # ψ(i,j,Δx) = (1/2) (i-j)² Δx² ln(|i-j|Δx)
                p = 0.5 * (i - j)**2 * dx**2 * np.log(np.abs(i - j) * dx)
            
            # Replace NaN values with 0.0
            p[np.isnan(p)] = 0.0
            
            return p
        
        # Default values are class properties
        if x is None:
            x = self.__x
        if disregistry is None:
            disregistry = self.__disregistry
        
        # Extract values
        d = disregistry
        dx = x[1] - x[0]
        cdiff = self.cdiffelastic
        Kij = self.__K_tensor
        
        rho = self.disldensity(x=x, disregistry=d, cdiff=cdiff)[1]
        j = np.arange(len(rho), dtype=int)
        
        energy = 0.0
        
        # Compute elastic energy (looping over i, vectorization over j)
        # 1/(4π) Σ_i Σ_j χ(i,j,Δx) K_lm ρ_l[i] ρ_m[j]
        for i in j:
            energy += np.sum(chi(i, j, dx) * np.inner(rho[i].dot(Kij), rho)) / (4 * np.pi)
        
        return energy
    
    def longrange_energy(self):
        """
        Computes the long-range elastic energy term for the dislocation using
        the K_tensor, Burgers vector and long-range cutoff.  This term is
        configuration-independent and thus the method takes no parameters.
        
            E_longrange = 1/(2π) K_lm b_l b_m ln(L)
        
        Returns
        -------
        float
            The long-range energy for the dislocation.
        """
        # Extract values
        Kij = self.__K_tensor
        b = self.__burgers
        L = self.__cutofflongrange
        
        # Compute long-range energy
        # 1/(2π) K_lm b_l b_m ln(L)
        return np.inner(b.dot(Kij), b) * np.log(L) / (2 * np.pi)
    
    def stress_energy(self, x=None, disregistry=None):
        """
        Computes the stress energy due to the applied stress, tau.
        If fullstress is True:
        
            E_stress = -1/2 Σ_i (x[i]² - x[i-1]²) ρ_l τ_2l
        
        If fullstress is False:
        
            E_stress = -1/2 Σ_i τ_2l (δ_l[i] + δ_l[i+1]) Δx
        
        Parameters
        ----------
        x : numpy.ndarray, optional
            x-coordinates.  Default value is the stored x-coordinates.
        disregistry : numpy.ndarray, optional
            (N, 3) shaped array of disregistry vectors at each x-coordinate.
            Default value is the stored disregistry values.
        
        Returns
        -------
        float
            The stress energy for the dislocation.
        """
        # Default values are class properties
        if x is None:
            x = self.__x
        if disregistry is None:
            disregistry = self.__disregistry
        
        # Extract values
        d = disregistry
        dx = x[1] - x[0]
        tau = self.__tau
        full = self.fullstress
        cdiff = self.cdiffstress
        
        if full is True:
            rho = self.disldensity(x=x, disregistry=d, cdiff=cdiff)[1]
            # -1/2 Σ_i (x[i]² - x[i-1]²) ρ_l τ_2l
            return -0.5 * np.sum((x[1:]**2 - x[:-1]**2) * np.inner(rho, tau[1,:]))
            
        else:
            # -1/2 Σ_i τ_2l (δ_l[i] + δ_l[i+1]) Δx
            return -0.5 * np.sum(np.inner(tau[1,:], (d[:-1] + d[1:]) * dx))
    
    def surface_energy(self, x=None, disregistry=None):
        """
        Computes the gradient surface energy correction using beta
        coefficients.
        
            E_surface = Σ_j β_lj / 4 Σ_i ρ_l[i]² Δx
        
        Parameters
        ----------
        x : numpy.ndarray, optional
            x-coordinates.  Default value is the stored x-coordinates.
        disregistry : numpy.ndarray, optional
            (N, 3) shaped array of disregistry vectors at each x-coordinate.
            Default value is the stored disregistry values.
        
        Returns
        -------
        float
            The surface energy for the dislocation.
        """
        # Default values are class properties
        if x is None:
            x = self.__x
        if disregistry is None:
            disregistry = self.__disregistry
            
        # Extract values
        d = disregistry
        dx = x[1] - x[0]
        beta = self.__beta
        cdiff = self.__cdiffsurface
        
        rho = self.disldensity(x=x, disregistry=d, cdiff=cdiff)[1]
        
        # Σ_j β_lj / 4 Σ_i ρ_l[i]² Δx
        return np.sum( np.inner(rho**2 * dx, beta) ) / 4
    
    def nonlocal_energy(self, x=None, disregistry=None):
        """
        Computes the nonlocal energy correction using alpha coefficient(s).
        
            E_nonlocal = Σ_m α_m Σ_i δ[i] (δ[i] - (δ[i+m] + δ[i-m]) / 2) Δx
        
        Parameters
        ----------
        x : numpy.ndarray, optional
            x-coordinates.  Default value is the stored x-coordinates.
        disregistry : numpy.ndarray, optional
            (N, 3) shaped array of disregistry vectors at each x-coordinate.
            Default value is the stored disregistry values.
        
        Returns
        -------
        float
            The nonlocal energy for the dislocation.
        """
        # Default values are class properties
        if x is None:
            x = self.__x
        if disregistry is None:
            disregistry = self.__disregistry
            
        # Extract values
        d = disregistry
        dx = x[1] - x[0]
        alpha = self.__alpha
        
        energy = 0.0
        
        # Σ_m α_m Σ_i δ[i] (δ[i] - (δ[i+m] + δ[i-m]) / 2) Δx
        for num, a in enumerate(alpha):
            m = num + 1
            dd = d[m:-m] - 0.5 * (d[2*m:] + d[:-2*m])
            energy += a * np.sum(d[m:-m] * dd * dx)
            
        return energy
    
    def total_energy(self, x=None, disregistry=None):
        """
        Computes the total energy for the dislocation.
        
            E_total = E_elastic + E_misfit + E_stress + E_surface + E_nonlocal
        
        Parameters
        ----------
        x : numpy.ndarray, optional
            x-coordinates.  Default value is the stored x-coordinates.
        disregistry : numpy.ndarray, optional
            (N, 3) shaped array of disregistry vectors at each x-coordinate.
            Default value is the stored disregistry values.
        
        Returns
        -------
        float
            The total energy for the dislocation.
        """
        # Default values are class properties
        if x is None:
            x = self.__x
        if disregistry is None:
            disregistry = self.__disregistry
            
        return (self.misfit_energy(x, disregistry)
                + self.elastic_energy(x, disregistry)
                + self.longrange_energy()
                + self.stress_energy(x, disregistry)
                + self.nonlocal_energy(x, disregistry)
                + self.surface_energy(x, disregistry))
    
    def load(self, model, gamma=None):
        """
        Load solution from a data model.
        
        Parameters
        ----------
        model : str or DataModelDict
            The semi-discrete-Peierls-Nabarro data model to load.
        gamma : atomman.defect.GammaSurface
            The gamma surface to use.
        """
        
        # Identify model element
        try:
            sdpn = DM(model).find('semidiscrete-variational-Peierls-Nabarro')
        except:
            sdpn = DM(model).find('semi-discrete-Peierls-Nabarro')
        # Load calculation parameters
        params = sdpn['parameter']
        self.__axes = uc.value_unit(params['axes'])
        self.__K_tensor = uc.value_unit(params['K_tensor'])
        self.__tau = uc.value_unit(params['tau'])
        self.__alpha = uc.value_unit(params['alpha'])
        if not isinstance(self.__alpha, list):
            self.__alpha = [self.__alpha]
        self.__beta = uc.value_unit(params['beta'])
        self.__cutofflongrange = uc.value_unit(params['cutofflongrange'])
        self.__burgers = uc.value_unit(params['burgers'])
        self.__fullstress = params['fullstress']
        self.__cdiffelastic = params['cdiffelastic']
        self.__cdiffsurface = params['cdiffsurface']
        self.__cdiffstress = params['cdiffstress']
        self.__min_method = params['min_method']
        self.__min_options = params['min_options']
        
        # Load gamma
        if gamma is None:
            try:
                gamma = GammaSurface(model)
            except:
                raise ValueError('No gamma surface in model or given!')
        elif not isinstance(gamma, GammaSurface):
            gamma = GammaSurface(gamma)
        self.__gamma = gamma
        
        # Load calculation solution
        solution = sdpn['solution']
        self.__x = uc.value_unit(solution['x'])
        self.__disregistry = uc.value_unit(solution['disregistry'])
    
    def model(self, length_unit='angstrom', energy_unit='eV',
              pressure_unit='GPa', include_gamma=False):
        """
        Generate a data model for the object.
        
        Parameters
        ----------
        length_unit : str
            The unit of length to save values as.  Default is 'angstrom'.
        energy_unit : str
            The unit of energy to save values as.  Default is 'eV'.
        pressure_unit : str
            The unit of pressure to save values as.  Default is 'GPa'.
        include_gamma : bool
            Flag indicating if the gamma surface data is to be included.
            Default value is False.
        """
        model = DM()
        model['semidiscrete-variational-Peierls-Nabarro'] = sdpn = DM()
        
        sdpn['parameter'] = params = DM()
        params['axes'] = uc.model(self.__axes, None)
        params['K_tensor'] = uc.model(self.__K_tensor, pressure_unit)
        params['tau'] = uc.model(self.__tau, pressure_unit)
        params['alpha'] = uc.model(self.__alpha, pressure_unit+'/'+length_unit)
        
        params['beta'] = uc.model(self.__beta, pressure_unit+'*'+length_unit)
        params['cdiffelastic'] = self.cdiffelastic
        params['cdiffsurface'] = self.cdiffsurface
        params['cdiffstress'] = self.cdiffstress
        params['cutofflongrange'] = uc.model(self.__cutofflongrange, length_unit)
        params['burgers'] = uc.model(self.__burgers, length_unit)
        params['fullstress'] = self.fullstress
        params['min_method'] = self.min_method
        params['min_options'] = self.min_options
        
        if include_gamma is True:
            sdpn['generalized-stacking-fault'] = self.__gamma.model(
                                                  length_unit=length_unit,
                                                  energy_unit=energy_unit,
                                                  pressure_unit=pressure_unit)
        
        sdpn['solution'] = solution = DM()
        solution['x'] = uc.model(self.__x, length_unit)
        solution['disregistry'] = uc.model(self.__disregistry, length_unit)
        
        return model