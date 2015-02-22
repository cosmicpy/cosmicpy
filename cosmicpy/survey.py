# Copyright (c) 2014-2015, CosmicPy Developers
# Licensed under CeCILL 2.1 - see LICENSE.rst
"""
:mod:`cosmicpy.survey` -- Galaxy survey representation
=====================================================

.. module:: cosmicpy.survey
    :synopsis: Contains surveys related objects and functions
.. moduleauthor:: Francois Lanusse <francois.lanusse@cea.fr>
.. moduleauthor:: Anais Rassat <anais.rassat@epfl.ch>
.. Created on Jun 13, 2013 by Francois Lanusse

.. autosummary::

    survey
    redshift_bin
    nuisance
    grid_nuisance

"""
from numpy import *
from scipy.integrate import romberg
from scipy.optimize import brentq
from scipy.special import erf
from scipy.interpolate import interp2d
from abc import ABCMeta, abstractmethod
from .utils import *


class nuisance(object):
    r"""Generic nuisance representation of a 2D function"""
    __metaclass__ = ABCMeta

    def __init__(self, name, function, Nx, Ny, xlim, ylim):
        self._name = name
        self._xlim = xlim
        self._ylim = ylim
        self._Ny = Ny
        self._Nx = Nx
        self._params_values = zeros((Nx, Ny))
        self._params_amplitude = 1.0
        self._model(function)

    def update(self, **kwargs):
        r"""Update the nuisance parameters"""
        if len(kwargs) is 0:
            return
        # Looks for keywords and update the survey parameters
        for kw in kwargs:
            if self._name in kw:
                # Extract the number of the parameter
                n = int(kw.split('_p')[1])
                if n == 0:
                    self._params_amplitude = kwargs[kw]
                else:
                    self._params_values.flat[n - 1] = kwargs[kw]

    @property
    def Np(self):
        r"""
        Number of nuisance parameters"""
        return self._Nx * self._Ny + 1

    def get_value(self, ind):
        if ind == 0:
            return self._params_amplitude
        else:
            return self._params_values.flat[ind - 1]

    @abstractmethod
    def __call__(self, x, y):
        r""" Evaluates the function for the given set of nuisance parameters"""
        pass

    @abstractmethod
    def _model(self, func):
        r""" Models the given function using the nuisance parameters"""
        pass


class grid_nuisance(nuisance):
    r""" Generic nuisance parametrisation on a grid
    """
    def _model(self, func):
        # Compute the position of the grid nodes
        self._x = logspace(log10(self._xlim[0]), log10(self._xlim[1]),
                           self._Nx + 2)
        self._y = logspace(log10(self._ylim[0]), log10(self._ylim[1]),
                           self._Ny + 2)
        self._params_values *= 0
        self._border = zeros((self._Nx + 2, self._Ny + 2))
        # for i in range(self._Nx + 2):
        #    for j in range(self._Ny + 2):
        #        self._border[i, j] = func(self._x[i], self._y[j])
        self._params_amplitude = 1.0
        self._func = func

    def __call__(self, x, y):
        x = atleast_1d(x)
        y = atleast_1d(y)
        # Compute the bilinear interpolation of the logarithmique grid at the
        # requested position
        self._border[1:-1, 1:-1] = self._params_values[:, :]
        lnQ = interp2d(log10(self._y), log10(self._x), self._border,
                       copy=False, bounds_error=False, fill_value=0)
        return self._params_amplitude * exp(lnQ(log10(y), log10(x)))*self._func(x, y)


class survey(object):
    r""" Contains all the parameters corresponding to a survey.

    Parameters
    ----------
    nzparams : dict, optional
        Dictionary containing the parameters of the desired redshift
        distribution.
    nzbins : int, optional
        Number of redshift bins in the survey (def: 1).
    bintype : str, optional
        Type of redshift binning. Either 'eq_dens' or 'eq_size'
        (def: 'eq_dens').
    ngal : float, optional
        Number of galaxies per square arcmins (def: 40).
    zmin : float, optional
        Minimum redshift considered in the survey (def: 0).
    zmax : float, optional
        Maximum redshift considered in the survey (def: 5).
    chicut : float, optional
        Level of the cut to apply to the selection function, relative to
        the maximum of the selection fuction. If 0, :math:`\chi_\max` is set to
        :math:`\chi(z_\max)` (def : 1e-5).
    zphot_sig : float, optional
        Standard variation of photometric redshift errors (def: 0.05).
    zphot_bias : float, optional
        Bias of photometric redshift errors (def: 0).
    fsky : float, optional
        Sky fraction covered by the survey (def: 0.4848).
    biastype : str, optional
        Type of galaxy bias. Either 'constant' or 'sqrt' (def: 'sqrt')

    Notes
    -----
    Several parametrisations for the redshifts can be specified:
        * 'smail'    : n(z) = z^a exp( - (z/z0)^b)
        * 'gaussian' : n(z) = exp(-(r/r0)^2) * r ^2 * dr/dz
    """

    def __init__(self, nuisances=[], **kwargs):
        self._nzparams = {'type': 'smail',
                          'a': 2.0,
                          'b': 1.5,
                          'z0': 0.9 / sqrt(2)}
        self._zmin = 0
        self._zmax = 5
        self._chicut = 1e-5
        self._nzbins = 1
        self._zphot_sig = 0.05
        self._zphot_bias = 0
        self._ngal = 40
        self._fsky = 0.4848
        self._bintype = 'eq_dens'
        self._biastype = 'sqrt'
        self.nuisances = {}

        self._norm = None
        self._zmean = None
        self._zmed = None

        self.update(**kwargs)

    def addNuisance(self, params):
        """ Adds nuisance parameters to the survey
            params is of the form
            {'name':'bias', 'type':'grid', 'Nx':3,'Ny':4,'xlim':[0.001,1],'ylim':[0.001,0.1]}
        """
        if params['name'] == 'bias':
            if params['type'] == 'grid':
                if self._biastype == 'sqrt':
                    def func(zp, k):
                        zp = atleast_1d(zp)
                        k = atleast_1d(k)
                        return einsum('i,j->ij', sqrt(zp), ones(len(k)))
                else:
                    def func(zp, k):
                        zp = atleast_1d(zp)
                        k = atleast_1d(k)
                        return ones((len(zp), len(k)))

                nuisance = grid_nuisance('bias', func, params['Nx'],
                                          params['Ny'], array(params['xlim'])+1, params['ylim'])
                self.nuisances['bias'] = nuisance

    def update(self, **kwargs):
        r""" Updates the survey after parameters have been modified.

        Parameters
        ----------
        nzparams : dict, optional
            Dictionary containing the parameters of the desired redshift
            distribution.
        nzbins : int, optional
            Number of redshift bins in the survey.
        bintype : str, optional
            Type of redshift binning. Either 'eq_dens' or 'eq_size'.
        ngal : float, optional
            Number of galaxies per square arcmins.
        zmin : float, optional
            Minimum redshift considered in the survey.
        zmax : float, optional
            Maximum redshift considered in the survey.
        chicut : float, optional
            Level of the cut to apply to the selection function, relative to
            the maximum of the selection fuction. If 0, :math:`\chi_\max` is
            set to :math:`\chi(z_\max)` (def : 1e-5).
        zphot_sig : float, optional
            Standard variation of photometric redshift errors.
        zphot_bias : float, optional
            Bias of photometric redshift errors.
        fsky : float, optional
            Sky fraction covered by the survey.
        biastype : str, optional
            Type of galaxy bias. Either 'constant' or 'sqrt'.
        """
        if len(kwargs) is 0:
            return
        # Looks for keywords and update the survey parameters
        for kw in kwargs:
            if kw == 'nzbins':
                self._nzbins = kwargs[kw]
            elif kw == 'nzparams':
                self._nzparams = kwargs[kw]
            elif kw == 'zmin':
                self._zmin = kwargs[kw]
            elif kw == 'zmax':
                self._zmax = kwargs[kw]
            elif kw == 'chicut':
                self._chicut = kwargs[kw]
            elif kw == 'zphot_sig':
                self._zphot_sig = kwargs[kw]
            elif kw == 'zphot_bias':
                self._zphot_bias = kwargs[kw]
            elif kw == 'ngal':
                self._ngal = kwargs[kw]
            elif kw == 'bintype':
                self._bintype = kwargs[kw]
            elif kw == 'fsky':
                self._fsky = kwargs[kw]
            elif kw == 'biastype':
                self._biastype = kwargs[kw]

        for n in self.nuisances.values():
            n.update(**kwargs)

        self._norm = None
        self._zmean = None
        self._zmed = None

        # Updates the redshift bins
        nztot = redshift_bin(self, zphot_min=self._zmin, zphot_max=self._zmax)
        self.zbins = nztot.subdivide(self._nzbins, bintype=self._bintype)

    def nz_unorm(self, z):
        r""" Computes the non normalized n(z)

        Parameters
        ----------
        z : array_like
            Redshift

        Returns
        -------
        n : ndarray, or float if input scalar
            Non normalised redshift distribution evaluated at the specified
            redshift
        """
        if self._nzparams['type'] == 'smail':
            a = self._nzparams['a']
            b = self._nzparams['b']
            z0 = self._nzparams['z0']
            p = z**a * exp(-(z / z0)**b)
        elif self._nzparams['type'] == 'gaussian':
            r0 = self._nzparams['r0']
            cosmo = self._nzparams['cosmology']
            a = z2a(z)
            r = cosmo.a2chi(a)
            p = exp(-(r / r0)**2) * r**2 * (cosmo.dchioverda(a) /
                                            cosmo.dzoverda(a))
        else:
            print ('ERROR : unknown redshift bin type ' + self.params['type'])
            p = 0
        return p

    def nz(self, z):
        r""" Computes the normalized n(z)

        Parameters
        ----------
        z : array_like
            Redshift

        Returns
        -------
        n : ndarray, or float if input scalar
            Normalised redshift distribution evaluated at the specified
            redshift

        See Also
        --------
        nz_unorm : Non normalised n(z)
        """
        return self.nz_unorm(z) / self.norm

    def phi(self, cosmo, chi):
        r""" Computes the normalised survey selection function at
        a given comoving distance.

        Parameters
        ----------
        cosmo : cosmology
            Cosmology object required to convert comoving distances.
        chi : array_like
            Comoving distance.

        Returns
        -------
        phi : ndarray, or float if input is scalar
            Normalised selection function of the survey.

        Notes
        -----
        The survey selection function :math:`\phi` verifies:

        .. math::

            \int d^3 r \phi (r) = V

        where :math:`V` is a characteristic volume of the survey chosen such
        that :math:`\phi \rightarrow 1` as :math:`V \rightarrow \infty`.

        .. warning::

            This function returns the selection function scaled by
            the survey volume i.e. :math:`\frac{\phi}{V}`
        """
        chi = atleast_1d(chi)
        a = cosmo.chi2a(chi)
        phi = self.nz(a2z(a)) / (4.*pi*chi**2) * (cosmo.dzoverda(a) /
                                                  cosmo.dchioverda(a))
        # Apply zmax cut
        phi[a2z(a) > self.zmax] = 0

        return phi

    def chimax(self, cosmo):
        r""" Computes the maximum comoving distance for the survey by
        corresponding to the point where the selection function reaches the cut
        level.

        Parameters
        ----------
        cosmo : cosmology
            Cosmology object required to convert comoving distances.
        cut : float, optional
            Level of the cut to apply to the selection function, relative to
            the maximum of the selection fuction (def : 1e-5)

        Returns
        -------
        chi : float
            Maximum survey comoving distance
        """
        # If chicut is set to 0, chimax is set to zmax
        if self._chicut <= 0:
            return cosmo.a2chi(z2a(self.zmax)).tolist()
        # First, find approximate maximum of the selection function
        chi = cosmo.a2chi(z2a(self.zmax))
        chi = linspace(1.0, chi)
        phitab = self.phi(cosmo, chi)
        m = log10(phitab.max()*self._chicut)
        phi = lambda x: log10(self.phi(cosmo, x)) - m
        chimax = brentq(phi, chi[argmax(phitab)],
                        cosmo.a2chi(cosmo._amin))

        if chimax >= cosmo.a2chi(cosmo._amin):
            print("Warning : reaching maximum tabulated comoving distance")
        return chimax

    @property
    def norm(self):
        r""" Normalisation factor of the redshift distribution"""
        if self._norm is None:
            if self._nzparams['type'] == 'smail':
                z0 = self._nzparams['z0']
            elif self._nzparams['type'] == 'gaussian':
                r0 = self._nzparams['r0']
                cosmo = self._nzparams['cosmology']
                z0 = a2z(cosmo.chi2a(r0))
            else:
                print('ERROR : unknown redshift bin type ' +
                      self.params['type'])
            # The normalisation factor is computed in 2 steps to avoid failure
            # of the integration scheme
            self._norm = romberg(self.nz_unorm, self.zmin, z0) + \
                romberg(self.nz_unorm, z0, self.zmax)
        return self._norm

    def bias(self, z, k):
        r""" Galaxy bias b(z,k)

        Parameters
        ----------
        z : array_like
            Redshift

        Returns
        -------
        b : ndarray, or float if input is scalar
            Galaxy bias evaluated at the specified redshift

        Notes
        -----
        The galaxy bias follows the
        """
        z = atleast_1d(z)
        k = atleast_1d(k)
        if 'bias' in self.nuisances:
            return self.nuisances['bias'](1+z, k)

        if self._biastype == 'sqrt':
            return einsum('i,j->ij', sqrt(1 + z), ones(len(k)))
        else:
            return ones((len(z), len(k)))

    @property
    def zmed(self):
        r""" Median of the redshift distribution.
        """
        if self._zmed is None:
            f = lambda x: romberg(self.nz, self.zmin, x) - 0.5
            self._zmed = brentq(f, self.zmin, self.zmax)
        return self._zmed

    @property
    def zmean(self):
        r""" Mean of the redshift distribution.
        """
        if self._zmean is None:
            self._zmean = romberg(lambda z: z * self.nz(z), self.zmin, self.zmax)
        return self._zmean

    @property
    def zmin(self):
        r""" Minimun Redshift
        """
        return self._zmin

    @zmin.setter
    def zmin(self, val):
        r""" Sets the minimum redshift
        """
        self.update(wa=val)

    @property
    def zmax(self):
        return self._zmax

    @property
    def zphot_sig(self):
        return self._zphot_sig

    @property
    def zphot_bias(self):
        return self._zphot_bias

    @property
    def ngal(self):
        return self._ngal

    @property
    def N(self):
        res = self._ngal
        res *= 60.0*60.0 # gal/deg^2
        res *= (180.0/pi)**2 # gal/stradian
        res *= 4.0 * pi * self.fsky
        return res

    @property
    def nzbins(self):
        return self._nzbins

    @property
    def fsky(self):
        return self._fsky


class redshift_bin(object):
    """ Represents a redshift bin

    """

    def __init__(self, surv, **kwargs):

        self.surv = surv
        self.zphot_min = surv.zmin
        self.zphot_max = surv.zmax

        for kw in kwargs:
            if kw == 'zphot_min':
                self.zphot_min = kwargs[kw]
            elif kw == 'zphot_max':
                self.zphot_max = kwargs[kw]
            else:
                print ('WARNING : Unknown keyword ' + kw)

        self.zmin = max([self.zphot_min - 10*self.surv.zphot_sig, 0])
        self.zmax = self.zphot_max + 10*self.surv.zphot_sig

        self._norm = None
        self._zmean = None
        self._zmed = None

    def nz_unorm(self, z):
        """ Computes the un-normalized n(z) """
        p = self.surv.nz(z)

        # Apply photo-z errors
        x = 1.0/(sqrt(2.0)*self.surv.zphot_sig*(1.0+z))
        res = 0.5* p *( erf((self.zphot_max - z + self.surv.zphot_bias)*x) - erf((self.zphot_min - z + self.surv.zphot_bias)*x))
        return res

    def nz(self, z):
        return self.nz_unorm(z)/self.norm

    def subdivide(self, nbins, bintype='eq_dens'):
        """ Divide this redshift bins into sub-bins 
            nbins : Number of bins to generate
            bintype : 'eq_dens' or 'eq_size'
        """
        # Compute the redshift boundaries for each bin generated
        zbounds = [self.zphot_min]
        bins = []
        n_per_bin = self.norm / nbins
        for i in range(nbins-1):
            if bintype == 'eq_dens':
                    zbound = brentq(lambda z: romberg(self.nz,self.zphot_min,z) - (i+1.0)*n_per_bin, zbounds[i], self.zmax)
            else:
                if bintype != 'eq_size':
                    print ('WARNING : unknown binning scheme ' + bintype + '. Assuming equal size bins')
                zbound = (i+1.) * (self.zphot_max - self.zphot_min)/nbins

            zbounds.append(zbound)

            new_bin = redshift_bin(self.surv, zphot_min=zbounds[i],zphot_max=zbounds[i+1])

            bins.append(new_bin)

        zbounds.append(self.zmax)
        new_bin = redshift_bin(self.surv, zphot_min=zbounds[nbins-1],zphot_max=zbounds[nbins])
        bins.append(new_bin)

        return bins

    def phi(self, cosmo, chi):
        r""" Computes the normalised selection function for the bin at
        a given comoving distance.

        Parameters
        ----------
        cosmo : cosmology
            Cosmology object required to convert comoving distances.
        chi : array_like
            Comoving distance.

        Returns
        -------
        phi : ndarray, or float if input is scalar
            Normalised selection function of the bin.

        Notes
        -----
        The normalised survey selection function :math:`\phi` verifies:

        .. math::

            \int d^3 r \phi (r) = 1
        """
        chi = atleast_1d(chi)
        a = cosmo.chi2a(chi)
        phi = self.nz(a2z(a)) / (4.*pi*chi**2) * (cosmo.dzoverda(a) /
                                                  cosmo.dchioverda(a))
        # Apply zmax cut
        phi[a2z(a) > self.zmax] = 0

        return phi

    @property
    def kmax_lin(self):
        """ Maximum linear scale at the median redshift in h^-1 Mpc"""
        return min(0.132*self.zmed, 0.25)

    @property
    def zmed(self):
        """ Median of the redshift distribution """
        if self._zmed is None:
            f = lambda x: romberg(self.nz, self.zmin, x) - 0.5
            self._zmed = brentq(f, self.zmin, self.zmax)

        return self._zmed

    @property
    def zmean(self):
        """ Mean of the distribution """
        if self._zmean is None:
            self._zmean = romberg(lambda z: z*self.nz(z), self.zmin, self.zmax)

        return self._zmean

    @property
    def norm(self):
        """ Normalisation of the distribution """
        if self._norm is None:
            self._norm = romberg(self.nz_unorm, self.zmin, self.zmax)
        return self._norm

    @property
    def ngal(self):
        """ Returns the number of galaxies in this photo-z bin in steradian"""
        res = self.surv.ngal*self.norm
        res *= 60.0*60.0 #gal/deg^2
        res *= (180.0/pi)**2 #gal/stradian
        return res