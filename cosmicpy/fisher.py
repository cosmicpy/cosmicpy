# Copyright (c) 2014-2015, CosmicPy Developers
# Licensed under CeCILL 2.1 - see LICENSE.rst
"""
:mod:`cosmicpy.fisher` -- Fisher forecasts
=========================================

.. module:: cosmicpy.fisher
    :synopsis: Computes Fisher Matrices
.. moduleauthor:: Francois Lanusse <francois.lanusse@cea.fr>
.. moduleauthor:: Anais Rassat <anais.rassat@epfl.ch>
.. Created on Jul 9, 2013 by Francois Lanusse

.. autosummary::

    fisher
    fisherTomo
    fisher3d

"""

from .cosmology import cosmology
from .spectra import spectra
from .utils import *
from copy import deepcopy
import itertools
from numpy import *
from scipy.integrate import *
from scipy.optimize import brentq
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.interpolate import interp1d
from abc import ABCMeta, abstractmethod


class fisher(object):
    """
    Base class to perform a Fisher Analysis from specified cosmology,
    survey and cosmological parameters.
    """
    __metaclass__ = ABCMeta

    def __init__(self, fid_cosmo, fid_survey, params, margin_params=[]):
        """
        Constructor
        """

        self.step = 0.003
        self.fid_cosmo = fid_cosmo
        self.fid_surv = fid_survey
        self.params = []

        # Check that the parameters provided are present in survey or cosmo
        for p in params:
            # First find the fiducial value for the parameter in question
            if p in dir(self.fid_surv) or p in dir(self.fid_cosmo):
                self.params.append(p)
            else:
                print("Warning, unknown parameter in derivative :" + p)

        # Checks that the marginalisation parameters are actually considered
        # in the Fisher analysis
        self.margin_params = []
        for p in margin_params:
            if self.fid_surv.nuisances:
                self.margin_params.append(p)
            else:
                print("Warning, requested marginalisation parameter " + p +
                      " is not included in the analysis")

        # Create a spectra object with a copy of the fiducial cosmology
        self.spectra = spectra(deepcopy(self.fid_cosmo),
                               deepcopy(self.fid_surv))

        # Precomputed Fisher matrix
        self._fullMat = None
        self._fullInvMat = None
        self._mat = None
        self._invmat = None

    @abstractmethod
    def _computeObservables(self):
        pass

    @abstractmethod
    def _computeFullMatrix(self):
        pass

    def Fij(self, param_i, param_j):
        """
            Returns the matrix element of the Fisher matrix for parameters
            param_i and param_j
        """
        i = self.params.index(param_i)
        j = self.params.index(param_j)

        return self.mat[i, j]

    def invFij(self, param_i, param_j):
        """
            Returns the matrix element of the inverse Fisher matrix for
            parameters param_i and param_j
        """
        i = self.params.index(param_i)
        j = self.params.index(param_j)

        return self.invmat[i, j]

    def sigma_fix(self, param):
        return 1.0 / sqrt(self.Fij(param, param))

    def sigma_marg(self, param):
        return sqrt(self.invFij(param, param))

    def sub_matrix(self, subparams):
        """
        Extracts a submatrix from the current fisher matrix using the
        parameters in params
        """
        params = []
        for p in subparams:
            # Checks that the parameter exists in the orignal matrix
            if p in self.params:
                params.append(p)
            else:
                print("Warning, parameter not present in original \
                    Fisher matrix, left ignored :" + p)
        newFisher = fisher(self.fid_cosmo, self.fid_surv, params)

        # Fill in the fisher matrix from the precomputed matrix
        newFisher._mat = zeros((len(params), len(params)))

        for i in range(len(params)):
            indi = self.params.index(params[i])
            for j in range(len(params)):
                indj = self.params.index(params[j])
                newFisher._mat[i, j] = self.mat[indi, indj]

        newFisher._invmat = linalg.inv(newFisher._mat)

        return newFisher

    def _marginalise(self, params):
        r""" Marginalises the Fisher matrix over unwanted parameters.

        Parameters
        ----------
        params: list
            List of parameters that should not be marginalised over.

        Returns
        -------
        (mat, invmat): ndarray
            Marginalised Fisher matrix and its invers
        """
        # Builds inverse matrix
        marg_inv = zeros((len(params), len(params)))
        for i in range(len(params)):
            indi = self.params.index(params[i])
            for j in range(len(params)):
                indj = self.params.index(params[j])
                marg_inv[i, j] = self.invmat[indi, indj]

        marg_mat = linalg.pinv(marg_inv)

        return (marg_mat, marg_inv)

    def corner_plot(self, nstd=2, labels=None, **kwargs):
        r""" Makes a corner plot including all the parameters in the Fisher analysis
        """

        if labels is None:
            labels = self.params

        for i in range(len(self.params)):
            for j in range(i):
                ax = plt.subplot(len(self.params)-1, len(self.params)-1 , (i - 1)*(len(self.params)-1) + (j+1))
                if i == len(self.params) - 1:
                    ax.set_xlabel(labels[j])
                else:
                    ax.set_xticklabels([])
                if j == 0:
                    ax.set_ylabel(labels[i])
                else:
                    ax.set_yticklabels([])

                self.plot(self.params[j], self.params[i], nstd=nstd, ax=ax, **kwargs)

        plt.subplots_adjust(wspace=0)
        plt.subplots_adjust(hspace=0)

    def plot(self, p1, p2, nstd=2, ax=None, **kwargs):
        r""" Plots confidence contours corresponding to the parameters
        provided.

        Parameters
        ----------
        """
        params = [p1, p2]

        def eigsorted(cov):
            vals, vecs = linalg.eigh(cov)
            order = vals.argsort()[::-1]
            return vals[order], vecs[:, order]

        mat, cov = self._marginalise(params)
        # First find the fiducial value for the parameter in question
        fid_param = None
        pos = [0, 0]
        for p in params:
            if p in dir(self.fid_surv):
                fid_param = getattr(self.fid_surv, p)
            else:
                fid_param = getattr(self.fid_cosmo, p)

            pos[params.index(p)] = fid_param

        if ax is None:
            ax = plt.gca()

        vals, vecs = eigsorted(cov)
        theta = degrees(arctan2(*vecs[:, 0][::-1]))

        # Width and height are "full" widths, not radius
        width, height = 2 * nstd * sqrt(vals)
        ellip = Ellipse(xy=pos, width=width,
                        height=height, angle=theta, **kwargs)

        ax.add_artist(ellip)
        sz = max(width, height)
        s1 = 1.5*nstd*self.sigma_marg(p1)
        s2 = 1.5*nstd*self.sigma_marg(p2)
        ax.set_xlim(pos[0] - s1, pos[0] + s1)
        ax.set_ylim(pos[1] - s2, pos[1] + s2)
        #ax.set_xlim(pos[0] - sz, pos[0] + sz)
        #ax.set_ylim(pos[1] - sz, pos[1] + sz)
        plt.draw()
        return ellip

    @property
    def FoM_DETF(self):
        """
            Computes the figure of merit from the Dark Energy Task Force
            Albrecht et al 2006
            FoM = 1/sqrt(det(F^-1_{w0,wa}))
        """
        det = (self.invFij('w0', 'w0') * self.invFij('wa', 'wa') -
               self.invFij('wa', 'w0') * self.invFij('w0', 'wa'))
        return 1.0 / sqrt(det)

    @property
    def FoM(self):
        """
            Total figure of merit : ln (1/det(F^{-1}))
        """
        return log(1.0 / abs(linalg.det(self.invmat)))

    @property
    def invmat(self):
        """
        Returns the inverse fisher matrix
        """
        if self._invmat is None:
            self._invmat = linalg.inv(self.mat)
        return self._invmat

    @property
    def mat(self):
        """
        Returns the fisher matrix marginalised over nuisance parameters
        """
        # If the matrix is not already computed, compute it
        if self._mat is None:
            self._fullMat = self._computeFullMatrix()
            self._fullInvMat = linalg.pinv(self._fullMat)

            # Apply marginalisation over nuisance parameters
            self._invmat = self._fullInvMat[0:len(self.params),
                                     0:len(self.params)]

            self._mat = linalg.pinv(self._invmat)
        return self._mat

    def _computeDerivatives(self):
        """ Computes all the derivatives of the specified observable with
        respect to the parameters and nuisance parameters in the analysis"""
        # List the derivatives with respect to all the parameters
        dcldp = []

        old_fid_param = None
        old_param = None
        # Computes all the derivatives with respect to the main parameters
        for p in self.params:
            print("varying :" + p)
            # First find the fiducial value for the parameter in question
            fid_param = None
            if p in dir(self.fid_surv):
                fid_param = getattr(self.fid_surv, p)
            else:
                fid_param = getattr(self.fid_cosmo, p)

            step = fid_param * self.step
            if fid_param == 0:
                step = self.step
            # Compute derivative using the 2 point formula
            if old_param is None:
                kw = {p: fid_param + step, 'makeFlat': True}
            else:
                kw = {p: fid_param + step, old_param: old_fid_param,
                      'makeFlat': True}

            self.spectra.update(**kw)

            clp = self._computeObservables()

            kw = {p: fid_param - step, 'makeFlat': True}
            self.spectra.update(**kw)

            clm = self._computeObservables()
            dcl = []
            for (pl, mi) in zip(clp, clm):
                dcl.append((pl - mi) / (2.0 * step))

            dcldp.append(dcl)
            old_fid_param = fid_param
            old_param = p

        # Reset everything to the fiducial value
        kw = {old_param: old_fid_param, 'makeFlat': True}
        self.spectra.update(**kw)

        # Now, computing derivatives with respect to the nuisance parameters
        for p in self.margin_params:
            print("Varying nuisance parameter :" + p)
            nuisance = self.fid_surv.nuisances[p]
            np = nuisance.Np

            # Start varying each nuisance parameter
            for ind_param in range(np):
                fid_param = nuisance.get_value(ind_param)

                step = fid_param * self.step
                if fid_param == 0:
                    step = self.step

                kw = {p + '_p' + str(ind_param): fid_param + step}
                self.spectra.update(**kw)

                clp = self._computeObservables()

                kw = {p + '_p' + str(ind_param): fid_param - step}
                self.spectra.update(**kw)

                clm = self._computeObservables()
                dcl = []
                for (pl, mi) in zip(clp, clm):
                    dcl.append((pl - mi) / (2.0 * step))

                dcldp.append(dcl)

                # Resetting nuisance parameter
                kw = {p + '_p' + str(ind_param): fid_param}
                self.spectra.update(**kw)

        return dcldp


class fisherTomo(fisher):
    """ Tomographic Fisher matrix """

    def __init__(self, fid_cosmo, fid_survey, params, probes, margin_params=[],
                 cutNonLinearScales=None, lmax=10000, nl=200, lmin=2, diagonal=False):
        """ Initializes a tomographic Fisher Matrix using the probes given
        in probes
        """
        self.probes = probes
        # Calls super class initialization
        super(fisherTomo, self).__init__(fid_cosmo, fid_survey, params,
                                         margin_params)

        self.diagonal = diagonal
        # Create a list of redshift bins
        if diagonal:
            self.bins = []
            for i in range(self.fid_surv.nzbins):
                self.bins.append((i,i))
        else:
            self.bins = list(itertools.combinations_with_replacement(
                    range(self.fid_surv.nzbins), r=2))

        # Create a list of spectra from the probes
        self.crossprobes = list(
            itertools.combinations_with_replacement(probes, r=2))

        # Create list of power spectra from probes and redshift bins
        self.cls = list(itertools.product(self.crossprobes, self.bins))

        # Compute l range for each redshift bin
        self.lmax = zeros(len(self.cls))

        for i in range(len(self.cls)):
            zmin = min(self.fid_surv.zbins[self.cls[i][1][0]].zmed,
                       self.fid_surv.zbins[self.cls[i][1][1]].zmed)
            if cutNonLinearScales is None:
                self.lmax[i] = lmax
            else:
                if cutNonLinearScales is 'optimistic':
                    kmax = 0.25
                else:
                    kmax = min(self.fid_surv.zbins[self.cls[i][1][0]].kmax_lin,
                               self.fid_surv.zbins[self.cls[i][1][1]].kmax_lin)
                self.lmax[i] = kmax * self.fid_cosmo.a2chi(z2a(zmin))

        lmax = self.lmax.max()

        # Generating an hybrid array
        self._l = logspace(0, log10(lmax), nl)
        for i in range(nl):
            if self._l[i] <= (i + lmin):
                self._l[i] = i + lmin

        self._nl = len(self._l)

    def _computeObservables(self, shotNoise=False):
        obs = []
        for cl in self.cls:
            obs.append(getattr(self.spectra, 'cl_' + cl[0][0] + cl[0][1])
                       (cl[1][0], cl[1][1], self._l, shotNoise=shotNoise))
        return obs

    def _computeFullMatrix(self):
        """
        Returns the full fisher matrix.
        """
        def find_index(a, b, i, j):
            if ((a, b),  (i, j)) in self.cls:
                return self.cls.index(((a, b), (i, j)))
            else:
                return self.cls.index(((b, a), (j, i)))

        # Prefactor
        geom = 1.0 / (2.0 * self._l + 1.0) / self.fid_surv.fsky

        print("Computing derivatives")
        self._dcldp = self._computeDerivatives()

        print("Computing covariance matrix")
        cl = self._computeObservables(shotNoise=True)

        # Precompute all the covariance matrices
        self._cov = zeros((len(self.cls), len(self.cls), self._nl))
        for ind1 in range(len(self.cls)):
            
            for ind2 in range(ind1 + 1):
                if self.diagonal and ind1 != ind2:
                    continue
                cl1 = self.cls[ind1]
                cl2 = self.cls[ind2]
                C02 = find_index(cl1[0][0], cl2[0][0],
                                 cl1[1][0], cl2[1][0])
                C13 = find_index(cl1[0][1], cl2[0][1],
                                 cl1[1][1], cl2[1][1])
                C03 = find_index(cl1[0][0], cl2[0][1],
                                 cl1[1][0], cl2[1][1])
                C12 = find_index(cl1[0][1], cl2[0][0],
                                 cl1[1][1], cl2[1][0])

                self._cov[ind1, ind2, :] = geom * (cl[C02] * cl[C13] +
                                                   cl[C03] * cl[C12])
                self._cov[ind2, ind1, :] = self._cov[ind1, ind2, :]
        # Precomputes the inverse of the covariance matrix for each l
        self._invcov = zeros_like(self._cov)
        for indl in range(len(self._l)):
            self._invcov[:, :, indl] = linalg.pinv(self._cov[:, :, indl])

        # Computes the fisher Matrix
        # The size of the Fisher matrix is given by the number of derivatives
        print("Computing full Fisher matrix")
        nparams = len(self._dcldp)
        mat = zeros((nparams, nparams))
        for i in range(nparams):
            for j in range(i + 1):
                res = zeros_like(self._l, dtype=double)
                for ind1 in range(len(self.cls)):
                    for ind2 in range(len(self.cls)):
                        temp = (self._dcldp[i][ind1] * self._dcldp[j][ind2] *
                                self._invcov[ind1, ind2, :])
                        lmax = min(self.lmax[ind1], self.lmax[ind2])
                        temp[self._l > lmax] = 0
                        res += temp
                mat[i, j] = simps(res, x=self._l.astype('float'))
                mat[j, i] = mat[i, j]
        return mat


class fisher3d(fisher):
    """ Full 3D fisher analysis using the Spherical Fourier-Bessel expansion
    """

    def __init__(self, fid_cosmo, fid_survey, params, margin_params=[],
                 cutNonLinearScales=True):
        super(fisher3d, self).__init__(fid_cosmo, fid_survey, params,
                                       margin_params)

        self.cutNonLinearScales = cutNonLinearScales

        # TODO : Adapt the l range
        self._l = arange(1, 100)
        self._l[29:] = around(logspace(log(30) / log(10),
                               log(900) / log(10),
                               70))
        self._l = self._l.astype('int')

        # Compute the mask of admissible values i.e. corresponding to the linear part of the power spectrum
        self._linearcut = ones_like(self._l) * 0.25
        if cutNonLinearScales:
            for i in range(len(self._l)):
                cut = lambda k: k * self.fid_cosmo.a2chi(z2a(k/0.132)) - self._l[i]
                self._linearcut[i] = brentq(cut, 0., 0.25)

    def _computeObservables(self, shotNoise=False):
        cl = self.spectra.cl_sfb(self._l, kmax=self._linearcut,
                                  shotNoise=shotNoise,
                                  fid_cosmo=self.fid_cosmo)
        obs = []
        for i in range(len(self._l)):
            obs.append(cl[i][1])
        return obs

    def _computeFullMatrix(self):

        # Prefactor
        geom = 0.5 * self.fid_surv.fsky * (2.0 * self._l + 1.0)

        print("Computing covariance matrix")
        # Computing eigen decomposition of the covariance matrix
        self._cl = self._computeObservables(shotNoise=True)
        # Invert the covariance matrix for every l
        self._invcov = []
        for i in range(len(self._l)):
            self._invcov.append(linalg.eigh(self._cl[i]))

        print("Computing derivatives")
        self._dcldp = self._computeDerivatives()
        nparams = len(self._dcldp)
        # Apply inverse covariance matrix
        for i in range(len(self._l)):
            v, q = self._invcov[i]
            for p in range(nparams):
                self._dcldp[p][i] = dot(dot(q.T, dot(self._dcldp[p][i], q)),
                                        diag(1.0/v))

        print("Computing Full Fisher matrix")
        # Computes the fisher Matrix
        self._integrand = zeros((nparams, nparams, len(self._l)))

        for i in range(nparams):
            for j in range(i+1):
                for n in range(len(self._l)):
                    self._integrand[i, j, n] = trace(dot(self._dcldp[i][n],
                                                         self._dcldp[j][n]))
                    self._integrand[j, i, n] = self._integrand[i, j, n]

        self._integrand *= geom
        mat = simps(self._integrand, x=self._l.astype('float'))
#        integrand = interp1d(self._l, self._integrand, kind='cubic')
#        l = arange(self._l.min(), self._l.max()+1)
#        self._mat = 0.5 * self.fid_surv.fsky * sum(integrand(l) *
#                                                   (2.0*l + 1.0),
#                                                   axis=2)
        return mat
