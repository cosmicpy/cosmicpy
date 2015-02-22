# Copyright (c) 2014-2015, CosmicPy Developers
# Licensed under CeCILL 2.1 - see LICENSE.rst
"""
:mod:`cosmicpy.spectra` -- Clustering and lensing power spectra
==============================================================

.. module:: cosmicpy.spectra
    :synopsis: Computes galaxy spectra
.. moduleauthor:: Francois Lanusse <francois.lanusse@cea.fr>
.. moduleauthor:: Anais Rassat <anais.rassat@epfl.ch>
.. Created on Jun 14, 2013 by Francois Lanusse

.. autosummary::

    spectra

"""

from scipy.integrate import *
from scipy.interpolate import interp1d
from numpy import *
from .tools import BesselWindow
from .utils import *
from . import constants as const
from scipy.special import jn_zeros


class spectra(object):
    """ Computes and stores all the window functions """

    def __init__(self, cosmo, surv, lmax=1000, **kwargs):

        self.cosmo = cosmo
        self.surv = surv
        self._n_a = 256 + 1
        self._lmax = lmax
        # Number of points in the dsbt
        self._chimax = self.surv.chimax(cosmo)

        self._kmax_bessel = 0.5

        bess_zeros = jn_zeros(2, 10000)
        self._nmax = where(bess_zeros > self._kmax_bessel * self._chimax)[0][0]

        self.besselWindow = BesselWindow(self._nmax,
                                         self._lmax,
                                         self._kmax_bessel,
                                         "qlnTable.dat")
        self.update(**kwargs)

    def update(self, **kwargs):

        if len(kwargs) > 0:
            self.cosmo.update(**kwargs)
            self.surv.update(**kwargs)

        # initializes the tabulation arrays
        (self._a, self._da) = linspace(z2a(self.surv.zmax),
                                       z2a(self.surv.zmin),
                                       self._n_a,
                                       retstep=True,
                                       endpoint=False)
        self._z = a2z(self._a)
        self._chi = self.cosmo.a2chi(self._a)

    def g(self, i, l, z):
        """ galaxy clustering window function """
        # Checks that the bin index is correct
        if (i < 0 or i >= self.surv.nzbins):
            print ("Error, bin number " + i + " is out of range.")
            return 0.0
        z = atleast_1d(z)
        l = atleast_1d(l)
        b = zeros((len(l), len(z)))
        r = self.cosmo.a2chi(z2a(z))

        for j in range(len(z)):
            k = (l + 0.5) / r[j]
            b[:, j] = self.surv.zbins[i].nz(z[j])*self.surv.bias(z[j], k)
        return b

    def cl_gg(self, i, j, l, shotNoise=False, linear=True, **kwargs):
        """ galaxy-galaxy angular power spectrum in the Limber approximation"""
        l = atleast_1d(l)
        (a, da) = linspace(z2a(max([self.surv.zbins[i].zmax,
                                    self.surv.zbins[j].zmax])),
                           z2a(min([self.surv.zbins[i].zmin,
                                    self.surv.zbins[j].zmin])),
                           self._n_a, retstep=True, endpoint=False)

        clintegrand = 1.0 / (self.cosmo.dchioverda(a) *
                             (self.cosmo.f_k(a)**2) * a**4)
        if linear:
            pl = self.cosmo.pl_lin(l, a, **kwargs)
        else:
            pl = self.cosmo.pl(l, a, **kwargs)
        res = romb(clintegrand * (self.g(i, l, a2z(a)) * self.g(j, l, a2z(a)) *
                                  pl), da)

        if shotNoise and i == j:
            res += 1.0 / self.surv.zbins[i].ngal

        return res

    def W(self, l, k1=None, k2=None, evol=True, fid_cosmo=None, kmax=0.25):
        r"""
        Computes the Spherical Fourier-Bessel window function for the survey

        Parameters
        ----------
        l : int
            Order of the Bessel functions.
        k1 : array_like, optional
            Value of the scale at which to evaluate the window function
            in h/Mpc. If None is provided, the window function is evaluated
            at discrete points defined by the zeroes of the Bessel
            functions up to the kmax parameter (def : None).
        k2 : array_like, optional
            Value of the scale at which to evaluate the window function
            in h/Mpc. If None is provided, a symmetric Bessel window is
            computed :math:`W_l(k1, k1)`  (def : None).
        evol : boolean, optional
            Flag to include the time dependent linear bias and growth in
            the computation of the window function (def : True).
        fid_cosmo : Cosmology, optional
            Fiducial cosmology to use for the computation of redshift to
            comoving distance in the case of a real survey. If None is
            provided, the true cosmology is used (def : None).
        kmax : float, optional
            Maximum scale at which to compute the window if k1 is not
            provided. This is usefull to avoid computing the
            window at non linear scales which would be cut afterwards
            (def : 0.25).

        Notes
        -----

        :math:: W_l (k1,k2) = \int k1 \phi(r)
                j_l(k1 r) j_l(k2 r) r^2 dr
        """

        # Computing the selection function for the requested modes l
        chi_l = self.besselWindow.rgrid(l)
        phi = self.surv.phi(self.cosmo, chi_l)

        # Include growth and bias if requested
        if evol:
            # Include growth
            a = self.cosmo.chi2a(chi_l)
            g = self.cosmo.G(a)
            phi *= g
            # Include scale dependent bias, either using the provided k or
            # using the default k_ln grid
            if k2 is None:
                k = self.besselWindow.kgrid(l, self._chimax)
                b = self.surv.bias(a2z(a), k)
            else:
                b = self.surv.bias(a2z(a), k2)
            phi = einsum('i,ij->ij', phi, b)
        else:
            if k2 is None:
                k = self.besselWindow.kgrid(l, self._chimax)
            else:
                k = k2
            phi = einsum('i,ij->ij', phi, ones((len(phi), len(k))))

        # Computing the comoving distance for a redshift space survey
        # according to a fiducial cosmology
        if fid_cosmo is None:
            s = chi_l
        else:
            s = fid_cosmo.a2chi(self.cosmo.chi2a(chi_l))

        # First case, if k1 is not specified, the window is computed for W_l(k_{l n}, k2)
        # with an adapted step in k2 to correctly sample the main lobe of the
        # window
        if k1 is None:
            return self.besselWindow.tabulated(l, self._chimax,
                                               kmax, s, phi)

        k1 = atleast_1d(k1)

        # Second case, k1 is provided, not k2. In this case,
        # an optimal window is computed
        if k2 is None:
            return self.besselWindow.optimal(l, self._chimax, k1, s, phi)

        k2 = atleast_1d(k2)

        # 3rd case, k1 and k2 are provided
        return self.besselWindow.custom(l, k1, k2, s, phi)

    def cl_sfb(self, l, k=None, shotNoise=False, onlyNoise=False,
               evol=True, fid_cosmo=None, kmax=0.25, **kwargs):

        l = atleast_1d(l)
        kmax = atleast_1d(kmax)
        if kmax.size != l.size:
            kmax = ones_like(l) * kmax[0]

        if k is None:
            # Computing linear scale cut
            res = []
            for i in range(len(l)):
                k1, k2, w = self.W(l[i], evol=evol, kmax=kmax[i],
                                   fid_cosmo=fid_cosmo)
                pk = self.cosmo.pk_lin(k2, **kwargs)
                mat = zeros((len(k1), len(k1)))
                integrand = w*k2**2*pk
                for j in range(len(k1)):
                    integ = (2.0 / pi)**2 * simps(integrand[j, :] * w[0:j+1, :], x=k2)
                    mat[j, 0:j + 1] = integ
                    mat[0:j + 1, j] = integ
                if shotNoise:
                    mat += sqrt(2) / pi * k1*self.W(l[i], k1, k1, evol=False) / self.surv.N

                res.append([k1, mat])

            if len(res) == 1:
                return res[0]
            else:
                return res
        else:
            k = atleast_1d(k)
            res = zeros((l.size, k.size))
            noise = zeros((l.size, k.size))
            if not onlyNoise:
                for i in range(len(l)):
                    kw, w = self.W(l[i], k, evol=evol)
                    pk = self.cosmo.pk_lin(kw[:, :], **kwargs)
                    integrand = square(kw * w) * pk
                    res[i] = (2.0 / pi)**2 * trapz(abs(integrand), x=kw)

            if shotNoise or onlyNoise:
                for i in range(l.size):
                    w = k * diag(self.W(l[i], k, k, evol=False)) / self.surv.N
                    noise[i, :] = sqrt(2) / pi * w

            if onlyNoise:
                return noise.squeeze()
            elif shotNoise:
                return (res + noise).squeeze()
            else:
                return res.squeeze()