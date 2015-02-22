# Copyright (c) 2014-2015, CosmicPy Developers
# Licensed under CeCILL 2.1 - see LICENSE.rst
"""
.. module:: theory
    :synopsis: Performs tests based on exact theoretical solutions
.. moduleauthor:: Francois Lanusse <francois.lanusse@cea.fr>
.. moduleauthor:: Anais Rassat <anais.rassat@epfl.ch>

.. Created on Jan 20, 2014 by Francois Lanusse
"""
from scipy.special import iv
import numpy as np
from math import *
import cosmicpy
from scipy.integrate import trapz


def gaussian_bessel_window(l, r0, k1, k2):
    return 0.25 * (np.pi * r0**2 * np.sqrt(k1 / k2) *
                   np.exp(-r0**2 * 0.25 * (k1**2 + k2**2)) *
                   iv(l+0.5, r0**2*k1*k2 / 2))


def check_Bessel_window(l, k=None, r0=100):
    r""" Checks the computation of the bessel window compared to the exact
    solution
    for Gaussian selection functions.

    Parameters
    ----------
    l : int
        Order of the multipole.
    r0 : real, optional
        Radius parameter for the survey Gaussian selection function in Mpc/h

    Returns
    -------
    (k, W, Wexact) : tuple of ndarrays
        Window compute the discrete spherical bessel transform and the exact
        window computed from theory
    """
    # Initialise survey and spectra
    cosmo = cosmicpy.cosmology()
    surv = cosmicpy.survey(nzparams={'type': 'gaussian',
                                    'cosmology': cosmo,
                                    'r0': r0},
                          zmax=0.5,
                          fsky=1.0)
    spectra = cosmicpy.spectra(cosmo, surv)
    survey_volume = pi**(3.0 / 2.0) * r0**3
    if k is None:
        k = np.logspace(log(0.001) / log(10), log(0.20) / log(10), 512,
                        endpoint=False)
    W = spectra.W(l, k, k, evol=False, kmax=0.5)

    W = W.squeeze()
    Wexact = np.zeros_like(W)
    for k1 in range(len(k)):
            Wexact[k1, :] = (gaussian_bessel_window(l, r0, k[k1], k) /
                             survey_volume)
    print("Maximum absolute error vs maximum value: " +
          str(abs(W - Wexact).max()) +
          " " +
          str(Wexact.max()))
    return (k, W, Wexact)


def check_Bessel_spectrum(l, k=None, r0=100):
      # Initialise survey and spectra
    cosmo = cosmicpy.cosmology()
    surv = cosmicpy.survey(nzparams={'type': 'gaussian',
                                    'cosmology': cosmo,
                                    'r0': r0},
                          zmax=5,
                          fsky=1.0)
    spectra = cosmicpy.spectra(cosmo, surv)
    survey_volume = pi**(3.0 / 2.0) * r0**3
    if k is None:
        k = np.logspace(log(0.0001) / log(10), log(0.25) / log(10), 512,
                        endpoint=False)

    kw, w = spectra.W(l, k, evol=False)
    kw = kw.squeeze()
    w = w.squeeze()
    w_exact = np.zeros_like(w)
    for k1 in range(len(k)):
        w_exact[k1, :] = gaussian_bessel_window(l,r0, k[k1], kw[k1,:])/ survey_volume

    pk = cosmo.pk_lin(kw)
    integrand = np.square(kw * w_exact) * pk
    res = (2.0 / pi)**2 * trapz(integrand, x=kw)

    cl = spectra.cl_sfb(l, k)

    return k, cl, res
