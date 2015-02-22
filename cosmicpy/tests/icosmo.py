# Copyright (c) 2014-2015, CosmicPy Developers
# Licensed under CeCILL 2.1 - see LICENSE.rst
import matplotlib.pyplot as plt
from scipy.io import readsav
from numpy import *
from scipy.integrate import romberg
from scipy.optimize import brentq
import cosmicpy


def _print_diff(name, param1, param2):
    r""" Helper function used to display differences between values.
    """
    div = param1
    if div == 0:
        div = 1
    print(name + ': {:15f} | {:15f} | {:15f} | {:.2%}'\
          .format(param1,
                  param2,
                  abs(param1 - param2),
                  abs((param1 - param2) / div)))


def import_cosmology(filename, structure_name="fid"):
    r""" Loads an icosmo cosmology from a fiducial structure stored in an
    idl save file into a cosmicpy cosmology.

    Parameters
    ----------
    filename : str
        Name of the idl save from which to load the cosmology.
    structure_name : str, optional
        Name of the icosmo fiducial structure stored in the save file.

    Returns
    -------
    cosmo : cosmology
        cosmicpy cosmology corresponding to the icosmo input.
    """
    icosmo_file = readsav(filename)
    icosmo = icosmo_file.get(structure_name)

    h = icosmo['cosmo'][0]['h'][0]
    Omega_m = icosmo['cosmo'][0]['omega_m'][0]
    Omega_de = icosmo['cosmo'][0]['omega_l'][0]
    Omega_b = icosmo['cosmo'][0]['omega_b'][0]
    w0 = icosmo['cosmo'][0]['w0'][0]
    wa = icosmo['cosmo'][0]['wa'][0]
    tau = icosmo['cosmo'][0]['tau'][0]
    n = icosmo['cosmo'][0]['n'][0]
    sigma8 = icosmo['cosmo'][0]['sigma8'][0]

    cosmo = cosmicpy.cosmology(h=h, Omega_m=Omega_m, Omega_de=Omega_de,
                              Omega_b=Omega_b, w0=w0, wa=wa, tau=tau,
                              n=n, sigma8=sigma8)
    return cosmo


def import_survey(filename, structure_name="fid", expt_name="sv1"):
    r""" Loads an icosmo survey from a fiducial structure stored in an
    idl save file into a cosmicpy survey.

    Parameters
    ----------
    filename : str
        Name of the idl save from which to load the cosmology.
    structure_name : str, optional
        Name of the icosmo fiducial structure stored in the save file.

    Returns
    -------
    surv : survey
        cosmicpy survey corresponding to the icosmo input.
    """
    icosmo_file = readsav(filename)
    icosmo = icosmo_file.get(structure_name)['expt'][0]

    nzbins = icosmo[expt_name][0]['n_zbin'][0]
    zerror = icosmo[expt_name][0]['zerror'][0]
    ng = icosmo[expt_name][0]['ng'][0]
    a_survey = icosmo[expt_name][0]['a_survey'][0]
    dndztype = icosmo[expt_name][0]['dndztype'][0]
    if dndztype != b"smail":
        print("unsupported galaxy distribution")
        return None
    a = icosmo[expt_name][0]['dndzp'][0][0]
    b = icosmo[expt_name][0]['dndzp'][0][1]
    zmed = icosmo[expt_name][0]['z_med'][0]

    biastype = icosmo[expt_name][0]['biastype'][0]
    if biastype == b'bias0':
        btype = 'constant'
    elif biastype == b'bias1':
        btype = 'sqrt'
    else:
        print("unsupported bias type")
        return None

    # Compute fsky
    fsky = a_survey/(180./pi)**2/(4.0*pi)

    # find z0 corresponding to zmedian
    def med_smail(z0):
        smail = lambda z: z**a * exp(-(z/z0)**b)
        smail_norm = romberg(smail, 0, 5)
        smailn = lambda z: (z**a * exp(-(z/z0)**b)) / smail_norm
        f = lambda x: romberg(smailn, 0, x) - 0.5
        return brentq(f, 0.01, 5) - zmed
    z0 = brentq(med_smail, 0.01, 3)

    surv = cosmicpy.survey(nzbins=nzbins,
                          ngal=ng,
                          zphot_sig=zerror,
                          fsky=fsky,
                          biastype=btype,
                          nzparams={'type': 'smail',
                                    'a': a,
                                    'b': b,
                                    'z0': z0})
    return surv

def check_cosmology(filename, structure_name="cosmo"):
    r""" Loads an icosmo cosmology from an idl save file and compares the
    content of the cosmology structure with results from cosmicpy.

    Parameters
    ----------
    filename : str
        Name of the idl save from which to load the cosmology.
    structure_name : str, optional
        Name of the icosmo cosmology structure stored in the save file.
    """
    icosmo_file = readsav(filename)
    icosmo = icosmo_file.get(structure_name)

    # Reading cosmological parameters and creating corresponding
    # cosmicpy object
    h = icosmo['const'][0]['h'][0]
    Omega_m = icosmo['const'][0]['omega_m'][0]
    Omega_de = icosmo['const'][0]['omega_l'][0]
    Omega_b = icosmo['const'][0]['omega_b'][0]
    w0 = icosmo['const'][0]['w0'][0]
    wa = icosmo['const'][0]['wa'][0]
    tau = icosmo['const'][0]['tau'][0]
    n = icosmo['const'][0]['n'][0]
    sigma8 = icosmo['const'][0]['sigma8'][0]

    cosmo = cosmicpy.cosmology(h=h, Omega_m=Omega_m, Omega_de=Omega_de,
                      Omega_b=Omega_b, w0=w0, wa=wa, tau=tau,
                      n=n, sigma8=sigma8)
    print(cosmo)

    # Comparing derived quantities
    Omega_k = icosmo['const'][0]['omega_k'][0]
    r0 = icosmo['const'][0]['r0'][0] * h
    gamma = icosmo['const'][0]['gamma'][0]
    sh = icosmo['const'][0]['sh'][0] * h

    print("\nComparing derived constants :")
    print("parameter:        [icosmo] |       [cosmicpy] |"
          "            diff | relative error in percents ")
    _print_diff("Omega_k  ", Omega_k, cosmo.Omega_k)
    _print_diff("gamma    ", gamma, cosmo.gamma)
    _print_diff("sh       ", sh, cosmo.sh_r)

    # Comparing evolved quantities
    rh = icosmo['const'][0]['rh'][0] * h
    z = icosmo['evol'][0]['z'][0]
    a = icosmo['evol'][0]['a'][0]
    chi = icosmo['evol'][0]['chi'][0] * rh
    hc = icosmo['evol'][0]['hc'][0] / h
    Omega_m = icosmo['evol'][0]['omega_m_a'][0]
    Omega_de_a = icosmo['evol'][0]['omega_l_a'][0]
    dzdr = icosmo['evol'][0]['dzdr'][0] / h
    da = icosmo['evol'][0]['da'][0] * h

    print("\nComparing evolved quantities, maximum relative error:")
    print("Radial Comoving Distance  :  " +
          str(abs((chi[1:] - cosmo.a2chi(a[1:])) / chi[1:]).max()))
    print("Angular Diameter Distance :  " +
          str(abs((da[1:] - a[1:] * cosmo.f_k(a[1:])) / da[1:]).max()))
    print("Hubble constant           :  " +
          str(abs((hc - cosmo.H(a)) / hc).max()))
    print("Omega_m(a)                :  " +
           str(abs((Omega_m - cosmo.Omega_m_a(a)) / Omega_m).max()))
    print("Omega_de(a)               :  " +
           str(abs((Omega_de_a - cosmo.Omega_de_a(a)) / Omega_de_a).max()))
    print("dzdr                      :  " +
           str(abs((dzdr - cosmo.dzoverda(a) /
                    cosmo.dchioverda(a)) / dzdr).max()))

    # Comparing Power spectrum
    k = icosmo['pk'][0]['k'][0]
    zpk = icosmo['pk'][0]['z'][0]
    pk = icosmo['pk'][0]['pk'][0]
    pk_l = icosmo['pk'][0]['pk_l'][0]

    print("\nComparing linear power spectra")
    pk2 = cosmo.pk_lin(k, cosmicpy.z2a(zpk))
    plt.figure()
    plt.contourf(k, zpk, abs(pk2.T - pk_l) / pk_l, 250)
    plt.xlabel(r"Scale in Mpc/h")
    plt.ylabel(r"Redshift")
    plt.xscale('log')
    plt.colorbar()
    plt.title('Relative Error in linear power specrum computation')
    plt.show()
    print("Maximum relative error on linear power specrum : " +
          str((abs(pk2.T - pk_l) / pk_l).max()))


