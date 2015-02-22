.. Copyright (c) 2014-2015, CosmicPy Developers
.. Licensed under CeCILL 2.1 - see LICENSE.rst
.. currentmodule:: cosmicpy.cosmology

Getting Started with CosmicPy
============================

First, start an :program:`ipython` sesion in pylab mode::

    $ ipython --pylab
    
This will automatically load :program:`numpy` and :program:`matplolib` and enable the
interactive plotting environment.

Once :program:`ipython` is running, you can import the cosmicpy package:

.. ipython::

    @suppress
    In [0]: import matplotlib; matplotlib.use('Agg'); from pylab import *;

    In [1]: from cosmicpy import *


Create and modify a cosmology
-----------------------------

The next step is to create a :py:class:`cosmology` object. Let us start with a default cosmology:

.. ipython::

    In [2]: cosmo = cosmology()

To get a summary of the cosmological parameters defined in this object, you can use the :py:func:`print` fucntion:

.. ipython::

    In [3]: print cosmo
    FLRW Cosmology with the following parameters:
        h:        0.7
        Omega_b:  0.045
        Omega_m:  0.25
        Omega_de: 0.75
        w0:       -0.95
        wa:       0.0
        n:        1.0
        tau:      0.09
        sigma8:   0.8

All the cosmological parameters are accessible as attributes of the :py:class:`cosmology` object. For instance, to check the value of the matter density today :math:`\Omega_m` run:

.. ipython::

    In [3]: cosmo.Omega_m
    Out[3]: 0.25

The values for the cosmological parameters can be modified which will automatically trigger an update of all the derived quantities.

To illustrate this, let's see the value of the sound horizon at recombination :math:`s_h` with the default cosmology, update :math:`\Omega_m` and check that the value of :math:`s_h` has changed:

.. ipython::

    In [4]: cosmo.sh_r
    Out[4]: 105.27381988350771

    In [5]: cosmo.Omega_m = 0.3

    In [6]: cosmo.sh_r
    Out[6]: 100.59308203228834

The cosmological parameters can also be set at the creation of the :py:class:`cosmology` object by specifying their values to the constructor.
For instance, to create a cosmology with :math:`\Omega_m = 0.2` and :math:`\Omega_{de} = 0.8`:

.. ipython::

    In [7]: cosmo2 = cosmology(Omega_m=0.2, Omega_de=0.8)

    In [8]: print cosmo2
    FLRW Cosmology with the following parameters:
        h:        0.7
        Omega_b:  0.045
        Omega_m:  0.2
        Omega_de: 0.8
        w0:       -0.95
        wa:       0.0
        n:        1.0
        tau:      0.09
        sigma8:   0.8


Compute distances
-----------------

Now that we have seen how to create a :py:class:`cosmology` object with the desired cosmological parameters, let's use it to compute some distances.

:py:class:`cosmology` offers a number of different methods to compute various cosmological distances:

.. autosummary::

    cosmology.a2chi
    cosmology.f_k
    cosmology.d_A

These methods all compute a distance as a function of the scale factor :math:`a`.

It is often usefull to work with redshifts instead of
scale factors. The conversion between the two quantities is performed by the 2 utility functions:

.. autosummary::

    cosmicpy.utils.z2a
    cosmicpy.utils.a2z

Here is an example to compute the radial comoving distance for an array of redshifts between 0 and 10:

.. ipython::

    In [7]: z = arange(0,10,0.01)

    In [8]: figure();

    # The z2a function is used to convert z to scale factors
    In [9]: plot(z,cosmo.a2chi(z2a(z))); grid();

    @savefig start-1.png
    In [10]: title('Radial comoving distance as a function of redshift');


Compute the matter power spectrum
---------------------------------

:py:class:`cosmology` implements several methods to compute the linear and non linear matter power spectrum:

.. autosummary::
    cosmology.pk_lin
    cosmology.pk

By default, these two functions will compute the power spectrum at :math:`a=1`
and using the Eisenstein & Hu transfer function with oscillations. The non-linear
correction implemented in :py:meth:`cosmology.pk` is computed based on :cite:`2003:Smith`.

Here is an example to compute these two power spectra:

.. ipython::

    # Creating a vector of wavenumbers
    In [11]: k = logspace(-3,0,100)

    @suppress
    In [11]: figure()

    In [12]: loglog(k,cosmo.pk(k), label='Smith et al. $P(k)$');

    In [13]: loglog(k,cosmo.pk_lin(k), label='Linear $P(k)$');

    @savefig start-2.png
    In [14]: legend(); xlabel('k [h/Mpc]');


Both power spectra can be computed with or without baryonic oscillations by setting the **type** keyword to
either:

* 'eisenhu_osc' : Include BAOs (default behavior)
* 'eisenhu' : Without BAOs

This example shows the matter power spectrum, with non-linear corrections, with and without BAO wiggles:

.. ipython::

    @suppress
    In [15]: figure()

    In [16]: loglog(k, cosmo.pk(k, type='eisenhu'), label='eisenhu');

    In [17]: loglog(k, cosmo.pk(k, type='eisenhu_osc'), label='eisenhu_osc');

    @savefig start-3.png
    In [18]: legend(); xlabel('k [h/Mpc]'); title('Matter power spectrum');


Finally, the power spectrum can be computed at different scale factors:

.. ipython::

    @suppress
    In [19]: figure()

    # Create an array of scale factors corresponding to z=0, 1 ,2 and 5
    In [20]: a = z2a(array([0, 1, 2, 5]))

    In [21]: loglog(k, cosmo.pk(k,a));

    @savefig start-4.png
    In [22]: xlabel('k [h/Mpc]'); title('Matter power spectrum at z=0, 1, 2, 5');


References
----------

.. bibliography:: biblio.bib
    :filter: docname in docnames
