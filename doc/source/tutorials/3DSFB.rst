.. Copyright (c) 2014-2015, CosmicPy Developers
.. Licensed under CeCILL 2.1 - see LICENSE.rst
3D SFB galaxy clustering
========================

In this tutorial we use cosmicpy to recover some of the results from :cite:`2013:Yoo`.

Setting up the cosmology and the surveys
----------------------------------------

For this tutorial, we are going to use galaxy surveys following a Gaussian selection function of the form:
    
.. math::    
        \Phi(r) = e^{-(r/r_0)^2}

Here, we are using :math:`r_0 = 2354 \quad h^{-1}` Mpc. The resulting volume of the survey is computed 
for a Gaussian selection function using :math:`V = \pi^{3/2} r_0^3`.

The average galaxy density is set to :math:`\bar{n} = 10^{-4} (h^{-1} \mathrm{Mpc})^{-3}` which corresponds
to a total number of galaxies in survey :math:`N` computed as:

.. math::
        N = \bar{n} V = 7,263,470

We create a survey following these parameters in cosmicpy using: 

.. ipython::
    :okwarning:

    @suppress
    In [0]: from pylab import *;
    
    In [1]: from cosmicpy import *

    In [2]: cosmo = cosmology()

    In [3]: r0 = 2354.0

    In [4]: survey_volume = pi**(3.0 / 2.0) * r0**3

    In [5]: N = survey_volume * 10**(-4)

    In [6]: ngal = N /(4.0* pi)/((180.0/pi)**2)/(60.0*60.0) # Conversion to gal/ arcmin^2

    In [7]: surv = survey(nzparams={'type':'gaussian','cosmology':cosmo, 'r0':2354.0}, zmax=7.0, fsky=1.0, ngal=ngal, biastype='constant')
    
    In [8]: surv.N     # Checking the number of galaxies in the survey                                                                 
    Out[8]: 7263470.626200338  
    
    
    In [11]: surv.ngal # Checking the mean galaxy density per square arcmin
    Out[11]: 0.048908749054432946

Below, we plot the selection function along with the normalized galaxy density distribution:

.. ipython::
    :okwarning:
      
    In [8]: z = arange(0.001,7.0,0.01)
    
    In [9]: pz = surv.nz(z)

    In [10]: r = cosmo.a2chi(z2a(z))     # Convert redshift array into comoving distance
    
    In [11]: phi = surv.phi(cosmo,r)     # Compute selection function for both survey
    
    In [12]: figure(); subplot(211);
    
    In [13]: xlabel('Redshift'); grid(); title('Galaxy density probability distributions');
    
    In [14]: plot(z,pz);
    
    In [16]: subplot(212); 
    
    In [17]: xlabel('Comoving distance [Mpc/h]'); grid();  title(r'Selection functions');
    
    @suppress
    In [18]: tight_layout()
    
    @savefig 3DSFB-1.png
    In [19]: plot(r,phi);

Computing 3D power spectra
--------------------------


.. ipython::
    :okwarning:

    In [20]: sp = spectra(cosmo,surv)
    
    In [21]: l = [2, 5, 10]             # Array of angular multipoles
    
    In [22]: k = logspace(-4,-0.65,100) # Array of wavenumbers
    
    In [23]: cl = sp.cl_sfb(l,k,shotNoise=False,evol=False)*(survey_volume**2)/r0*(2*sqrt(2*pi))
    
    In [24]: cl_evol = sp.cl_sfb(l,k,shotNoise=False,evol=True)*(survey_volume**2)/r0*(2*sqrt(2*pi))
    
    In [25]: cln = sp.cl_sfb(l,k,onlyNoise=True)*(survey_volume**2)/r0*(2*sqrt(2*pi))
    
    In [26]: figure(); subplot(321);
    
    In [27]: loglog(k,cl[0,:]); loglog(k,cln[0,:]); loglog(k,cosmo.pk_lin(k));
    
    In [30]: ylim(1e2,6e4); xlim(0.0001,0.23);
    
    In [31]: subplot(323);
    
    In [27]: loglog(k,cl[1,:]);  loglog(k,cln[1,:]); loglog(k,cosmo.pk_lin(k));

    In [30]: ylim(1e2,6e4); xlim(0.0001,0.23);
    
    In [30]: subplot(325);
    
    In [27]: loglog(k,cl[2,:]); loglog(k,cln[2,:]); loglog(k,cosmo.pk_lin(k));
    
    In [30]: ylim(1e2,6e4); xlim(0.0001,0.23);
    
    In [30]: subplot(322);
    
    In [27]: loglog(k,cl_evol[0,:]); loglog(k,cln[0,:]); loglog(k,cosmo.pk_lin(k));
    
    In [30]: ylim(1e2,6e4); xlim(0.0001,0.23);
    
    In [30]: subplot(324);
    
    In [27]: loglog(k,cl_evol[1,:]); loglog(k,cln[1,:]); loglog(k,cosmo.pk_lin(k));
    
    In [30]: ylim(1e2,6e4); xlim(0.0001,0.23);
    
    In [30]: subplot(326);
    
    In [27]: loglog(k,cl_evol[2,:]); loglog(k,cln[2,:]); loglog(k,cosmo.pk_lin(k));

    In [30]: ylim(1e2,6e4); xlim(0.0001,0.23);

    @savefig 3DSFB-2.png
    In [31]: subplots_adjust(hspace=0);


References
----------

.. bibliography:: biblio.bib
    :filter: docname in docnames