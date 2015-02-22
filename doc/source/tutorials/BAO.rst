.. Copyright (c) 2014-2015, CosmicPy Developers
.. Licensed under CeCILL 2.1 - see LICENSE.rst
.. currentmodule:: cosmicpy.cosmology

3D Spherical Analysis of BAOs
=============================

In this tutorial, we recover the results of :cite:`2012:Rassat` using cosmicpy.


Setting up the cosmology and the surveys
----------------------------------------

For this tutorial, we are going to use galaxy surveys following a Gaussian selection function of the form:

.. math::    
    
    \Phi(r) = e^{-(r/r_0)^2}

We define 2 such surveys, one shallow with :math:`r_0 = 100 \quad h^{-1}` Mpc and one deep with :math:`r_0 = 1400 \quad h^{-1}` Mpc:

.. ipython::
    :okwarning:

    @suppress
    In [0]: from pylab import *;

    In [1]: from cosmicpy import *

    In [2]: cosmo = cosmology()
          
    In [3]: surv_deep = survey(nzparams={'type':'gaussian', 'cosmology':cosmo, 'r0':1400}, zmax=3.0, fsky=1.0)

    In [4]: surv_shallow = survey(nzparams={'type':'gaussian', 'cosmology':cosmo, 'r0':100}, zmax=0.5, fsky=1.0, chicut=0)
      
Now, we plot the selection functions along with the galaxy density probability distribution for these 2 surveys normalized by their maximum values:

.. ipython::
    :okwarning:
    
    # Create an array of sampling redshifts between z=0 and z=2
    In [5]: z = arange(0.0001,2.0,0.001)
    
    # Convert redshift array into comoving distance
    In [8]: r = cosmo.a2chi(z2a(z)) 
    
    # Compute n(z) for both survey
    In [6]: pzd = surv_deep.nz(z)
    
    In [7]: pzs = surv_shallow.nz(z)
    
    # Compute selection function for both survey
    In [9]: phid = surv_deep.phi(cosmo,r)
        
    In [10]: phis = surv_shallow.phi(cosmo,r)

Now, let us plot the selection functions and galaxy distributions for both surveys:
    
.. ipython::
    :okwarning:

    In [10]: figure();
    
    In [11]: subplot(211); xlabel('Redshift'); grid();
    
    In [12]: plot(z,pzd/pzd.max(),label='deep');
    
    In [13]: plot(z,pzs/pzs.max(),label='shallow');
    
    In [14]: legend(); title('Galaxy density probability distributions');
    
    In [15]: subplot(212); xlabel('Comoving distance [Mpc/h]'); grid();
    
    In [16]: plot(r,phid/phid.max(),label='deep');
    
    In [17]: plot(r,phis/phis.max(),label='shallow');
    
    @suppress
    In [19]: tight_layout()
    
    @savefig BAO-1.png
    In [18]: legend(); title(r'Selection functions');
    
Impact of the depth of the survey on the SFB window
---------------------------------------------------

The depth of the survey directly impacts the SFB window function :math:`W_\ell(k, k')` : the deeper the survey, the closer to a Dirac  :math:`\delta(k - k')` it becomes. 
We illustrate this by computing the window function for :math:`\ell = 3` for the 2 surveys we consider using :py:meth:`spectra.W` :

.. ipython::
    :okwarning:

    In [4]: k = logspace(-3,log(0.25)/log(10),512)
    
    In [5]: sp_deep = spectra(cosmo,surv_deep)
        
    In [6]: sp_shallow = spectra(cosmo,surv_shallow)
                
    In [7]: subplot(121); xscale('log') ; yscale('log')
        
    In [8]: contourf(k, k, sp_shallow.W(3, k, k), 256);
    
    In [9]: title('Shallow survey'); xlabel('k [h/Mpc]'); ylabel('k [h/Mpc]');
        
    In [9]: subplot(122); xscale('log') ; yscale('log')
        
    In [10]: contourf(k, k, sp_deep.W(3, k, k), 256);
    
    @savefig BAO-2.png
    In [11]: title(r'Deep survey'); xlabel('k [h/Mpc]'); ylabel('k [h/Mpc]');


Impact of BAOs on the power spectrum
------------------------------------

To emphasize the BAO features on the power spectrum, we can plot the ratio of linear power spectra today with and without BAOs. The fitting formulae used here are
taken from :cite:`1998:EisensteinHu` . The 2 kinds of power spectra are easily computed using the :py:meth:`cosmology.pk_lin` method:

.. ipython::
    :okwarning:

    In [12]: pk_osc = cosmo.pk_lin(k,1.0,type='eisenhu_osc')
    
    In [13]: pk = cosmo.pk_lin(k,1.0,type='eisenhu')
    
    In [14]: figure(); grid();
    
    In [14]: semilogx(k,pk_osc/pk);
    
    @savefig BAO-3.png
    In [15]:  xlabel('k [h/Mpc]'); title('Ratio of power spectra with and without BAO');
    

Impact of BAOs on the SFB power spectrum
----------------------------------------

In Spherical Fourier-Bessel space, we consider the ratio :math:`R_\ell^C (k)` giben by:

.. math::    
        R_\ell^C(k) = \frac{C_\ell^\mathrm{b}(k)}{C_\ell^{\mathrm{nob}}(k)}

where :math:`C_\ell^\mathrm{b}(k)` is the diagonal SFB power spectrum including the physical effects of baryons and :math:`C_\ell^\mathrm{nob}(k)` is the 
smooth part of the SFB power spectrum. Here we use the :py:meth:`spectra.cl_sfb` method:

.. ipython::
    :okwarning:

    In [11]: k = logspace(-3,log(0.25)/log(10),100)

    In [15]: l = arange(1,65)

    In [16]: R_deep = sp_deep.cl_sfb(l,k,shotNoise=False,type='eisenhu_osc')/sp_deep.cl_sfb(l,k,shotNoise=False,type='eisenhu')

    In [17]: R_shallow = sp_shallow.cl_sfb(l,k,shotNoise=False,type='eisenhu_osc')/sp_shallow.cl_sfb(l,k,shotNoise=False,type='eisenhu')

    In [18]: figure(1); subplot(121); xscale('log') ; yscale('log'); xlim([1,65]);
    
    In [19]: contourf(l,k,R_shallow.T,256);
    
    In [20]: xlabel('l'); ylabel('k [h/Mpc]'); title(r'Wide and Shallow survey');
    
    In [21]: subplot(122); xscale('log') ; yscale('log'); xlim([1,65]);
    
    In [22]: contourf(l,k,R_deep.T,256);

    @savefig BAO-4.png
    In [23]: xlabel('l'); ylabel('k [h/Mpc]'); title(r'Wide and Deep survey');
   

    