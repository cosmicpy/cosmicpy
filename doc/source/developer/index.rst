.. Copyright (c) 2014-2015, CosmicPy Developers
.. Licensed under CeCILL 2.1 - see LICENSE.rst
=======================
Developer Documentation
=======================


Building the documentation
==========================

Requirements
------------

To be able to build the documentation, the following softwares are required:
    
* Sphinx
* `sphinxcontrib-napoleon <https://pypi.python.org/pypi/sphinxcontrib-naopleon>`_
* `sphinxcontrib-bibtex <https://pypi.python.org/pypi/sphinxcontrib-bibtex/>`_
* `sphinx-bootstrap-theme <https://pypi.python.org/pypi/sphinx-bootstrap-theme/>`_
    
On a `Mac`, Sphinx and Matplotlib can be installed through Macports::
    
    $ sudo port install py27-sphinx py27-matplotlib
    
To be able to add bibtex references, you need the sphinxcontrib-bibtex software. The numpydoc module allows the use of human readable
documentation markings that are still able to be processed by sphinx. Installing both of these extensions can be easily done using the following command::
    
    $ sudo pip install sphinx-bootstrap-theme sphinxcontrib-napoleon sphinxcontrib-bibtex
    
    
And that's it, your should now be able to build the documentation.


Making the build
----------------

To build the documentation, go to the `doc` folder and execute::
    
    $ make html
    
This should execute sphinx and the documentation should be generated in the doc/build directory. To open the documentation::
    
    $ open buid/html/index.html
    
    
Documentation guidelines
------------------------

The documentation style followed in cosmicpy is taken from the Numpy documentation. 
Guidelines for the documentation are discussed `here <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_. 
An example is provided by Numpy in `example.py <http://github.com/numpy/numpy/blob/master/doc/example.py>`_.



Testing cosmicpy
===============

A test module is included in cosmicpy to check the results of the computations.

Testing against iCosmo
----------------------

A testing suite is implemented to compare cosmicpy against
`iCosmo <http://icosmo.ethz.ch/Initiative_Web/Initiative.html>`_.
If you have :program:`iCosmo` installed, you can save a cosmology structure into
an :program:`IDL` save file using the following commands::
    
    $ idl
    IDL> FID = SET_FIDUCIAL()
    IDL> COSMO = MK_COSMO(FID)
    IDL> SAVE, COSMO, FILENAME='icosmo.xdr'
    
In cosmicpy, the function to test the cosmology in implemented in
:py:meth:`cosmicpy.tests.icosmo.check_cosmology`:
    
.. sourcecode:: ipython

    In [1]: import cosmicpy.tests.icosmo as test
    
    In [2]: test.check_cosmology('icosmo.xdr')
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
        
    Comparing derived constants :
        parameter:        [icosmo] |       [cosmicpy] |            diff | relative error in percents 
        Omega_k  :        0.000000 |        0.000000 |        0.000000 | 0.00%
        gamma    :        0.135207 |        0.135207 |        0.000000 | 0.00%
        sh       :      105.273820 |      105.273820 |        0.000000 | 0.00%
        
    Comparing evolved quantities, maximum relative error:
        Radial Comoving Distance  :  2.27991896578e-08
        Angular Diameter Distance :  2.27991896578e-08
        Hubble constant           :  2.31725750094e-16
        Omega_m(a)                :  4.23296149845e-16
        Omega_de(a)               :  3.50323656321e-16
        dzdr                      :  4.0784692469e-16
        
    Comparing linear power spectra
    Maximum relative error on linear power specrum : 0.000379879249849
            