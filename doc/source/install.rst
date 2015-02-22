Installation
============

Requirements
------------

CosmicPy requires the following softwares and libraries:

* `CMake <http://www.cmake.org>`_ version 2.6 or above
* `Python <http://www.python.org/>`_ version 2.7
* `Boost <http://www.boost.org>`_ version 1.50 or above, boost.math and boost.python packages required 
* `GSL <http://www.gnu.org/software/gsl>`_ version 1.16 or above recommended
* `Numpy <http://www.numpy.org>`_ version 1.8 or above recommended
* `SciPy <http://www.scipy.org/scipylib/index.html>`_ version 0.14 or above recommended
* `Matplotlib <http://matplotlib.org/>`_ version 1.4.2 or above recommended

It is also extremely recommended to install the `ipython <http://ipython.org>`_ shell. The different requirements can
easily be installed on your system using a package manager:

* Setting up requirements on **Linux**:

  On `Linux`, just use your favorite package manager to install the dependencies.
  For instance, on Ubuntu Linux::
      
      $ sudo apt-get install cmake boost gsl ipython numpy scipy matplotlib


* Setting up requirements on **Mac OS X**:

  The recommended way of installing the dependencies for CosmicPy on a Mac is through `MacPorts <http://www.macports.org/>`_.
  Provided that MacPorts is installed on your system, all the dependencies can be installed with the following command::
      
      $ sudo port install cmake boost gsl py27-ipython py27-numpy py27-scipy py27-matplotlib pkgconfig

  .. Note::
      Due to the Apple's removal of gcc in OSX Mavericks and Yosemite, the default macports compiler is now clang.
      Because **clang lacks OpenMP support**, it is not currently easily possible to compile CosmicPy with OpenMP
      parallelisation under Mac OS X. It is possible but requires compiling boost with gcc.


Download from GitHub
--------------------

CosmicPy is hosted on a public GitHub repository at this address: https://github.com/cosmicpy/cosmicpy.git

The package can be retrieved with the following command::
    
    $ git clone https://github.com/cosmicpy/cosmicpy.git

This will clone the latest release of the code into a local folder **cosmicpy**.


Install
-------

Once the requirements are installed and the CosmicPy package is downloaded, go to the cosmicpy source folder::

    $ cd cosmicpy
    
and run the following command::
    
    $ sudo python setup.py install
    
This should compile and install CosmicPy on your computer and it should be ready to use.


You can check that CosmicPy has correctly been installed by starting :program:`ipython` and running

.. sourcecode:: ipython

    In [1]: import cosmicpy
    
If no error message appears then the package is correctly installed.