.. Copyright (c) 2014-2015, CosmicPy Developers
.. Licensed under CeCILL 2.1 - see LICENSE.rst
Installation
============

Requirements
------------

CosmicPy requires the following softwares and libraries:

* `CMake <http://www.cmake.org>`_ version 2.6 or above
* `Python <http://www.python.org/>`_ version 2.7
* `Numpy <http://www.numpy.org>`_ version 1.8 or above recommended
* `SciPy <http://www.scipy.org/scipylib/index.html>`_ version 0.14 or above recommended
* `Matplotlib <http://matplotlib.org/>`_ version 1.4.2 or above recommended
* `GSL <http://www.gnu.org/software/gsl>`_ version 1.16 (Optional)

It is also extremely recommended to install the `ipython <http://ipython.org>`_ shell.

All dependencies can be resolved using PyPI except for CMake which can be installed on your system
using a package manager:

* Setting up requirements on **Linux**:

  On `Linux`, just use your favorite package manager to install the dependencies.
  For instance, on Ubuntu Linux::

      $ sudo apt-get install cmake libgsl0-dev


* Setting up requirements on **Mac OS X**:

  The recommended way of installing the dependencies for CosmicPy on a Mac is through `MacPorts <http://www.macports.org/>`_.
  Provided that MacPorts is installed on your system, all the dependencies can be installed with the following command::

      $ sudo port install cmake libgsl0-dev pkgconfig


Install CosmicPy from PyPI
--------------------------

To install CosmicPy, the simplest option is to use the following command

      $ pip install cosmicpy


You can check that CosmicPy has correctly been installed by starting :program:`ipython` and running

.. sourcecode:: ipython

    In [1]: import cosmicpy

If no error message appears then the package is correctly installed.


Install CosmicPy from GitHub
----------------------------

CosmicPy is hosted on a public GitHub repository at this address: https://github.com/cosmicpy/cosmicpy.git

The package can be retrieved with the following command::

    $ git clone https://github.com/cosmicpy/cosmicpy.git

This will clone the latest release of the code into a local folder **cosmicpy**.

Once the requirements are installed and the CosmicPy package is downloaded, go to the cosmicpy source folder::

    $ cd cosmicpy

and run the following command::

    $ sudo python setup.py install

This should compile and install CosmicPy on your computer and it should be ready to use.
