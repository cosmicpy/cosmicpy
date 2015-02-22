#!/usr/bin/env python
# Copyright (c) 2014-2015, CosmicPy Developers
# Licensed under CeCILL 2.1 - see LICENSE.rst

from distutils.command.build_py import build_py as _build_py
from distutils.core import setup
from subprocess import call
import os


class build_py(_build_py):
    """Specialized Python source builder."""
    def run(self):
        call(["mkdir", "-p", "build"])
        call(["cmake", "-H.", "-Bbuild"])
        call(["make", "-Cbuild", "install"])
        _build_py.run(self)

setup(name='cosmicpy',
      version='0.1',
      requires=['numpy', 'scipy', 'matplotlib'],
      description='An interactive python package for cosmology and parameter forecasts',
      author='Francois Lanusse, Anais Rassat',
      author_email='francois.lanusse@cea.fr',
      url='http://cosmicpy.github.io',
      provides=['cosmicpy'],
      packages=['cosmicpy','cosmicpy.tests'],
      package_data={'cosmicpy': ['tools.so']},
      cmdclass={'build_py': build_py},
      license='CeCILL',
      classifiers=[
	  'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: POSIX :: Linux',
          'Programming Language :: C++',
          'Programming Language :: Python :: 2.7',
          'Topic :: Scientific/Engineering :: Physics'
      ]
      )
