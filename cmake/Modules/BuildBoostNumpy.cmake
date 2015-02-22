# Copyright (c) 2014-2015, CosmicPy Developers
# Licensed under CeCILL 2.1 - see LICENSE.rst
#==============================================================#
# Build the Boost.NumPy dependencies for the project using a specific version of python  #
#==============================================================#
 
 # Downloads and compiles the Boost.NumPy package
# The library and include files are located in the build/extern directory
ExternalProject_Add(BoostNumpy
    PREFIX BoostNumpy
    GIT_REPOSITORY https://github.com/ndarray/Boost.NumPy
    CMAKE_ARGS     -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/extern
                                     -DLIBRARY_TYPE=STATIC
                                     -DBUILD_EXAMPLES=OFF
                                     -DBUILD_TEST=OFF
                                     -DBOOST_ROOT=${CMAKE_BINARY_DIR}/extern
)
set(BoostNumpy_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/extern/include)
set(BoostNumpy_LIBRARY_DIRS ${CMAKE_BINARY_DIR}/extern/lib)
set(BoostNumpy_LIBRARIES -lboost_numpy)
add_dependencies(BoostNumpy Boost)
