# Copyright (c) 2014-2015, CosmicPy Developers
# Licensed under CeCILL 2.1 - see LICENSE.rst
#========================================================#
# Build the Boost dependencies for the project using a specific version of python #
#========================================================#

set(BoostVersion 1.57.0)
set(BoostMD5 1be49befbdd9a5ce9def2983ba3e7b76)

string(REGEX REPLACE "beta\\.([0-9])$" "beta\\1" BoostFolderName ${BoostVersion})
string(REPLACE "." "_" BoostFolderName ${BoostFolderName})
set(BoostFolderName boost_${BoostFolderName})

ExternalProject_Add(Boost
    PREFIX Boost
    URL  http://sourceforge.net/projects/boost/files/boost/${BoostVersion}/${BoostFolderName}.tar.bz2/download
    URL_MD5 ${BoostMD5}
    CONFIGURE_COMMAND ./bootstrap.sh
                                                        --with-libraries=python,math
                                                        --with-python=${PYTHON_EXECUTABLE}
    BUILD_COMMAND           ./b2 install 
                                                        variant=release
                                                        link=static
                                                        cxxflags='-fPIC'
                                                        --prefix=${CMAKE_BINARY_DIR}/extern
                                                        -d 0
                                                        -j8
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
    )

set(Boost_LIBRARY_DIR ${CMAKE_BINARY_DIR}/extern/lib/ )
set(Boost_INCLUDE_DIR ${CMAKE_BINARY_DIR}/extern/include/boost/ )

if(${PYTHON_VERSION_STRING} GREATER 3.0)
  message(STATUS "Using Python3")
  set(Boost_LIBRARIES -lboost_math_tr1 -lboost_python3)
else()
  message(STATUS "Using Python2")
  set(Boost_LIBRARIES -lboost_math_tr1 -lboost_python)
endif()
