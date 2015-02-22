# Copyright (c) 2014-2015, CosmicPy Developers
# Licensed under CeCILL 2.1 - see LICENSE.rst
cmake_minimum_required(VERSION 2.8)

set(GSLVersion 1.16)
set(GSLMD5 e49a664db13d81c968415cd53f62bc8b)

# Adding customized cmake module
list(APPEND CMAKE_MODULE_PATH  "${CMAKE_SOURCE_DIR}/cmake/Modules/")

include(ExternalProject)
include(FindPkgConfig)

project(cosmicpy)

# Find required packages
find_package(PythonInterp REQUIRED)
find_package(PythonLibsNew REQUIRED)

# Try to find the GSL dependency
pkg_check_modules(GSL gsl)

# If GSL is not installed, lets go ahead and compile it
if(NOT GSL_FOUND)
   message(STATUS "GSL not found, downloading and compiling from source")
    ExternalProject_Add(GSL
        PREFIX GSL
        URL http://ftp.igh.cnrs.fr/pub/gnu/gsl/gsl-${GSLVersion}.tar.gz
        URL_MD5 ${GSLMD5}
        CONFIGURE_COMMAND ./configure
                                                            --prefix=${CMAKE_BINARY_DIR}/extern
                                                            --enable-shared=no
                                                            --with-pic=yes
        BUILD_COMMAND             make -j8
        INSTALL_COMMAND         make install
        BUILD_IN_SOURCE 1)
        set(GSL_LIBRARY_DIR ${CMAKE_BINARY_DIR}/extern/lib/ )
        set(GSL_INCLUDE_DIR ${CMAKE_BINARY_DIR}/extern/include/)
        set(GSL_LIBRARIES -lgsl -lgslcblas -lm)
endif()


# Compile boost packages to match the specific system environment
include(BuildBoost)
include(BuildBoostNumpy)

# Define include and library directories
include_directories(
  ${Boost_INCLUDE_DIRS}
  ${BoostNumpy_INCLUDE_DIRS}
  ${GSL_INCLUDE_DIRS}
  ${PYTHON_INCLUDE_DIRS}
)
link_directories(${BoostNumpy_LIBRARY_DIRS} ${Boost_LIBRARY_DIRS}  ${GSL_LIBRARY_DIRS})

# Compilation flags
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang, disabling OpenMP support
  set(CMAKE_CXX_FLAGS "-O3 -fomit-frame-pointer -fno-common -fPIC")
else()
  set(CMAKE_CXX_FLAGS "-O3 -fomit-frame-pointer -fno-common -fPIC -fopenmp")
endif()

# Build the tools module
add_library(tools SHARED cosmicpy/cxx/tools.cpp cosmicpy/cxx/besselwindow.cpp)
if(NOT GSL_FOUND)
add_dependencies(tools GSL)
endif()
add_dependencies(tools BoostNumpy Boost)
target_link_libraries(tools ${BoostNumpy_LIBRARIES} ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${PYTHON_LIBRARIES})
set_target_properties(tools PROPERTIES SUFFIX .so)
set_target_properties(tools PROPERTIES PREFIX "")

# Installs the module inside the python package 
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})


install(TARGETS tools DESTINATION cosmicpy)

