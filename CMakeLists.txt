# Copyright (c) 2014-2015, CosmicPy Developers
# Licensed under CeCILL 2.1 - see LICENSE.rst
cmake_minimum_required(VERSION 2.8)

include(ExternalProject)
include(FindPkgConfig)

# Adding customized cmake module
list(APPEND CMAKE_MODULE_PATH  "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(cosmicpy)

# Find required packages
find_package(Boost REQUIRED COMPONENTS math_tr1 python)
find_package(PythonInterp REQUIRED)
find_package(PythonLibsNew REQUIRED)

pkg_check_modules(GSL REQUIRED gsl)

# Downloads and compiles the Boost.NumPy package
# The library and include files are located in the build/extern directory
ExternalProject_Add(BoostNumpy
    GIT_REPOSITORY https://github.com/ndarray/Boost.NumPy
    CMAKE_ARGS     -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/extern
                   -DLIBRARY_TYPE=STATIC
)
set(BoostNumpy_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/extern/include)
set(BoostNumpy_LIBRARY_DIRS ${CMAKE_BINARY_DIR}/extern/lib)
set(BoostNumpy_LIBRARIES -lboost_numpy)


# Define include and library directories
include_directories(
  ${GSL_INCLUDE_DIRS}
  ${PYTHON_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${BoostNumpy_INCLUDE_DIRS}
)
link_directories(${BoostNumpy_LIBRARY_DIRS} ${GSL_LIBRARY_DIRS})

# Compilation flags
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang, disabling OpenMP support
  set(CMAKE_CXX_FLAGS "-O3 -fomit-frame-pointer -fno-common -fPIC")
else()
  set(CMAKE_CXX_FLAGS "-O3 -fomit-frame-pointer -fno-common -fPIC -fopenmp")
endif()


# Build the tools module
add_library(tools SHARED cosmicpy/cxx/tools.cpp cosmicpy/cxx/besselwindow.cpp)
add_dependencies(tools BoostNumpy)
target_link_libraries(tools ${BoostNumpy_LIBRARIES} ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${PYTHON_LIBRARIES})
set_target_properties(tools PROPERTIES SUFFIX .so)
set_target_properties(tools PROPERTIES PREFIX "")

# Installs the module inside the python package 
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})


install(TARGETS tools DESTINATION cosmicpy)

