# CMake Project file for the basic research code 
# Changxi Zheng (cxzheng@cs.cornell.edu)
# Copyright @ Jul 2019
project(sploosh)

cmake_minimum_required(VERSION 3.0)

# General Configuration: find all the required libraries.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(UpdateCacheVariable)

# compile option
option(BUILD_DEBUG "Turn on the debug mode" OFF)
option(BUILD_TOOLS "Build tools" ON)
option(BUILD_DEMO "Build demos" ON)

include(sploosh_utils)
include(CMakeDependentOption)

set(BIN_DIR ${CMAKE_CURRENT_BINARY_DIR}/bin)

#===================================================================
## Compiler
# set compiler flags for debug/release
if ( BUILD_DEBUG )
    add_definitions(-DDEBUG)
    set(CMAKE_BUILD_TYPE Debug)
else ()
    set(CMAKE_BUILD_TYPE Release)
endif ()

add_definitions(-Wall)

# only a fix on Mac to disable boost-caused warnings
if ( APPLE )
    add_definitions(-Wno-unused-local-typedef)
endif ()

config_compiler_and_linker()

#===================================================================
## Libraries
find_package(Threads 1.0 REQUIRED)

# check Intel MKL
find_package(MKL REQUIRED)
set(MKL_LIBS mkl_intel_lp64 mkl_intel_thread mkl_core ${OPENMP_LIB})
add_definitions(-DUSE_MKL)
add_definitions(-DEIGEN_USE_MKL_ALL)

# Boost
set(Boost_ADDITIONAL_VERSIONS "1.51" "1.51.0")
find_package(Boost 1.50 REQUIRED 
    COMPONENTS "filesystem" "system" "program_options")

# GSL
find_package(GSL REQUIRED)

# Eigen
find_package(Eigen REQUIRED)

if (EIGEN_FOUND)
  if (EIGEN_VERSION VERSION_LESS 3.1.0)
    message(FATAL_ERROR "-- Ceres requires Eigen version >= 3.1.0 in order "
      "that Eigen/SparseCore be available, detected version of Eigen is: "
      "${EIGEN_VERSION}")
  endif (EIGEN_VERSION VERSION_LESS 3.1.0)

  message("-- Found Eigen version ${EIGEN_VERSION}: ${EIGEN_INCLUDE_DIRS}")
endif (EIGEN_FOUND)

include(sploosh_external)

#===================================================================

if ( BUILD_TOOLS OR BUILD_DEMO )
    # package for opengl and glut
    find_package(GLUT REQUIRED)
    find_package(OpenGL REQUIRED)

    # check the existence of libQGLViewer
    find_package(Qt5 COMPONENTS Core Widgets Xml OpenGL Multimedia)
    find_package(QGLViewer 2.0 REQUIRED)
endif ()

add_subdirectory(src)
