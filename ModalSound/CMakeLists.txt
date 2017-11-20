# CMake Project file for the basic research code 
# Changxi Zheng (cxzheng@cs.cornell.edu)
# Copyright @ Jul 2019
project(sploosh)

cmake_minimum_required(VERSION 3.0)

# General Configuration: find all the required libraries.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(UpdateCacheVariable)

list(APPEND CMAKE_PREFIX_PATH "/home/dli/Qt5.7.0/5.7/gcc_64/lib/cmake")

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
set(MKL_LIBS ${IOMP5_LIBRARY} ${MKL_LIBRARY_RT} ${MKL_LIBRARY_THREAD}  ${MKL_LIBRARY_ILP64} ${MKL_LIBRARY} ${OPENMP_LIB})
message(${MKL_LIBS})
add_definitions(-DUSE_MKL)
add_definitions(-DEIGEN_USE_MKL_ALL)
add_definitions(-Wl,--no-as-needed)
add_definitions(-lpthread)
add_definitions(-lm)

# Boost
set(Boost_ADDITIONAL_VERSIONS "1.51" "1.51.0")
find_package(Boost 1.50 REQUIRED 
    COMPONENTS "filesystem" "system" "program_options")

# GSL
find_package(GSL REQUIRED)

# Eigen
#find_package(Eigen REQUIRED)

if (EIGEN_FOUND)
  message("-- Found Eigen version ${EIGEN_VERSION}: ${Eigen_INCLUDE_DIR}")
endif (EIGEN_FOUND)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include/eigen
)

include(sploosh_external)

#===================================================================

if ( BUILD_TOOLS OR BUILD_DEMO )
    # package for opengl and glut
    find_package(GLUT REQUIRED)
    find_package(OpenGL REQUIRED)


    # check the existence of libQGLViewer
    if(APPLE AND EXISTS /usr/local/opt/qt5)
	    # Homebrew installs Qt5 (up to at least 5.9.1) in
	    # /usr/local/qt5, ensure it can be found by CMake since
	    # it is not in the default /usr/local prefix.
	    list(APPEND CMAKE_PREFIX_PATH "/usr/local/opt/qt5")
    endif()
    find_package(Qt5 COMPONENTS Core Widgets Xml OpenGL Multimedia)
    find_package(QGLViewer REQUIRED)
endif ()

add_subdirectory(src)
