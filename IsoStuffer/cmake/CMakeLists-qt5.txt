# CMake Project file for isosurface stuffing  
# Changxi Zheng (cxzheng@cs.cornell.edu)
# Copyright @ Jan 2009
project(isostuffer)

cmake_minimum_required(VERSION 2.6)

# General Configuration: find all the required libraries.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# compile option
option(USE_DEBUG "Turn on the debug mode" OFF)
option(USE_OPENMP "Turn on the OpenMP feature when compiling" OFF)
option(USE_GUI   "Compile the program with GUI" ON)
option(USE_64BIT_ARCH "Compile the 64bit executable" OFF)

#===================================================================
## Compiler
# set compiler flags for debug/release
if ( USE_DEBUG )
    add_definitions(-DDEBUG)
    set(CMAKE_BUILD_TYPE Debug)
else ( USE_DEBUG )
    set(CMAKE_BUILD_TYPE Release)
endif ( USE_DEBUG )

if ( USE_64BIT_ARCH )
    add_definitions(-m64)
    set(LINK_FLAGS -m64)
endif ( USE_64BIT_ARCH )

#===================================================================
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()
#===================================================================
## Libraries
# check boost libraries
find_package(Boost 1.48 REQUIRED COMPONENTS
    "filesystem" "system" "program_options")

if ( USE_GUI )
    # package for opengl and glut
    find_package(OpenGL REQUIRED)
    find_package(GLUT REQUIRED)
    find_package(Qt5 COMPONENTS Core Widgets Xml OpenGL)

    # check the existence of libQGLViewer
    find_package(QGLViewer 2.0 REQUIRED)
endif ( USE_GUI )

add_definitions(-Wall)

if ( USE_OPENMP )
    add_definitions(-DUSE_OPENMP)
endif ( USE_OPENMP )

set(LINK_FLAGS " ")
# check if we're using Intel's compiler
if ( CMAKE_CXX_COMPILER MATCHES ".*icpc$" )
    add_definitions(-wd981 -wd383 -wd444 -wd1224 -wd1572)
    if ( NOT USE_DEBUG )
        set(CMAKE_CXX_FLAGS_RELEASE "-O3 -no-prec-div -xHost -DNDEBUG")
    endif ( NOT USE_DEBUG )

    if ( USE_OPENMP )
        find_package(Threads 1.0 REQUIRED)
        set(ICC_LIBS iomp5 ${CMAKE_THREAD_LIBS_INIT})
        add_definitions(-openmp)
    endif ( USE_OPENMP )
else ( CMAKE_CXX_COMPILER MATCHES ".*icpc$" )
    if ( USE_OPENMP )
        add_definitions(-fopenmp)
    endif ( USE_OPENMP )
endif ( CMAKE_CXX_COMPILER MATCHES ".*icpc$" )

add_subdirectory(src)

