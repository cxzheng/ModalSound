# IsoStuffer

_IsoStuffer_ constructs a tetrahedral mesh from a water-tight triangle mesh. It is an implementation of the following [paper](http://www.cs.berkeley.edu/~jrs/papers/stuffing.pdf). 

```
Isosurface Stuffing: Fast Tetrahedral Meshes with Good Dihedral Angles
François Labelle, Jonathan Shewchuk
SIGGRAPH 2007
```

We have a command-line only version and a GUI version. The screenshots illustrate some typical use cases.

## Input & Output

The input is a water-tight triangular surface mesh in obj format. The output is the tetrahedronized version. By default the export format is *.tet which can be opened and visualized by tetviewer in ModalSound project. It is also possible to export to TetGen format.

##  Dependencies

Isostuffer requires Boost and several openGL-related dependencies:

It is also highly recommended to install Qt5 and libQGLViewer(Qt5) to enable the GUI option. 

## Build instructions

The build has been tested with Intel compiler and GNU compiler. For this instruction I will use GNU compiler.

### Ubuntu 16.04 LTS

1. Install GNU compiler gcc-6 and g++-6:
```bash
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install gcc-6 g++-6
```

2. Install the dependencies:
```bash
sudo apt-get install libboost-all-dev xorg-dev libgl u1-mesa-dev freeglut3-dev qtmultimedia5-dev  qt5-default mesa-utils libqglviewer-dev  libxmu-dev libxi-dev libconfig++-dev libprotobuf-dev libgsl-dev protobuf-compiler
```

3. Clone this repository and change into the project root: 
```bash
git clone https://github.com/dingzeyuli/ModalSound.git
cd ModalSound
```

4. Create a build directory under the project root and change into this directory:
```bash
mkdir gcc-build
cd gcc-build
```

5. Run CMake to create the build system using gcc-6 and g++-6:
```bash
CC=gcc-6 CXX=g++-6 cmake ..
```

In case you do not have Qt5 or libqglviewer installed, run cmake with the GUI option disabled:
```bash
CC=gcc-6 CXX=g++-6 cmake -USE_GUI=OFF ..
```
6. Build Isostuffer:

```bash
make -j
```
7. Run Isostuffer:

Without GUI:
```bash
 ./src/isostuffer
```
With GUI:
```bash
 ./src/isostuffer-gui
```

###  Other Linux distros & macOS

For macOS, I use homebrew to download most of the required packages. If certain packages are missing from brew, install them from source code. As for other Linux distros, similar steps follow.



### Windows

I have not tried compiling Isostuffer on Windows machines. Please let me know if you successfully compile and run this code in Windows. 

### Travis CI

I've prepared an automatic continuous integration script for [Travis CI](https://travis-ci.org/dingzeyuli/ModalSound). Checkout [the script](https://github.com/dingzeyuli/ModalSound/blob/master/.travis.yml) for the full automatic process.

## Issues

If you run into problems during compiling or running, please checkout [the common issues page](../Issues.md). 