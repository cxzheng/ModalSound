# ModalSound

_ModalSound_ synthesizes rigid body sounds using the modal sound model. The implementation is based on the following SIGGRAPH papers.

```
Toward High-Quality Modal Contact Sound.
Changxi Zheng, Doug L. James
SIGGRAPH 2011
```
```
Rigid-Body Fracture Sound with Precomputed Soundbanks
Changxi Zheng,  Doug L. James
SIGGRAPH 2010
```
```
Precomputed Acoustic Transfer: Output-sensitive, accurate sound generation for geometrically complex vibration sources
Doug L. James, Jernej Barbić, Dinesh K. Pai
SIGGRAPH 2006
```

## Input & Output

To be updated. 

##  Dependencies

Isostuffer requires Boost, Intel MKL, eigen3, GSL, Qt5, and libQGLViewer(Qt5) dependencies.

## Build instructions

The build has been tested with Intel compiler (Ubuntu), GNU compiler (Ubuntu), and Clang (macOS) . For this instruction I will use GNU compiler.

### Ubuntu 16.04 LTS

1. Install GNU compiler gcc-6 and g++-6:
```bash
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install gcc-6 g++-6
```

2. Install the APT dependencies:
```bash
sudo apt-get install git libboost-all-dev xorg-dev libglu1-mesa-dev libgoogle-glog-dev freeglut3-dev qtmultimedia5-dev  qt5-default mesa-utils libqglviewer-dev  libxmu-dev libxi-dev libconfig++-dev libprotobuf-dev libgsl-dev protobuf-compiler
```

3. Install Intel MKL library.

For Ubuntu, the follow script works:
```bash
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
sudo sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
sudo apt-get update
sudo apt-get install intel-mkl-64bit-2017.4-061
```
For other platforms, please visit the [official page for Intel MKL](https://software.intel.com/en-us/mkl) for installation hints. 

4. Install the CMake 3.9.6 -- I tried with an earlier version and it didn't work with Qt5 UIC, 3.9.6 is known to work.
```bash
cd path-to-download-and-compile-cmake
version=3.9.6
wget https://cmake.org/files/v3.9/cmake-$version.tar.gz
tar -xzvf cmake-$version.tar.gz 1> /dev/null
cd cmake-$version/
./bootstrap 1> /dev/null
make -j4 1>/dev/null
sudo make install
```

5. Assuming the repo has been cloned, now change into the subfolder "ModalSound/ModalSound": 
```bash
cd path-to/ModalSound/ModalSound/
```

6. Download eigen3 into local folder:
```bash
chmod 755 get_eigen.sh
./get_eigen.sh
```

7. Create a build directory under the project root and change into this directory:
```bash
mkdir gcc-build
cd gcc-build
```

8. Initialize the Intel MKL variables:
```bash
source /opt/intel/mkl/bin/mklvars.sh intel64
source /opt/intel/bin/compilervars.sh intel64
```

9. Run CMake to create the build system using gcc-6 and g++-6:
```bash
CC=gcc-6 CXX=g++-6 cmake ..
```

10. Build ModalSound:

```bash
make -j
```

11. Run ModalSound:

```bash
 ./bin/tetviewer
 ./bin/click_synth
```

###  Other Linux distros & macOS

For macOS, I use homebrew to download most of the required packages. If certain packages are missing from brew, install them from source code. As for other Linux distros, similar steps follow.

### Windows

I have not tried compiling Isostuffer on Windows machines. Please let me know if you successfully compile and run this code in Windows. 

### Travis CI

I've prepared an automatic continuous integration script for [Travis CI](https://travis-ci.org/dingzeyuli/ModalSound). Checkout [the script](https://github.com/dingzeyuli/ModalSound/blob/master/.travis.yml) for the full automatic process.

## Issues

If you run into problems during compiling or running, please checkout [the common issues page](../Issues.md). 