# ModalSound 
[![Build Status](https://travis-ci.org/dingzeyuli/ModalSound.svg?branch=master)](https://travis-ci.org/dingzeyuli/ModalSound)

This project releases the educational C++ code used at the following SIGGRAPH 2016 Course. 

Physically Based Sound for Computer Animation and Virtual Environments
Doug L. James, Timothy R. Langlois, Ravish Mehra, and Changxi Zheng

The slides and lecture notes can be found on the [course webpage](http://graphics.stanford.edu/courses/sound/).  The conference video recording is [here](https://youtu.be/4OGeAfyDa4Y).

There are two sub-projects, _IsoStuffer_ and _ModalSound_.

1. Given a water-tight triangle surface mesh  _IsoStuffer_ construct a volumetric tetrahedral mesh.

2. _ModalSound_ is a demostration of synthesizing rigid body sounds using the modal sound model. It includes the code of computing stiffness and mass matrices from a volumetric tetrahedral mesh and use the precomputed acoustic multipole coefficients to evaluate the acoustic transfer values at runtime.

Build instructions are given in the [respective](IsoStuffer/) [folders](ModalSound/). 

