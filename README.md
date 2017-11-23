# ModalSound 
[![Build Status](https://travis-ci.org/dingzeyuli/ModalSound.svg?branch=master)](https://travis-ci.org/dingzeyuli/ModalSound)

This project releases the educational C++ code used at SIGGRAPH 2016 Course "Physically Based Sound for Computer Animation and Virtual Environments".

There are two sub-projects, _IsoStuffer_ and _ModalSound_.

_IsoStuffer_ is an implementation of [[Labelle and Shewchuk 2007]](http://www.cs.berkeley.edu/~jrs/papers/stuffing.pdf). It provides the functionality to construct a tetrahedral mesh from a water-tight triangle mesh.

_ModalSound_ is a demostration of synthesizing rigid body sounds using the modal sound model. It includes the code of computing stiffness and mass matrices from a volumetric mesh and use the precomputed acoustic multipole coefficients to evaluate the acoustic transfer values at runtime.
