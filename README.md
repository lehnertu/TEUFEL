# TEUFEL
THz Emission from Undulators and Free-Electron Lasers
=======================================

This is a C++ project to support tracking of charged particles in
arbitrary external fields and to compute the electromagnetic radiation
emitted by these particles using the Liénard-Wiechert formula.

Functionality
=============
- particle tracking using different pusher algorithms (Euler, Vay)<br>
  see \ref ChargedParticle
- radiation emission towards a single observation point
- several external field objects (homogeneous dipole, planar undulator ...)

Testcases
=========
- circular trajectory of an electron in a homogeneous dipole magnet - tests/teufel.magnet.cpp
- single-electron trajectory in a planar undulator - tests/teufel.undulator.cpp
