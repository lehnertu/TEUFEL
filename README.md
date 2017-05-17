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

A number of test cases is provided which serve both for code benchmarking
against known results and a coding examples. All tests can be built from
the main folder  by

```make tests```

For running all the checks in
sequence, a script is provided in the main folder:

```./run_tests```

- circular trajectory of an electron in a homogeneous dipole magnet - tests/teufel.magnet.cpp
- single-electron trajectory in a planar undulator - tests/teufel.undulator.cpp
- field generated by an electron in circular motion at the center of the trajectory - tests/teufel.loop.cpp
- straight trajectory in crossed electric and magnetic fields - tests/teufel.EcrossB.cpp

To check for memory leaks the tool [Valgrind](http://valgrind.org) is recommended.

```valgrind --tool=memcheck tests/teufel.xxx```