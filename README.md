# TEUFEL
THz Emission from Undulators and Free-Electron Lasers
=======================================

This is a C++ project to support tracking of charged particles in
arbitrary external fields and to compute the electromagnetic radiation
emitted by these particles using the Liénard-Wiechert formula.

The particles are distributed over several compute nodes which comminicate
using OpenMPI. In this computing model a Bunch() is an ensemble of particles
residing on a single compute node. The Beam() may contain a number of bunches
with different properties (like mass and charge of the particles).
The Beam() is distributed over all available compute nodes which
each track their own particles independently. 

Every particle stores its full trajectory data. All field computation is done
after the tracking is finished. Every compute node computes the fields of
its own particles which are then gatherd onto the root node for file output.

Functionality
-------------
- particle tracking using different pusher algorithms (Euler, Vay)<br>
  see \ref ChargedParticle
- radiation emission towards a single observation point or observation screen
- several external field objects (homogeneous dipole, planar undulator ...)
- the trackig problem is defined by a human-readable XML file
- a built-in calculator allows easy computation of input parameters within the input file

Build and Installation
----------------------

We are using the [cmake](https://cmake.org/) build system to allow an easy build on a variety of platforms.
Here we describe a typical out-of-source build. When cmake ist installed it usually defines a variable
CMAKE_MODULE_PATH which set the directories where cmake tries to find its module definitions. TEUFEL adds
(in CMakeLists.txt) it directory teufel/lib to this search path. This way some custom FindXXX.cmake files
can be found.

First one should obtain the sources by cloning the repository from Github.

<pre>
git clone http://github.com/lehnertu/teufel.git
cd teufel
</pre>

A few libraries are required to build the TEUFEL executable.

- For implementing parallel execution on multi-core hardware we use the
  [OpenMPI](https://www.open-mpi.org/) message paasing interface. It should be installed from the
  system distribution. Other MPI implementations may work as well and may even be found by cmake
  but are not tested. For development currently OpenMPI-1.8 ist used.

- To simplify building TEUFEL version 3.4 of the the [SDDS toolkit](https://ops.aps.anl.gov/SDDSInfo.shtml)
  is contained in the repository.

- On most systems the [HDF5 library](https://support.hdfgroup.org/HDF5/) can be installed
  from the distribution repositories.
  It is necessary to set an environment variable for the library to be found by the build system

  ```export HDF5_ROOT=_path_to_library```

- For linear algebra calculations we use the [Eigen3](http://eigen.tuxfamily.org) library.
  It should preferably be installed through the systems package management system.
  In case this is not possible (i.e. no root access) it can be used from the git repository.
  
  ```cd lib/```
  
  ```git clone https://gitlab.com/libeigen/eigen.git```
  
  No build process is required for this library.
  
  A cmake script is provided which will find the library in either case.

- We use [PugiXML](https://pugixml.org/) to parse the XML input files for TEUFEL.
  On most systems this can be installed from the distribution repositories. For Linux-X86 systems we
  provide part of the library in teufel/lib/pugixml.
  The file teufel/lib/FindPUGIXML.cmake is provided for cmake to find the library.
  One has to set an environment variable pointing to the installation directory
  (if in teufel/lib or any other less common installation position -
  typical Linux installation positions are found without this hint) for the library to be found.

  ```export PUGIXML_ROOT=_path_to_library```

- We use [muParser](https://github.com/beltoforion/muparser) to provide
  an inline scientific calculator that allows calculations to be performed
  inside the XML input file. If it is not installed on your system it should
  be cloned an built under teufel/lib/muparser/.
  <pre>
  cd lib
  git clone https://github.com/beltoforion/muparser.git
  cd muparser/
  mkdir build
  cd build/
  cmake ..
  make
  </pre>
  In both cases the script
  teufel/lib/FindMUPARSER.cmake will find and include the library.

Then we create a build directory in the downloaded source directory.

<pre>
mkdir build
cd build
</pre>

Then we build the makefile from CMakeLists.txt contained in the root directory.

```cmake ..```

One can check the libraries and tools found and change the make options.
This can be usefull if it is desired to build the documentation by default
or to skip the build of the test executables (enabled by default).

```ccmake ..```

After that 

```make```

creates the executable in the build directory.

Documentation
-------------

We aim at fully documenting the code for easy reuse and maintenance.
The documentation can be built using doxygen.

```make docs```

The documentation can then be accessed with a browser starting from `doc/html/index.html`.

Testcases
---------

A number of test cases is provided which serve both for code benchmarking
against known results and as coding examples. All tests are built by default
in the build directory. For running all the checks in sequence,
right away from the build directory, a script is provided in the main folder:

```./run_tests```

- circular trajectory of an electron in a homogeneous dipole magnet - tests/teufel.magnet.cpp
- single-electron trajectory in a planar undulator - tests/teufel.undulator.cpp
- field generated by an electron in circular motion at the center of the trajectory - tests/teufel.loop.cpp
- straight trajectory in crossed electric and magnetic fields - tests/teufel.EcrossB.cpp

To check for memory leaks the tool [Valgrind](http://valgrind.org) is recommended.

```valgrind --tool=memcheck tests/teufel.xxx```

For some testcases and examples python scripts for visualizing the data are
provided in the scrips/ directory. For reading HDF5 files these scripts use the
[h5py](http://www.h5py.org/) library which can be installed from the "python-h5py" package on most Linux systems.
