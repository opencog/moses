
MOSES -- Meta-Optimizing Semantic Evolutionary Search
=====================================================
                         August, 2012


INTRODUCTION
------------
MOSES is an evolutionary program learner. It is mostly based on Moshe
Looks' thesis, "Competent Porgram Evolution", 2006 (Washington 
University, Missouri) http://metacog.org/main.pdf. Moshe is also one of
the primary authors of this code.


LICENSE
-------
MOSES is under double license, Apache 2.0 and GNU AGPL 3.


DOCUMENTATION
-------------
Documentation can be found in the docs directory, which includes a 
"QuickStart.pdf" that reviews the algorithms and data structures
used withie MOSES.  There is alo a considerable amount of information
int the OpenCog wiki: 
http://wiki.opencog.org/w/Meta-Optimizing_Semantic_Evolutionary_Search

Prerequisites
-------------
To build and run MOSES, the packages listed below are required. With a
few exceptions, most Linux distributions will provide these packages.

###### boost
> C++ utilities package
> http://www.boost.org/ | libboost-dev

###### cmake
> Build management tool; v2.8 or higher recommended.
> http://www.cmake.org/ | cmake

###### cxxtest
> Test framework
> http://cxxtest.sourceforge.net/ |
> https://launchpad.net/~opencog-dev/+archive/ppa
> Currently, opencog requires cxxtest version 3, and is not compatible
  with version 4.

###### cogutil
> Common OpenCog C++ utilities
> http://github/opencog/cogutils
> It uses exactly the same build proceedure as this pakcage. Be sure
  to `sudo make install` at the end.

Optional Prerequisites
----------------------
The following packages are optional. If they are not installed, some
optional parts of OpenCog will not be built.  The CMake command, during
the build, will be more precise as to which parts will not be built.

###### MPI
> Message Passing Interface
> Required for compute-cluster version of MOSES
> Use either MPICHV2 or OpenMPI

Building MOSES
--------------
Peform the following steps at the shell prompt:
```
    cd to project root dir
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make
```
Libraries will be built into subdirectories within build, mirroring the
structure of the source directory root. The flag
-DCMAKE_BUILD_TYPE=Release
results in binaries that are optimized for for performance; ommitting
this flag will result in faster builds, but slower executables.

Unit tests
----------
To build and run the unit tests, from the ./build directory enter (after
building moses as above):
```
    make test
```

INSTALLATION
------------
Just say `sudo make install`  after finishing the build.


EXAMPLES DIRECTORY
------------------
To build the examples, say:
```
    make examples
```
> example-ant: Santa Fe Trail ant example
> example-data: Simple data sets on which moses can be run.
> example-progs: Other example programs.
