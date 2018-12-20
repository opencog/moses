
MOSES -- Meta-Optimizing Semantic Evolutionary Search
=====================================================

[![CircleCI](https://img.shields.io/circleci/project/github/singnet/moses/master.svg)](https://circleci.com/gh/singnet/moses/tree/master)

MOSES is a machine-learning tool; it is an "evolutionary program
learner". It is capable of learning short programs that capture
patterns in input datasets.  These programs can be output in either
the `combo` programming language, or in python.  For a given data
input, the programs will roughly recreate the dataset on which they
were trained.

MOSES has been used in several commercial applications, including
the analysis of medical patient and physician clinical data, and
in several different financial systems.  It is also used by OpenCog
to learn automated behaviors, movements and actions in response to
perceptual stimulus of artificial-life virtual agents (i.e. pet-dog
game avatars). Future plans including using it to learn behavioral
programs that control real-world robots, via the OpenPsi implementation
of Psi-theory and ROS nodes running on the OpenCog AtomSpace.

The term "evolutionary" means that MOSES uses genetic programming
techniques to "evolve" new programs. Each program can be thought
of as a tree (similar to a "decision tree", but allowing intermediate
nodes to be any programming-language construct).  Evolution proceeds
by selecting one exemplar tree from a collection of reasonably fit
individuals, and then making random alterations to the program tree,
in an attempt to find an even fitter (more accurate) program.

It is derived from the ideas forumlated in Moshe Looks' thesis,
"Competent Program Evolution", 2006 (Washington University, Missouri)
http://metacog.org/main.pdf.  Moshe is also one of the primary authors
of this code.


LICENSE
-------
MOSES is under double license, Apache 2.0 and GNU AGPL 3.


DOCUMENTATION
-------------
Documentation can be found in the `/docs` directory, which includes a
"QuickStart.pdf" that reviews the algorithms and data structures
used within MOSES.  A detailed man-page can be found in
`/moses/moses/man/moses.1`.  There is also a considerable amount of
information in the OpenCog wiki:
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
> https://launchpad.net/~opencog-dev/+archive/ppa | cxxtest

###### cogutil
> Common OpenCog C++ utilities
> http://github.com/opencog/cogutil
> It uses exactly the same build procedure as this package. Be sure
  to `sudo make install` at the end.

Optional Prerequisites
----------------------
The following packages are optional. If they are not installed, some
optional parts of MOSES will not be built.  The CMake command, during
the build, will be more precise as to which parts will not be built.

###### MPI
> Message Passing Interface
> Required for compute-cluster version of MOSES
> Use either MPICHV2 or OpenMPI |
> http://www.open-mpi.org/ | libopenmpi-dev

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
`-DCMAKE_BUILD_TYPE=Release`
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
MOSES can be used in one of two ways: either directly from the command
line, or by embedding its low-level API into C++ programs. For almost
all users, the command-line interface is strongly recommended.

For those who absolutely must used the low-level C++ programming API,
there is the `/examples` directory.  To build the examples, say:
```
    make examples
```
* example-ant: Santa Fe Trail ant example
* example-data: Simple data sets on which moses can be run.
* example-progs: Other example programs.
