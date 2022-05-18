# README

UltraCold is a modular and extensible collection of C++ libraries for the study of ultra-cold atomic systems in the 
context of Gross-Pitaevskii theory. 

## What UltraCold can do

UltraCold contains several solver classes for different flavors of Gross-Pitaevskii and Bogolyubov equations, allowing
for the description of ultra-cold systems of bosons at the mean-field level, studying their ground-state properties,
the dynamics, and elementary excitations.

Right now, all the solver classes take advantage of OpenMP parallelization.

The complete documentation is available either in [html](https://smroccuzzo.github.io/UltraCold) or
[pdf](./docs/manual.pdf).   

## Prerequisites and platforms

UltraCold is built on top of Intel's Math Kernel Library, and relies
upon [arpack-ng](https://github.com/opencollab/arpack-ng) (which is provided as a bundled package) for the solution of 
Bogolyubov equations. Hence, in order to use UltraCold, you first need to download and install a distribution of Intel's
software.

Right now, the package has been only tested with Intel oneAPI, although it should also work with previous versions 
of Intel Parallel Studio. The Intel oneAPI package can be downloaded **for free** from
[here](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html). 
In particular, UltraCold relies on the
[Intel oneAPI Base Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html) 
and on the
[Intel oneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html).

The package has been tested only on Linux machines, including the High Performance Computing cluster
[Galileo100](https://www.hpc.cineca.it/hardware/galileo100) from the italian supercomputing 
consortium [CINECA](https://www.cineca.it/). 

## Installation

To get UltraCold, clone it into your machine

    git clone https://github.com/smroccuzzo/UltraCold.git

Then, enter the directory UltraCold, and follow the usual steps required to build a project using
[cmake](https://cmake.org/). By default, the build type is *Release*. So, all you have to do is

    cd UltraCold
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=</your/install/path> ..
    make
    make install

##  Usage and examples

UltraCold comes packed with several solver classes for different flavors of Gross-Pitaevskii and Bogolyubov
equations. The complete list of available solvers, as well as other useful classes (e.g. for data output) is
available under the namespace [UltraCold](https://smroccuzzo.github.io/UltraCold/html/namespace_ultra_cold.html). 

The [examples](./examples) folder contains several folders called *example-n*, each containing

- a source file called *example-n.cpp*, containing an example on how to write an executable that uses UltraCold 
  libraries,
- a *CMakeLists.txt* containing instructions on how to configure and build an executable based on
  UltraCold,
- an eventual parameter file, called *example-n.prm*,

To run the examples, follow the usual steps required to build a project using *cmake*, namely, open a
terminal in the folder containing the example you are interested in and type

    mkdir build
    cd build
    cmake -DULTRACOLD_DIR=/path/to/the/directory/where/you/installed/UltraCold ..
    make

This will create an executable called *example-n*, which (if everything went fine) should be ready to be
executed. Now you can just copy the eventual file *example-n.prm* from the parent folder and run the example

    cp ../example-n.prm .
    ./example-n

The output of course will depend on the particular example, and is documented fully for each of them.

Although a prior knowledge of C++ is highly recommended, to use these examples also very basic knowledge is more
than sufficient. The documentation tries to be as pedantic as possible, so that extending these examples for user's
need shouldn't be too difficult.

Here is the complete list of all examples and a brief description of what each of them does. Refer to the detailed
description available in the documentation.

- [example-1](./examples/example-1) Defined in file [example-1.cpp](./examples/example-1/example-1.cpp) Ground state 
 and simple dynamics of a three-dimensional, harmonically trapped Bose gas.
- [example-2](./examples/example-2) Defined in file [example-2.cpp](./examples/example-2/example-2.cpp) Elementary 
  excitations via Bogolyubov equations in a two-dimensional, harmonically trapped Bose gas.
- [example-3](./examples/example-3) Defined in file [example-3.cpp](./examples/example-3/example-3.cpp) Supersolid 
  ground state of a trapped dipolar Bose-Einstein condensate.
- [example-4](./examples/example-4) Defined in file [example-4.cpp](./examples/example-4/example-4.cpp) Excitation 
spectrum of a trapped dipolar Bose-Einstein condensate across the superfluid-supersolid phase transition.

## Contributing

UltraCold is developed using [Git](https://git-scm.com/) as a version control tool, and [GitHub](https://github.com) as
the central host of the source code.

If you find some issue, and/or have suggestions for additions and/or improvements, please open a 
[GitHub issue](https://docs.github.com/en/issues/tracking-your-work-with-issues).

If you want to actively contribute, after opening an issue, follow this (pretty standard) workflow, based on the 
[fork and pull model](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/getting-started/about-collaborative-development-models):

- If you already have a GitHub account, [sign in](https://github.com/login). Otherwise, 
[create one](https://github.com/signup?ref_cta=Sign+up&ref_loc=header+logged+out&ref_page=%2F&source=header-home)
(it's free, and we are pretty confident that it will always be).
- Fork the UltraCold project.
- Make a local clone of your fork to your own computer.
- Create a new branch on which you will be making changes. This marks the point from which your copy of the project
starts to differ from that of the main development branch. 
- Start to make your changes, for example by modifying existing files and/or creating new ones. Once you are
satisfied with your changes, you can commit each change, in such a way that Git can keep track of them. With each commit,
write a short message describing what your particular set of changes does.
- When you're finished committing all of your changes to your local repository, you can push them all upstream to your
GitHub repository.
- Finally, open a pull request on GitHub to the main development repository (i.e., by now,
[this](https://github.com/smroccuzzo/UltraCold) one). 

## License 

UltraCold is distributed as free software under the GPL3. See the file [LICENSE.md](./LICENSE.md).