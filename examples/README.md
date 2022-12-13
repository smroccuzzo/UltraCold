This folder contains several folders called *example-n*, each containing

- a source file called *example-n.cpp*, containing an example on how to write an executable that uses UltraCold
  libraries,
- a *CMakeLists.txt* containing instructions on how to configure and build an executable based on
  UltraCold,
- an eventual parameter file, called *example-n.prm*,

To run the examples, follow the usual steps required to build a project using *cmake*, namely, open a
terminal in the folder containing the example you are interested in and type

    mkdir build
    cd build
    cmake -DULTRACOLD_DIR=/path/to/the/directory/where/you/installed/UltraCold \
          -DCMAKE_BUILD_TYPE=Release .. \
          -DCMAKE_CXX_COMPILER=icpc
    make

This will create an executable called *example-n*, which (if everything went fine) should be ready to be
executed. Notice that we here pick explicitly the Intel C++ compiler. Sometimes, in some system, the gnu 
compiler would be picked by default, but then the program wouldn't work. 
After successful compilation, you can just copy the eventual file *example-n.prm* from the parent folder 
and run the example

    cp ../example-n.prm .
    ./example-n

The output of course will depend on the particular example, and is documented fully for each of them.

Although a prior knowledge of C++ is highly recommended, to use these examples also very basic knowledge is more
than sufficient. The documentation tries to be as pedantic as possible, so that extending these examples for user's
need shouldn't be too difficult.

Here is the complete list of all examples and a brief description of what each of them does. Refer to the detailed
description available in the documentation.

- [example-1](./example-1) Defined in file [example-1.cpp](./example-1/example-1.cpp) Ground state
  and simple dynamics of a three-dimensional, harmonically trapped Bose gas.
- [example-2](./example-2) Defined in file [example-2.cpp](./example-2/example-2.cpp) Elementary
  excitations via Bogolyubov equations in a two-dimensional, harmonically trapped Bose gas.
- [example-3](./example-3) Defined in file [example-3.cpp](./example-3/example-3.cpp) Supersolid
  ground state of a trapped dipolar Bose-Einstein condensate.
- [example-4](./example-4) Defined in file [example-4.cpp](./example-4/example-4.cpp) Excitation
  spectrum of a trapped dipolar Bose-Einstein condensate across the superfluid-supersolid phase transition.
- [example-5](./example-5) Defined in file [example-5.cpp](./example-5/example-5.cpp) Simple dynamics of two vortices in 
  a two-dimensional dipolar Bose gas, **using GPU acceleration**. 
- [example-6](./example-6) Defined in file [example-6.cpp](./example-6/example-6.cpp) Far from equilibrium dynamics in a
  two-dimensional dipolar Bose gas, **using GPU acceleration**. 
- [example-7](./example-7) Defined in file [example-7.cpp](./example-7/example-7.cpp) Far from equilibrium dynamics in a
  three-dimensional dipolar Bose gas, **using GPU acceleration**. 
  
