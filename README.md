# README

## What it is

UltraCold is a modular and extensible collection of C++ libraries for the study of ultra-cold atomic systems in the 
context of Gross-Pitaevskii theory. 

## What UltraCold can do

UltraCold contains several solver classes for different flavors of Gross-Pitaevskii-like equations, allowing to 
describe ultra-cold systems of Bosons at the mean-field levels and to study both their ground-state configuration and 
the dynamics. UltraCold also offers several solver classes for Bogolyubov equations, allowing to study the elementary 
excitations on top of certain ground states. Most of the solver classes can take advantage of several parallelization 
strategies, including the use of OpenMP, MPI, and CUDA. The complete list of available solver classes, as well as 
detailed examples of usage, can be found in the [documentation](./docs/html/index.html).   

## Prerequisites and platforms

UltraCold has been tested only on Linux machines, including the High Performance Computing cluster 
[Galileo100](https://www.hpc.cineca.it/hardware/galileo100) from the italian supercomputing center 
[CINECA](https://www.cineca.it/). 

UltraCold is entirely based on the Math Kernel Library from 

## Installation



## Examples

Check the [examples](./examples) folder. Always source the setvarsh.sh script 

## License 

UltraCold is distributed as free software under the GPL3. See also the file LICENSE.md