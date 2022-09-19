# zernike3d

## Content
zernike3D contains routines to compute 3D Zernike polynomials and 3D Zernike moments on point clouds or triangular mesh surfaces. It can also compute basic rotational invariants.

It also contains 3 standalone programs:

## Shape2Zernike
Shape2Zernike computes Zernike moments for triangular mesh surfaces given in OFF format.

## Zernike2Shape
Zernike2Shape takes a set of 3D Zernike moments as computed by Shape2Zernike and builds the corresponding shape as a triangular mesh surface in OFF format.

## MakeShape
MakeShape is a little tool to build shapes in OFF format.

## Installation
This is made for linux, it seems to work on MacOS. It may or may not work on other operating systems including Windows.
It is written in standard c++-11 and uses std::threads for parallelization.

1. Clone the repository.
2. Call `make` in the root directory. The Makefile assumes a g++ compiler change the first lines of Makefile to adapt. You need a c++11 compiler for this.
3. You can also use `make clean` to start anew.
4. You can use `make doc` to build the documentation (to be found in refman.pdf and index.html). You need to have doxygen installed for this. 
5. You can finally use `make bin` to build and copy the binaries to your ~/bin directory.
