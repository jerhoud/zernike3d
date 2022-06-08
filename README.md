# zernike3d

## Content
zernike3D contains routines to compute 3D Zernike polynomials and 3D Zernike moments on point clouds or triangular mesh surfaces. It can also compute basic rotational invariants.

It also contains 3 standalone programs:

## zm
zm computes Zernike moments for triangular mesh surfaces given in OFF format.

## rzm
rzm takes a set of 3D Zernike moments as computed by zm and builds the corresponding shape as a triangular mesh surface in OFF format.

## makeOFF
makeOFF is a little tool to build shapes in OFF format.

## Installation
This is made for linux, I have not the knowledge for other operating systems.

1. Clone the repository.
2. Call `make` in the root directory. The Makefile assumes a g++ compiler change the first lines of Makefile to adapt. You need a c++11 compiler for this.
3. You can also use `make clean` to start anew.
4. You can use `make doc` to build the documentation (to be found in refman.pdf and index.html). You need to have doxygen installed for this. 
5. You can finally use `make bin` to build and copy the binaries to your ~/bin directory.
