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
To install the programs you need the cmake tool (see cmake.org) and a c++-11 compliant compiler to build it. If found it will use the c++ thread library for parallelization.

To install, first dowload the project and go to its root directory then execute the following commands:

1. To prepare the installation

    cmake -S src -B build
2. To build the project

    cmake --build build
3. To install the executables in directory myinstallpath (on linux or macOS you typically use $HOME or /usr/local as myinstallpath so the executables are in $HOME/bin or /usr/local/bin)

    cmake --install build --prefix myinstallpath

If you need to restart the installation simply remove the build directory and go to step 1
