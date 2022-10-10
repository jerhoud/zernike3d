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
To install the programs you need the git program to clone the project (or some other mean to get a copy) and the cmake tool and a c++-11 compliant compiler to build it. If found it will use the c++ thread library for parallelization.

To install execute the following commands in the directory where you want to place your copy of the project

1. To get your local copy of the project
git clone https://github.com/jerhoud/zernike3d.git
2. To set the current directory in the project
cd zernike3d
3. To prepare the installation
cmake -S src -B build
4. To build the project
cmake --build build
5. To install the executables in directory myinstallpath (on linux you typically use $HOME or /usr/local as myinstallpath so the executables are in $HOME/bin or /usr/local/bin)
cmake --install build --prefix myinstallpath

If you need to restart the installation simply remove the build directory and restart from step 3
