/** \file zernike_data.cpp
  Implementation of zernike_dat.hpp.
  It includes raw data from the root_data directory.
  \author J. Houdayer
*/
#include "zernike_data.hpp"

const int zri_n_roots[ZER_MAX_MOMENTS] = {
  #ifdef ZER_N50
  #include "root_data/zernike_n_roots_50.dat"
  #endif
  #ifdef ZER_N100
  #include "root_data/zernike_n_roots_100.dat"
  #endif
  #ifdef ZER_N150
  #include "root_data/zernike_n_roots_150.dat"
  #endif
  #ifdef ZER_N200
  #include "root_data/zernike_n_roots_200.dat"
  #endif
};

const double zri_roots[ZER_MAX_ROOTS] = {
  #ifdef ZER_N50
  #include "root_data/zernike_roots_50.dat"
  #endif
  #ifdef ZER_N100
  #include "root_data/zernike_roots_100.dat"
  #endif
  #ifdef ZER_N150
  #include "root_data/zernike_roots_150.dat"
  #endif
  #ifdef ZER_N200
  #include "root_data/zernike_roots_200.dat"
  #endif
};
