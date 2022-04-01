/** \file zernike_data.hpp
  Defines global arrays containing data to compute
  integrated Zernike radial part.
  \author J. Houdayer
*/
#ifndef ZERNIKE_DATA_HPP
#define ZERNIKE_DATA_HPP

#include "zernike_def.hpp"

/** Number of real roots of the integrated Zernike radial part. */
extern const int zri_n_roots[ZER_MAX_MOMENTS];
/** Roots of the integrated Zernike radial part. */
extern const double zri_roots[ZER_MAX_ROOTS];

#endif
