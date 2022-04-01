/** \file zernike_def.hpp
  Defines macros to select the amount of data to include.
  Defining ZER_N50 (the default) includes data up to N = 50.
  Likewise for ZER_N100, ZER_N150 and ZER_N200. The size of the data is cubic with N.
  Size of zernike_data.o respectively: 56K, 378K, 1.2M, 2.8M
  Size of invariants_data.o respectively: 27K, 185K, 596K, 1.4M
  \author J. Houdayer
*/
#ifndef ZERNIKE_DEF_HPP
#define ZERNIKE_DEF_HPP

#ifdef ZER_N200
#define ZER_N150
#define ZER_MAX_N 201
#endif

#ifdef ZER_N150
#define ZER_N100
#ifndef ZER_N200
#define ZER_MAX_N 151
#endif
#endif

#ifdef ZER_N100
#define ZER_N50
#ifndef ZER_N150
#define ZER_MAX_N 101
#endif
#endif

#ifndef ZER_N50
#define ZER_N50
#endif

#ifndef ZER_N100
#define ZER_MAX_N 51
#endif

#define ZER_MAX_MOMENTS ((ZER_MAX_N / 2 + 1) * (ZER_MAX_N / 2 + 2))
#define ZER_MAX_ROOTS ((ZER_MAX_N / 2 + 1) * (ZER_MAX_N / 2 + 2) * (ZER_MAX_N / 2 + 3) / 3)

#endif
