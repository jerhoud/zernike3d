/** \file s_root_data.hpp
  This files defines the Snl root data.
  \author J. Houdayer
*/

#ifndef S_ROOT_DATA_HPP
#define S_ROOT_DATA_HPP

#ifndef DEF_S_50
#define DEF_S_50
#endif

#ifdef DEF_S_100
#define S_MAX 100
#else
#define S_MAX 50
#endif

#define S_SIZE ((S_MAX + 1) * (S_MAX + 2) * (S_MAX + 3) / 6)

extern const double s_root_data[S_SIZE];

#endif
