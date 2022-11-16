/** \file s_root_data.cpp
  imlementation of s_root_data.hpp
  \author J. Houdayer
*/

#include "s_root_data.hpp"

extern const double s_root_data[S_SIZE] = {
#ifdef DEF_S_50
#include "s_root_data_50.cpp"
#endif
#ifdef DEF_S_100
#include "s_root_data_100.cpp"
#endif
};