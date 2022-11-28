/** \file brent.hpp
  Minimization in 1D without derivative.
  \author J. Houdayer
*/

#ifndef BRENT_HPP
#define BRENT_HPP

#include <functional>

double brent(std::function<double(double)> f, double start, double start_scale, double thresh);

#endif