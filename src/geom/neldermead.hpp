/** \file neldermead.hpp
  Minimization without derivative in 3D space.
  \author J. Houdayer
*/

#ifndef NELDERMEAD_HPP
#define NELDERMEAD_HPP

#include <functional>
#include "vec.hpp"

vec minimize(std::function<double(const vec &)> f, const vec &start, double start_scale, double thresh, int itermax);

#endif