/** \file invariants.hpp
  A class to compute spherical bessel functions
  \author J. Houdayer
*/

#ifndef SPHERICALB_HPP
#define SPHERICALB_HPP

#include <vector>

/** A class for computing spherical bessel functions up to a given order. */
class spherical_bessel
{
public:
  const int N;
  spherical_bessel(int n);
  void eval(double x);
  const std::vector<double> get_bsl() const
  { return bsl; }
private:
  std::vector<double> bsl;

  void asccending(double ix);
  void descending(double ix, int lmax);
};

#endif