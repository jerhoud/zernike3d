/** \file invariants.hpp
  A class to compute spherical bessel functions
  \author J. Houdayer
*/

#ifndef SPHERICALB_HPP
#define SPHERICALB_HPP

#include <vector>

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

  void asccending(int lmin, int lmax, double ix);
  void descending(std::vector<double> &b1, std::vector<double> &b2, int lmax, int lmin, double ix);
};

#endif