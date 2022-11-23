/** \file k_distance.cpp
  Implementation of k_distance.hpp.
  \author J. Houdayer
*/

#include "k_distance.hpp"

double delta(int n0, const inv_k3 &k1, const inv_k3 &k2, double alpha)
{
  const std::vector<double> &c1 = k1.get_d();
  const std::vector<double> c2 = k2.resized(n0, alpha);
  double sum = 0;
  for (int n = 1 ; n <= n0 ; n++) {
    const double diff = c1[n] - c2[n];
    sum += diff * diff / (2 * n + 3);
  }
  return sum;
}

double dep_distance(int n, const inv_k3 &k1, const inv_k3 &k2)
{

}

double simple_undep_distance(int n, const inv_k3 &k1, const inv_k3 &k2);
double min_undep_distance(int n, const inv_k3 &k1, const inv_k3 &k2);
