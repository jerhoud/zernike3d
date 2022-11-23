/** \file k_distance.hpp
  To compute distances between invariant sets.
  \author J. Houdayer
*/

#ifndef K_DISTANCE_HPP
#define K_DISTANCE_HPP

#include "invariants.hpp"

double dep_distance(int n, const inv_k3 &k1, const inv_k3 &k2);
double simple_undep_distance(int n, const inv_k3 &k1, const inv_k3 &k2);
double min_undep_distance(int n, const inv_k3 &k1, const inv_k3 &k2);

#endif