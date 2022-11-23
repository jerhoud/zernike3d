/** \file brent.cpp
  Implementation of brent.hpp.
  \author J. Houdayer
*/

#include <cmath>
#include <utility>
#include "brent.hpp"

class xval {
public:
  double x, fx;
  xval(double v, std::function<double(double)> f): x(v), fx(f(v)) {}
  void set(double v, std::function<double(double)> f)
  { x = v; fx = f(v); }
};

constexpr double g = 0.381966011250105; // (3 - sqrt(5)) / 2

double minimize(std::function<double(double)> f, double start, double scale, double thresh)
{
  xval va(start, f), vb(start + scale, f);
  xval *v1 = NULL, *v2 = NULL;
  
  double d;
  if (va.fx > vb.fx) {
    v1 = &va;
    v2 = &vb;
    d = 2 * scale;
  }
  else {
    v1 = &vb;
    v2 = &va;
    d = -scale;
  }
  // descending
  xval vc(start + d, f);
  xval *v3 = &vc;
  while (v3->fx < v2->fx && fabs(v1->x - v3->x) > thresh) {
    d *= 2;
    v1->set(v2->x + d, f);
    xval *tmp = v1;
    v1 = v2;
    v2 = v3;
    v3 = tmp;
  }

  if (v1->x > v3->x)
    std::swap(v1, v3);

  double e1 = v3->x - v1->x;
  double e0 = e1;

  while (e1 > thresh) {
    // compute quadratic interpolation
    const double d1 = v1->x - v2->x;
    double df1 = v1->fx - v2->fx;
    const double d3 = v3->x - v2->x;
    double df3 = v3->fx - v2->fx;
    df1 *= d3;
    df3 *= d1;
    const double p = d1 * df3 - d3 * df1;
    const double q = 2 * (df3 - df1);
    double d = p / q;; // quadratic interpolation step
    const double e2 = e1;
    e1 = e0;
    e0 = fabs(d);
    if (!(e0 < e2 / 2 && d >= d1 && d <= d3)) { // change too large, golden ratio step instead
      if (-d1 > d3) {
        d = g * d1;
        e1 = -d1;
        e0 = -d; 
      }
      else {
        d = g * d3;
        e1 = d3;
        e0 = d;
      }
    }
    xval vc(v2->x + d, f);
    if (vc.fx <= v2->fx) { //new point is better, put it in the middle
      if (vc.x < v2->x) {
        *v3 = vc;
        std::swap(v2, v3);
      }
      else {
        *v1 = vc;
        std::swap(v2, v1);
      }
    }
    else { // new point is not better, put it on the side
      if (vc.x < v2->x)
        *v1 = vc;
      else
        *v3 = vc;
    }
  }
  return v2->x;
}
