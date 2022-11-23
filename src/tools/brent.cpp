/** \file brent.cpp
  Implementation of brent.hpp.
  \author J. Houdayer
*/

#include <iostream>
#include <iomanip>
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


void descending(xval *&v1, xval *&v2, xval *&v3, std::function<double(double)> f)
{
  std::cout << "descending: {" << v1->x << ", " << v2->x << ", " << v3->x << "} -> " << v2->fx << "\n";
  double x = 3 * v3->x - 2 * v2->x;
  v1->set(x, f);
  xval *tmp = v1;
  v1 = v2;
  v2 = v3;
  v3 = tmp;
}

bool insert(xval *&v1, xval *&v2, xval *&v3, double d, std::function<double(double)> f)
{
  xval vc(v2->x + d, f);

  if (vc.fx < v2->fx) { //new point is better, put it in the middle
    std::cout << "better\n";
    if (vc.x < v2->x) {
      *v3 = vc;
      std::swap(v2, v3);
    }
    else {
      *v1 = vc;
      std::swap(v2, v1);
    }
    return true;
  }
  else { // new point is not better, put it on the side
    if (vc.x < v2->x)
      *v1 = vc;
    else
      *v3 = vc;
    return false;
  }
}

double minimize(std::function<double(double)> f, double start, double scale, double thresh)
{
  const double g = 0.381966011250105; // (3 - sqrt(5)) / 2
  
  std::cout << std::setprecision(12);

  xval va(start, f), vb(start + scale, f), vc(vb);
  xval *v1 = NULL, *v2 = NULL, *v3 = &vc;
  
  if (va.fx > vb.fx) {
    v1 = &va;
    v2 = &vb;
    v3->set(start + 3 * scale, f);
  }
  else {
    v1 = &vb;
    v2 = &va;
    v3->set(start - 2 * scale, f);
  }
  // descending
  while (v3->fx < v2->fx)
    descending(v1, v2, v3, f);

  if (v1->x > v3->x)
    std::swap(v1, v3);

  int m = 3;
  double e1 = v3->x - v1->x, e0 = e1;

  while (e1 > thresh) {
    std::cout << "reducing: {" << v1->x << ", " << v2->x << ", " << v3->x << "} -> " << v2->fx << "\n";
    // compute quadratic interpolation
    const double d1 = v1->x - v2->x;
    double df1 = v1->fx - v2->fx;
    const double d3 = v3->x - v2->x;
    double df3 = v3->fx - v2->fx;
    df1 *= d3;
    df3 *= d1;
    double d = (d1 * df3 - d3 * df1) / (2 * (df3 - df1)); // quadratic interpolation
    double emax = e1 / 2;
    e1 = e0;
    e0 = fabs(d);
    bool inter = m-- != 0 && e0 < emax && d >= d1 && d <= d3;
    if (!inter) { // golden ratio step
      if (-d1 > d3) {
        d = g * d1;
        e0 = -d;
        std::cout << "right ratio:\n";
      }
      else {
        d = g * d3;
        e0 = d;
        std::cout << "left ratio\n";
      }
    }
    else
      std::cout << "interpolation\n";
    if (fabs(d) < thresh)
      d = (d < 0) ? -thresh : thresh;

    bool better = insert(v1, v2, v3, d, f);
    if (!inter) {
      if (better)
        m = 0;
      else
        m = 3;
    }
  }
  return v2->x;
}
