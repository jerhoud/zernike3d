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
};

double brent(std::function<double(double)> f, double start, double scale, double thresh)
{ 
  xval vb(start, f), v1(start + scale, f);
  
  if (vb.fx > v1.fx)
    std::swap(vb, v1);

  xval v2(vb.x + 2 * (vb.x - v1.x), f);

  // descending
  while (v2.fx < vb.fx) {
    double x = v2.x + 2 * (v2.x - vb.x);
    v1 = vb;
    vb = v2;
    v2 = xval(x, f);
  }

  if (v1.x > v2.x)
    std::swap(v1, v2);

  int m = 2;

  while (v2.x - v1.x > 2 * thresh) {
    double d;
    bool split = true;
    const double d1 = v1.x - vb.x;
    const double d2 = v2.x - vb.x;
    if (m-- != 0) {
      // compute quadratic interpolation
      const double df1 = d2 * (v1.fx - vb.fx);
      const double df2 = d1 * (v2.fx - vb.fx);
      d = (d1 * df2 - d2 * df1) / (2 * (df2 - df1)); // quadratic interpolation
      split = !(d >= d1 && d <= d2);
    }
    if (split) { // perfect ratio step
      if (d2 + d1 > 0)
        d = d1 + sqrt(d1 * (d1 - d2));
      else
        d = d2 - sqrt(d2 * (d2 - d1));
    }
    
    if (fabs(d) < thresh)
      d = 0.5 * ((d < 0) ? -thresh : thresh);
    if (d - d1 < thresh)
      d = d1 + thresh;
    else if (d2 - d < thresh)
      d = d2 - thresh;

    xval vn(vb.x + d, f);

    if (vn.fx < vb.fx) { //new point is better, put it in the middle
      if (vn.x < vb.x)
        v2 = vb;
      else
        v1 = vb;
      vb = vn;
      if (split)
        m = 0;
    }
    else { // new point is not better, put it on the side
      if (vn.x < vb.x)
        v1 = vn;
      else
        v2 = vn;
      if (split)
        m = 2;
      else if (m == 1)
        m = 0;
    }
  }
  return vb.x;
}
