/** \file neldermead.cpp
  Implementation of neldermead.hpp.
  \author J. Houdayer
*/

#include <utility>
#include "neldermead.hpp"

using std::swap;

class vecval
{
public:
  vec v;
  double val;
};

void sort4(vecval &v1, vecval &v2, vecval &v3, vecval &v4)
{
  if(v1.val>v2.val) swap(v1, v2);
  if(v3.val>v4.val) swap(v3, v4);
  if(v1.val>v3.val) swap(v1, v3);
  if(v2.val>v4.val) swap(v2, v4);
  if(v2.val>v3.val) swap(v2, v3);
}

vec minimize(std::function<double(const vec &)> f, const vec &start, double scale, double fthresh, double vthresh, int itermax)
{
  vthresh *= vthresh;
  
  vecval v1{start, f(start)};
  vecval v2{start + vec{scale, 0, 0}, f(start + vec{scale, 0, 0})};
  vecval v3{start + vec{0, scale, 0}, f(start + vec{0, scale, 0})};
  vecval v4{start + vec{0, 0, scale}, f(start + vec{0, 0, scale})};

  sort4(v1, v2, v3, v4);

  while (itermax-- != 0 && (v1.v - v4.v).length_square() > vthresh && fabs(v1.val - v4.val) > fthresh) {
    // new points
    vec b = (v1.v + v2.v + v3.v) / 3;
    vec d = b - v4.v;
    vec r = b + d;
    vecval vr{r, f(r)};
    if (vr.val < v1.val) { // expansion
      v4 = v3;
      v3 = v2;
      v2 = v1;
      vec e = r + d;
      vecval ve{e, f(e)};
      if (ve.val < v1.val)
        v1 = ve;
      else
        v1 = vr;
      continue;
    }
    if (vr.val < v3.val) { // reflection
      v4 = v3;
      if (vr.val < v2.val) {
        v3 = v2;
        v2 = vr;
      }
      else
        v3 = vr;
      continue;
    }
    d /= 2;
    vec c;
    if (vr.val < v4.val) { // outer contraction
      c = b + d;
      v4 = vr;
    }
    else // inner contraction
      c = b - d;
    vecval vc{c, f(c)};
    if (vc.val < v4.val) { // contraction
      if (vc.val >= v3.val)
        v4 = vc;
      else {
        v4 = v3;
        if (vc.val >= v2.val)
          v3 = vc;
        else {
          v3 = v2;
          if (vc.val >= v1.val)
            v2 = vc;
          else {
            v2 = v1;
            v1 = vc;
          }
        }
      }
      continue;
    }
    else { // shrink
      v2.v += (v2.v - v1.v) / 2;
      v2.val = f(v2.v);
      v3.v += (v3.v - v1.v) / 2;
      v3.val = f(v3.v);
      v4.v += (v4.v - v1.v) / 2;
      v4.val = f(v4.v);
    }
  }
  return v1.v;
}