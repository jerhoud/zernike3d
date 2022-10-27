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
  vecval(const vec &vc, std::function<double(const vec &)> f):
  v(vc), val(f(vc)) {}
  void set(const vec &vc, std::function<double(const vec &)> f)
  { v = vc; val = f(vc); }
};

void sort4(vecval *&v1, vecval *&v2, vecval *&v3, vecval *&v4)
{
  if(v1->val > v2->val) swap(v1, v2);
  if(v3->val > v4->val) swap(v3, v4);
  if(v1->val > v3->val) swap(v1, v3);
  if(v2->val > v4->val) swap(v2, v4);
  if(v2->val > v3->val) swap(v2, v3);
}

int insert4(vecval *&v1, vecval *&v2, vecval *&v3, vecval *&v4, const vecval &nv)
{
  const double v = nv.val;
  if (v >= v4->val)
    return 5;
  *v4 = nv;
  if (v >= v3->val)
    return 4;
  vecval *tmp = v4;
  v4 = v3;
  if (v >= v2->val) {
    v3 = tmp;
    return 3;
  }
  v3 = v2;
  if (v >= v1->val) {
    v2 = tmp;
    return 2;
  }
  v2 = v1;
  v1 = tmp;
  return 1;
}

vec minimize(std::function<double(const vec &)> f, const vec &start, double scale, double fthresh, double vthresh, int itermax)
{
  vthresh *= vthresh;
  
  vecval va(start, f);
  vecval vb(start + vec{scale, 0, 0}, f);
  vecval vc(start + vec{0, scale, 0}, f);
  vecval vd(start + vec{0, 0, scale}, f);

  vecval *v1 = &va, *v2 = &vb, *v3 = &vc, *v4 = &vd;
  sort4(v1, v2, v3, v4);

  while (itermax-- != 0 && (v1->v - v4->v).length_square() > vthresh && fabs(v1->val - v4->val) > fthresh) {
    vec b = (v1->v + v2->v + v3->v) / 3;
    vec d = b - v4->v;

    const int p = insert4(v1, v2, v3, v4, vecval(b + d, f)); // reflection
    if (p == 1) { // expansion
      vecval ve(v1->v + d, f);
      if (ve.val < v1->val)
        *v1 = ve;
    }
    else if (p >= 4) { // contraction
      d /= 2;
      if (p == 4)
        b += d;
      else
        b -= d;
      const int p2 = insert4(v1, v2, v3, v4, vecval(b, f));
      if (p2 == 5) { // shrink
        v2->set((v2->v + v1->v) / 2, f);
        v3->set((v3->v + v1->v) / 2, f);
        v4->set((v4->v + v1->v) / 2, f);
        sort4(v1, v2, v3, v4);
      }
    }
  }
  return v1->v;
}