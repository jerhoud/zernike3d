/** \file triangle.hpp
  Classes for manipulating triangles and integration.
  \author J. Houdayer
*/
#ifndef TRIANGLE_HPP
#define TRIANGLE_HPP
#include "vec.hpp"

/** A class representing a triangle. */
class triangle
{
public:
  vec p1, p2, p3;

  /** Area of the triangle. */
  double area() const
  { return cross(p2 - p1, p3 - p1).length() / 2; }

  /** Volume of the tetrahedron defined by the origin and the triangle. */
  double volume() const
  { return dot(p1, cross(p2, p3)) / 6; }
};

/** A point in an integration scheme. */
class triquad_point
{
public:
  const double weight;
  const double c1, c2, c3;

  /** Computes the weighted point corresponding to the given triangle and overall weight. */
  w_vec point(const triangle &t, double w) const
  { return {weight * w, c1 * t.p1 + c2 * t.p2 + c3 * t.p3}; }
};

/** An integration scheme over the triangle. */
class triquad_scheme
{
public:
  const int order; /**< order of integration. */

  triquad_scheme(int o, const std::vector<triquad_point> &d):
  order(o), data(d) {}

  double check_unity() const;
  double check_weights() const;
  double check_poly(const int n1, const int n2) const;
  double check_poly_auto() const;
  int check_sorted() const;

  /** Integrates over a triangle.
    @param s The integration scheme.
    @param v The integrator object. It should have an "add" member that takes a w_vec.
    @param w An overall weight.
  */
  template<class T>
  void integrate(const triangle &t, T &v, double w, int n=0) const
  {
    if (n <= 0)
      for (auto &i: data)
        v.add(i.point(t, w));
    else {
      n--;
      w /= 4;
      const vec p12 = (t.p1 + t.p2) / 2;
      const vec p23 = (t.p2 + t.p3) / 2;
      const vec p31 = (t.p3 + t.p1) / 2;
      integrate({t.p1, p12, p31}, v, w, n);
      integrate({t.p2, p23, p12}, v, w, n);
      integrate({t.p3, p31, p23}, v, w, n);
      integrate({p12, p23, p31}, v, w, n);
    }
  }

  const std::vector<triquad_point> data; /**< The integration points. */
};

std::ostream &operator <<(std::ostream &os, const triquad_scheme &s);

class triquad_selector
{
public:
  triquad_selector();
  const triquad_scheme &get_scheme(int n) const;
  int max_order() const;

  std::vector<triquad_scheme> schemes;
};

#endif
