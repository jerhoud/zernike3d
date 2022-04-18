/** \file gauss.hpp
  Classes for integrating on a segment.
  \author J. Houdayer
*/

#include <vector>
#include <iostream>


class gauss_point
{
public:
    const double x, weight;
};

class gauss_scheme
{
public:
  const int order; /**< order of integration. */

  gauss_scheme(int o, const std::vector<gauss_point> &d):
  order(o), data(d) {}

  double check_weights() const;
  double check_poly() const;
  int check_sorted() const;

  /** Integrates over a triangle.
    @param s The integration scheme.
    @param v The integrator object. It should have an "add" member that takes a w_vec.
    @param w An overall weight.
  */
  template<class T>
  void integrate(T &v, double a, double b, double w) const
  {
    const double s = b - a;
    for (auto &i: data)
      v.add(a + s * i.x, w * i.weight);
  }

  const std::vector<gauss_point> data; /**< The integration points. */
};

std::ostream &operator <<(std::ostream &os, const gauss_scheme &s);

class gauss_selector
{
public:
  gauss_selector();
  const gauss_scheme &get_scheme(int n) const;
  int max_order() const;

  std::vector<gauss_scheme> schemes;
};

