/** \file sphericalbessel.cpp
  imlementation of sphericalbessel.hpp
  \author J. Houdayer
*/

#include <cmath>
#include "sphericalbessel.hpp"

spherical_bessel::spherical_bessel(int n):
N(n), bsl(N + 1)
{}

void spherical_bessel::eval(double x)
{
  if (x == 0) {
    for (auto &b: bsl)
      b = 0;
    bsl[0] = 1;
    return;
  }

  const double ax = fabs(x);
  const double x2 = x * x;
  
  // computing j0 = sin x / x;
  if (ax < 0.5) {
    const double a1 = 1. / 6;
    const double a2 = 1. / 120;
    const double a3 = 1. / 5040;
    const double a4 = 1. / 362880;
    const double a5 = 1. / 39916800;
    const double a6 = 1. / 6227020800;
    bsl[0] = 1. - x2 * (a1 - x2 * (a2 - x2 * (a3 - x2 * (a4 - x2 * (a5 - x2 * a6)))));
  }
  else
    bsl[0] = sin(x) / x;

  if (N == 0)
    return;

  // computing j1 = sin x / x^2 - cos x / x
  if (ax < 0.25) {
    const double a1 = 1. / 10;
    const double a2 = 1. / 280;
    const double a3 = 1. / 15120;
    const double a4 = 1. / 1330560;
    const double a5 = 1. / 172972800;
    bsl[1] = x * (1 - x2 * (a1 - x2 * (a2 - x2 * (a3 - x2 * (a4 - x2 * a5))))) / 3;
  }
  else
    bsl[1] = (sin(x) / x - cos(x)) / x;
  
  if (N == 1)
    return;
  
  const double ix = 1 / x;
  if (x >= N)
    asccending(ix);
  else
    descending(ix, N + 30);
}

void spherical_bessel::asccending(double ix)
{
  double jpp = bsl[0];
  double jp = bsl[1];
  for (int l = 2 ; l <= N ; l++) {
    double j = -jpp + (2 * l - 1) * ix * jp;
    jpp = jp;
    jp = j;
    bsl[l] = j;
  }
}

void spherical_bessel::descending(double ix, int lmax)
{
  const double j0 = bsl[0];
  const double j1 = bsl[1];

  for(int l = 2 ; l <= N ; l++)
    bsl[l] = 0;

  if (lmax <= 1)
    return;

  const double vs = 1e-290;
  for (int l = lmax + 1 ; l <= N ; l++)
    bsl[l] = 0;

  double jpp = 0;
  double jp = vs;
  int l = lmax;
  for (; l > N ; l--) {
    double j = -jpp + (2 * l + 3) * ix * jp;
    jpp = jp;
    jp = j;
  }

  for(; l >= 0 ; l--) {
    double j = -jpp + (2 * l + 3) * ix * jp;
    jpp = jp;
    jp = j;
    bsl[l] = j;
  }

  if (!std::isfinite(bsl[0])) {
    int i = 1;
    for (; i <= N ; i++)
      if (std::isfinite(bsl[i]))
        break;
    bsl[0] = j0;
    bsl[1] = j1;
    return descending(ix, lmax - i - 5);
  }
  
  double c, err;
  if (fabs(j0) >= fabs(j1)) {
    c = j0 / bsl[0];
    err = fabs(j1 - c * bsl[1]);
  }
  else {
    c = j1 / bsl[1];
    err = fabs(j0 - c * bsl[0]);
  }

  bsl[0] = j0;
  bsl[1] = j1;

  if (err > 1e-14)
    return descending(ix, lmax + 50);
  
  for (int l = 2 ; l<= N ; l++)
    bsl[l] *= c;
}
