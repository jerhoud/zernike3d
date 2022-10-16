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
  const double ix = 1 / x;
  if (N <= 0)
    return;
  
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
  if (N <= 1)
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
  
  int la = (N > x + 1) ? x + 1 : N;
  if (la == N - 1)
    la = N;
  asccending(2, la, ix);
  if (N == la)
    return;
  
  std::vector<double> b1(N + 1), b2(N + 1);
  descending(b1, b2, N, la, ix);
}

void spherical_bessel::asccending(int lmin, int lmax, double ix)
{
  double jpp = bsl[lmin - 2];
  double jp = bsl[lmin - 1];
  for (int l = lmin ; l <= lmax ; l++) {
    double j = -jpp + (2 * l - 1) * ix * jp;
    jpp = jp;
    jp = j;
    bsl[l] = j;
  }
}

void spherical_bessel::descending(std::vector<double> &bsl1, std::vector<double> &bsl2, int lmax, int lmin, double ix)
{
  if (bsl[lmin - 1] == 0 && bsl[lmin] == 0) {
    for (int l = lmin + 1 ; l <= lmax ; l++)
      bsl[l] = 0;
    return;
  }

  const double vs = 1e-300;
  double jpp1 = bsl1[lmax] = 0;
  double jp1 = bsl1[lmax - 1] = vs;
  double jpp2 = bsl2[lmax] = vs;
  double jp2 = bsl2[lmax - 1] = 0;

  for (int l = lmax - 2 ; l >= lmin - 1 ; l--) {
    const double cf = (2 * l + 3) * ix;
    double j1 = -jpp1 + cf * jp1;
    jpp1 = jp1;
    jp1 = j1;
    bsl1[l] = j1;
    double j2 = -jpp2 + cf * jp2;
    jpp2 = jp2;
    jp2 = j2;
    bsl2[l] = j2;
  }
  if (std::isinf(bsl1[lmin - 1]) || std::isinf(bsl2[lmin - 1])) {
    int lmid = (lmin + lmax) / 2;
    descending(bsl1, bsl2, lmid, lmin, ix);
    descending(bsl1, bsl2, lmax, lmid, ix);
    return;
  }

  const double a = bsl[lmin - 1];
  const double b = bsl[lmin];
  double a1 = bsl1[lmin - 1];
  double b1 = bsl1[lmin];
  double a2 = bsl2[lmin - 1];
  double b2 = bsl2[lmin];

  // normalization to prevent overflow / underflow
  const double sz = 1. / std::max(std::max(fabs(a1), fabs(a2)), std::max(fabs(b1), fabs(b2)));
  a1 *= sz;
  a2 *= sz;
  b1 *= sz;
  b2 *= sz;

  const double det = (a1 * b2 - a2 * b1);
  const double c1 = (a * b2 - b * a2) / det * sz;
  const double c2 = (a * b1 - b * a1) / det * sz;

  for (int l = lmin + 1 ; l <= lmax ; l++)
    bsl[l] = c1 * bsl1[l] + c2 * bsl2[l];
}

