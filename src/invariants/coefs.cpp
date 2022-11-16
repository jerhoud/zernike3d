/** \file coefs.hpp
  Classes to compute sets of large integers using gmp
  \author J. Houdayer
*/

#include <gmp.h>
#include "coefs.hpp"

mpz_class mul_2exp(const mpz_class &x, u_int n)
{
  mpz_class r;
  mpz_mul_2exp(r.get_mpz_t(), x.get_mpz_t(), n);
  return r;
}

factorials::factorials(int N):
facs(N + 1, 1)
{
  for (int n = 2 ; n <= N ; n++)
    facs[n] = n * facs[n - 1];
}

double_factorials::double_factorials(int N):
dfacs(N + 1, 1)
{
  for (int n = 1 ; n <= N ; n++)
    dfacs[n] = (2 * n + 1) * dfacs[n - 1];
}

binomials::binomials(int N):
bins((N / 2 + 1) * (N / 2 + 2), 1)
{
  for (int n = 2 ; n <= N ; n++)
    for (int l = 1 ; l <= n / 2 ; l++)
      bins[index(n, l)] = bins[(index(n - 1, l - 1))] + bins[(index(n - 1, l))];
}

void unl::make_d()
{
  for (size_t i = 0 ; i < u.size() ; i++)
    ud[i] = u[i].get_d();
}

std::vector<mpq_class> unl::apply(const std::vector<mpq_class> &v) const
{
  const int nmax = v.size() - 1;
  std::vector<mpq_class> r(nmax + 1);
  for (int n = 0, idx = 0 ; n <= nmax ; n++) {
    mpq_class sum = 0;
    for (int l = 0 ; l <= n ; l++, idx++)
      sum += u[idx] * v[l];
    r[n] = sum;
  }
  return r;
}

std::vector<double> unl::apply(const std::vector<double> &v) const
{
  const int nmax = v.size() - 1;
  std::vector<double> r(nmax + 1);
  for (int n = 0, idx = 0 ; n <= nmax ; n++) {
    double sum = 0;
    for (int l = 0 ; l <= n ; l++, idx++)
      sum += ud[idx] * v[l];
    r[n] = sum;
  }
  return r;
}

unl3::unl3(int n, const binomials &bins):
unl(n)
{
  for (int n = 0, idx = 0 ; n <= N ; n++)
    for (int l = 0 ; l <= n ; l++, idx++)
      u[idx] = ((l & 1) ? -1 : 1) * bins.get(n, l) * bins.get(n + l + 1, l + 1) *
               ((2 * n + 3) * (2 * l + 3) * (n + l + 2)) / 18_mpq;
  make_d();
}

vnl3::vnl3(int n, const factorials &f, const double_factorials &df, const binomials &b):
unl(n)
{
  for (int n = 0, idx = 0 ; n <= N ; n++)
    for (int l = 0 ; l <= n ; l++, idx++) {
      u[idx] = mpq_class(((l & 1) ? -9 : 9) * f.get(n) * b.get(2 * n + 3, n - l),
                          mul_2exp((2 * n + 3) * df.get(n + 1), n));
      u[idx].canonicalize();
    }
  make_d();
}


ucompose::ucompose(const unl &a, const unl &b):
unl(a.N)
{
  for (int n = 0, idx = 0 ; n <= N ; n++)
    for (int l = 0 ; l <= n ; l++, idx++) {
      mpq_class sum = 0;
      for (int i = l ; i <= n ; i++)
        sum += a.get(n, i) * b.get(i, l);
      u[idx] = sum;
    }
  make_d();
}

void coefs::make_d()
{
  for (size_t i = 0 ; i < c.size() ; i++)
    cd[i] = c[i].get_d();
}

std::vector<mpq_class> coefs::apply(const std::vector<mpq_class> &f) const
{
  std::vector<mpq_class> r(N + 1);
  for (int l = 0, idx = 0 ; l <= N ; l++) {
    const size_t sz = (l + 1) * (l + 2) / 2;
    if (sz > f.size()) {
      r.resize(l);
      break;
    }
    mpq_class sum = 0;
    for (size_t i = 0 ; i < sz ; i++, idx++)
      sum += c[idx] * f[i];
    r[l] = sum;
  }
  return r;
}

std::vector<double> coefs::apply(const std::vector<double> &f) const
{
  std::vector<double> r(N + 1);
  for (int l = 0, idx = 0 ; l <= N ; l++) {
    const size_t sz = (l + 1) * (l + 2) / 2;
    if (sz > f.size()) {
      r.resize(l);
      break;
    }
    double sum = 0;
    for (size_t i = 0 ; i < sz ; i++, idx++)
      sum += cd[idx] * f[i];
    r[l] = sum;
  }
  return r;
}

theta::theta(int N, const factorials &f, const double_factorials &df, const binomials &b):
coefs(N)
{
  for (int l = 0, idx = 0 ; l <= N ; l++) {
    for (int n = 0 ; n <= l ; n++)
      for (int k = 0 ; k <= n ; k++, idx++) {
        c[idx] = mpq_class(((n & 1) ? -4 : 4) * b.get(2 * l + 3, l - n) * f.get(2 * l + 1),
                             mul_2exp(df.get(l + k + 1) * df.get(l - k + 1), 2 * l));
        c[idx].canonicalize();
      }
  }
  make_d();
}

omega::omega(const unl &u, const coefs &th):
coefs(u.N)
{
  for (int m = 0, idx = 0 ; m <= N ; m++)
    for (int n = 0 ; n <= m ; n++)
      for (int k = 0 ; k <= n ; k++, idx++) {
        mpq_class sum = 0;
        for (int l = n ; l <= m ; l++)
          sum += u.get(m, l) * th.get(l, n, k);
        c[idx] = sum;
      }
  make_d();
}

inv_coefs::inv_coefs(int n):
N(n), facs(2 * n + 1), dfacs(2 * n + 1), bins(2 * n + 3),
u3(n, bins), v3(n, facs, dfacs, bins), o3()
{}