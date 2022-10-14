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

unl0::unl0(int N, const binomials &bins):
unl(N)
{
  for (int n = 0, idx = 0 ; n <= N ; n++)
    for (int l = 0 ; l <= n ; l++, idx++)
      u[idx] = ((l & 1) ? -1 : 1) * bins.get(n, l) * bins.get(n + l, l);
}

unl3::unl3(int n, const binomials &bins):
unl(n)
{
  for (int n = 0, idx = 0 ; n <= N ; n++)
    for (int l = 0 ; l <= n ; l++, idx++)
      u[idx] = ((l & 1) ? -1 : 1) * bins.get(n, l) * bins.get(n + l + 1, l + 1) *
               ((2 * n + 3) * (2 * l + 3) * (n + l + 2)) / 18_mpq;
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
}

HK_coefs::HK_coefs(int n)
{
  factorials f(2 * n + 1);
  double_factorials df(2 * n + 2);
  binomials b(2 * n + 3);
  unl0 u0(n, b);
  unl3 u3(n, b);
  theta t(n, f, df, b);
  omega o0(u0, t);
  omega o3(u3, t);

  h_coefs = t.extract();
  k0_coefs = o0.extract();
  k3_coefs = o3.extract();
}

