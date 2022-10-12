/** \file largeint.hpp
  Classes to compute sets of large integers using gmp
  \author J. Houdayer
*/

#include "largeint.hpp"

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

unl::unl(int N, const factorials & facs, const binomials &bins):
us((N + 1) * (N + 2) / 2)
{
  for (int n = 0, idx = 0 ; n <= N ; n++)
    for (int l = 0 ; l <= n ; l++, idx++)
      us[idx] = ((l & 1) ? -4 : 4) * facs.get(2 * l + 1) * bins.get(n, l) * bins.get(n + l, l);
}

theta::theta(int N, const double_factorials &facs, const binomials &bins):
th((N + 1) * (N + 2) * (N + 3) / 6)
{
  for (int l = 0, idx = 0 ; l <= N ; l++)
    for (int n = 0 ; n <= l ; n++)
      for (int k = 0 ; k <= n ; k++, idx++) {
        th[idx] = mpq_class(((n & 1) ? -1 : 1) * bins.get(2 * l + 1, l - n),
                             facs.get(l + k) * facs.get(l - k));
        th[idx].canonicalize();
      }
}

omega::omega(int N):
omg((N + 1) * (N + 2) * (N + 3) / 6)
{
  factorials f(2 * N + 1);
  double_factorials df(2 * (N + 1));
  binomials b(2 * N + 3);
  unl u(N, f, b);
  theta t(N + 1, df, b);

  for (int m = 0, idx = 0 ; m <= N ; m++)
    for (int n = 0 ; n <= m ; n++)
      for (int k = 0 ; k <= n ; k++, idx++) {
        mpq_class sum = 0;
        for (int l = n ; l <= m ; l++)
          sum -= u.get(m, l) * t.get(l + 1, n + 1, k);
        omg[idx] = sum;
      }
}
