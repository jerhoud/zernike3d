/** \file largeint.hpp
  Classes to compute sets of large integers using gmp
  \author J. Houdayer
*/

#ifndef LARGEINT_HPP
#define LARGEINT_HPP

#include <gmpxx.h>
#include <vector>

class factorials
{
public:
  factorials(int n);
  const mpz_class &get(int n) const
  { return facs[n]; }

private:
  std::vector<mpz_class> facs;
};

class double_factorials
{
public:
  double_factorials(int N);
  const mpz_class &get(int n) const
  { return dfacs[n]; }
private:
  std::vector<mpz_class> dfacs;
};

class binomials
{
public:
  binomials(int N);
  int index(int n, int l) const
  {
    if (n - l < l)
      l = n - l;
    int idx = (n / 2) * (n / 2 + 1);
    if (n & 1)
      idx += n / 2 + 1;
    return idx + l;
  }

  const mpz_class &get(int n, int l) const
  {
    return bins[index(n, l)];
  }

private:
  std::vector<mpz_class> bins;
  
};

inline mpz_class binomial(int n, int l)
{
  mpz_class r;
  mpz_bin_uiui(r.get_mpz_t(), n, l);
  return r;
}

class unl
{
public:
  unl(int N, const factorials &, const binomials &);
  const mpz_class &get(int n, int l) const
  { return us[n * (n + 1) / 2 + l]; }
private:
  std::vector<mpz_class> us;
};

class theta
{
public:
  theta(int N, const double_factorials &, const binomials &);
  const mpq_class &get(int l, int n, int k) const
  { return th[l * (l + 1) * (l + 2) / 6 + n * (n + 1) / 2 + k]; }
private:
  std::vector<mpq_class> th;
};

class omega
{
public:
  omega(int N);
  const mpq_class get(int m, int n, int l)
  { return omg[m * (m + 1) * (m + 2) / 6 + n * (n + 1) / 2 + l]; }
private:
  std::vector<mpq_class> omg;
};

#endif