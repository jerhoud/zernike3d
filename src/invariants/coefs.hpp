/** \file coefs.hpp
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
  const int N;
  unl(int n): N(n), u((N + 1) * (N + 2) / 2) {}
  const mpq_class &get(int n, int l) const
  { return u[n * (n + 1) / 2 + l]; }
protected:
  std::vector<mpq_class> u;
};

class unl0: public unl
{
public:
  unl0(int n, const binomials &);
};

class unl3: public unl
{
public:
  unl3(int n, const binomials &);
};

class coefs
{
public:
  const int N;
  coefs(int n): N(n), c((N + 1) * (N + 2) * (N + 3) / 6) {}
  const mpq_class &get(int l, int n, int k) const
  { return c[l * (l + 1) * (l + 2) / 6 + n * (n + 1) / 2 + k]; }
  std::vector<mpq_class> extract() const
  { return c; }
protected:
  std::vector<mpq_class> c;
};

class theta: public coefs
{
public:
  theta(int N, const factorials &, const double_factorials &, const binomials &);
};

class omega: public coefs
{
public:
  omega(const unl &u, const coefs &th);
};

class HK_coefs
{
public:
  HK_coefs(int n);
  std::vector<mpq_class> h_coefs, k0_coefs, k3_coefs;
};

#endif