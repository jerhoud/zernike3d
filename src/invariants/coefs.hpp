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
  unl(int n): N(n), u((N + 1) * (N + 2) / 2), ud(u.size()) {}
  const mpq_class &get(int n, int l) const
  { return u[n * (n + 1) / 2 + l]; }
  const std::vector<mpq_class> &get_mpq() const
  { return u; }
  const std::vector<double> &get_d() const
  { return ud; }
  std::vector<mpq_class> apply(const std::vector<mpq_class> &v) const;
  std::vector<double> apply(const std::vector<double> &v) const;
protected:
  std::vector<mpq_class> u;
  std::vector<double> ud;
  void make_d();
};

// unl0 needs binomials up to 2n
class unl0: public unl
{
public:
  unl0(int n, const binomials &);
};

// unl3 needs binomials up to 2n+1
class unl3: public unl
{
public:
  unl3(int n, const binomials &);
};

// vnl0 needs factorials up to 2n+1
class vnl0: public unl
{
public:
  vnl0(int n, const factorials &f);
};

// vnl3 needs factorials up to n
//            double_fac up to n+1
//            binomials  up to 2n+3
class vnl3: public unl
{
public:
  vnl3(int n, const factorials &f, const double_factorials &df, const binomials &b);
};

class coefs
{
public:
  const int N;
  coefs(int n): N(n), c((N + 1) * (N + 2) * (N + 3) / 6), cd(c.size()) {}
  const mpq_class &get(int l, int n, int k) const
  { return c[l * (l + 1) * (l + 2) / 6 + n * (n + 1) / 2 + k]; }
  const std::vector<mpq_class> &get_mpq() const
  { return c; }
  const std::vector<double> &get_d() const
  { return cd; }
  std::vector<mpq_class> apply(const std::vector<mpq_class> &f) const;
  std::vector<double> apply(const std::vector<double> &f) const;
protected:
  std::vector<mpq_class> c;
  std::vector<double> cd;
  void make_d();
};

class ucompose: public unl
{
public:
  ucompose(const unl &u, const unl &v);
};

// theta needs factorials up to 2n+1
//             double_fac up to 2n+1
//             binomials  up to 2n+3
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

#endif