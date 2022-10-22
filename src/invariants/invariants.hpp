/** \file invariants.hpp
  A class to compute rotational and translational invariants
  \author J. Houdayer
*/


#ifndef INVARIANTS_HPP
#define INVARIANTS_HPP

#include "zernike.hpp"
#include "coefs.hpp"

/** Class for computing rotational invariants from Zernike moments.
  The normalization of the result corresponds to the one of the z given.
  First use z.orthonormalize() if needed.
*/
class rotational_invariants
{
public:
  const int N;

  rotational_invariants(int n = 0);

  void eval(const zernike &zm);

  /** Index of a given invariant in the storage.
    @param n1 Should be positive.
    @param n2 Should be positive, not larger than n1 and with the same parity.
    @param l Should be between 0 and n2 and with the same parity.
    @return The index of invariant n1, n2, l in ri.
  */
  int index(int n1, int n2, int l) const
  {
    const int n1_2 = n1 / 2;
    const int n2_2 = n2 / 2;
    return l + n2_2 * (n2_2 + 1) + n1_2 * (n1_2 + 1) * (n1_2 + 2) / 3;
  }

  /** Value of a given invariant.
    @param n1 Should be positive and no larger than maximum order.
    @param n2 Should be positive, not larger than n1 and with the same parity.
    @param l Should be between 0 and n2 and with the same parity.
    @return The value of invariant n1, n2, l in ri.
  */
  double get(int n1, int n2, int l) const
  { return ri[index(n1, n2, l)]; }

  /** Direct access to data. */
  const std::vector<double> &get_ri() const
  { return ri; }

protected:
  std::vector<double> ri; /**< The storage for the result. */
};


class fnk
{
public:
  const int N;
  fnk(int n = 0);

  void eval(const rotational_invariants &ri);

  int index(int n, int k) const
  { return n * (n + 1) / 2 + k; }
  
  double get(int n, int k) const
  { return f[index(n, k)]; }

  const std::vector<double> &get_f() const
  { return f; }
protected:
  std::vector<double> f;
};


class invariants
{
public:
  const int N;
  invariants(int n);
  std::vector<double> fk0(const fnk &f)
  { return o0.apply(f.get_f()); }
  std::vector<double> fk3(const fnk &f)
  { return o3.apply(f.get_f()); }
  std::vector<mpq_class> hk0(const std::vector<mpq_class> &h)
  { return u0.apply(h); }
  std::vector<mpq_class> hk3(const std::vector<mpq_class> &h)
  { return u3.apply(h); }
  std::vector<double> k0k3(const std::vector<double> &k0)
  { return m03.apply(k0); }
  std::vector<double> k3k0(const std::vector<double> &k3)
  { return m30.apply(k3); }
private:
  factorials facs;
  double_factorials dfacs;
  binomials bins;
  unl0 u0;
  unl3 u3;
  vnl0 v0;
  vnl3 v3;
  theta t;
  ucompose m03, m30;
  omega o0, o3;
};

#endif
