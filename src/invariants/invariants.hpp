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

class inv
{
public:
  inv(inv_coefs &ic): cfs(ic) {}
  bool isexact() const
  { return exact; }
  const std::vector<mpq_class> & get_q() const
  { return cq; }
  const std::vector<double> & get_d() const
  { return cd; }
  double get_D() const
  { return D; }
  void set(double sz, const std::vector<mpq_class> q);
  void set(double sz, const std::vector<double> d);
  void set_scale(double d)
  { D = d; }
  void noexact()
  { exact = false; }
  void normalize();
protected:
  inv_coefs &cfs;
  bool exact;
  double D;
  std::vector<mpq_class> cq;
  std::vector<double> cd;
};

std::ostream &operator <<(std::ostream &os, const inv &i);
smart_input &operator >>(smart_input &is, inv &i);

class inv_h;

class inv_k3: public inv
{
public:
  inv_k3(inv_coefs &ic): inv(ic) {}
  void eval(double sz, const fnk &f);
  void eval(const inv_h &h);
  void resize(double alpha);
};

inline std::ostream &operator <<(std::ostream &os, const inv_k3 &k3)
{ return os << "K3\n" << static_cast<inv>(k3); }

class inv_h: public inv
{
public:
  inv_h(inv_coefs &ic): inv(ic) {}
  void eval(double sz, const fnk &f);
  void eval(const inv_k3 &k3);
};

inline std::ostream &operator <<(std::ostream &os, const inv_h &h)
{ return os << "H\n" << static_cast<inv>(h); }

class hball: public inv_h
{
public:
  hball(inv_coefs &ic, double sz);
};

class hcube: public inv_h
{
public:
  hcube(inv_coefs &ic, double sz);
};

#endif