/** \file zernike.hpp
  Classes to compute spherical harmonics and Zernike polynomials.
  \author J. Houdayer
*/

#ifndef ZERNIKE_HPP
#define ZERNIKE_HPP

#include "iotools.hpp"
#include "vec.hpp"

/** A pair of double.

  It stores coefficients for computing spherical harmonics and integrated zernike radial part.
  Used by spherical_harmonics and zernike_int0.
*/
class help2
{
public:
  double c1, c2;
  void set_sh(int l, int m);
  void set_int0(int n, int l);
};

/** A class for computing spherical harmonics.

  It allows the computation of all spherical harmonics
  up to order N given at creation.

  Usage:
    1. create one instance with the maximum order needed.
    2. use spherical_harmonics::eval with chosen parameters.
    3. get results with spherical_harmonics::get.
    4. go to step 2.
*/
class spherical_harmonics
{
public:
  const int N; /**< Maximum order. */

  spherical_harmonics(int n);
  void eval_sh(double theta, double phi);

  /**
    Index of element l, m in the storage.

    @param l Must be between 0 and N (included).
    @param m Must be between -l and l (included).
    @return Index of element l, m in storage. Between 0 and (N+1)^2 (excluded).
  */
  int index(int l, int m) const
  { return m + l * (l + 1); }

  /**
    Value of element l, m as computed by eval.

    @param l Must be between 0 and N (included).
    @param m Must be between -l and l (included).
    @return Value of spherical harmonic l, m computed by spherical_harmonics::eval.
  */
  double get(int l, int m) const
  { return sh[index(l, m)]; }

  /** A direct access to data. */
  const std::vector<double> &get_sh() const
  { return sh;}

protected:
  std::vector<double> sh; /**< Storage for the result. */
private:
  std::vector<help2> help; /**< Fixed coefficients needed by computation. */
};

/** Base class for computing radial part of zernike polynomials.
*/
class zernike_radial
{
public:
  const int N; /**< Maximum order. Positive. */

  zernike_radial(int n);

  /** Index of element n, l in the storage.
    Elements are ordered like this:
    (0, 0), (1, 1), (2, 0), (3, 1), (2, 2), (3, 3),
    (4, 0), (5, 1), (4, 2), (5, 3), (4, 4), (5, 5) ...

    @param n Must be between 0 and N (included).
    @param l Must be between 0 and N (included), with the same parity as n.
    @return Index of element n, l in storage. Between 0 and (N / 2 + 1) * (N / 2 + 2) (excluded).
  */
  int index(int n, int l) const
  { return l + (n / 2) * (n / 2 + 1); }

    /** Value of element n, l.

    @param n Must be between 0 and N (included).
    @param l Must be between 0 and N (included), with the same parity as n.
    @return Value of element n, l.
  */
  double get(int n, int l) const
  { return zr[index(n, l)]; }

  /** A direct access to data. */
  const std::vector<double> &get_zr() const
  { return zr;}

  /** Reset all elements to zero. */
  void reset_zr();

protected:
  std::vector<double> zr; /**< Storage for the result. */
};

/** A truple of coefficients used by zernike_r and zernike_int2. */
class help3
{
public:
  double c1, c2, c3;
  void set_r(int n, int l);
  void set_int2(int n, int l);
};

/** A class to compute the radial part of the Zernike polynomials.

  It allows the computation of all radial parts
  up to order N given at creation.

  Normalization does not include the usual \f$\sqrt{2n+3}\f$ term. So that
  \f[ \int Z_{nlm}^2 = \frac 1{2n+3}. \f]

  Usage:
    1. create one instance with the maximum order needed.
    2. use zernike_r::eval_zr with chosen parameters.
    3. get results with zernike_r::get_zr.
    4. go to step 2.
*/
class zernike_r: public zernike_radial
{
public:
  zernike_r(int n);
  void eval_zr(double r, double weight = 1);
private:
  std::vector<help3> help; /**< Fixed coefficients used in the computation. */
};

/** A class to compute the integrated radial part of the Zernike polynomials.

  It allows the computation of all integrated radial parts
  up to order N given at creation.

  More precisely it computes:
  \f[ \int_0^r R_{n,l}(x)\mathrm dx, \f]
  where the normalization of \f$R\f$ is the same as in zernike_r.

  Usage:
    1. create one instance with the maximum order needed.
    2. use zernike_int0::eval_zr with chosen parameters.
    3. get results zernike_int0::get_zr.
    4. go to step 2.
*/
class zernike_int0: public zernike_radial
{
public:
  zernike_int0(int n);
  void eval_zr(double r, double weight = 1);
private:
  std::vector<help2> help; /**< Fixed coefficients used in the computation. */
  zernike_r base_r;
};

/** A class to compute the integrated radial part of the Zernike polynomials.

  It allows the computation of all integrated radial parts
  up to order N given at creation.

  More precisely it computes:
  \f[ \int_0^r x^2 R_{n,l}(x)\mathrm dx, \f]
  where the normalization of \f$R\f$ is the same as in zernike_r.

  Usage:
    1. create one instance with the maximum order needed.
    2. use zernike_int2::eval_zr with chosen parameters.
    3. get results zernike_int2::get_zr.
    4. go to step 2.
*/
class zernike_int2: public zernike_radial
{
public:
  zernike_int2(int n);
  void eval_zr(double r, double weight = 1);
private:
  std::vector<help3> help; /**< Fixed coefficients used in the computation. */
  zernike_int0 base_0;
};

/** Enumeration to represent possible Zernike moments normalizations.
 raw is ortho / sqrt{2n+3}, it is used for evaluation of the moments
 ortho is the orthonormal normalization
 dual is ortho * sqrt{2n+3}, it is used for reconstruction of the original density

 raw_n, ortho_n, dual_n are the same with an additional 3/(4pi) factor.
 All zernike moments are created raw. Use zernike::normalize to change this.
*/
enum class zm_norm {raw, ortho, dual, raw_n, ortho_n, dual_n};

/** Creates zm_norm object. 
 By defaut (all args false) it gives ortho, the args sets the different possibilities.
*/
zm_norm make_norm(bool raw, bool dual, bool norm);

/** Enumeration to represent types of output for Zernike moments.*/
enum class zm_output {real, complex, real_p, complex_p};

zm_output make_output(bool cplx, bool phase);

/** A base class to compute Zernike moments. */
class zernike
{
public:
  zernike(int n=0);
  zernike(int n, const zernike &source);

  /** Maximum order available .*/
  int order() const
  { return N; }

  /** Index of element n, l, m in storage.

    The order for n and l is the same as zernike_radial
    and for each value of n and l, m change from -l to l, namely:
    (0, 0, 0), (1, 1, -1), (1, 1, 0), (1, 1, 1),
    (2, 0, 0), (3, 1, -1), (3, 1, 0), (3, 1, 1),
    (2, 2, -2) ... (2, 2, 2), (3, 3, -3) ... (3, 3, 3), (4, 0, 0) ...

    @param n Between 0 and N (included).
    @param l Between 0 and n (included), should have the same parity as n.
    @param m Between -l and l (included).
    @return The index of element n, l, m in storage,
    between 0 and 2 * (n / 2 + 1) * (n / 2 + 2) * (2 * (n / 2) + 3) / 3) (excluded).
  */
  int index(int n, int l, int m) const
  {
    const int n2 = n / 2;
    return m + l * (l + 1) + (2 * n2 * (n2 + 1) * (2 * n2 + 1)) / 3;
  }

  /** Value of element n, l, m.
    @param n Between 0 and N (included).
    @param l Between 0 and n (included), should have the same parity as n.
    @param m Between -l and l (included).
    @return The value of element n, l, m.
  */
  double get(int n, int l, int m) const
  { return zm[index(n, l, m)]; }

  /** A direct access to data. */
  const std::vector<double> &get_zm() const
  { return zm;}

  /** The current norm used. */
  zm_norm get_norm() const
  { return norm; }

  /** An estimation of the error made during evaluation.*/
  double get_error() const
  { return sqrt(variance); }

  void reset_zm();
  void normalize(zm_norm new_norm);
  double operator()(const vec &v) const;
  void finish();
  double distance(const zernike &z) const;
  zernike &operator +=(const zernike &z);

  friend smart_input &operator >>(smart_input &, zernike &);
  friend zernike operator -(const zernike &z1, const zernike &z2);
 
  double variance;
  zm_output output;
protected:
  int N;
  zm_norm norm;
  bool odd_clean;
  std::vector<double> zm; /**< Storage for the results. */

  void add_core(const std::vector<double> &z, const std::vector<double> &sh,
                double weight);
};

zernike operator -(const zernike &z1, const zernike &z2);

std::ostream &operator <<(std::ostream &, const zernike &);
smart_input &operator >>(smart_input &, zernike &);

/** Class for computing weighted sums of zernike polynomials.

  Normalization is the one from zernike_r.

  Usage :
    1. create one instance with the maximum order needed.
    2. Use zernike::reset_zm to start from 0.
    3. Repeatedly call zernike_m_r::add to add the
    corresponding polynomials with the given weights.
    4. normalize if needed with zernike_m::normalize to fix element 0,0,0.
    5. use result
    6. go to 2.
 */
class zernike_m_r:
public zernike_r, public spherical_harmonics, public zernike
{
public:
  zernike_m_r(int n);
  void add(const w_vec &p);
};

/** Class for computing weighted sums of integrated zernike polynomials.

  Normalization is the one from zernike_int2.

  Usage :
    1. create one instance with the maximum order needed.
    2. Use zernike::reset_zm to start from 0.
    3. Repeatedly call zernike_m_int::add to add the
    corresponding integrated polynomials with the given weights.
    4. normalize if needed with zernike_m::normalize to fix element 0,0,0.
    5. use result
    6. go to 2.
 */
class zernike_m_int:
public zernike_int2, public spherical_harmonics, public zernike
{
public:
  zernike_m_int(int n);
  void add(const w_vec &p);
};

/** Class for computing rotational invariants from Zernike moments.
  The normalization of the result corresponds to the one of the z given.
  First use z.orthonormalize() if needed.
*/
class rotational_invariants
{
public:
  rotational_invariants(int n = 0);

  void eval_ri(const zernike &z);

  /** Maximum order available .*/
  int order() const
  { return N; }

  /** The current norm used. */
  zm_norm get_norm() const
  { return norm; }

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

  friend smart_input &operator >>(smart_input&, rotational_invariants &);
  friend rotational_invariants operator -(const rotational_invariants &, const rotational_invariants &);

private:
  int N;
  zm_norm norm;
  std::vector<double> ri; /**< The storage for the result. */
};

rotational_invariants operator -(const rotational_invariants &r1, const rotational_invariants &r2);

std::ostream &operator <<(std::ostream &, const rotational_invariants &);
smart_input &operator >>(smart_input &, rotational_invariants &);

class signature_invariants
{
public:
  signature_invariants(int n = 0);

  void eval_si(const zernike &z);

   /** Maximum order available .*/
  int order() const
  { return N; }

  /** The current norm used. */
  zm_norm get_norm() const
  { return norm; }

  /** Value of a given invariant.
    @param n Should be positive and no larger than maximum order.
    @return The value of invariant n.
  */
  double get(int n) const
  { return si[n]; }

  /** Direct access to data. */
  const std::vector<double> &get_si() const
  { return si; }

  friend smart_input &operator >>(smart_input &, signature_invariants &);
  friend signature_invariants operator -(const signature_invariants &s1, const signature_invariants &s2);
private:
  int N;
  zm_norm norm;
  std::vector<double> si; /**< The storage for the result. */
};

signature_invariants operator -(const signature_invariants &s1, const signature_invariants &s2);
std::ostream &operator <<(std::ostream &, const signature_invariants &);
smart_input &operator >>(smart_input &, signature_invariants &);

#endif
