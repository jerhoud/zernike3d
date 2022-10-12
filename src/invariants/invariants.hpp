/** \file invariants.hpp
  A class to compute rotational and translational invariants
  \author J. Houdayer
*/

#include "zernike.hpp"

#ifndef INVARIANTS_HPP
#define INVARIANTS_HPP

/** Class for computing rotational invariants from Zernike moments.
  The normalization of the result corresponds to the one of the z given.
  First use z.orthonormalize() if needed.
*/
class rotational_invariants: public zernike
{
public:
  rotational_invariants(int n = 0);

  void eval_ri();

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


/** A class for computing K rotational and translational invariants.
*/
class invariants_K: public rotational_invariants
{
public:
  invariants_K(int n);

  void eval_K();

/** The value of the l-th invariant.
  @param l Must be between 0 and the maximum order / 2.
*/
  double get(int l) const
  { return K0[l]; }

  const std::vector<double> get_K0() const
  { return K0; }
private:
  void eval_fnk();

  std::vector<double> K0;
  std::vector<double> fnk;
  std::vector<double> coefs;
};

#endif
