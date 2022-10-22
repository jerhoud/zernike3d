/** \file invariants.cpp
  Implementation of invariants.hpp
  \author J. Houdayer
*/

#include "invariants.hpp"


/** Dummy constructor for operator >>.
 @param n The maximum order available. 
*/
rotational_invariants::rotational_invariants(int n):
N(n), ri((N / 2 + 1) * (N  / 2 + 2) * (N / 2 + 3) / 3, 0)
{}

/** Constructor.
 @param zm The moments to use to compute the rotational invariants.
 You can normalize it first if needed.
*/
void rotational_invariants::eval(const zernike &zm)
{
  const std::vector<double> &z = zm.get_zm();
  for (int n1_2 = 0, i = 0 ; n1_2 <= N / 2 ; n1_2++)
    for (int n2_2 = 0 ; n2_2 <= n1_2 ; n2_2++)
      for (int l = 0 ; l <= 2 * n2_2 + 1 ; l++, i++) {
        const int idx1 = zm.index(2 * n1_2 + (l & 1), l, 0);
        const int idx2 = zm.index(2 * n2_2 + (l & 1), l, 0);
        double sum = 0;
        for (int m = -l ; m <= l ; m++)
          sum += z[idx1 + m] * z[idx2 + m];
        ri[i] = sum;
      }
}

fnk::fnk(int n):
N(n), f((n + 1) * (n + 2) / 2, 0)
{}

void fnk::eval(const rotational_invariants &ri)
{
  for (int n = 0, i = 0 ; n <= N / 2 ; n++) {
    int sgn = 1, epsilon = 1;
    for (int k = 0 ; k <= n ; k++, i++, sgn = -sgn, epsilon = 2) {
      double sum = 0;
      const int idx = ri.index(n + k, n - k, 0);
      for (int l = (n - k) & 1 ; l <= n - k ; l += 2)
        sum += ri.get_ri()[idx + l];
      f[i] = sgn * epsilon * (2 * (n + k) + 3) * (2 * (n - k) + 3) * sum;
    }
  }
}

invariants::invariants(int n):
N(n), facs(2 * n + 1), dfacs(2 * n + 1), bins(2 * n + 3),
u0(N, bins), u3(n, bins), v0(n, facs), v3(n, facs, dfacs, bins),
t(n, facs, dfacs, bins), m03(u3, v0), m30(u3, v0), o0(u0, t), o3(u3, t)
{}
