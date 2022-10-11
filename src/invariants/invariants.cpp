/** \file invariants.cpp
  Implementation of invariants.hpp
  \author J. Houdayer
*/

#include "invariants.hpp"


/** Dummy constructor for operator >>.
 @param n The maximum order available. 
*/
rotational_invariants::rotational_invariants(int n):
zernike(n), ri((N / 2 + 1) * (N  / 2 + 2) * (N / 2 + 3) / 3, 0)
{}

/** Constructor.
 @param zm The moments to use to compute the rotational invariants.
 You can normalize it first if needed.
*/
void rotational_invariants::eval_ri()
{
  for (int n1_2 = 0, i = 0 ; n1_2 <= N / 2 ; n1_2++)
    for (int n2_2 = 0 ; n2_2 <= n1_2 ; n2_2++)
      for (int l = 0 ; l <= 2 * n2_2 + 1 ; l++, i++) {
        const int idx1 = zernike::index(2 * n1_2 + (l & 1), l, 0);
        const int idx2 = zernike::index(2 * n2_2 + (l & 1), l, 0);
        double sum = 0;
        for (int m = -l ; m <= l ; m++)
          sum += zm[idx1 + m] * zm[idx2 + m];
        ri[i] = sum;
      }
}

invariants_K::invariants_K(int n):
rotational_invariants(2*n), fnk((n + 1) * (n + 2) / 2, 0)
{}

void invariants_K::eval_fnk()
{
  eval_ri();
  for (int n = 0, i = 0 ; n <= N / 2 ; n++) {
    int sgn = 1, epsilon = 1;
    for (int k = 0 ; k <= n ; k++, i++, sgn = -sgn, epsilon = 2) {
      double sum = 0;
      const int idx = rotational_invariants::index(n + k, n - k, 0);
      for (int l = (n - k) & 1 ; l <= n - k ; l += 2)
        sum += ri[idx + l];
      fnk[i] = sgn * epsilon * (2 * (n + k) + 3) * (2 * (n - k) + 3) * sum;
    }
  }
}

// invariants_K::invariants_K(zernike &z):
// rotational_invariants(z), k0(z.order() / 2 + 1, 0)
// {
//   // compute fnk
//   const std::vector<double> &ri = get_ri();
//   std::vector<double> fnk((z.order() / 2 + 1) * (z.order() / 2 + 2) / 2, 0);
//   for (int n = 0, i = 0 ; n <= z.order() / 2 ; n++) {
//     int sgn = 1, epsilon = 1;
//     for (int k = 0 ; k <= n ; k++, i++, sgn = -sgn, epsilon = 2) {
//       double sum = 0;
//       int idx = rotational_invariants::index(n + k, n - k, 0);
//       for (int l = (n - k) & 1 ; l <= n - k ; l += 2)
//         sum += ri[idx + l];
//       fnk[i] = sgn * epsilon * (2 * (n + k) + 3) * (2 * (n - k) + 3) * sum;
//     }
//   }

//   // compute k0
//   for (int l = 0, i = 0 ; l <= z.order() / 2 ; l++) {
//     double sum = 0;
//     for (int n = 0, j = 0 ; n <= l ; n++)
//       for (int k = 0 ; k <= n ; k++, i++, j++)
//         sum += fnk_omega[i] * fnk[j];
//     k0[l] = sum;
//   }
// }
