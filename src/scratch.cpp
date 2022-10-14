#include <iostream>

#include "coefs.hpp"

int main()
{
  const int N = 51;

  HK_coefs hk(N);
  for (int m = 0, idx = 0 ; m <= N ; m++)
    for (int n = 0 ; n <= m ; n++)
      for (int k = 0 ; k <= n ; k++, idx++)
        std::cout << m << ", " << n << ", " << k << " = " << 9 * hk.k3_coefs[idx].get_d() / 4 << "\n"; 
}