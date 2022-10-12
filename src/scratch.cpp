#include <iostream>

#include "largeint.hpp"

int main()
{
  const int N = 5;
  // factorials f(N);

  // for (int n = 0 ; n <= N ; n++)
  //   std::cout << n << "! = " << f.get(n) << "\n";
  
  // double_factorials f(N);

  // for (int n = 0 ; n <= N ; n++)
  //   std::cout << 2 * n + 1 << "!! = " << f.get(n) << "\n";

  // binomials b(N);
  // for (int n = 0 ; n <= N ; n++)
  //   for (int l = 0 ; l <= n ; l++)
  //     std::cout << n << ", " << l << " = " << b.get(n, l) << "\n";

  // factorials f(2 * N + 1);
  // binomials b(2 * N);
  // unl u(N, f, b);
  // for (int n = 0 ; n <= N ; n++)
  //   for (int l = 0 ; l <= n ; l++)
  //     std::cout << n << ", " << l << " = " << u.get(n, l) << "\n";

  // double_factorials df(2 * N + 1);
  // theta th(N, df, b);

  omega omg(N);
  for (int m = 0 ; m <= N ; m++)
    for (int n = 0 ; n <= m ; n++)
      for (int k = 0 ; k <= n ; k++)
        std::cout << m << ", " << n << ", " << k << " = " << omg.get(m, n, k) << "\n"; 
}