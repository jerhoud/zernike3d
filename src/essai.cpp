#include <iostream>
#include <cmath>
#include "brent.hpp"

double f(double x)
{
  return x * x - 3.135645 * x + 7;
}

double g(double x)
{
  x *= x;
  return (x - 20) * x;
}

double h(double x)
{
  if (x<0)
    return -x;
  else
    return 0.0001 * x;
}

int main() {
  std::cout << minimize(h, 3, 1, 1e-8) << "\n";

  return 0;
}