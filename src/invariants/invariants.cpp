/** \file invariants.cpp
  Implementation of invariants.hpp
  \author J. Houdayer
*/

#include "invariants.hpp"
#include "s_root_data.hpp"


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
  for (int n = 0, i = 0 ; n <= N ; n++) {
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

void inv::set(double sz, const std::vector<mpq_class> q)
{
  exact = true;
  D = sz;
  cq = q;
  cd = std::vector<double>(cq.size());
  for (size_t i = 0 ; i < q.size() ; i++)
    cd[i] = cq[i].get_d();
}

void inv::set(double sz, const std::vector<double> d)
{
  exact = false;
  D = sz;
  cd = d;
}

void inv::normalize()
{
  if (exact) {
    if (cq.size() == 0)
      return;
    const mpq_class d = cq[0];
    if (d == 0)
      return;
    for (auto &c : cq)
      c /= d;
    double dd = d.get_d();
    for (auto &c : cd)
      c /= dd;
  }
  else {
    if (cd.size() == 0)
      return;
    const double d = cd[0];
    if (d == 0)
      return;
    for (auto &c : cd)
      c /= d;
  }
}

std::ostream &operator <<(std::ostream &os, const inv &i)
{
  if (i.isexact()) {
    const std::vector<mpq_class> &v = i.get_q();
    os << v.size() - 1 << " " << i.get_D() << "\n";
    for (auto &x: v)
      if (x.get_den()==1)
        os << x << "/1\n";
      else
        os << x << "\n";
  }
  else {
    const std::vector<double> &v = i.get_d();
    os << v.size() - 1 << " " << i.get_D() << "\n";
    for (auto &x: v)
      os << x << "\n";
  }
  return os;
}

smart_input &operator >>(smart_input &is, inv &i)
{
  std::istringstream s;
  if (!is.next_line(s)) //remove first line containing "H", "K0" or "K3"
    return is.failed();
  if (!is.next_line(s))
    return is.failed();
  int n0;
  double sz;
  s >> n0 >> sz;
  if (!s || n0 < 0)
    return is.failed();
  is.peek_line(s);
  std::string testslash;
  s >> testslash;
  if (!s)
    return is.failed();
  if (testslash.find('/')==std::string::npos) { // read reals
    std::vector<mpq_class> qs(n0 + 1);
    for (auto &q : qs) {
      is.next_line(s);
      s >> q;
      if (!s)
        return is.failed();
    }
    if (!is)
      return is.failed();
    i.set(sz, qs);
  }
  else { // read rationals
    std::vector<double> ds(n0 + 1);
    for (auto &d : ds) {
      is.next_line(s);
      s >> d;
      if (!s)
        return is.failed();
    }
    if (!is)
      return is.failed();
    i.set(sz, ds);
  }
  return is;
}

void inv_k3::eval(double sz, const fnk &f)
{ set(sz, cfs.get_o3().apply(f.get_f())); }

void inv_k3::eval(const inv_h &h)
{
  if (h.isexact())
    set(h.get_D(), cfs.u3.apply(h.get_q()));
}

std::vector<double> inv_k3::resized(int n0, double alpha) const
{
  const double x = 1. / (alpha * alpha);
  std::vector<double> k(n0 + 1);
  int idx = 0;
  for (int n = 0 ; n <= n0 ; n++) {
    double sum = 0;
    double xl = 1;
    for (int l = 0 ; l <= n ; l++, xl *= x) {
      double a = xl * s_root_data[idx++];
      for (int m = 0 ; m < n - l ; m++)
        a *= x - s_root_data[idx++];
      sum += a;
    }
    k[n] = sum;
  }
  return k;
}

void inv_h::eval(double sz, const fnk &f)
{ set(sz, cfs.get_t().apply(f.get_f())); }

void inv_h::eval(const inv_k3 &k3)
{
  if (k3.isexact())
    set(k3.get_D(), cfs.v3.apply(k3.get_q()));
  else
    set(k3.get_D(), cfs.v3.apply(k3.get_d()));
}

hball::hball(inv_coefs &ic, double sz):
inv_h(ic)
{
  const int N = cfs.N;
  std::vector<mpq_class> h(N + 1);
  for (int n = 0 ; n <= N ; n++)
    h[n] = 18_mpq / ((n + 2) * (n + 3) * (2 * n + 3));
  set(sz, h);
}

hcube::hcube(inv_coefs &ic, double sz):
inv_h(ic)
{
  const int N = cfs.N;
  std::vector<mpq_class> h(N + 1);
  mpz_class pow3n = 1;
  for (int n = 0 ; n <= N ; n++, pow3n *= 3) {
    mpq_class sum = 0;
    for (int n1 = 0 ; n1 <= n ; n1++) {
      const mpz_class t1 = (n1 + 1) * (2 * n1 + 1);
      for (int n2 = 0 ; n2 <= n - n1 ; n2++) {
        const mpz_class t2 = (n2 + 1) * (2 * n2 + 1) * t1;
        const int n3 = n - n1 - n2;
        mpq_class tmp = mpq_class(cfs.bins.get(n, n1) * cfs.bins.get(n - n1, n2),
                                  (n3 + 1) * (2 * n3 + 1) * t2);
        tmp.canonicalize();
        sum += tmp;
      }
    }
    h[n] = sum / pow3n;
  }
  set(sz, h);
}
