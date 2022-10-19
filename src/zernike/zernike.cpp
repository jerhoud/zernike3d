/** \file zernike.cpp
  Implementation of zernike.hpp.
  \author J. Houdayer
*/

#include "zernike.hpp"
#include <sstream>
#include <iomanip>
#include <numeric>

#ifndef M_PI
#define M_PI 3.141592653589793238
#endif

void help2::set_sh(int l, int m)
{
  if (m < l - 1) {
    const double a = (2 * l + 1) / (double) ((l + m) * (l - m));
    c1 = sqrt(a * (2 * l - 1));
    c2 = sqrt(a * (l + m - 1) * (l - m - 1) / (double) (2 * l - 3));
  }
  else if (m == l - 1)
    c1 = sqrt(2 * l + 1);
  else
    c1 = sqrt(1 + 0.5 / l);
}

void help2::set_int0(int n, int l)
{
  c1 = (2 * l + 3) / (double) ((2 * n + 3) * (l + 1));
  c2 = (l + 2) / (double) (l + 1);
}

/** Constructor.
  @param n Maximum order N for the computation. Should be positive.
*/
spherical_harmonics::spherical_harmonics(int n):
N(n), sh((N + 1) * (N + 1), 0), help((N + 1) * (N + 2)  / 2, {0, 0})
{
  for (int l = 1, i = 1 ; l <= N ; l++)
    for (int m = 0 ; m <= l ; m++, i++)
      help[i].set_sh(l, m);
}

/** Runs the computation.

  Evaluates all spherical harmonics up to order N. For the given parameters.
  Results are accessed by spherical_harmonics::get.
  The normalization is the usual one:
  spherical harmonics form an orthonormal base.

  @param theta The colatitude. Must be between 0 and pi.
  @param phi The longitude. Should be between -pi and pi (but works in any case).
*/
void spherical_harmonics::eval_sh(double theta, double phi)
{
  const double x = cos(theta);
  const double sx = - sin(theta);
  double mm = 1 / sqrt(4 * M_PI);

  sh[0] = mm;
  mm *= sqrt(2);
  for (int l = 1, i = 2, j = 1 ; l <= N ; l++, i += 2 * l) {
    for (int m = 0 ; m < l - 1 ; m++, j++) {
      const double xc1 = x * help[j].c1, c2 = help[j].c2;
      sh[i + m] = xc1 * sh[i + m - 2 * l] - c2 * sh[i + m - 4 * l + 2];
      sh[i - m] = xc1 * sh[i - m - 2 * l] - c2 * sh[i - m - 4 * l + 2];
    }
    const double xc1 = x * help[j++].c1;
    sh[i + l - 1] = xc1 * sh[i - 1 - l];
    sh[i - l + 1] = xc1 * sh[i - 3 * l + 1];
    mm *= sx * help[j++].c1;
    const double lphi = l * phi;
    sh[i + l] = mm * cos(lphi);
    sh[i - l] = mm * sin(lphi);
  }
}

/** Constructor.
  @param n Maximum order needed. Should be positive.
*/
zernike_radial::zernike_radial(int n):
N(n), zr((N / 2 + 1) * (N / 2 + 2), 0)
{}

void zernike_radial::reset_zr()
{
  for (auto &v: zr)
    v = 0;
}

void help3::set_r(int n, int l)
{
  if (n > l) {
    const double np1 = 2 * n + 1;
    const double nm1 = np1 - 2;
    const double nm3 = np1 - 4;
    const double lp1 = 2 * l + 1;
    c1 = nm1 * np1 / (double) ((n - l) * (n + l + 1));
    c2 = lp1 * lp1 / (2 * nm3 * np1) + 0.5;
    c3 = (n - l - 2) * (n + l - 1) / (nm3 * nm1);
  }
}

void help3::set_int2(int n, int l)
{
  if (n > l) {
    const double np1 = 2 * n + 1;
    const double np3 = np1 + 2;
    const double np5 = np3 + 2;
    c1 = (n - l + 2) * (n + l + 3) / (np3 * np5);
    c2 = 0.5 * (1 + (2 * l + 1) * (2 * l + 1) / (np5 * np1));  
    c3 = (n - l) * (n + l + 1) / (np3 * np1);
  }
}

/** Constructor.
  @param n Maximum order needed. Should be positive.
*/
zernike_r::zernike_r(int n):
zernike_radial(n), help((N / 2 + 1) * (N / 2 + 2), {0, 0, 0})
{
  int i = 0;
  for (int n2 = 0 ; n2 <= N / 2 ; n2++)
    for (int l2 = 0 ; l2 <= n2 ; l2++) {
      help[i++].set_r(2 * n2, 2 * l2);
      help[i++].set_r(2 * n2 + 1, 2 * l2 + 1);
    }
}

/** Runs the computation.
  @param r The radial parameter. Between 0 and 1.
*/
void zernike_r::eval_zr(double r, double weight)
{
  const double r2 = r * r;
  const help3 *h;
  double rn = r2 * weight;

  zr[0] = weight;
  zr[1] = r * weight;
  for (int n2 = 1, i = 2 ; n2 <= N / 2 ; n2++, rn *= r2, i++) {
    for (int l = 0 ; l < 2 * n2 - 2 ; l++, i++) {
      h = &(help[i]);
      zr[i] = h->c1 * ((r2 - h->c2) * zr[i - 2 * n2]
                      - h->c3 * zr[i - 4 * n2 + 2]);
    }
    h = &(help[i]);
    zr[i] = h->c1 * (r2 - h->c2) * zr[i - 2 * n2];
    h = &(help[++i]);
    zr[i] = h->c1 * (r2 - h->c2) * zr[i - 2 * n2];
    zr[++i] = rn;
    zr[++i] = r * rn;
  }
}

/** Constructor.
  @param n Maximum order needed. Should be positive.
*/
zernike_int0::zernike_int0(int n):
zernike_radial(n), help((N / 2 + 1) * (N / 2 + 2), {0, 0}), base_r(n+1)
{
  int i = 0;
  for (int n2 = 0 ; n2 <= N / 2 ; n2++)
    for (int l2 = 0 ; l2 <= n2 ; l2++) {
      help[i++].set_int0(2 * n2, 2 * l2);
      help[i++].set_int0(2 * n2 + 1, 2 * l2 + 1);
    }
}

/** Runs the computation.
  @param r The radial parameter. Between 0 and 1.
  @param weight An optional multiplicative weight.
*/
void zernike_int0::eval_zr(double r, double weight)
{
  base_r.eval_zr(r, weight);
  const std::vector<double> zr0 = base_r.get_zr();
  const help2 *h;
  const double r2 = r * r;
  double rn1 = r * weight;
  zr[0] = rn1;
  zr[1] = 0.5 * r * rn1;
  rn1 *= r2;
  for (int n2 = 1, i = 6 ; n2 <= N / 2 ; n2++, rn1 *= r2, i += 4 * n2 + 2) {
    double todd = zr[--i] = r * rn1 / (double) (2 * n2 + 2);
    double teven = zr[--i] = rn1 / (double) (2 * n2 + 1);
    for (int l2 = n2 - 1 ; l2 >=0 ; l2--) {
      h = &(help[--i]);
      todd = zr[i] = h->c1 * (zr0[i + 2 * n2 + 3] - zr0[i + 1]) - h->c2 * todd;
      h = &(help[--i]);
      teven = zr[i] = h->c1 * (zr0[i + 1] - zr0[i - 2 * n2 + 1]) - h->c2 * teven;
    }
  }
}

/** Constructor.
  @param n Maximum order needed. Should be positive.
*/
zernike_int2::zernike_int2(int n):
zernike_radial(n), help((N / 2 + 1) * (N / 2 + 2), {0, 0, 0}), base_0(n+2)
{
  int i = 0;
  for (int n2 = 0 ; n2 <= N / 2 ; n2++)
    for (int l2 = 0 ; l2 <= n2 ; l2++) {
      help[i++].set_int2(2 * n2, 2 * l2);
      help[i++].set_int2(2 * n2 + 1, 2 * l2 + 1);
    }
}

/** Runs the computation.
  @param r The radial parameter. Between 0 and 1.
  @param weight An optional multiplicative weight.
*/
void zernike_int2::eval_zr(double r, double weight)
{
  base_0.eval_zr(r, weight);
  const std::vector<double> zr0 = base_0.get_zr();
  const double r2 = r * r;
  const help3 *h;
  double rn = r2 * r * weight;

  zr[0] = rn / 3;
  zr[1] = r * rn / 4;
  rn *= r2;
  for (int n2 = 1, i = 2 ; n2 <= N / 2 ; n2++, rn *= r2, i++) {
    for (int l = 0 ; l <= 2 * n2 - 2 ; l++, i++) {
      h = &(help[i]);
      zr[i] = h->c1 * zr0[i + 2 * n2 + 2] + h->c2 * zr0[i] + h->c3 * zr0[i - 2 * n2];
    }
    zr[++i] = rn / (double) (2 * n2 + 3);
    zr[++i] = r * rn / (double) (2 * n2 + 4);
  }
}

/** Output operator for \a zm_norm. */
std::ostream &operator <<(std::ostream &os, zm_norm norm)
{
  switch (norm) {
    case zm_norm::raw:
      os << "RAW";
      break;
    case zm_norm::ortho:
      os << "ORTHO";
      break;
    case zm_norm::dual:
      os << "DUAL";
      break;
     case zm_norm::raw_n:
      os << "RAW_N";
      break;
    case zm_norm::ortho_n:
      os << "ORTHO_N";
      break;
    case zm_norm::dual_n:
      os << "DUAL_N";
      break;
  }
  return os;
}

/** Input operator for \a zm_norm. */
std::istream &operator >>(std::istream &is, zm_norm &norm)
{
  std::string s;
  is >> s;
  if (s=="RAW") {
    norm = zm_norm::raw;
    return is;
  }
  if (s=="ORTHO") {
    norm = zm_norm::ortho;
    return is;
  }
  if (s=="DUAL") {
    norm = zm_norm::dual;
    return is;
  }
  if (s=="RAW_N") {
    norm = zm_norm::raw_n;
    return is;
  }
  if (s=="ORTHO_N") {
    norm = zm_norm::ortho_n;
    return is;
  }
  if (s=="DUAL_N") {
    norm = zm_norm::dual_n;
    return is;
  }
  return failed(is);
}

zm_norm make_norm(bool raw, bool dual, bool norm)
{
  if (norm) {
    if (raw) return zm_norm::raw_n;
    if (dual) return zm_norm::dual_n;
    return zm_norm::ortho_n;
  }
  else {
    if (raw) return zm_norm::raw;
    if (dual) return zm_norm::dual;
    return zm_norm::ortho;   
  }
}

/** Output operator for \a zm_output. */
std::ostream &operator <<(std::ostream &os, zm_output output)
{
  switch (output) {
    case zm_output::real:
      os << "REAL";
      break;
    case zm_output::complex:
      os << "COMPLEX";
      break;
    case zm_output::real_p:
      os << "REAL_P";
      break;
    case zm_output::complex_p:
      os << "COMPLEX_P";
      break;
  }
  return os;
}

/** Input operator for \a zm_output. */
std::istream &operator >>(std::istream &is, zm_output &output)
{
  std::string s;
  is >> s;
  if (s=="REAL") {
    output = zm_output::real;
    return is;
  }
  if (s=="COMPLEX") {
    output = zm_output::complex;
    return is;
  }
  if (s=="REAL_P") {
    output = zm_output::real_p;
    return is;
  }
  if (s=="COMPLEX_P") {
    output = zm_output::complex_p;
    return is;
  }
  return failed(is);
}

zm_output make_output(bool cplx, bool phase)
{
  return (zm_output) ((int) cplx + 2 * (int) phase);
}

bool flip_out(zm_output output)
{
  return output == zm_output::real_p || output == zm_output::complex_p;
}

bool real_out(zm_output output)
{
  return output == zm_output::real || output == zm_output::real_p;
}

/** Constructor.
  @param n Maximum order needed. Should be positive.
*/
zernike::zernike(int n):
variance(0), output(zm_output::real),
N(n), norm(zm_norm::raw), odd_clean(false),
zm(2 * (n / 2 + 1) * (n / 2 + 2) * (2 * (n / 2) + 3) / 3, 0)
{}

/** Constructor.
  @param n Maximum order needed. Should be positive.
  @param source Moments to copy from (truncated at N = n.).
*/
zernike::zernike(int n, const zernike &source):
variance(0), output(source.output),
N(n), norm(source.norm), odd_clean(false),
zm(2 * (n / 2 + 1) * (n / 2 + 2) * (2 * (n / 2) + 3) / 3, 0)
{
  std::copy(source.get_zm().begin(), source.get_zm().begin() + zm.size(),
  zm.begin());
  finish();
}

/** The core of the computation.
  Nearly all computation time is concentrated here.

  This multiplies the radial part with the spherical harmonics
  and adds it with the corresponding weight.

  @param z The radial part. Its maximum order should be the same (or larger).
  @param sh The spherical harmonics. Its maximum order should be the same (or larger).
  @param weight The weight.
*/
void zernike::add_core(const std::vector<double> &z,
                       const std::vector<double> &sh,
                       double weight)
{
  int idzr = 0;
  int idz = 0;
  for (int n2 = 0 ; n2 <= N / 2 ; n2++) {
    int idsh = 0;
    for (int l = 0 ; l <= 2 * n2 + 1 ; l++, idzr++) {
      double r = weight * z[idzr];
      for (int m = -l ; m <= l ; m++, idsh++, idz++)
        zm[idz] += r * sh[idsh];
    }
  }
}

/** Reset the computation to 0. */
void zernike::reset_zm()
{
  norm = zm_norm::raw;
  odd_clean = false;
  for (auto &i: zm)
    i = 0;
}

/** To use to get the correct normalization for orthonormal Zernike polynomials.
  Do not use this before the computation of the moments is over.
*/
void zernike::normalize(zm_norm new_norm)
{
  finish();
  int dn = ((int)new_norm) % 3 - ((int)norm) % 3;
  int fn = ((int)new_norm) / 3 - ((int)norm) / 3;
  if (dn == 0 && fn == 0) return;
  double fact = 1;
  if (fn == 1)
    fact = sqrt(3 / (4 * M_PI));
  if (fn == -1)
    fact = sqrt(4 * M_PI / 3);
  bool r = dn < 0;
  int f = abs(dn);
  double sn[2];
  int idz = 0;
  for (int n2 = 0 ; n2 <= N / 2 ; n2++) {
    if (f==1) {
      sn[0] = sqrt(4 * n2 + 3);
      sn[1] = sqrt(4 * n2 + 5);
    }
    else if (f==2) {
      sn[0] = 4 * n2 + 3;
      sn[1] = 4 * n2 + 5;
    }
    else
      sn[0] = sn[1] = 1;
    if (r) {
      sn[0] = 1 / sn[0];
      sn[1] = 1 / sn[1];
    }
    sn[0] *= fact;
    sn[1] *= fact;
    for (int l = 0 ; l <= 2 * n2 + 1 ; l++) {
      double s = sn[l & 1]; // s = sqrt(3/(4pi))^fn * sqrt(2n+3)^dn
      for (int m = -l ; m <= l ; m++, idz++)
        zm[idz] *= s;
    }
  }
  norm = new_norm;
}

/** Evaluates the function corresponding to this moments at the given point.
 important: this suppose that the normalization is dual.
 use normalize to set this before calling.
*/
double zernike::operator()(const vec &v) const
{
  if (v.length_square() > 1)
    return 0;

  zernike_r r(N);
  spherical_harmonics s(N);
  
  s_vec sp = v.spherical();
  r.eval_zr(sp.r);
  s.eval_sh(sp.theta, sp.phi);

  const std::vector<double> &z = r.get_zr();
  const std::vector<double> &sh = s.get_sh();

  int idzr = 0;
  int idz = 0;
  double sum = 0;
  for (int n2 = 0 ; n2 <= N / 2 ; n2++) {
    int idsh = 0;
    for (int l = 0 ; l <= 2 * n2 + 1 ; l++, idzr++) {
      double sum_m = 0;
      for (int m = -l ; m <= l ; m++, idsh++, idz++)
        sum_m += sh[idsh] * zm[idz];
      sum += z[idzr] * sum_m;
    }
  }
  return sum;
}

/** Call this after the moments have been computed and before using them.
 \a normalize do it for you.
 It is ok but useless to do it more than once.
*/
void zernike::finish()
{
  if (odd_clean || (N & 1) == 1 || N == 0) return;
  int n = N - 1;
  int idx = 2 * (n / 2 + 1) * (n / 2 + 2) * (2 * (n / 2) + 3) / 3;
    for (int l = 0 ; l <= N + 1 ; l++) {
      if ((l & 1) == 0)
        idx += 2 * l + 1;
      else
        for (int m = -l ; m <= l ; m++, idx++)
          zm[idx] = 0;
    }
  odd_clean = true;
}

double zernike::distance(const zernike &z) const
{
  if (norm != z.get_norm() || order() != z.order())
    return 1;
  double d = 0;
  for (size_t i = 0 ; i != zm.size(); i++) {
    const double v = fabs(zm[i] - z.get_zm()[i]);
    if (v > d)
      d = v;
  }
  return d;
}
  
zernike &zernike::operator +=(const zernike &z)
{
  if (norm != z.get_norm())
    return *this;
  const std::vector<double> &z2 = z.get_zm();
  const size_t n = std::min(zm.size(), z2.size());
  for (size_t i = 0 ; i != n ; i++)
    zm[i] += z2[i];
  variance += z.variance;
  return *this;
}

zernike operator -(const zernike &z1, const zernike &z2)
{
  if (z1.norm != z2.norm)
    return zernike();
  int n1 = z1.order();
  int n2 = z2.order();
  int n = (n1 < n2) ? n1 : n2;
  zernike z(n);
  z.norm = z1.norm;
  for (size_t i = 0 ; i<z.zm.size() ; i++)
    z.zm[i] = z1.zm[i] - z2.zm[i];
  z.finish();
  return z;
}

/** Writes a zernike in ZM format.
 Format:
 ZM     <-- just the two letters
 norm N
 n1 m1 l1 z1
 n2 m2 l2 z2
 ...
 zero entries are not output

*/
std::ostream &operator <<(std::ostream &os, const zernike &zm)
{
  os << "ZM" << std::endl;
  os << zm.get_norm() << " " << zm.order() << " " << zm.output << std::endl;
  const bool flip = flip_out(zm.output);
  const bool real = real_out(zm.output);
  for (int n = 0 ; n <= zm.order() ; n++)
    for (int l = n & 1 ; l <= n ; l+=2) {
      if (real)
        for (int m = -l ; m <= l ; m++) {
          double z = zm.get(n, l, m);
          if (flip && (m & 1))
            z = -z;
          if (z != 0)
            os << n << " " << l << " " << m << " " << z << std::endl;
        }
      else {
        const double z0 = zm.get(n, l, 0);
        if (z0 != 0)
          os << n << " " << l << " 0 " << z0 << std::endl;
        for (int m = 1 ; m <= l ; m++) {
          double r = sqrt(0.5) * zm.get(n, l, m);
          double i = - sqrt(0.5) * zm.get(n, l, -m);
          if (flip && (m & 1)) {
            r = -r;
            i = -i;
          }
          if (r != 0 || i != 0)
            os << n << " " << l << " " << m << " " << r << " " << i << std::endl;
        }
      }
    }
  return os;
}

/** Reads a zernike in ZM format.
 See format in operator >> doc.
*/
smart_input &operator >>(smart_input &is, zernike &z)
{
  std::istringstream s;
  if (!is.next_line(s)) //remove first line containing "ZM"
    return is.failed();
  if (!is.next_line(s))
    return is.failed();
  int n0;
  zm_norm norm;
  zm_output output;
  s >> norm >> n0 >> output;
  if (!s || n0 < 0)
    return is.failed();
  const bool flip = flip_out(output);
  const bool real = real_out(output);
  zernike z0(n0);
  z0.norm = norm;
  z0.output = output;
  z0.odd_clean = true;
  while(is.next_line(s)) {
    int n, l, m;
    double r, i;
    s >> n >> l >> m >> r;
    if (!s || n < 0 || n > n0 || l < 0 || l > n || (l ^ n) == 1
        || m < -l || m > l)
      return is.failed();
    if (flip && (m & 1))
      r = -r;
    if (m == 0 || real)
      z0.zm[z0.index(n, l, m)] = r;
    else if (m < 0)
      return is.failed();
    else {
      s >> i;
      if (!s)
        return is.failed();
      if (flip && (m & 1))
        i = -i;
      z0.zm[z0.index(n, l, m)] = sqrt(2) * r;
      z0.zm[z0.index(n, l, -m)] = - sqrt(2) * i;
    }
  }
  if (is.eof()) {
    z = z0;
    is.clear();
  }
  return is;
}

/** Constructor.
  @param n Maximum order needed. Should be positive.
*/
zernike_m_r::zernike_m_r(int n) :
zernike_r(n), spherical_harmonics(n), zernike(n)
{}

/** Add Zernike polynomials for the given point and weight.
  @param p The weight point to use.
*/
void zernike_m_r::add(const w_vec &p)
{
  s_vec sp = p.v.spherical();
  eval_zr(sp.r);
  eval_sh(sp.theta, sp.phi);
  add_core(zr, sh, p.weight);
}

/** Constructor.
  @param n Maximum order needed. Should be positive.
*/
zernike_m_int::zernike_m_int(int n):
zernike_int2(n), spherical_harmonics(n), zernike(n)
{}

/** Add integrated Zernike polynomials for the given point and weight.
  @param p The weight point to use.
*/
void zernike_m_int::add(const w_vec &p)
{
  s_vec sp = p.v.spherical();
  if (sp.r != 0) {
    eval_zr(sp.r, 1 / (sp.r * sp.r * sp.r));
    eval_sh(sp.theta, sp.phi);
    add_core(zr, sh, p.weight);
  }
}

/** Dummy constructor for operator >>.
 @param n The maximum order available. 
*/
rotational_invariants::rotational_invariants(int n):
N(n), norm(zm_norm::raw),
ri((N / 2 + 1) * (N  / 2 + 2) * (N / 2 + 3) / 3, 0)
{}

/** Constructor.
 @param zm The moments to use to compute the rotational invariants.
 You can normalize it first if needed.
*/
void rotational_invariants::eval_ri(const zernike &zm)
{
  norm = zm.get_norm();
  const std::vector<double> &z = zm.get_zm();
  for (int n1_2 = 0, i = 0 ; n1_2 <= N / 2 ; n1_2++)
    for (int n2_2 = 0 ; n2_2 <= n1_2 ; n2_2++)
      for (int l = 0 ; l <= 2 * n2_2 + 1 ; l++, i++) {
        int idx1 = zm.index(2 * n1_2 + (l & 1), l, 0);
        int idx2 = zm.index(2 * n2_2 + (l & 1), l, 0);
        double sum = 0;
        for (int m = -l ; m <= l ; m++)
          sum += z[idx1 + m] * z[idx2 + m];
        ri[i] = sum;
      }
}

rotational_invariants operator -(const rotational_invariants &r1, const rotational_invariants &r2)
{
  if (r1.norm != r2.norm)
    return rotational_invariants();
  int n1 = r1.order();
  int n2 = r2.order();
  int n = (n1 < n2) ? n1 : n2;
  rotational_invariants r(n);
  r.norm = r1.norm;
  for (size_t i = 0 ; i<r.ri.size() ; i++)
    r.ri[i] = r1.ri[i] - r2.ri[i];
  return r;
}

std::ostream &operator <<(std::ostream &os, const rotational_invariants &ri)
{
  os << "ZRI" << std::endl;
  os << ri.get_norm() << " " << ri.order() << std::endl;
  for (int n1 = 0 ; n1 <= ri.order() ; n1++)
    for (int n2 = (n1 & 1) ; n2 <= n1 ; n2+=2)
      for (int l = n1 & 1 ; l <= n2 ; l+=2) {
        double z = ri.get(n1, n2, l);
        if (z != 0)
          os << n1 << " " << n2 << " " << l << " " << z << std::endl;
      }
  return os;
}

smart_input &operator >>(smart_input &is, rotational_invariants &ri)
{
  std::istringstream s;
  if (!is.next_line(s)) //remove first line containing "ZRI"
    return is.failed();
  if (!is.next_line(s))
    return is.failed();
  int n0;
  zm_norm norm;
  s >> norm >> n0;
  if (!s || n0 < 0)
    return is.failed();
  rotational_invariants ri0(n0);
  ri.norm = norm;
  while(is.next_line(s)) {
    int n1, n2, l;
    double z;
    s >> n1 >> n2 >> l >> z;
    if (!s || n1 < 0 || n1 > n0 || n2 < 0 || n2 > n1 || l < 0 || l > n2
        || (n1 | n2) == 1 || (l | n2) ==1)
      return is.failed();
    ri0.ri[ri0.index(n1, n2, l)] = z;
  }
  if (is.eof()) {
    ri = ri0;
    is.clear();
  }
  return is;
}

/** Dummy constructor for operator >>.
 @param n The maximum order available. 
*/
signature_invariants::signature_invariants(int n):
N(n), norm(zm_norm::raw), si(n + 1, 0)
{}

/** Constructor.
 @param zm The moments to use to compute the signature invariants.
 You can normalize it first if needed.
*/
void signature_invariants::eval_si(const zernike &zm)
{
  rotational_invariants ri(N);
  ri.eval_ri(zm);
  norm = ri.get_norm();
  for (int n = 0 ; n <= N ; n++) {
    double sum = 0;
    for (int l = n % 2 ; l <= n ; l += 2)
      sum += ri.get(n, n, l);
    si[n] = sum;
  }
}

signature_invariants operator -(const signature_invariants &s1, const signature_invariants &s2)
{
  if (s1.norm != s2.norm)
    return signature_invariants();
  int n1 = s1.order();
  int n2 = s2.order();
  int n = (n1 < n2) ? n1 : n2;
  signature_invariants s(n);
  s.norm = s1.norm;
  for (size_t i = 0 ; i<s.si.size() ; i++)
    s.si[i] = s1.si[i] - s2.si[i];
  return s;
}

std::ostream &operator <<(std::ostream &os, const signature_invariants &si)
{
  os << "ZSI" << std::endl;
  os << si.get_norm() << " " << si.order() << std::endl;
  for (int n = 0 ; n <= si.order() ; n++) {
    double z =si.get(n);
    if (z != 0)
      os << n << " " << z << std::endl;
  }
  return os;
}

smart_input &operator >>(smart_input &is, signature_invariants &si)
{
  std::istringstream s;
  if (!is.next_line(s)) //remove first line containing "ZSI"
    return is.failed();
  if (!is.next_line(s))
    return is.failed();
  int n0;
  zm_norm norm;
  s >> norm >> n0;
  if (!s || n0 < 0)
    return is.failed();
  signature_invariants si0(n0);
  si.norm = norm;
  while(is.next_line(s)) {
    int n;
    double z;
    s >> n >> z;
    if (!s || n < 0 || n > n0)
      return is.failed();
    si0.si[n] = z;
  }
  if (is.eof()) {
    si = si0;
    is.clear();
  }
  return is;
}
