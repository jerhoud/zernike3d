/** \file moments.cpp
  Implementation of moments.hpp
  \author J. Houdayer
*/

#include "moments.hpp"
#include "parallel.hpp"

class cloud_sumer:
public zernike_m_r
{
public:
  cloud_sumer(int n): zernike_m_r(n) {}
  std::string collect(const vec &v) {
    add({1, v});
    variance += 1e-30;
    return "";
  }
  std::string collect(const w_vec &v) {
    add(v);
    variance += 1e-30;
    return "";
  }
  void collect(const cloud_sumer &cs) {
    *this += cs;
  }
};

/** Computes the Zernike moments for a cloud.
  Use z.orthonormalize() afterwards if needed.
*/
zernike cloud_integrate(const cloud &c, int n, int nt, bool verbose)
{
  if (n < 0)
    return zernike();
  cloud_sumer sumer(n);
  return parallel_collect(nt, c.points, sumer, verbose);
}

/** Compute the Zernike moments for a weighted cloud.
  Use z.orthonormalize() afterwards if needed.
*/
zernike cloud_integrate(const w_cloud &c, int n, int nt, bool verbose)
{
  if (n < 0)
    return zernike();
  cloud_sumer sumer(n);
  return parallel_collect(nt, c.points, sumer, verbose);
}

class mesh_exact_sumer:
public zernike_m_int
{
public:
  const mesh &msh;
  const triquad_scheme &sch;
  
  mesh_exact_sumer(int n, const mesh &m, const triquad_scheme &s): zernike_m_int(n), msh(m), sch(s) {}
  std::string collect(const t_mesh &i)
  {
    const triangle t = i.get_triangle(msh);
    sch.integrate(t, *this, 3 * t.volume());
    variance += 1e-28;
    return "";
  }
  void collect(const mesh_exact_sumer &ms)
  {
    *this += ms;
  }
};

/** Computes the Zernike moments of a mesh.
  This suppose that the order sought is not larger than the order of the integration scheme.
*/
zernike mesh_exact_integrate(const mesh &m, int n, const triquad_selector &ts, int nt, bool verbose)
{
  if (n < 0)
    return zernike();
  
  mesh_exact_sumer sumer(n, m, ts.get_scheme(n));
  return parallel_collect(nt, m.triangles, sumer, verbose);
}

/** Maximum number of splits in four of a facet in facet_approx_integrate
  (each split multiplies the cost by four). */
const int max_split = 4;

/** Integrates over one facet with increasing quadrature orders, then splits, until two successive results differ by less than error.
  @return The order of the rule used, -n after n splits, 0 if the precision was not reached.
*/
int facet_approx_integrate(const triangle &t, double error, const triquad_selector &ts, zernike_m_int *&za, zernike_m_int *&zb)
{
  double err = 0;
  za->reset_zm();
  const double w = 3 * t.volume();
  for (auto &s: ts.schemes) {
    zb->reset_zm();
    s.integrate(t, *zb, w);
    zb->finish();
    err = za->distance(*zb);
    std::swap(za, zb);
    if (za->order() <= s.order) {
      za->variance = 1e-28;
      return s.order;
    }
    else if (err < error) {
      za->variance = err * err;
      return s.order;
    }
  }

  const triquad_scheme &s = ts.schemes.back(); 
  for (int n = 1 ; n <= max_split ; n++) {
    zb->reset_zm();
    s.integrate(t, *zb, w, n);
    zb->finish();
    err = za->distance(*zb);
    std::swap(za, zb);
    if (err < error) {
      za->variance = err * err;
      return -n;
    }
  }
  za->variance = err * err;
  return 0;
}

/** Error allowed on each facet by mesh_approx_integrate to reach the global error.
  If every facet error is below it, the total stays below error even if the facet errors are correlated.
*/
double facet_tolerance(double error, const mesh &m)
{ return error / m.triangles.size(); }

class mesh_approx_sumer:
public zernike
{
public:
  const mesh &msh;
  const triquad_selector &sel;
  const double err;
  zernike_m_int z1, z2;
  size_t failures = 0; /**< Number of facets where the precision was not reached. */

  mesh_approx_sumer(int n, const mesh &m, const triquad_selector &s, double e):
  zernike(n), msh(m), sel(s), err(facet_tolerance(e, m)), z1(n), z2(n) {}
  std::string collect(const t_mesh &i)
  {
    const triangle t = i.get_triangle(msh);
    zernike_m_int *za = &z1, *zb = &z2;
    const int res = facet_approx_integrate(t, err, sel, za, zb);
    if (res == 0)
      failures++;
    *this += *za;
    return " order: " + std::to_string(res);
  }
  void collect(const mesh_approx_sumer &ms)
  {
    *this += ms;
    failures += ms.failures;
  }
};

/** Computes the Zernike moments of a mesh with the adaptive method.
  @param failures Set to the number of facets where the precision was not reached.
*/
zernike mesh_approx_integrate(const mesh &m, int n, double error, const triquad_selector &ts, size_t &failures, int nt, bool verbose)
{
  failures = 0;
  if (n < 0)
    return zernike();
  mesh_approx_sumer sumer(n,m, ts, error);
  const mesh_approx_sumer result = parallel_collect(nt, m.triangles, sumer, verbose);
  failures = result.failures;
  return result;
}
