/** \file moments.cpp
  Implementation of moments.hpp
*/

#include "moments.hpp"
#include "iotools.hpp"
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
  if (n <= 0)
    return zernike();
  cloud_sumer sumer(n);
  return parallel_collect(nt, c.points, sumer, verbose);
}

/** Compute the Zernike moments for a weighted cloud.
  Use z.orthonormalize() afterwards if needed.
*/
zernike cloud_integrate(const w_cloud &c, int n, int nt, bool verbose)
{
  if (n <= 0)
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
  if (n <= 0)
    return zernike();
  
  mesh_exact_sumer sumer(n, m, ts.get_scheme(n));
  return parallel_collect(nt, m.triangles, sumer, verbose);
}

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
  for(int n = 1 ; ; n++) {
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
}

class mesh_approx_sumer:
public zernike
{
public:
  const mesh &msh;
  const triquad_selector &sel;
  const double err;
  zernike_m_int z1, z2;
  
  mesh_approx_sumer(int n, const mesh &m, const triquad_selector &s, double e):
  zernike(n), msh(m), sel(s), err(e / m.triangles.size()), z1(n), z2(n) {}
  std::string collect(const t_mesh &i)
  {
    const triangle t = i.get_triangle(msh);
    zernike_m_int *za = &z1, *zb = &z2;
    const int res = facet_approx_integrate(t, err, sel, za, zb);
    *this += *za;
    return " order: " + std::to_string(res);
  }
  void collect(const mesh_approx_sumer &ms)
  {
    *this += ms;
  }
};


zernike mesh_approx_integrate(const mesh &m, int n, double error, const triquad_selector &ts, int nt, bool verbose)
{
  if (n <= 0)
    return zernike();
  mesh_approx_sumer sumer(n,m, ts, error);
  return parallel_collect(nt, m.triangles, sumer, verbose);
}
