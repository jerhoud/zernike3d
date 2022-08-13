/** \file moments.cpp
  Implementation of moments.hpp
*/

#include "moments.hpp"
#include "iotools.hpp"

/** Computes the Zernike moments for a cloud.
  Use z.orthonormalize() afterwards if needed.
*/
void cloud_integrate(const cloud &c, zernike_m_r &z, bool verbose)
{
  z.reset_zm();
  progression prog(c.points.size(), verbose);
  for (auto &pt: c.points) {
    z.add({1, pt});
    prog.progress();
  }
}

/** Compute the Zernike moments for a weighted cloud.
  Use z.orthonormalize() afterwards if needed.
*/
void cloud_integrate(const w_cloud &c, zernike_m_r &z, bool verbose)
{
  z.reset_zm();
  progression prog(c.points.size(), verbose);
  for (auto &pt: c.points) {
    z.add(pt);
    prog.progress();
  }
}

/** Computes the Zernike moments of a mesh.
  This suppose that the order sought is not larger than the order of the integration scheme.
*/
zernike mesh_exact_integrate(const mesh &m, int n, const triquad_selector &ts, bool verbose)
{
  if (n <= 0)
    return zernike();
  
  zernike_m_int zm(n);
  const triquad_scheme &s = ts.get_scheme(n);
  progression prog(m.triangles.size(), verbose);
  for (auto &i: m.triangles) {
    const triangle t =  i.get_triangle(m);
    s.integrate(t, zm, 3 * t.volume());
    prog.progress();
  }
  zm.error = 1e-13;
  return zm;
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
      za->error = 1e-14;
      return s.order;
    }
    else if (err < error) {
      za->error = err;
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
      za->error = err;
      return -n;
    }
  }
}

zernike mesh_approx_integrate(const mesh &m, int n, double error, const triquad_selector &ts, bool verbose)
{
  if (n <= 0)
    return zernike();
  zernike z(n);
  zernike_m_int z1(n), z2(n);
  zernike_m_int *za = &z1, *zb = &z2;
  error /= m.triangles.size();
  progression prog(m.triangles.size(), verbose);
  for (auto &t: m.triangles) {
    const int res = facet_approx_integrate(t.get_triangle(m), error, ts, za, zb);
    z += *za;
    prog.progress(" order: " + std::to_string(res));
  }
  return z;
}
