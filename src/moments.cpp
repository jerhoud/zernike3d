/** \file moments.cpp
  Implementation of moments.hpp
*/

#include "moments.hpp"

/** Computes the Zernike moments for a cloud.
  Use z.orthonormalize() afterwards if needed.
*/
void cloud_integrate(const cloud &c, zernike_m_r &z)
{
  z.reset_zm();
  for (auto &pt: c.points)
    z.add({1, pt});
}

/** Compute the Zernike moments for a weighted cloud.
  Use z.orthonormalize() afterwards if needed.
*/
void cloud_integrate(const w_cloud &c, zernike_m_r &z)
{
  z.reset_zm();
  for (auto &pt: c.points)
    z.add(pt);
}

/** Computes the Zernike moments of a mesh.
  This suppose that the order sought is not larger than the order of the integration scheme.
*/
zernike mesh_exact_integrate(const mesh &m, int n, const triquad_selector &ts, const gauss_selector &gs)
{
  if (n <= 0)
    return zernike();
  
  zernike_m_int_alt zm(n, gs);
  //zernike_m_int zm(n);
  const triquad_scheme &s = ts.get_scheme(n);
  for (auto &i: m.triangles) {
    const triangle t =  i.get_triangle(m);
    s.integrate(t, zm, 3 * t.volume());
  }
  return zm;
}

void facet_approx_integrate(const triangle &t, double error, const triquad_selector &ts, zernike_m_int_alt *&za, zernike_m_int_alt *&zb)
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
    if (err < error) {
      std::cerr << "scheme order: " << s.order << ", error: " << err << "\n";
      break;
    }
  }

  const triquad_scheme &s = ts.schemes.back(); 
  for(int n = 1 ; err >= error ; n++) {
    zb->reset_zm();
    s.integrate(t, *zb, w, n);
    zb->finish();
    err = za->distance(*zb);
    std::swap(za, zb);
    if (err < error)
      std::cerr << "splitting depth: " << n << ", error: " << err << "\n";
  }

  za->error = err;
}

zernike mesh_approx_integrate(const mesh &m, int n, double error, const triquad_selector &ts, const gauss_selector &gs)
{
  if (n <= 0)
    return zernike();
  zernike z(n);
  zernike_m_int_alt z1(n, gs), z2(n, gs);
  zernike_m_int_alt *za = &z1, *zb = &z2;
  error /= m.triangles.size();
  for (auto &t: m.triangles) {
    facet_approx_integrate(t.get_triangle(m), error, ts, za, zb);
    z += *za;
  }
  return z;
}
