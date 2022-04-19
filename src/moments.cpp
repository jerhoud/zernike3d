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
  Use z.orthonormalize() afterwards if needed.
*/
void mesh_exact_integrate(const mesh &m, const scheme &s, zernike_m_int &z)
{
  z.reset_zm();
  for (auto &i: m.triangles) {
    const triangle t =  i.get_triangle(m);
    s.integrate(t, z, 3 * t.volume());
  }
}

double mesh_approx_integrate(const mesh &m, const scheme_selector &ss, zernike_m_int &z, double error)
{
  const int N = z.zernike_radial::N;
  double err = 1. / 0.;
  double error_t = 0;
  zernike_m_int z1(N), z2(N);
  zernike_m_int *za = &z1, *zb = &z2;
  z.reset_zm();
  error /= m.triangles.size();
  for (auto &i: m.triangles) {
    za->reset_zm();
    const triangle t =  i.get_triangle(m);
    const double w = 3 * t.volume();  
    for (auto &s: ss.schemes) {
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

    const scheme &s = ss.schemes.back(); 
    for(int n = 1 ; err >= error ; n++) {
      zb->reset_zm();
      s.integrate(t, *zb, w, n);
      zb->finish();
      err = za->distance(*zb);
      std::swap(za, zb);
      if (err < error)
        std::cerr << "splitting depth: " << n << ", error: " << err << "\n";
    }
    error_t += err;
    z += *za;
  }
  return error_t;
}