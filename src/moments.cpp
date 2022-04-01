/** \file moments.cpp
  Implementation of moments.hpp
*/
#include "moments.hpp"

/** Computes the Zernike moments for a cloud.
  Use z.orthonormalize() afterwards if needed.
*/
void cloud_integrate(const cloud &c, zernike_m_r &z)
{
  z.reset();
  for (auto &pt: c.points)
    z.add({1, pt});
}

/** Compute the Zernike moments for a weighted cloud.
  Use z.orthonormalize() afterwards if needed.
*/
void cloud_integrate(const w_cloud &c, zernike_m_r &z)
{
  z.reset();
  for (auto &pt: c.points)
    z.add(pt);
}

/** Computes the Zernike moments of a mesh.
  This suppose that the order sought is not larger than the order of the integration scheme.
  Use z.orthonormalize() afterwards if needed.
*/
void mesh_simple_integrate(const mesh &m, const scheme &s, zernike_m_int &z)
{
  z.reset();
  for (auto &i: m.triangles) {
    triangle t =  i.get_triangle(m);
    s.integrate(t, z, 3 * t.volume());
  }
}
