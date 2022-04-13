/** \file moments.hpp
  Functions to integrate Zernike polynomials on clouds and meshes.
*/
#ifndef MOMENTS_HPP
#define MOMENTS_HPP

#include "mesh.hpp"
#include "zernike.hpp"

void cloud_integrate(const cloud &c, zernike_m_r &zm);
void cloud_integrate(const w_cloud &c, zernike_m_r &zm);
void mesh_exact_integrate(const mesh &m, const scheme &s, zernike_m_int &zm);
double mesh_approx_integrate(const mesh &m, const scheme_selector &s, zernike_m_int &zm, double error);

#endif
