/** \file moments.hpp
  Functions to integrate Zernike polynomials on clouds and meshes.
*/
#ifndef MOMENTS_HPP
#define MOMENTS_HPP

#include "mesh.hpp"
#include "zernike.hpp"

void cloud_integrate(const cloud &c, zernike_m_r &zm, bool verbose = false);
void cloud_integrate(const w_cloud &c, zernike_m_r &zm, bool verbose = false);
zernike mesh_exact_integrate(const mesh &m, int n, const triquad_selector &ts, const gauss_selector &gs, bool verbose = false);
zernike mesh_approx_integrate(const mesh &m, int n, double error, const triquad_selector &ts, const gauss_selector &gs, bool verbose = false);

#endif
