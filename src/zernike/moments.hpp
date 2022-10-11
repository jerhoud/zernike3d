/** \file moments.hpp
  Functions to integrate Zernike polynomials on clouds and meshes.
  \author J. Houdayer
*/

#ifndef MOMENTS_HPP
#define MOMENTS_HPP

#include "mesh.hpp"
#include "zernike.hpp"

zernike cloud_integrate(const cloud &c, int n, int nt = 1, bool verbose = false);
zernike cloud_integrate(const w_cloud &c, int n, int nt = 1, bool verbose = false);
zernike mesh_exact_integrate(const mesh &m, int n, const triquad_selector &ts, int nt = 1, bool verbose = false);
zernike mesh_approx_integrate(const mesh &m, int n, double error, const triquad_selector &ts, int nt = 1, bool verbose = false);

#endif
