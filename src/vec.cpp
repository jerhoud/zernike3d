/** \file vec.cpp
  Implementation of vec.hpp.
  \author J. Houdayer
*/
#include <iomanip>
#include <algorithm>
#include "iotools.hpp"
#include "vec.hpp"

/** The cross product. */
vec cross(const vec &v1, const vec &v2)
{ return {v1.y * v2.z - v1.z * v2.y,
          v1.z * v2.x - v1.x * v2.z,
          v1.x * v2.y - v1.y * v2.x};
}

/** Normalizes the vector to unit length.
  @return The normalized vector.
*/
vec &vec::normalize()
{
  double l = length();
  if (l != 0)
    *this /= l;
  return *this;
}

mat operator * (const mat &m1, const mat &m2)
{
  vec vx = {m2.mx.x, m2.my.x, m2.mz.x};
  vec vy = {m2.mx.y, m2.my.y, m2.mz.y};
  vec vz = {m2.mx.z, m2.my.z, m2.mz.z};
  return {
    {dot(m1.mx, vx), dot(m1.mx, vy), dot(m1.mx, vz)},
    {dot(m1.my, vx), dot(m1.my, vy), dot(m1.my, vz)},
    {dot(m1.mz, vx), dot(m1.mz, vy), dot(m1.mz, vz)},
  };
}

const mat mat_id = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

mat diag_mat(const vec &v)
{ return {{ v.x, 0, 0}, {0, v.y, 0}, {0, 0, v.z}}; }

mat cross_mat(const vec &v)
{ return {{ 0, -v.z, v.y }, { v.z, 0, -v.x }, { -v.y, v.x, 0 }}; }

mat dot_mat(const vec &v1, const vec& v2)
{ return { v1.x * v2, v1.y * v2, v1.z * v2}; }

mat rotation_mat(const vec &v, double angle)
{
  double s = sin(angle);
  double c = cos(angle);
  return c * mat_id + (1 - c) * dot_mat(v, v) + s * cross_mat(v);
}

/** Spherical coordinates representation. */
s_vec vec::spherical() const
{
  double r = length();
  return {r, (r == 0) ? 0 : acos(z / r), atan2(y, x)};
}

/** output a vector.
  Looks like this "x y z".
*/
std::ostream& operator<<(std::ostream& os, const vec &v)
{
  os << v.x << " " << v.y << " " << v.z;
  return os;
}

/** input a vector.
  Should look like this "x y z".
*/
std::istream& operator>>(std::istream& is, vec &v)
{
  double x, y, z;
  is >> x >> y >> z;
  if (is)
    v = {x, y, z};
  return is;
}

/** Cartesian coordinates representation. */
vec s_vec::cartesian() const
{
  double rs = r * sin(theta);
  return {rs * cos(phi), rs * sin(phi), r * cos(theta)};
}

std::ostream& operator<<(std::ostream& os, const s_vec &v)
{
  int p = os.precision();
  os << std::setw(p + 9) << v.r
     << std::setw(p + 9) << v.theta
     << std::setw(p + 9) << v.phi;
  return os;
}

std::istream& operator>>(std::istream& is, s_vec &v)
{
  double r, theta, phi;
  is >> r >> theta >> phi;
  if (!is || r < 0)
    return failed(is);
  v = {r, theta, phi};
  return is;
}
