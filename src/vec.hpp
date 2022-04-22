/** \file vec.hpp
  Classes for 3D-vector, Cartesian and spherical coordinates.
  \author J. Houdayer
*/
#ifndef VEC_HPP
#define VEC_HPP

#include <vector>
#include <iostream>
#include <cmath>

class mat;
class s_vec;

/** 3D-points / vectors in Cartesian coordinates. */
class vec
{
public:
  double x, y, z;

  vec(): x(0), y(0), z(0) {}
  vec(double x0, double y0, double z0): x(x0), y(y0), z(z0) {}

  vec &operator += (const vec &v)
  { return *this = {x + v.x, y + v.y, z + v.z};  }

  vec &operator -= (const vec &v)
  { return *this = {x - v.x, y - v.y, z - v.z}; }

  vec &operator *= (double scalar)
  { return *this = {x * scalar, y * scalar, z * scalar}; }

  vec &operator /= (double scalar)
  { return *this = {x / scalar, y / scalar, z / scalar}; }

  vec &operator - ()
  { return *this = {-x, -y, -z}; }

  double length_square() const
  { return x * x + y * y + z * z; }

  double length() const
  { return sqrt(length_square()); }

  vec &normalize();
  s_vec spherical() const;
};

inline vec operator + (const vec &v1, const vec &v2)
{ return {v1.x + v2.x, v1.y + v2.y, v1.z + v2.z}; }

inline vec operator - (const vec &v1, const vec &v2)
{ return {v1.x - v2.x, v1.y - v2.y, v1.z - v2.z}; }

inline vec operator * (double scalar,  const vec &v)
{ return {v.x * scalar, v.y * scalar, v.z * scalar}; }

inline vec operator / (const vec &v, double scalar)
{ return {v.x / scalar, v.y / scalar, v.z / scalar}; }

inline double dot(const vec &v1, const vec &v2)
{ return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z ; }

vec cross(const vec &v1, const vec &v2);

std::ostream& operator<<(std::ostream& os, const vec &v);
std::istream& operator>>(std::istream& is, vec &v);

class mat
{
public:
  vec mx, my, mz;
};

inline vec operator * (const mat &m, const vec &v)
{ return {dot(m.mx, v), dot(m.my, v), dot(m.mz, v)}; }

inline mat operator + (const mat &m1, const mat &m2)
{ return { m1.mx + m2.mx, m1.my + m2.my, m1.mz + m2.mz}; }

inline mat operator - (const mat &m1, const mat &m2)
{ return { m1.mx - m2.mx, m1.my - m2.my, m1.mz - m2.mz}; }

inline mat operator * (double scalar, const mat &m)
{ return { scalar * m.mx, scalar * m.my, scalar * m.mz}; }

extern const mat mat_id;

mat diag_mat(const vec &v);
mat cross_mat(const vec &v);
mat dot_mat(const vec &v1, const vec& v2);
mat rotation_mat(const vec &v, double angle);

/** a class for representing weighted points. */
class w_vec {
public:
  double weight;
  vec v;
};

std::ostream& operator<<(std::ostream& os, const w_vec &v);
std::istream& operator>>(std::istream& is, w_vec &v);

/** 3D-points / vectors in spherical coordinates.
  The usual physical convention is used.
  theta is the colatitude and phi the longitude.
*/
class s_vec
{
public:
  double r, theta, phi;

  vec cartesian() const;
};

std::ostream& operator<<(std::ostream& os, const s_vec &v);
std::istream& operator>>(std::istream& is, s_vec &v);

#endif
