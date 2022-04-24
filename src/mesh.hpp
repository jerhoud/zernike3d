/** \file mesh.hpp
  Classes for clouds and meshes.
  \author J. Houdayer
*/
#ifndef MESH_HPP
#define MESH_HPP

#include <functional>
#include "iotools.hpp"
#include "triangle.hpp"

/** A cloud of points. */
class cloud
{
public:
  std::vector<vec> points; /**< The points. */

  /** returns true if there are no points in the cloud.*/
  bool empty() const
  { return points.empty(); }
 
  /** Adds a point to the cloud. */
  size_t add_point(const vec &pt)
  {
    points.push_back(pt);
    return points.size() - 1;
  }

  void read_point(smart_input &is);
  cloud &operator += (const vec &v);
  cloud &operator -= (const vec &v);
  cloud &operator *= (double scalar);
  cloud &operator /= (double scalar);
  cloud &apply(const mat &m);
  vec mass_center() const;
  double radius() const;
  void sphere_project();
  void torus_project(double r);
};

/** A cloud of weighted points. */
class w_cloud
{
public:
  std::vector<w_vec> points; /**< The points. */

  /** Adds a point to the cloud. */
  size_t add_point(const w_vec &pt)
  { 
    points.push_back(pt);
    return points.size() - 1;
  }

  w_cloud &operator += (const vec &v);
  w_cloud &operator -= (const vec &v);
  w_cloud &operator *= (double scalar);
  w_cloud &operator /= (double scalar);
  w_cloud &reweight(double scalar);
  w_cloud &apply(const mat &m);
  w_vec mass_center() const;
  double radius() const;
};

/** A truple of ints to represent a triangle by indices.
  Used by mesh.
*/
class t_mesh
{
public:
  size_t i1, i2, i3;

  void move(size_t offset)
  { i1 += offset; i2 += offset; i3 += offset; }

  triangle get_triangle(const cloud &cld) const
  { return {cld.points[i1], cld.points[i2], cld.points[i3]}; }

  bool collapsed() const
  { return (i1 == i2) || (i2 == i3) || (i3 == i1); }
};

std::istream &operator >>(std::istream &, t_mesh &t);
std::ostream &operator <<(std::ostream &, const t_mesh &t);

/** a class to store information about a mesh.*/
class edge_report {
public:
  size_t count, border, strange;
};

/** A triangular mesh. */
class mesh: public cloud
{
public:
  std::vector<t_mesh> triangles; /**< the triangles. */

  double volume() const;
  double area() const;

  /** Adds a triangle to the mesh, by giving the indices of the vertices.
   @param t the triangle to add
   @param rev whether to reverse the orientation of the triangle
  */
  void add_triangle(const t_mesh &t, bool rev=false)
  { 
    if (t.collapsed())
      return;
    if (rev)
      triangles.push_back({t.i1, t.i3, t.i2});
    else
      triangles.push_back(t);
  }
  void add_polygon(const std::vector<size_t> &p);
  void add_strip(const std::vector<size_t> &l1, const std::vector<size_t> &l2, bool rev);

  void read_triangle(smart_input &is);
  void add(const mesh &m);
  mesh split() const;
  edge_report edges() const;
};

smart_input &operator >>(smart_input &is, mesh &m);
std::ostream &operator <<(std::ostream &os, const mesh &m);

mesh make_cube();
mesh make_tetrahedron();
mesh make_icosahedron();
mesh make_octahedron();
mesh make_dodecahedron();
mesh make_torus(double r);

/** a class to store the resolution on a axis for marching_tetrahedra.*/
class mt_coord
{
public:
  double min, max;
  int N;

  double step() const
  { return (max - min) / N; }

  double pos(int n, int shift) const
  { return min + (0.5 * shift + (n - 1)) * step(); }

  int maxN() const
  { return N + 3;}
};

mesh marching_tetrahedra(const mt_coord &sx, const mt_coord &sy, const mt_coord &sz,
                      std::function<double(const vec &)> f, double thresh, bool regularized, bool verbose = false);

#endif
