/** \file mesh.hpp
  Classes for clouds and meshes.
  \author J. Houdayer
*/
#ifndef MESH_HPP
#define MESH_HPP

#include "triangle.hpp"

/** A cloud of points. */
class cloud
{
public:
  std::vector<vec> points; /**< The points. */

  /** Adds a point to the cloud. */
  void add_point(const vec &pt)
  { points.push_back(pt); }

  void read_point(std::istream &is);
  cloud &operator += (const vec &v);
  cloud &operator -= (const vec &v);
  cloud &operator *= (double scalar);
  cloud &operator /= (double scalar);
  cloud &apply(const mat &m);
  vec mass_center() const;
  double radius() const;
  void sphere_project(double r);
};

/** A cloud of weighted points. */
class w_cloud
{
public:
  std::vector<w_vec> points; /**< The points. */

  /** Adds a point to the cloud. */
  void add_point(const w_vec &pt)
  { points.push_back(pt); }

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
  int i1, i2, i3;

  void move(int offset)
  { i1 += offset; i2 += offset; i3 += offset; }

  triangle get_triangle(const cloud &cld) const
  { return {cld.points[i1], cld.points[i2], cld.points[i3]}; }

  bool collapsed() const
  { return (i1 == i2) || (i2 == i3) || (i3 == i1); }
};

std::istream &operator >>(std::istream &, t_mesh &t);
std::ostream &operator <<(std::ostream &, const t_mesh &t);

class edge_report {
public:
  int count, border, strange;
};

/** A triangular mesh. */
class mesh: public cloud
{
public:
  std::vector<t_mesh> triangles; /**< the triangles. */

  double volume() const;
  double area() const;

  /** Adds a triangle to the mesh, by giving the indices of the vertices. */
  void add_triangle(const t_mesh &t)
  { 
    if (!t.collapsed())
      triangles.push_back(t);
  }

  void read_triangle(std::istream &is);
  void add(const mesh &m);
  mesh split() const;
  edge_report edges() const;
};

std::istream &operator >>(std::istream &is, mesh &m);
std::ostream &operator <<(std::ostream &os, const mesh &m);

mesh make_cube(double r);
mesh make_tetrahedron(double r);
mesh make_icosahedron(double r);
mesh make_octahedron(double r);
mesh make_dodecahedron(double r);

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

class mt_node
{
public:
  vec pos;
  double val;
  bool inside, collapsed;
  int signature;
  int vertex;

  int operator()(int n) const;
};

class marching_tetrahedra
{
public:
  const mt_coord szx, szy, szz;
  const double threshold;

  template<class T> marching_tetrahedra(const mt_coord &sx, const mt_coord &sy, const mt_coord &sz, T &f, double thresh):
  szx(sx), szy(sy), szz(sz), threshold(thresh), node(2*szx.maxN()*szy.maxN()*szz.maxN())
  {
    int idx = 0;
    int s = 0;
    for (int nz = 0 ; nz < szz.maxN() * 2 ; nz++, s = 1 - s) {
      for (int ny = 0 ; ny < szy.maxN() - s ; ny++) // remove one line every two layers 
        for (int nx = 0 ; nx < szx.maxN() ; nx++, idx++) {
          mt_node &nd = node[idx];
          nd.pos = { szx.pos(nx, s), szy.pos(ny, s), szz.pos(nz / 2, s) };
          nd.val = f(nd.pos);
          nd.inside = nd.val > thresh 
                        && nx > 0 && nx < szx.maxN() - 1
                        && ny > 0 && ny < szy.maxN() - 1;
          if (nd.inside) // select all inside points
            in_node.push_back(idx);
        }
      if (s == 1) // remove one point every two layers
        idx--;
    }
  }

  mesh build();

protected:
  std::vector<mt_node> node;
  std::vector<int> in_node;
};


#endif
