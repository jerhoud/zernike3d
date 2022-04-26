/** \file mesh.cpp
  Implementation of mesh.hpp.
  \author J. Houdayer
*/

#include <unordered_map>
#include <algorithm>
#include "iotools.hpp"
#include "mesh.hpp"

#ifndef M_PI
#define M_PI 3.141592653589793238
#endif

/** Displaces the cloud by adding a vector. */
cloud &cloud::operator += (const vec &v)
{
  for (auto &pt: points)
    pt += v;
  return *this;
}

/** Displaces the cloud by subtracting a vector. */
cloud &cloud::operator -= (const vec &v)
{
  for (auto &pt: points)
    pt -= v;
  return *this;
}

/** Rescales the cloud by multiplying by a scalar. */
cloud &cloud::operator *= (double scalar)
{
  for (auto &pt: points)
    pt *= scalar;
  return *this;
}

/** Rescales the cloud by divinding by a scalar. */
cloud &cloud::operator /= (double scalar)
{
  for (auto &pt: points)
    pt /= scalar;
  return *this;
}

/** Applies a matrix to the cloud.
  @param mat The matrix to apply.
  @return The transformed cloud.
*/
cloud &cloud::apply(const mat &m)
{
  for (auto &pt: points)
    pt = m * pt;
  return *this;
}

/** Returns the center of mass. */
vec cloud::mass_center() const
{
  vec r;
  for (auto &pt: points)
    r += pt;
  return r / points.size();

}

/** returns the radius of the cloud. */
double cloud::radius() const
{
  double max = 0;
  for (auto &pt: points) {
    const double l2 = pt.length_square();
    if (max < l2)
      max = l2;
  }
  return sqrt(max);
}

/** projects all points on the unit sphere.*/
void cloud::sphere_project()
{
  for (auto &pt: points) {
    const double l = pt.length();
    if (l != 0)
      pt /= l;
  }
}

/** projects all points on the torus of outer radius one and inner radius r.*/
void cloud::torus_project(double r)
{
  const double r0 = (1 + r) / 2;
  const double r1 = (1 - r) / 2;
  for (auto &pt: points) {
    vec v = {pt.x, pt.y, 0};
    const double vl = v.length();
    if (vl != 0) {
      v *= r0 / vl;
      vec w = pt - v;
      double const wl = w.length();
      if (wl != 0) {
        w *= r1 / wl;
        pt = v + w;
      }
    }
  }
}

/** reads a point and adds it to the cloud. */
void cloud::read_point(smart_input &is)
{
  std::istringstream s;
  if (is.next_line(s)) {
    vec v;
    s >> v;
    if (s)
      add_point(v);
    else
      is.failed();
  }
}

/** Displaces the cloud by adding a vector. */
w_cloud &w_cloud::operator += (const vec &v)
{
  for (auto &pt: points)
    pt.v += v;
  return *this;
}

/** Displaces the cloud by subtracting a vector. */
w_cloud &w_cloud::operator -= (const vec &v)
{
  for (auto &pt: points)
    pt.v -= v;
  return *this;
}

/** Rescales the cloud by multiplying by a scalar. */
w_cloud &w_cloud::operator *= (double scalar)
{
  for (auto &pt: points)
    pt.v *= scalar;
  return *this;
}

/** Rescales the cloud by dividing by a scalar. */
w_cloud &w_cloud::operator /= (double scalar)
{
  for (auto &pt: points)
    pt.v /= scalar;
  return *this;
}

/** Applies a matrix to the cloud.
  @param mat The matrix to apply.
  @return The transformed cloud.
*/
w_cloud &w_cloud::apply(const mat &m)
{
  for (auto &pt: points)
    pt.v = m * pt.v;
  return *this;
}

/** Rescales the weights of the cloud. */
w_cloud &w_cloud::reweight(double scalar)
{
  for (auto &pt: points)
    pt.weight *= scalar;
  return *this;
}

/** Weighted center of mass. */
w_vec w_cloud::mass_center() const
{
  double weight = 0;
  vec r;
  for (auto &pt: points) {
    weight += pt.weight;
    r += pt.weight * pt.v;
  }
  return {weight, r / weight};
}

/** The radius of the cloud. */
double w_cloud::radius() const
{
  double max = 0;
  for (auto &pt: points) {
    const double l2 = pt.v.length_square();
    if (max < l2)
      max = l2;
  }
  return sqrt(max);
}

std::istream &operator >>(std::istream &is, t_mesh &t)
{
  size_t a, b, c;
  is >> a >> b >> c;
  if (is)
    t = { a, b, c };
  return is;  
}

std::ostream &operator <<(std::ostream &os, const t_mesh &t)
{
  os << t.i1 << " " << t.i2 << " " << t.i3;
  return os;
}

/** Volume of the mesh. */
double mesh::volume() const
{
  double sum = 0;
  for (auto &i: triangles) {
    sum += i.get_triangle(*this).volume();
  }
  return sum;
}

/** Area of the mesh. */
double mesh::area() const
{
  double sum = 0;
  for (auto &i: triangles)
    sum += i.get_triangle(*this).area();
  return sum;
}

/** Reads one triangle and adds it to the mesh. */
void mesh::read_triangle(smart_input &is)
{
  std::istringstream s;
  if (is.next_line(s)) {
    int dummy;
    t_mesh t;
    s >> dummy >> t;
    if (s)
      add_triangle(t);
    else
      is.failed();
  }
}

/** concats two meshes.*/
void mesh::add(const mesh &m)
{
  const size_t offset = points.size();
  const size_t n = triangles.size();
  points.insert(points.end(), m.points.begin(), m.points.end());
  triangles.insert(triangles.end(), m.triangles.begin(), m.triangles.end());
  for (auto i = triangles.begin() + n ; i != triangles.end() ; i++)
    i->move(offset);
}

/** adds a polygon with the given points.
 It create an additional point at the center.
*/
void mesh::add_polygon(const std::vector<size_t> &p)
{
  if (p.size() <= 2)
    return;
  vec mid;
  for (auto i: p)
    mid += points[i];
  
  mid /= p.size();
  const size_t i0 = add_point(mid);
  add_triangle({p.back(), p.front(), i0});
  for (size_t i = 0 ; i < p.size() - 1 ; i++)
    add_triangle({p[i], p[i+1], i0});
}

void mesh::add_strip(const std::vector<size_t> &l1, const std::vector<size_t> &l2, bool rev)
{
  const int s1 = l1.size() - 1;
  const int s2 = l2.size() - 1;

  int i1 = 0;
  int i2 = 0;
  for (; i1 < s1 ; i1++) {
    for (; i2 < s2 && i2 * s1 <= (i1 + 0.5) * s2 ; i2++)
      add_triangle({l2[i2], l2[i2+1], l1[i1]}, rev);
    add_triangle({l1[i1], l2[i2], l1[i1+1]}, rev);
  }
  for (; i2 < s2 && i2 * s1 <= (i1 + 0.5) * s2 ; i2++)
      add_triangle({l2[i2], l2[i2+1], l1[i1]}, rev);
}

class edge {
public:
  const size_t i1, i2;

  bool operator == (const edge &e) const
  { return (i1 == e.i1 && i2 == e.i2) || (i1 == e.i2 && i2 == e.i1); }
};

class hash_edge {
public:
  size_t operator()(const edge &e) const
  { return (std::hash<size_t>{}(e.i1) * std::hash<size_t>{}(e.i2)); }
};

typedef std::unordered_map<edge, size_t, hash_edge> middle_map;

size_t get_middle(cloud &c, middle_map &mdls, const edge &e)
{
  const auto i = mdls.find(e);
  if (i == mdls.end()) {
    size_t mi = c.add_point((c.points[e.i1] + c.points[e.i2]) / 2);
    mdls.insert(std::make_pair(e, mi));
    return mi;
  }
  else
    return i->second;
}

/** Creates a new mesh where each triangle is split in four. */
mesh mesh::split() const
{
  mesh m;
  m.points = points;
  middle_map middles;
  for (auto &t: triangles) {
    const size_t m12 = get_middle(m, middles, {t.i1, t.i2});
    const size_t m23 = get_middle(m, middles, {t.i2, t.i3});
    const size_t m31 = get_middle(m, middles, {t.i3, t.i1});
    m.add_triangle({t.i1, m12, m31});
    m.add_triangle({t.i2, m23, m12});
    m.add_triangle({t.i3, m31, m23});
    m.add_triangle({m12, m23, m31});
  }
  return m;
}

/** a class to gather data about an edge in a mesh.*/
class edge_info {
public:
  size_t count, order;
};

typedef std::unordered_map<edge, edge_info, hash_edge> edge_map;

void edge_register(edge_map &m, const edge &e)
{
  auto i = m.find(e);
  if (i == m.end())
    i = m.insert({e, {0, 0}}).first;
  i->second.count++;
  i->second.order += (e.i1 <= e.i2) ? 1 : -1;
}

/** builds a report about the mesh, to check for consistency.*/
edge_report mesh::edges() const
{
  edge_map m;
  for (auto &t: triangles) {
    edge_register(m, {t.i1, t.i2});
    edge_register(m, {t.i2, t.i3});
    edge_register(m, {t.i3, t.i1});
  }
  edge_report r = {m.size(), 0, 0};
  for (auto &i: m) {
    if (i.second.count > 2)
      r.strange++;
    else if (i.second.order != 0)
      r.border++;
  }
  return r;
}

/** Reads a mesh in OFF format. */
smart_input &operator >>(smart_input &is, mesh &m)
{
  mesh m0;
  std::istringstream s;
  if (!is.next_line(s)) //remove first line containing "OFF"
    return is;
  if (!is.next_line(s))
    return is;
  size_t n_points, n_faces, dummy;
  s >> n_points >> n_faces >> dummy;
  if (!s)
    return is.failed();
  for (size_t i = 0 ; i < n_points ; i++)
    m0.read_point(is);
  for (size_t i = 0 ; i < n_faces ; i++)
    m0.read_triangle(is);
  if (is)
    m = m0;
  return is;
}

/** Writes a mesh in OFF format. */
std::ostream &operator <<(std::ostream &os, const mesh &m)
{
  os << "OFF" << std::endl;
  os << m.points.size() << " " << m.triangles.size() << " " << 0 << std::endl;
  for (auto &pt: m.points)
    os << pt << std::endl;
  for (auto &t: m.triangles)
    os << 3 << " " << t << std::endl;
  return os;
}

/** builds a mesh representing a cube with 12 facets.*/
mesh make_cube()
{
  mesh m;
  m.points =
    {{1, 1, 1}, {-1, 1, 1}, {-1, -1, 1}, {1, -1, 1},
     {1, 1, -1}, {-1, 1, -1}, {-1, -1, -1}, {1, -1, -1}};
  m.triangles =
    {{0, 1, 2}, {2, 3, 0}, {0, 4, 1}, {1, 4, 5},
     {1, 5, 2}, {2, 5, 6}, {2, 6, 3}, {3, 6, 7},
     {3, 7, 0}, {0, 7, 4}, {4, 6, 5}, {7, 6, 4}};
  m *= 1 / sqrt(3);
  return m;
}

/** builds a mesh representing a regular tetrahedron with 4 facets.*/
mesh make_tetrahedron()
{
  mesh m;
  m.points =
    {{1, 1, 1}, {-1, -1, 1}, {1, -1, -1}, {-1, 1, -1}};
  m.triangles =
    {{0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {1, 3, 2}};
  m *= 1 / sqrt(3);
  return m;
}

/** builds a mesh representing a regular icosahedron with 20 facets.*/
mesh make_icosahedron()
{
  const double g = (1 + sqrt(5)) / 2;
  mesh m;
  m.points =
    {{1, 0, g}, {-1, 0, -g}, {-1, 0, g}, {1, 0, -g},
     {g, 1, 0}, {-g, -1, 0}, {g, -1, 0}, {-g, 1, 0},
     {0, g, 1}, {0, -g, -1}, {0, g, -1}, {0, -g, 1}};
  m.triangles =
    {{0, 2, 11}, {1, 10, 3}, {0, 8, 2}, {1, 3, 9},
     {0, 6, 4}, {1, 5, 7}, {3, 4, 6}, {2, 7, 5},
     {4, 10, 8}, {5, 9, 11}, {7, 8, 10}, {6, 11, 9},
     {0, 4, 8}, {1, 9, 5}, {0, 11, 6}, {1, 7, 10},
     {2, 8, 7}, {3, 6, 9}, {2, 5, 11}, {3, 10, 4}};
  m *= 1 / sqrt(2 + g);
  return m;
}

/** builds a mesh representing a regular octahedron with 8 facets.*/
mesh make_octahedron()
{
  mesh m;
  m.points =
    {{0, 0, 1}, {1, 0, 0}, {0, 1, 0},
     {-1, 0, 0}, {0, -1, 0}, {0, 0, -1}};
  m.triangles =
    {{0, 1, 2}, {0, 2, 3}, {0, 3, 4}, {0, 4, 1},
     {5, 1, 4}, {5, 2, 1}, {5, 3, 2}, {5, 4, 3}};
  return m;
}

std::vector<size_t> add_circle(mesh &m, vec v, double a)
{
  std::vector<size_t> l;
  double r = hypot(v.x, v.y);
  size_t n = ceil(2 * M_PI * r / a);
  if (r > 0 && n < 4)
    n = 4;
  if (n == 0)
    l.push_back(m.add_point(v));
  else {
    v = rotation_mat({0, 0, 1}, - M_PI / n) * v;
    mat rot = rotation_mat({0, 0, 1}, 2 * M_PI / n);
    while (n--) {
      l.push_back(m.add_point(v));
      v = rot * v;
    }
    l.push_back(l.front());
  }
  return l;
}

/** builds a mesh representing a torus.
  The outer radius is 1 and the inner radius is r.
  The mesh has at least 123 facets.*/
mesh make_torus(double r)
{
  mesh m;
  const double a = (1 - r) / 2;
  const double h = sqrt(3)/2 * a;
  const vec u = {0, 0, h};
  const std::vector<size_t> l1 = add_circle(m, {r, 0, 0}, a);
  const vec p2 = (r == 0) ? vec(0.25, 0, 0) : (r + a / 2) / r * m.points[l1[0]];
  const std::vector<size_t> l2u = add_circle(m, p2 + u, a);
  m.add_strip(l1, l2u, false);
  const std::vector<size_t> l2d = add_circle(m, p2 - u, a);
  m.add_strip(l1, l2d, true);
  const vec p3 = (r + 3 * a / 2) / (r + a / 2) * (m.points[l2u[0]] - u);
  const std::vector<size_t> l3u = add_circle(m, p3 + u, a);
  m.add_strip(l2u, l3u, false);
  const std::vector<size_t> l3d = add_circle(m, p3 - u, a);
  m.add_strip(l2d, l3d, true);
  const vec p4 = (m.points[l3u[0]] - u) / (r + 3 * a / 2);
  const std::vector<size_t> l4 = add_circle(m, p4, a);
  m.add_strip(l3u, l4, false);
  m.add_strip(l3d, l4, true);

  return m;
}

/** builds a mesh representing a regular dodecahedron with 60 facets.*/
mesh make_dodecahedron()
{
  const double g = (1 + sqrt(5)) / 2;
  const double h = g - 1;
  mesh m;
  m.points =
    {{1, 1, 1}, {-1, 1, 1}, {-1, -1, 1}, {1, -1, 1},
     {1, 1, -1}, {-1, 1, -1}, {-1, -1, -1}, {1, -1, -1},
     {0, h, g}, {0, -h, g}, {0, h, -g}, {0, -h, -g},
     {g, 0, h}, {g, 0, -h}, {-g, 0, h}, {-g, 0, -h},
     {h, g, 0}, {-h, g, 0}, {h, -g, 0}, {-h, -g, 0}};
  const std::vector<size_t> facets[12] =
    {{8, 9, 3, 12, 0},   {8, 1, 14, 2, 9},
     {10, 11, 6, 15, 5}, {10, 4, 13, 7, 11},
     {12, 13, 4, 16, 0}, {12, 3, 18, 7, 13},
     {14, 15, 6, 19, 2}, {14, 1, 17, 5, 15},
     {16, 17, 1, 8, 0},  {16, 4, 10, 5, 17},
     {18, 19, 6, 11, 7}, {18, 3, 9, 2, 19}
    };
  for (auto &f: facets)
    m.add_polygon(f);
  m *= 1 / sqrt(3);
  return m;
}

/** a class to represent a tetrahedron for marching_tetrahedra.*/
class tetrahedron
{
public:
  int n1, n2, n3;
};

#define TETRA(a, b, c) {a, b, c}, {a ^ 1, c ^ 1, b ^ 1}

const tetrahedron tetras[24] = { // the 24 tetrahedron around a point of the lattice
  TETRA(0, 2, 6), TETRA(0, 6, 8), TETRA(0, 8, 4), TETRA(0, 4, 2),
  TETRA(10, 2, 4), TETRA(10, 4, 7), TETRA(10, 7, 9), TETRA(10, 9, 2),
  TETRA(12, 2, 9), TETRA(12, 9, 5), TETRA(12, 5, 6), TETRA(12, 6, 2)
};



typedef u_int16_t sig_t;

inline bool out(sig_t sig, int n)
{ return sig & (1 << n); }

int count_bit(sig_t sig)
{
  int n = 0;
  for (; sig; sig &= sig - 1)
    n++;
  return n;
}

class mt_node {
public:
  size_t pos;
  sig_t signature, rest;
  std::vector<sig_t> groups;
  vec vertex[14];

  size_t operator()(int n) const
  {
    const sig_t mask = 1 << n;
    size_t vpos = pos;
    for (auto g: groups) {
      if (g & mask)
        return vpos;
      ++vpos;
    }
    vpos += count_bit(rest & (mask - 1));
    return vpos;
  }
};

class neighbors
{
public:
  std::vector<int> ngh;
  std::vector<std::vector<int>> nngh;

  neighbors(int dx, int dy, int dz);

  int operator[](int n) const
  { return ngh[n]; }

  int operator()(int n1, int n2) const
  {
    int r = nngh[n1][n2];
    return r;
  }

  bool connected(int n1, int n2) const
  { return nngh[n1][n2] != -1; }

  const std::vector<sig_t> component(sig_t sig, int &n_in) const;
};

#define NEIGH(x) x, -(x)

neighbors::neighbors(int dx, int dy, int dz)
{
  ngh = {
    NEIGH(-2*dz + dy + dx),
    NEIGH(-dz), NEIGH(dx-dz), NEIGH(dy-dz), NEIGH(dx+dy-dz),
    NEIGH(-dy), NEIGH(-dx)
  };
  nngh.assign(14, std::vector<int>(14, -1));
  for (int n1 = 1 ; n1 < 14 ; n1++)
    for (int n2 = 0 ; n2 < n1 ; n2++) {
      int d = ngh[n2] - ngh[n1];
      for (int n = 0 ; n < 14 ; n++)
        if (d == ngh[n]) {
          nngh[n1][n2] = n;
          nngh[n2][n1] = n ^ 1;  
          break;
        }
    }
}

const std::vector<sig_t> neighbors::component(sig_t sig, int &n_in) const
{
  int num = 0;
  std::vector<int> mark(14, 0);
  std::vector<int> stack;
  std::vector<sig_t> groups;
  n_in = 0;

  for (int i = 0 ; i < 14 ; i++) {
    if (mark[i] != 0)
      continue;
    bool t = out(sig, i);
    mark[i] = ++num;
    sig_t c = 1 << i;
    stack.push_back(i);
    while (!stack.empty()) {
      int j = stack.back();
      stack.pop_back();
      for (int k = 0 ; k < 14 ; k++)
        if (mark[k] == 0 && out(sig, k) == t && connected(j, k)) {
          mark[k] = num;
          c |= 1 << k;
          stack.push_back(k);
        }
    }
    if (!t)
      n_in++;
    else
      groups.push_back(c);
  }
  return groups;
}

vec pos_vertex(const mt_coord &sx, const mt_coord &sy, const mt_coord &sz, size_t idx)
{
  const size_t layer_sz = sy.maxN() * sx.maxN();
  const size_t dbl_layer_sz = 2 * layer_sz - sx.maxN() - 1;
  int s = 0;
  const size_t nz = idx / dbl_layer_sz;
  idx -= nz * dbl_layer_sz;
  if (idx > layer_sz) {
    s = 1;
    idx -= layer_sz;
  }
  const size_t ny = idx / sx.maxN();
  const size_t nx = idx - ny * sx.maxN();
  return vec {sx.pos(nx, s), sy.pos(ny, s), sz.pos(nz, s)}; 
}

mesh marching_tetrahedra(const mt_coord &sx, const mt_coord &sy, const mt_coord &sz,
                      std::function<double(const vec &)> f, double thresh, bool regularized, bool verbose)
{
  const size_t N = sz.maxN()*((2 * sy.maxN() - 1) * sx.maxN() - 1);
  std::vector<float> val(N+1);
  std::vector<size_t> in_node;
  { // Phase 1
    if (verbose)
      std::cerr << "Phase 1/4, function evaluation" << std::endl;
    progression prog(N + sz.maxN(), verbose);
    int s = 0;
    size_t idx = 0;
    for (int nz = 0 ; nz < sz.maxN() * 2 ; nz++, s = 1 - s) {
      for (int ny = 0 ; ny < sy.maxN() - s ; ny++) // remove one line every two layers 
        for (int nx = 0 ; nx < sx.maxN() ; nx++, idx++) {
          const vec pos { sx.pos(nx, s), sy.pos(ny, s), sz.pos(nz / 2, s) };
          double v = f(pos) - thresh;
          if (v > 0 && (nx <= 0 || nx >= sx.maxN() - 1
                    || ny <= 0 || ny >= sy.maxN() - 1))
            v = -1;
          val[idx] = v;
          if (v > 0)
            in_node.push_back(idx);
          prog.progress();
        }
      if (s == 1) // remove one point every two layers
        idx--;
    }
  }

  if (in_node.empty())
    return mesh();
  
  std::unordered_map<size_t, mt_node> surface;
  neighbors ngh(1, sx.maxN(), sx.maxN() * sy.maxN());
  { // Phase 2
    if (verbose)
      std::cerr << "Phase 2/4, computing vertices" << std::endl;
    progression prog(in_node.size(), verbose);
    std::unordered_map<sig_t, std::vector<sig_t>> components;
    std::unordered_map<sig_t, bool> collapsable;

    for (auto idx: in_node) { // go through all points
      mt_node node;
      prog.progress();
      const vec pos = pos_vertex(sx, sy, sz, idx);
      const double v = val[idx];
      sig_t sig = 0;
      for (int i = 0 ; i < 14 ; i++) {// check all neighbors
        const size_t idx_ngh = idx + ngh[i];
        const double v_ngh = val[idx_ngh];
        if (v_ngh <= 0) {
          sig |= 1 << i;
          const vec pos_ngh = pos_vertex(sx, sy, sz, idx_ngh);
          node.vertex[i] = pos - v / (v_ngh - v) * (pos_ngh - pos);
        }
      }

      if (sig == 0)
        continue;

      node.signature = sig;
      int n = count_bit(sig);
      if (regularized && n >= 2 && n <=12) {
        auto it = components.find(sig);
        if (it != components.end())
          node.groups = it->second;
        else {
          int n_in;
          std::vector<sig_t> compo = ngh.component(sig, n_in);
          for (auto c:compo) {
            bool collapse;
            auto it = collapsable.find(c);
            if (it != collapsable.end())
              collapse = it->second;
            else {
              ngh.component(c, n_in);
              collapse = n_in == 1;
              collapsable[c] = collapse;
            }
            if (collapse) {
              node.groups.push_back(c);
              sig ^= c;
            }
          }
          components[node.signature] = node.groups;
        }
      }
      node.rest = sig;
      surface.insert({idx, node});
    }
  }

  val.clear();
  in_node.clear();

  std::vector<size_t> reject(1, -1);
  int count = 0;
  while (true) {
    count++;
    mesh m;
    { // Phase 3
      if (verbose)
        std::cerr << "Phase 3/4, grouping vertices" << std::endl;
      progression prog(surface.size(), verbose);
      size_t idx = 0, r_idx = 0;
      for (auto &v: surface) { // go through all inside points on surface
        mt_node &node = v.second;
        node.pos = m.points.size();
        sig_t dispatch = 0;
        for (size_t i = 0 ; i < node.groups.size() ; i++) {
          sig_t g = node.groups[i];
          if (idx == reject[r_idx]) {
            r_idx++;
            dispatch |= g;
            node.groups.erase(node.groups.begin() + i--);
          }
          else {
            vec sum {0, 0, 0};
            for (int j = 0 ; j < 14 ; j++)
              if (out(g, j))
                sum += node.vertex[j];
            sum /= count_bit(g);
            m.points.push_back(sum);
          }
          idx++;
        }
        if ((node.rest | dispatch) != 0)
          for (int i = 0 ; i < 14 ; i++) {
            if (out(node.rest, i)) {
              if (idx == reject[r_idx])
                r_idx++;
              m.points.push_back(node.vertex[i]);
              idx++;
            }
            else if(out(dispatch, i))
              m.points.push_back(node.vertex[i]);
          }
        node.rest |= dispatch;
        prog.progress();
      }
    }

    { // Phase 4
      if (verbose)
        std::cerr << "Phase 4/4, computing triangles" << std::endl;
      progression prog(surface.size(), verbose);
      for (auto &v: surface) {// go through all inside points on surface
        prog.progress();
        const size_t idx = v.first;
        const mt_node &node = v.second;
        for (auto &tetra: tetras) { // go through all possible tetra at this point
          int n1 = tetra.n1;
          int n2 = tetra.n2;
          int n3 = tetra.n3;
          bool i1 = !out(node.signature, n1);
          bool i2 = !out(node.signature, n2);
          bool i3 = !out(node.signature, n3);
          int n = i1 + i2 + i3; // numbers of neighbors inside
          if (n == 3) // tetra is completly inside, do nothing
            continue;
          if (n == 0) // only current point inside
            m.add_triangle({node(n1), node(n2), node(n3)});
          else {      // 2 or 3 points inside
            while (!i1 || i3) { // rotate points, after that inside points are n1 and possibly n2
              bool ti = i1;
              i1 = i2;
              i2 = i3;
              i3 = ti;
              int tn = n1;
              n1 = n2;
              n2 = n3;
              n3 = tn;
            }
            if (ngh[n1] < 0) // stop if tetra has already been done by n1
              continue;
            const mt_node &node1 = surface[idx + ngh[n1]];
            if (n == 1) { // tetra has 2 points in and 2 out
              size_t p02 = node(n2);
              size_t p03 = node(n3);
              size_t p12 = node1(ngh(n1, n2));
              size_t p13 = node1(ngh(n1, n3));
              if ((m.points[p02] - m.points[p13]).length_square()
                  > (m.points[p12] - m.points[p03]).length_square()) { // cut through the smallest diagonal
                m.add_triangle({p03, p13, p12});
                m.add_triangle({p03, p12, p02});
              }
              else {
                m.add_triangle({p13, p12, p02});
                m.add_triangle({p13, p02, p03});
              }
            }

            else { // tetra has 3 points in and 1 out
              if (ngh[n2] < 0) // stop if tetra has been done by n2
                continue;
              const mt_node &node2 = surface[idx + ngh[n2]];        
              m.add_triangle({node(n3),
                              node1(ngh(n1, n3)),
                              node2(ngh(n2, n3))});
            }
          }
        }
      }
    }

    if (count == 1) {
      if (verbose)
        std::cerr << "Checking for bad bonds" << std::endl;
      edge_map em;
      for (auto &t: m.triangles) {
        edge_register(em, {t.i1, t.i2});
        edge_register(em, {t.i2, t.i3});
        edge_register(em, {t.i3, t.i1});
      }
      reject.clear();
      for (auto &i: em) {
        if (i.second.count > 2)
          reject.push_back(i.first.i1);
      }
      if (verbose)
        std::cerr << "Found: " << reject.size()  << std::endl;
      
      reject.push_back(-1);
      std::sort(reject.begin(), reject.end());
      std::unique(reject.begin(), reject.end());
    }
    else
      return m;
  }
}
