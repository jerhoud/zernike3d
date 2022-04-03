/** \file mesh.cpp
  Implementation of mesh.hpp.
  \author J. Houdayer
*/
#include <unordered_map>
#include "iotools.hpp"
#include "mesh.hpp"


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
    double l2 = pt.length_square();
    if (max < l2)
      max = l2;
  }
  return sqrt(max);
}

/** projects all points on the centered sphere of given radius.*/
void cloud::sphere_project(double r)
{
  for (auto &pt: points)
    pt *= r / pt.length();
}

/** reads a point and adds it to the cloud. */
void cloud::read_point(std::istream &is)
{
  std::istringstream s;
  if (next_line(is, s)) {
    vec v;
    s >> v;
    if (s)
      add_point(v);
    else
      failed(is);
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
    double l2 = pt.v.length_square();
    if (max < l2)
      max = l2;
  }
  return sqrt(max);
}

std::istream &operator >>(std::istream &is, t_mesh &t)
{
  int a, b, c;
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
void mesh::read_triangle(std::istream &is)
{
  std::istringstream s;
  if (next_line(is, s)) {
    int dummy;
    t_mesh t;
    s >> dummy >> t;
    if (s)
      add_triangle(t);
    else
      failed(is);
  }
}

void mesh::add(const mesh &m)
{
  int offset = points.size();
  int n = triangles.size();
  points.insert(points.end(), m.points.begin(), m.points.end());
  triangles.insert(triangles.end(), m.triangles.begin(), m.triangles.end());
  for (auto i = triangles.begin() + n ; i != triangles.end() ; i++)
    i->move(offset);
}

void mesh::add_polygon(const std::vector<int> &p)
{
  if (p.size() <= 2)
    return;
  vec mid;
  for (auto i: p)
    mid += points[i];
  
  mid /= p.size();
  int i0 = points.size();
  add_point(mid);
  add_triangle({p.back(), p.front(), i0});
  for (size_t i = 0 ; i < p.size() - 1 ; i++)
    add_triangle({p[i], p[i+1], i0});
}

class edge {
public:
  int i1, i2;

  bool operator == (const edge &e) const
  { return (i1 == e.i1 && i2 == e.i2) || (i1 == e.i2 && i2 == e.i1); }
};

class hash_edge {
public:
  size_t operator()(const edge &e) const
  { return (std::hash<int>{}(e.i1) * std::hash<int>{}(e.i2)); }
};

typedef std::unordered_map<edge, int, hash_edge> middle_map;

int get_middle(cloud &c, middle_map &mdls, const edge &e)
{
  auto i = mdls.find(e);
  if (i == mdls.end()) {
    int mi = c.points.size();
    mdls.insert(std::make_pair(e, mi));
    c.add_point((c.points[e.i1] + c.points[e.i2]) / 2);
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
    int m12 = get_middle(m, middles, {t.i1, t.i2});
    int m23 = get_middle(m, middles, {t.i2, t.i3});
    int m31 = get_middle(m, middles, {t.i3, t.i1});
    m.add_triangle({t.i1, m12, m31});
    m.add_triangle({t.i2, m23, m12});
    m.add_triangle({t.i3, m31, m23});
    m.add_triangle({m12, m23, m31});
  }
  return m;
}

class edge_info {
public:
  int count, order;
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

edge_report mesh::edges() const
{
  edge_map m;
  for (auto &t: triangles) {
    edge_register(m, {t.i1, t.i2});
    edge_register(m, {t.i2, t.i3});
    edge_register(m, {t.i3, t.i1});
  }
  edge_report r = {(int)m.size(), 0, 0};
  for (auto &i: m) {
    if (i.second.count % 2)
      r.border++;
    else if (i.second.order != 0)
      r.strange++;
  }
  return r;
}

/** Reads a mesh in OFF format.
  Adds the shape to the given mesh. */
std::istream &operator >>(std::istream &is, mesh &m)
{
  mesh m0;
  std::istringstream s;
  if (!next_line(is, s)) //remove first line containing "OFF"
    return is;
  if (!next_line(is, s))
    return is;
  int n_points, n_faces, dummy;
  s >> n_points >> n_faces >> dummy;
  if (!s || n_points < 0 || n_faces < 0)
    return failed(is);
  for (int i = 0 ; i < n_points ; i++)
    m0.read_point(is);
  for (int i = 0 ; i < n_faces ; i++)
    m0.read_triangle(is);
  if (is)
    m.add(m0);
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

mesh make_cube(double r)
{
  mesh m;
  m.points =
    {{1, 1, 1}, {-1, 1, 1}, {-1, -1, 1}, {1, -1, 1},
     {1, 1, -1}, {-1, 1, -1}, {-1, -1, -1}, {1, -1, -1}};
  m.triangles =
    {{0, 1, 2}, {2, 3, 0}, {0, 4, 1}, {1, 4, 5},
     {1, 5, 2}, {2, 5, 6}, {2, 6, 3}, {3, 6, 7},
     {3, 7, 0}, {0, 7, 4}, {4, 6, 5}, {7, 6, 4}};
  m *= r / sqrt(3);
  return m;
}

mesh make_tetrahedron(double r)
{
  mesh m;
  m.points =
    {{1, 1, 1}, {-1, -1, 1}, {1, -1, -1}, {-1, 1, -1}};
  m.triangles =
    {{0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {1, 3, 2}};
  m *= r / sqrt(3);
  return m;
}

mesh make_icosahedron(double r)
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
  m *= r / sqrt(2 + g);
  return m;
}

mesh make_octahedron(double r)
{
  mesh m;
  m.points =
    {{0, 0, 1}, {1, 0, 0}, {0, 1, 0},
     {-1, 0, 0}, {0, -1, 0}, {0, 0, -1}};
  m.triangles =
    {{0, 1, 2}, {0, 2, 3}, {0, 3, 4}, {0, 4, 1},
     {5, 1, 4}, {5, 2, 1}, {5, 3, 2}, {5, 4, 3}};
  m *= r;
  return m;
}

#define PT(a, b, c, d, e) {a, b, c}, {a, c, d}, {a, d, e}

mesh make_dodecahedron(double r)
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
  const std::vector<int> facets[12] =
    {{8, 9, 3, 12, 0},   {8, 1, 14, 2, 9},
     {10, 11, 6, 15, 5}, {10, 4, 13, 7, 11},
     {12, 13, 4, 16, 0}, {12, 3, 18, 7, 13},
     {14, 15, 6, 19, 2}, {14, 1, 17, 5, 15},
     {16, 17, 1, 8, 0},  {16, 4, 10, 5, 17},
     {18, 19, 6, 11, 7}, {18, 3, 9, 2, 19}
    };
  for (auto &f: facets)
    m.add_polygon(f);
  m *= r / sqrt(3);
  return m;
}

class tetrahedron
{
public:
  int n1, n2, n3;
};

#define TETRA(a, b, c) {a, b, c}, {a ^ 1, c ^ 1, b ^ 1}

const tetrahedron tetras[24] = {
  TETRA(0, 2, 6), TETRA(0, 6, 8), TETRA(0, 8, 4), TETRA(0, 4, 2),
  TETRA(10, 2, 4), TETRA(10, 4, 7), TETRA(10, 7, 9), TETRA(10, 9, 2),
  TETRA(12, 2, 9), TETRA(12, 9, 5), TETRA(12, 5, 6), TETRA(12, 6, 2)
};

#define NEIGH(x) x, -(x)

inline bool out(int sig, int n)
{ return sig & (1 << n); }

inline int count_bit(int sig)
{
  int n = 0;
  for (; sig; sig &= sig - 1)
    n++;
  return n;
}

int mt_node::operator()(int n) const
{
  if (collapsed)
    return vertex;
  else
    return vertex + count_bit(signature & ((1 << n) - 1));
}

class neighbors
{
public:
  std::vector<int> ngh;
  std::vector<std::vector<int>> nngh;
  std::vector<bool> collapse;

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

  int component(int sig) const;
};

#define CORNER(x) x, ((x) & 0x1555) << 1 | ((x) & 0x2aaa) >> 1

int corners[14] = {
  CORNER(0x155), CORNER(0x1405), CORNER(0x2411), CORNER(0x1841),
  CORNER(0x2901), CORNER(0x694), CORNER(0x1264)
};

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
  collapse.assign(16384, false);
  for (int sig = 0 ; sig < 16384 ; sig++) {
    int n = count_bit(sig);
    if (n < 2 || n > 10)
      continue;
    bool ok = false;
    for (auto c: corners)
      if ((c & sig) == 0) {
        ok = true;
        break;
      }
    if (ok && component(sig) == 2)
      collapse[sig] = true;
  }
}

int neighbors::component(int sig) const
{
  int num = 0;
  std::vector<int> mark(14, 0);
  std::vector<int> stack;

  for (int i = 0 ; i < 14 ; i++) {
    if (mark[i] != 0)
      continue;
    bool t = out(sig, i);
    mark[i] = num++;
    stack.push_back(i);
    while (!stack.empty()) {
      int j = stack.back();
      stack.pop_back();
      for (int k = 0 ; k < 14 ; k++)
        if (mark[k] == 0 && out(sig, k) == t && connected(j, k)) {
          mark[k] = num;
          stack.push_back(k);
        }
    }
  }
  return num;
}

marching_tetrahedra::marching_tetrahedra(const mt_coord &sx, const mt_coord &sy, const mt_coord &sz,
                      std::function<double(const vec &)> f, double thresh):
szx(sx), szy(sy), szz(sz), threshold(thresh), func(f), node(2*szx.maxN()*szy.maxN()*szz.maxN())
{
  int idx = 0;
  int s = 0;
  for (int nz = 0 ; nz < szz.maxN() * 2 ; nz++, s = 1 - s) {
    for (int ny = 0 ; ny < szy.maxN() - s ; ny++) // remove one line every two layers 
      for (int nx = 0 ; nx < szx.maxN() ; nx++, idx++) {
        mt_node &nd = node[idx];
        nd.pos = { szx.pos(nx, s), szy.pos(ny, s), szz.pos(nz / 2, s) };
        nd.val = func(nd.pos);
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


mesh marching_tetrahedra::build()
{
  neighbors ngh(1, szx.maxN(), szx.maxN() * szy.maxN());
  std::vector<int> surface;
  mesh m;

  for (auto idx: in_node) { // go through all points
    mt_node &nd0 = node[idx];
    nd0.vertex = m.points.size();

    vec sum;
    for (int n = 0 ; n < 14 ; n++) {// check all neighbors
      mt_node &nd1 = node[idx + ngh[n]];
      if(!nd1.inside) {
        nd0.signature |= 1 << n;
        vec v = nd0.pos
          + (threshold - nd0.val) / (nd1.val - nd0.val)
          * (nd1.pos - nd0.pos);
        sum += v;
        m.points.push_back(v);
      }
    }

    if (nd0.signature == 0)
      continue;

    surface.push_back(idx);
    nd0.collapsed = ngh.collapse[nd0.signature];

    if (nd0.collapsed) {
      int n = m.points.size() - nd0.vertex;
      sum /= n;
      m.points.resize(nd0.vertex);
      m.points.push_back(sum);
    }
  }

  for (auto idx: surface) {// go through all inside points on surface
    mt_node &nd0 = node[idx];
    for (auto &tetra: tetras) { // go through all possible tetra at this point
      int n1 = tetra.n1;
      int n2 = tetra.n2;
      int n3 = tetra.n3;
      bool i1 = !out(nd0.signature, n1);
      bool i2 = !out(nd0.signature, n2);
      bool i3 = !out(nd0.signature, n3);
      int n = i1 + i2 + i3; // numbers of neighbors inside
      if (n == 3) // tetra is completly inside, do nothing
        continue;
      if (n == 0) // only current point inside
        m.add_triangle({nd0(n1), nd0(n2), nd0(n3)});
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
        mt_node &nd1 = node[idx + ngh[n1]];
        if (n == 1) { // tetra has 2 points in and 2 out
          int p02 = nd0(n2);
          int p03 = nd0(n3);
          int p12 = nd1(ngh(n1, n2));
          int p13 = nd1(ngh(n1, n3));
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
          mt_node &nd2 = node[idx + ngh[n2]];        
          m.add_triangle({nd0(n3),
                          nd1(ngh(n1, n3)),
                          nd2(ngh(n2, n3))});
        }
      }
    }
  }
  return m;
}
