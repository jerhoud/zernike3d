/** \file makeOFF.cpp
  A standalone program to create simple shapes in OFF format.
*/
#include "arg_parse.hpp"
#include "iotools.hpp"
#include "mesh.hpp"

using namespace std;
using namespace argparse;

string sh =
  "Creates or modifies a shape in OFF format.\n"
  "Outputs it on standard output in OFF format.";
string eh =
  "Operations are executed in the order presented above.\n"
  "One can read or create multiple shapes into one.\n"
  "If no file is given and no shapes created, then reads standard input.\n"
  "The shapes created by makeOFF have originally a radius equal to one.\n\n"
  "Examples:\n"
  "makeOFF --cube: creates a cube.\n"
  "makeOFF -cr 1 file: reads shape from file, centers it and fixes its radius to 1.\n"
  "makeOFF -v \"1 0 0\" file: translates the given form with the given vector.\n"
  "makeOFF file1 file2 file3: creates a form containing the 3 forms given in files.\n"
  "makeOFF --sphere -s4: creates a sphere with 5120 facets.\n";
string FILES_help =  "reads files FILES in OFF format";
string c_help = "centers the shape around the center of mass";
string r_help = "rescales the shape to fix the radius to R";
string s_help = "multiplies the number of facets by four N times,\n"
                "projects the new points for the sphere and the torus";
string v_help =
  "translates the shape along the given vector (3 numbers inside a quote)\n"
  "when used with -a, sets the rotation axis instead";
string a_help =
  "rotates the shape with the given angle in degrees\n"
  "the axis is given with -v (default is vertical axis)";
string d_help =
  "applies a diagonal matrix (i.e. expands along the axis)\n"
  "the three expansion factors are given inside a quote";
string f_help = "numerical precision in fixed notation";
string e_help = "numerical precision in scientific notation";
string i_help = "shows informations about the shape and stops";
string cub_help = "creates a cube with 12 facets";
string ico_help = "creates a regular icosahedron with 20 facets";
string oct_help = "creates a regular octahedron with 8 facets";
string tet_help = "creates a regular tetrhedron with 4 facets";
string dod_help = "creates a regular dodecahedron with 60 facets";
string sph_help = "creates a sphere with 24 facets";
string tor_help = "creates a torus with the given inner radius\n"
                  "the number of facets increases with radius starting at 123";

int main (int argc, char *argv[])
{
  int digit = 6;
  int s = 0;
  double r = 0;
  vec v = {0, 0, 1};
  vec d;
  double a = 0;
  double tor_radius;
  vector<string> fs;
  bool project = true;

  parser p(sh, eh);
  p.prog_name = "makeOFF";
  p.rest_arg("FILES", fs, FILES_help);
  p.flag("", "cube", cub_help);
  p.flag("", "icosahedron", ico_help);
  p.flag("", "octahedron", oct_help);
  p.flag("", "tetrahedron", tet_help);
  p.flag("", "dodecahedron", dod_help);
  p.flag("", "sphere", sph_help);
  p.option("", "torus", "RADIUS", tor_radius, tor_help);
  p.option("s", "", "N", s, s_help);
  p.flag("c", "", c_help);
  p.option("r", "", "R", r, r_help);
  p.option("d", "", "FACTORS", d, d_help);
  p.option("v", "", "VEC", v, v_help);
  p.option("a", "", "ANG", a, a_help);
  p.option("f", "", "DIGITS", digit, f_help);
  p.option("e", "", "DIGITS", digit, e_help);
  p.flag("i", "", i_help);

  p.run(argc, argv);

  if (fs.empty() &&
      p.none({"cube", "icosahedron", "octahedron", "tetrahedron",
              "dodecahedron", "sphere", "torus"}))
    fs.push_back("-");

  mesh m;
  for (auto &filename: fs) {
    string err = read_file(filename, m);
    if (!err.empty())
      p.die(err);
  }
 
  if (p("cube"))
    m.add(make_cube());
  if (p("icosahedron"))
    m.add(make_icosahedron());
  if (p("octahedron"))
    m.add(make_octahedron());
  if (p("tetrahedron"))
    m.add(make_tetrahedron());
  if (p("dodecahedron"))
    m.add(make_dodecahedron());
  if (p("sphere")) {
    project = m.empty();
    m.add(make_icosahedron());
  }
  if (p("torus")) {
    project = m.empty();
    m.add(make_torus(tor_radius));
  }

  for (int i = 0 ; i < s ; i++) {
    m = m.split();
    if (project) {
      if (p("sphere"))
        m.sphere_project();
      else if (p("torus"))
        m.torus_project(tor_radius);
    }
  }

  if (p("c"))
    m -= m.mass_center();

  if (p("r"))
    m *= r / m.radius();

  if (p("d"))
    m.apply(diag_mat(d));
  
  if (p("a"))
    m.apply(rotation_mat(v.normalize(), a * 3.141592653589793238 / 180));
  else if (p("v"))
    m += v;

  if (p("f"))
    cout << fixed << setprecision(digit);
  else if (p("e"))
    cout << scientific << setprecision(digit);

  if (p("i")) {
    vec mc = m.mass_center();
    m -= mc;
    edge_report r = m.edges();
    cout << "Number of vertices: " << m.points.size() << endl;
    cout << "Number of facets: " << m.triangles.size() << endl;
    cout << "Number of edges: " << r.count << endl;
    cout << "V - E + F = "
         << int(m.points.size() - r.count + m.triangles.size()) << endl;
    if (r.border != 0)
      cout << "There are " << r.border << " border edges." << endl;
    if (r.strange != 0)
      cout << "There are " << r.strange << " ill oriented edges." << endl;
    cout << endl;
    cout << "Center of mass: " << mc << endl;
    cout << "Radius from center of mass: " << m.radius() << endl;
    cout << "Area: " << m.area() << endl;
    cout << "Volume: " << m.volume() << endl;
  }
  else
    cout << m;
}
