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
  "Operations are executed in the order of the command line.\n"
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
string t_help =
  "translates the shape along the given vector (3 numbers inside a quote)";
string a_help =
  "rotates the shape with the given angle in degrees and axis,\n"
  "4 numbers inside a quote (\"x y z angle\")";
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
  vector<double> tor_radius;
  vector<int> s;
  vector<double> r;
  vector<vec> t;
  vector<vec> d;
  vector<w_vec> a;
  vector<string> fs;
  double project = 0;
  vector<opt_recorder> rec;

  parser p(sh, eh);
  p.prog_name = "makeOFF";
  p.rest_arg("FILES", fs, FILES_help);
  p.rec_flag("", "cube", rec, cub_help);
  p.rec_flag("", "icosahedron", rec, ico_help);
  p.rec_flag("", "octahedron", rec, oct_help);
  p.rec_flag("", "tetrahedron", rec, tet_help);
  p.rec_flag("", "dodecahedron", rec, dod_help);
  p.rec_flag("", "sphere", rec, sph_help);
  p.rec_list_option("", "torus", "RADIUS", tor_radius, rec, tor_help);
  p.rec_list_option("s", "", "N", s, rec, s_help);
  p.rec_flag("c", "", rec, c_help);
  p.rec_list_option("r", "", "R", r, rec, r_help);
  p.rec_list_option("d", "", "FACTORS", d, rec, d_help);
  p.rec_list_option("t", "", "VEC", t, rec, t_help);
  p.rec_list_option("a", "", "VEC_ANGLE", a, rec, a_help);
  p.option("f", "", "DIGITS", digit, f_help);
  p.option("e", "", "DIGITS", digit, e_help);
  p.flag("i", "", i_help);

  p.run(argc, argv);

  if (p("f"))
    cout << fixed << setprecision(digit);
  else if (p("e"))
    cout << scientific << setprecision(digit);


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
 
  for (auto opt: rec) {
    const string &n = opt.name;
    if (n == "sphere") {
      project = (m.empty()) ? -1 : -2;
      m.add(make_icosahedron());
      continue;
    }
    if (n == "torus") {
      project = (m.empty()) ? tor_radius[opt.pos] : -2;
      m.add(make_torus(tor_radius[opt.pos]));
      continue;
    }
    if (n == "s") {
      for (int i = 0 ; i < s[opt.pos] ; i++) {
        m = m.split();
        if (project == -1)
          m.sphere_project();
        else if (project >=0)
          m.torus_project(project);
        }
      continue;
    }
    project = -2; 
    if (n == "cube")
      m.add(make_cube());
    else if (n == "icosahedron")
      m.add(make_icosahedron());
    else if (n == "octahedron")
      m.add(make_octahedron());
    else if (n == "tetrahedron")
      m.add(make_tetrahedron());
    else if (n == "dodecahedron")
      m.add(make_dodecahedron());
    else if (n == "c")
      m -= m.mass_center();
    else if (n == "r")
      m *= r[opt.pos] / m.radius();
    else if (n == "d")
      m.apply(diag_mat(d[opt.pos]));
    else if (n == "a")
      m.apply(rotation_mat(a[opt.pos].v.normalize(),
              a[opt.pos].weight * 3.141592653589793238 / 180));
    else if (n == "t")
      m += t[opt.pos]; 
  }

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
