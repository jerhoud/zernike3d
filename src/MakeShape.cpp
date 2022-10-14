/** \file MakeShape.cpp
  A standalone program to create simple shapes in OFF format.
  \author J. Houdayer
*/

#include <map>
#include "version.hpp"
#include "arg_parse.hpp"
#include "mesh.hpp"

using namespace std;
using namespace argparse;

string sh =
  "Creates or modifies shapes in OFF format.\n"
  "Outputs it on standard output in OFF format.";
string eh =
  "Operations (except global options) are executed in the order of the command line.\n"
  "One can read or create multiple shapes into one.";
string ex =
  "MakeShape --cube                        Creates a cube.\n"
  "MakeShape --sphere -s5 -o sphere.off    Creates a sphere with 20480 facets (20 * 4^5) and saves it to file.\n"
  "MakeShape -l shape.off -cr 1            Reads file shape.off, centers it and sets its outer radius to 1.\n"
  "MakeShape -l shape1.off -l shape2.off   Combines the two given shape into one.\n"
  "MakeShape -l shape.off -i               Gives information on shape.off";
string l_help = "adds file FILE in OFF format to the current shape";
string o_help = "save current shape in OFF format to file FILE";
string c_help = "centers the shape around the center of mass";
string r_help = "rescales the shape to set its outer radius to R";
string s_help = "multiplies the number of facets by four N times,\n"
                "projects the new points for the sphere and the torus";
string t_help =
  "translates the shape along the given vector: -t \"dx dy zy\"";
string a_help =
  "rotates the shape with the given angle in degrees and axis: -a \"x y z angle\"";
string e_help =
  "applies a diagonal matrix (i.e. expands along the axis): -e \"fx fy fz\"";
string d_help = "number of significant digits printed in the output (default is 6)";
string i_help = "shows informations about the shape and stops";
string cub_help = "adds a cube with 12 facets";
string ico_help = "adds a regular icosahedron with 20 facets";
string oct_help = "adds a regular octahedron with 8 facets";
string tet_help = "adds a regular tetrhedron with 4 facets";
string dod_help = "adds a regular dodecahedron with 60 facets";
string sph_help = "adds a sphere with 20 facets";
string tor_help = "adds a torus with the given inner radius\n"
                  "the number of facets increases with the inner radius starting at 123";
string mem_help = "memorizes the current shape under the given NAME";
string rec_help = "adds the shape memorized under NAME to the current shape";
string clear_help = "starts with a fresh empty shape";
string bad_output_msg = "Cannot open output file: ";

int main (int argc, char *argv[])
{
  int digit = 6;
  vector<string> string_dat;
  vector<double> double_dat;
  vector<int> int_dat;
  vector<vec> vec_dat;
  vector<w_vec> wvec_dat;
  double project = 0;
  vector<opt_recorder> rec;

  parser p(sh, eh, ex);
  p.prog_name = "MakeShape";

  p.group("Global options");
  p.option("d", "digits", "DIGITS", digit, d_help);

  p.group("Load / save options");
  p.rec_list_option("o", "save", "FILE", string_dat, rec, o_help);
  p.rec_list_option("l", "load", "FILE", string_dat, rec, l_help);
  p.rec_list_option("", "memorize", "NAME", string_dat, rec, mem_help);
  p.rec_list_option("", "recall", "NAME", string_dat, rec, rec_help);
  p.rec_flag("", "clear", rec, clear_help);

  p.group("Shape options (shapes are created with radius 1)");
  p.rec_flag("", "cube", rec, cub_help);
  p.rec_flag("", "icosahedron", rec, ico_help);
  p.rec_flag("", "octahedron", rec, oct_help);
  p.rec_flag("", "tetrahedron", rec, tet_help);
  p.rec_flag("", "dodecahedron", rec, dod_help);
  p.rec_flag("", "sphere", rec, sph_help);
  p.rec_list_option("", "torus", "RADIUS", double_dat, rec, tor_help);

  p.group("Transformation options");
  p.rec_list_option("s", "", "N", int_dat, rec, s_help);
  p.rec_flag("c", "", rec, c_help);
  p.rec_list_option("r", "", "R", double_dat, rec, r_help);
  p.rec_list_option("e", "expand", "FACTORS", vec_dat, rec, e_help);
  p.rec_list_option("t", "", "VEC", vec_dat, rec, t_help);
  p.rec_list_option("a", "", "VEC_ANGLE", wvec_dat, rec, a_help);

  p.group("Miscellaneous");
  p.rec_flag("i", "", rec, i_help);

  p.run(argc, argv);

  if (digit <= 0)
    digit = 1;

  map<string, mesh> memo;
  mesh m;
  
  if (rec.empty() || rec.back().name!="o") {
    int pos = string_dat.size();
    string_dat.push_back("-");
    rec.push_back(opt_recorder("o", pos));
  }

  for (auto opt: rec) {
    const string &n = opt.name;

    if (n == "l") {
      mesh m0;
      string err = read_file(string_dat[opt.pos], m0);
      if (!err.empty())
        p.die(err);
      m.add(m0);
    }
    else if (n == "o") {
      smart_output out(string_dat[opt.pos]);
      if (!out)
        p.die(bad_output_msg + n + " (" + strerror(errno) + ")");
      out << setprecision(digit) << m;
    }
    else if (n == "sphere") {
      project = (m.empty()) ? -1 : -2;
      m.add(make_icosahedron());
      continue;
    }
    else if (n == "torus") {
      project = (m.empty()) ? double_dat[opt.pos] : -2;
      m.add(make_torus(double_dat[opt.pos]));
      continue;
    }
    else if (n == "s") {
      for (int i = 0 ; i < int_dat[opt.pos] ; i++) {
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
    
    if (n == "c")
      m -= m.mass_center();
    else if (n == "r")
      m *= double_dat[opt.pos] / m.radius();
    else if (n == "e")
      m.apply(diag_mat(vec_dat[opt.pos]));
    else if (n == "a")
      m.apply(rotation_mat(wvec_dat[opt.pos].v.normalize(),
              wvec_dat[opt.pos].weight * 3.141592653589793238 / 180));
    else if (n == "t")
      m += vec_dat[opt.pos];
    
    if (n == "memorize")
      memo[string_dat[opt.pos]] = m;
    else if (n == "recall")
      m.add(memo[string_dat[opt.pos]]);
    else if (n == "clear")
      m = mesh();
    else if (n == "i") {
      const vec mc = m.mass_center();
      m -= mc;
      edge_report r = m.edges();
      cout << setprecision(digit);
      cout << "Number of vertices: " << m.points.size() << endl;
      cout << "Number of facets: " << m.triangles.size() << endl;
      cout << "Number of edges: " << r.count << endl;
      cout << "V - E + F = "
          << int(m.points.size() - r.count + m.triangles.size()) << endl;
      if (r.border != 0)
        cout << "There are " << r.border << " border edges." << endl;
      if (r.bad_orient != 0)
        cout << "There are " << r.bad_orient << " ill oriented edges." << endl;
      if (r.strange != 0)
        cout << "There are " << r.strange << " edges connected to more than 2 facets." << endl;
      cout << endl;
      cout << "Center of mass: " << mc << endl;
      cout << "Radius from center of mass: " << m.radius() << endl;
      cout << "Area: " << m.area() << endl;
      cout << "Volume: " << m.volume() << endl;
      exit(0);
    }
  }
}
