/** \file rzm.cpp.
  A standalone program to compute a mesh from zernike moments.

  It reads Zernike moments from a ZM file (as produced by zm)
  and build a mesh (in OFF format) using marching tetrahedra.
*/

#include "iotools.hpp"
#include "arg_parse.hpp"
#include "mesh.hpp"
#include "zernike.hpp"

using namespace std;
using namespace argparse;

string sh =
  "Computes a mesh from Zernike moments.\n"
  "Input should be in ZM format as produced by zm.";
string eh = "";
string v_help = "Outputs additional informations including progression bars";
string f_help = "Numerical precision in fixed notation";
string e_help = "Numerical precision in scientific notation";
string t_help = "Threshold value which separates the inside\n"
                "from the outside (default is 1/2)";
string N_help = "The maximum order of Zernike moments to use (if available)";
string RES_help = "Résolution of the mesh (i.e. number of intervals between -1 and 1)";
string FILE_help = "Reads FILE in ZM format (default is standard input)";
string die_N_msg = "N must be positive.";
string warn_N_msg = "N larger than maximum moment available. Adapting.";

parser p(sh, eh);
int N = 0;
int digit = 6;
int res;
double thresh = 0.5;

int main (int argc, char *argv[])
{
  string filename = "-";
  p.prog_name = "rzm";

  p.flag("v", "verbose", v_help);
  p.option("f", "", "DIGITS", digit, f_help);
  p.option("e", "", "DIGITS", digit, e_help);
  p.option("t", "threshold", "THRESH", thresh, t_help);
  p.arg("N", N, N_help);
  p.arg("RES", res, RES_help);
  p.opt_arg("FILE", filename, FILE_help);

  p.run(argc, argv);

  if (digit < 0)
    digit = 6;
  if (p("f"))
    cout << fixed << setprecision(digit);
  else if (p("e"))
    cout << scientific << setprecision(digit);

  if (N < 0)
    p.die(die_N_msg);

  cout << "# Produced by rzm (" << p.version_text << ") from file: " << filename << endl;
  cout << "# Date: " << now() << endl;

  zernike zm;
  string err = read_file(filename, zm, p("v"));
  if (!err.empty())
    p.die(err);
  if (N > zm.order()) {
    p.warn(warn_N_msg);
    N = zm.order();
  }
  zm.normalize(zm_norm::dual);

  zernike_build f(zernike(N, zm));
  marching_tetrahedra mt({-1, 1, res}, {-1, 1, res}, {-1, 1, res}, f, thresh, p("v"));
  cout << mt.build();
}