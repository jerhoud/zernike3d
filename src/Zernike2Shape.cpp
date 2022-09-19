/** \file Zernike2Shape.cpp.
  A standalone program to compute a mesh from zernike moments.

  It reads Zernike moments from a ZM file (as produced by Shape2Zernike)
  and build a mesh (in OFF format) using marching tetrahedra.
*/

#include "parallel.hpp"
#include "arg_parse.hpp"
#include "mesh.hpp"
#include "zernike.hpp"

using namespace std;
using namespace argparse;

string sh =
  "Computes a mesh from Zernike moments.\n"
  "Input should be in ZM format as produced by Shape2Zernike.";
string eh = "";
string v_help = "Outputs additional informations including progression bars";
string t_help = "number of threads to use in parallel, use 0 to adapt to the machine";
string f_help = "Numerical precision in fixed notation";
string e_help = "Numerical precision in scientific notation";
string thresh_help = "Threshold value which separates the inside\n"
                "from the outside (default is 1/2)";
string r_help = "Does not regularized the mesh";
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
  elapsed timer;
  string filename = "-";
  int nt = 1;
  p.prog_name = "Zernike2Shape";

  p.flag("v", "verbose", v_help);
  p.option("t", "threads", "N_THREAD", nt, t_help);
  p.option("f", "", "DIGITS", digit, f_help);
  p.option("e", "", "DIGITS", digit, e_help);
  p.flag("r", "raw", r_help);
  p.option("", "threshold", "THRESH", thresh, thresh_help);
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

  // number of threads
  if (nt < 0)
    nt = 1;
  if (nt == 0) {
    nt = thread::hardware_concurrency();
    if (nt == 0)
      nt = 1;
    if (p("v"))
      cerr << "Choosing to run on " << nt << " threads" << endl;
  }

  cout << "# Produced by " << p.prog_name << " (" << p.version_text << ") from file: " << filename << endl;
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

  mesh m = marching_tetrahedra({-1, 1, res}, {-1, 1, res}, {-1, 1, res}, zm, thresh, !p("r"), nt, p("v"));
  
  cout << m;
  
  if (p("v"))
    cerr << p.prog_name << " used " << (int) (timer.seconds() * 100) / 100. << " seconds to run.\n"; 
}