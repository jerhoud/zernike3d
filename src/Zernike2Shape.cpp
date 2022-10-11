/** \file Zernike2Shape.cpp.
  A standalone program to compute a mesh from zernike moments.

  It reads Zernike moments from a ZM file (as produced by Shape2Zernike)
  and build a mesh (in OFF format) using marching tetrahedra.
  \author J. Houdayer
*/

#include "version.hpp"
#include "parallel.hpp"
#include "arg_parse.hpp"
#include "mesh.hpp"
#include "zernike.hpp"

using namespace std;
using namespace argparse;

string sh =
  "Computes a shape from Zernike moments.\n"
  "Input should be in ZM format as produced by Shape2Zernike.\n"
  "Output is in OFF format.";
string ex = "Zernike2Shape 50 100 mom.zm                     Builds an OFF shape from the given moments up to order 50 on a 100^3 lattice\n"
            "Zernike2Shape -vt4 -o shape.off 50 100 mom.zm   Same running on 4 threads with progression bar and output saved to file";
string eh = "";
string v_help = "Outputs additional informations including progression bars";
string o_help = "Save output to the given file instead of standard output";
string t_help = "number of threads to use in parallel, use 0 to adapt to the machine";
string d_help = "Number of significant digits printed in the output (default is 6)";
string thresh_help = "Threshold value which separates the inside from the outside (default is 1/2)";
string N_help = "The maximum order of Zernike moments to use (if available)";
string RES_help = "Resolution of the mesh (i.e. number of intervals between -1 and 1)";
string FILE_help = "Reads FILE in ZM format (default is standard input)";
string die_N_msg = "N must be positive.";
string warn_N_msg = "N larger than maximum moment available. Adapting.";
string bad_output_msg = "Cannot open output file: ";

int main (int argc, char *argv[])
{
  int N = 0;
  int digit = 6;
  int res;
  double thresh = 0.5;
  elapsed timer;
  string filename = "-";
  string output = "-";
  int nt = 1;
  
  parser p(sh, eh, ex);
  p.prog_name = "Zernike2Shape";
  p.flag("v", "verbose", v_help);
  p.option("t", "threads", "THREAD", nt, t_help);
  p.option("o", "output", "FILE", output, o_help);
  p.option("d", "digits", "DIGITS", digit, d_help);
  p.option("", "threshold", "THRESH", thresh, thresh_help);
  p.arg("N", N, N_help);
  p.arg("RES", res, RES_help);
  p.opt_arg("FILE", filename, FILE_help);

  p.run(argc, argv);

  if (N < 0)
    p.die(die_N_msg);

  smart_output out(output);
  if (!out)
    p.die(bad_output_msg + output + " (" + strerror(errno) + ")");

  if (digit <= 0)
    digit = 1;
  out << setprecision(digit);

  // number of threads
  #ifdef NO_THREADS
    if (nt != 1)
      p.warn("Threads are not available in this build, running on one thread");
  #else
  if (nt < 0)
    nt = 1;
  if (nt == 0) {
    nt = max_threads();
    if (nt == 0)
      nt = 1;
    if (p("v"))
      cerr << "Choosing to run on " << nt << " threads" << endl;
  }
  #endif

  zernike zm;
  string err = read_file(filename, zm, p("v"));
  if (!err.empty())
    p.die(err);
  if (N > zm.order()) {
    p.warn(warn_N_msg);
    N = zm.order();
  }
  zm.normalize(zm_norm::dual);

  mesh m = marching_tetrahedra({-1, 1, res}, {-1, 1, res}, {-1, 1, res}, zm, thresh, true, nt, p("v"));

  out << "# Produced by " << p.prog_name << " (" << p.version_text << ") from file: " << filename << "\n";
  out << "# Date: " << now() << "\n";
  out << m;

  if (p("v"))
    cerr << p.prog_name << " used " << (int) (timer.seconds() * 100) / 100. << " seconds to run.\n"; 
}