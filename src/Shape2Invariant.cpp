/** \file Shape2Invariant.cpp
  A standalone program to compute invariants from an OFF file or a simple shape.
  \author J. Houdayer
*/

#include "version.hpp"
#include "arg_parse.hpp"
#include "parallel.hpp"
#include "moments.hpp"
#include "invariants.hpp"
#include "coefs.hpp"

using namespace std;
using namespace argparse;

string sh =
  "Computes shape invariants, input should be in OFF format.";
string eh = "";
string ex = "Shape2Invariant 20 shape.off                   Computes the invariants of shape.off up to order 20\n"
            "Shape2Zernike -t4 -o result.si 20 shape.off    Same running on four threads with results written to file\n"
            "Shape2Zernike -e --cube 20                     Exact invariants of a cube";
string v_help = "outputs more informations, including progression bars";
string q_help = "represses all warnings and error messages";
string o_help = "save output to the given file instead of standard output";
string t_help = "number of threads to use in parallel, use 0 to adapt to the machine";
string d_help = "number of significant digits printed in the output (default is 8)\n"
                "the precision of computation is set accordingly";

string FILE_help = "reads FILE in OFF or ZM format (default is standard input)";
string N_help = "the maximum order of invariant computed";
string die_N_msg ="N must be positive";
string approx_warning = "Warning; requested precision is very small, program may not halt. Allowed error by facet: ";
string die_unknown_format = "Unknown file format (should be OFF or ZM): ";
string warn_N_msg = "N too large for the maximum moment available. Adapting.";
string bad_output_msg = "Cannot open output file: ";

int main(int argc, char *argv[])
{
  // Initialization

  elapsed timer;
  
  int N = 0;
  int digit = 8;
  int nt = 1;

  string filename = "-";
  string output = "-";
  string zm_filename;

  const triquad_selector triquad_schemes;


  // Set command line options 

  parser p(sh, eh, ex);
  p.prog_name = "Shape2Invariant";

  p.flag("v", "verbose", v_help);
  p.flag("q", "quiet", q_help);
  p.option("o", "output", "FILE", output, o_help);
  p.option("t", "threads", "THREAD", nt, t_help);
  p.option("d", "digits", "DIGITS", digit, d_help);

  p.arg("N", N, N_help);
  p.opt_arg("FILE", filename, FILE_help);

  p.quiet("q");
  p.exclusion({"v", "q"});


// Parse command line

  p.run(argc, argv);

  smart_output out(output);
  if (!out)
    p.die(bad_output_msg + output + " (" + strerror(errno) + ")");

  if (digit <= 0)
    digit = 1;
  out << setprecision(digit);
  const double approx_err = pow(0.1, digit);


  // check N
    
  if (N < 0)
    p.die(die_N_msg);
  
  
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


  // Prepare input stream 
  
  smart_input is(filename);
  if (!is)
    p.die(cannot_open_msg + is.name + " (" + strerror(errno) + ")");


  // Write ouput header

  out << "# Produced by " << p.prog_name << " (" << p.version_text << ") from file: " << is.name << "\n";
  out << "# Date: " << now() << "\n";


  // Identify type of input file
  
  zernike zm;

  istringstream iss;
  is.peek_line(iss);
  string filetype;
  iss >> filetype;

  // it is a ZM file
  if (filetype == "ZM" || filetype == "zm") {
    zernike zm2;
    string err = read_object(is, zm2, p("v"));
    if (!err.empty())
      p.die(err);
    if (2 * N > zm.order()) {
      p.warn(warn_N_msg);
      N = zm.order() / 2;
    }
    zm = zernike(2 * N, zm2);
  }
  // it is an OFF file
  else if (filetype == "OFF" || filetype == "off") {
    mesh m;
    string err = read_object(is, m, p("v"));
    if (!err.empty())
      p.die(err);

    double rad = m.radius();
    out << "# Mesh: " << m.points.size() << " vertices, "
        << m.triangles.size() << " facets, "
        << "radius: " << rad << "\n";

    // compute moments
    const double facet_error = approx_err / sqrt(m.triangles.size());
    if (facet_error < 1e-13) {
      ostringstream out;
      out << scientific << facet_error;
      p.warn(approx_warning + out.str());
    }
    zm = mesh_approx_integrate(m, 2 * N, approx_err, triquad_schemes, nt, p("v"));
    out << "# estimation of approximation error on the moments: " << zm.get_error() << "\n";
  }
  else
    p.die(die_unknown_format + is.name);

  
  // compute rotational invariants
  
  rotational_invariants ri(N);
  ri.eval(zm);
  fnk f(N);
  f.eval(ri);


  // compute coefficients

  factorials facs(2 * N + 1);
  double_factorials dfacs(2 * N + 1);
  binomials bins(2 * N + 3);
  unl0 u0(N, bins);
  theta t(N, facs, dfacs, bins);
  omega o0(u0, t);

  

  // ciao !

  if (p("v"))
    cerr << p.prog_name << " used " << (int) (timer.seconds() * 100) / 100. << " seconds to run.\n"; 
}