/** \file Shape2Zernike.cpp
  A standalone program to compute zernike moments from an OFF file.
  \author J. Houdayer
*/

#include "version.hpp"
#include "arg_parse.hpp"
#include "parallel.hpp"
#include "moments.hpp"

using namespace std;
using namespace argparse;

const triquad_selector triquad_schemes;
const int N_exact = triquad_schemes.max_order();
const string n_exact = to_string(N_exact);

string sh =
  "Computes Zernike moments, input should be in OFF format.";
string eh = "Currently works up to N = " + n_exact
            + " for the exact computation of the moments.\n"
            "No limit for N when using -a.\n"
            "The shape must fit into the unit ball (no implicit centering or rescaling, use MakeShape to do this).";
string ex = "Shape2Zernike 50 shape.off                     Computes the Zernike moments of shape.off up to order 50\n"
            "Shape2Zernike -a 8 -o result.zm 50 shape.off   Same using approximate algorithm with 8 digit precision and results written to file\n"
            "Shape2Zernike -vt 4 50 shape.off               Same running on four threads, with progression bar";
string v_help = "outputs more informations, including progression bars";
string q_help = "represses all warnings and error messages";
string o_help = "save output to the given file instead of standard output";
string t_help = "number of threads to use in parallel, use 0 to adapt to the machine";
string tests_help = "runs internal sanity checks and exits";
string n_help = "multiplies the moments by sqrt(3/4pi)";
string a_help = "computes the moments using approximate methods to get the required correct DIGITS";
string r_help = "the Zernike moments are output in real form instead of complex";
string p_help = "multiplies the moments by the phase factor (-1)^m";
string diff_help = "reads Zernike moments in ZM format and substract them from the computed moments";
string d_help = "number of significant digits printed in the output (default is 8)";

string FILE_help = "reads FILE in OFF or ZM format (default is standard input)";
string N_help = "the maximum order of Zernike moments computed";
string die_N_msg ="N must be positive and no more than "
                       + n_exact + " for exact computation of the moments";
string radius_warning =
  "Warning: shape radius is larger than one. Risks of imprecisions.";
string approx_warning = "Warning; requested precision is very small, program may not halt. Allowed error by facet: ";
string die_unknown_format = "Unknown file format (should be OFF or ZM): ";
string bad_output_msg = "Cannot open output file: ";

int main (int argc, char *argv[])
{
  // Initialization

  elapsed timer;
  int N = 0;
  int digit = 8;
  int approx = 13;
  int nt = 1;

  string filename = "-";
  string output = "-";
  string zm_filename;


  // Set command line options 

  parser p(sh, eh, ex);
  p.prog_name = "Shape2Zernike";

  p.flag("v", "verbose", v_help);
  p.flag("q", "quiet", q_help);
  p.option("o", "output", "FILE", output, o_help);
  p.option("t", "threads", "THREAD", nt, t_help);
  p.option("a", "approximate", "DIGITS", approx, a_help);
  p.option("d", "digits", "DIGITS", digit, d_help);

  p.hidden(true);
  p.flag("r", "real", r_help);
  p.flag("n", "normalize", n_help);
  p.flag("p", "phase", p_help);
  p.option("", "diff", "ZMFILE", zm_filename, diff_help);
  p.flag("", "tests", tests_help);

  p.arg("N", N, N_help);
  p.opt_arg("FILE", filename, FILE_help);

  p.quiet("q");
  p.exclusion({"v", "q"});

  // Parse command line

  p.run(argc, argv);

  smart_output out(output);
  if (!out)
    p.die(bad_output_msg + output + " (" + strerror(errno) + ")");

  // Apply options -a and -d

  const double approx_err = pow(0.1, approx);
  if (p("a") && !p("d"))
    digit = approx + 1;
  if (digit <= 0)
    digit = 1;
  out << setprecision(digit);


  // Option --tests auto tests and exits

  if (p("tests")) {
    out << "checking primary quadratures on the triangle\n";
    for (auto &s: triquad_schemes.schemes)
      out << s;
    out << "checking secondary quadratures on the triangle\n";
    for (auto &s: triquad_schemes.secondary_schemes)
      out << s;
    return 0;
  }


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

  // it is a ZM file, read it
  if (filetype == "ZM" || filetype == "zm") {
    zernike zm2;
    string err = read_object(is, zm2, p("v"));
    if (!err.empty())
      p.die(err);
    zm = zernike(N, zm2);
  }
  // it is an OFF file, compute moments
  else if (filetype == "OFF" || filetype == "off") {
    // check max bound on N
    if ((!p("a") && N > N_exact))
      p.die(die_N_msg);
    // read OFF file
    mesh m;
    string err = read_object(is, m, p("v"));
    if (!err.empty())
      p.die(err);
    // output some infos
    double rad = m.radius();
    out << "# Mesh: " << m.points.size() << " vertices, "
        << m.triangles.size() << " facets, "
        << "radius: " << rad << "\n";
    if (rad > 1.001)
      p.warn(radius_warning);

    // compute moments
    if (p("a")) {
      const double facet_error = approx_err / sqrt(m.triangles.size());
      if (facet_error < 1e-13) {
        ostringstream out;
        out << scientific << facet_error;
        p.warn(approx_warning + out.str());
      }
      zm = mesh_approx_integrate(m, N, approx_err, triquad_schemes, nt, p("v"));
      out << "# approximation error estimate: " << zm.get_error() << "\n";
    }
    else {
      zm = mesh_exact_integrate(m, N, triquad_schemes, nt, p("v"));
      out << "# error estimate: " << zm.get_error() << "\n";
    }
  }
  // unknown file type
  else
    p.die(die_unknown_format + is.name);
  

  // Select normalization and apply output options -r -p and -n

  zm.normalize(make_norm(false, false, p("n")));
  zm.output = make_output(!p("r"), p("p"));

  // option -d read a secondary zm file

  if (p("diff")) {
    zernike zm2;
    string err = read_file(zm_filename, zm2, p("v"));
    if (!err.empty())
      p.die(err);
    zm2.normalize(zm_norm::ortho);
    out << "# Substracted data from file " << zm_filename << "\n";
    zm = zm - zm2;
  }


  // write result

  out << zm;

  // ciao !

  if (p("v"))
    cerr << p.prog_name << " used " << (int) (timer.seconds() * 100) / 100. << " seconds to run.\n"; 

  return 0;
}
