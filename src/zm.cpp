/** \file zm.cpp
  A standalone program to compute zernike moments from an OFF file.
  It can also computes rotational invariants.
*/
#include "arg_parse.hpp"
#include "iotools.hpp"
#include "moments.hpp"

using namespace std;
using namespace argparse;

const triquad_selector triquad_schemes;
const gauss_selector gauss_schemes;
const int N_exact = triquad_schemes.max_order();
const string n_exact = to_string(N_exact);

string sh =
  "Computes Zernike moments or invariants derived from them.\n"
  "Input should be in OFF or ZM format.";
string eh = "Currently works up to N = " + n_exact
            + " for the exact computation of the moments.\n";
string v_help = "outputs more informations, including progression bars";
string t_help = "runs internal sanity checks and exits";
string m_help = "Computes Zernike moments";
string i_help = "Computes Zernike rotational invariants";
string s_help = "Computes signature invariants";
string n_help = "Alternative normalization (first polynomial is 1): moments * sqrt(3/4pi)";
string a_help = "Computes the moments using approximate methods to get the required correct DIGITS";
string c_help = "The Zernike moments are output in complex form";
string p_help = "Multiplies the moments by the phase factor (-1)^m";
string chop_help = "Chops to 0 very small Zenike moments";
string raw_help = "Divides Zernike moments by sqrt(2n+3)";
string dual_help = "Multiplies Zernike moments by sqrt(2n+3)";
string d_help = "Reads Zernike moments in ZM format and substract them from the computed moments";
string f_help = "Numerical precision in fixed notation";
string e_help = "Numerical precision in scientific notation";

string FILE_help = "Reads FILE in OFF or ZM format (default is standard input)";
string N_help = "The maximum order of Zernike moments computed";
string die_N_msg ="N must be positive and no more than "
                       + n_exact + " for exact compututation of the moments.";
string radius_warning =
  "Warning: shape radius is larger than one. Risks of imprecisions.";
string die_approx_msg = "-a option: ERROR must be positive";
string approx_warning = "Warning; requested precision is very small, program may not halt. Allowed error by facet: ";
string die_unknown_format = "Unknown file format (should be OFF or ZM): ";


int main (int argc, char *argv[])
{
  elapsed timer;
  parser p(sh, eh);
  int N = 0;
  int digit = 6;
  int approx = 13;

  string filename = "-";
  string zm_filename;
  p.prog_name = "zm";

  // Set command line options 

  p.flag("v", "verbose", v_help);
  p.flag("t", "tests", t_help);
  p.flag("m", "moments", m_help);
  p.flag("i", "invariants", i_help);
  p.flag("s", "signatures", s_help);
  p.flag("n", "normalize", n_help);
  p.flag("", "chop", chop_help);
  p.flag("", "raw", raw_help);
  p.flag("", "dual", dual_help);
  p.flag("c", "complex", c_help);
  p.flag("p", "phase", p_help);
  p.option("a", "approximate", "DIGITS", approx, a_help);
  p.option("d", "diff", "ZMFILE", zm_filename, d_help);
  p.option("f", "", "DIGITS", digit, f_help);
  p.option("e", "", "DIGITS", digit, e_help);
  p.arg("N", N, N_help);
  p.opt_arg("FILE", filename, FILE_help);

  p.selection({"t", "m", "i", "s"});
  p.exclusion({"raw", "dual"});

  // Parse command line

  p.run(argc, argv);

  // Apply options -e and -f

  const double approx_err = pow(0.1, approx);
  if (p("a") && !p("e") && !p("f"))
    digit = approx + 1;
  if (digit < 0)
    digit = 6;
  if (p("f"))
    cout << fixed;
  else if (p("e"))
    cout << scientific;
  cout << setprecision(digit);

  // Option -t auto tests and exits

  if (p("t")) {
    cout << "checking gauss quadratures on the segment" << endl;
    for (auto &s: gauss_schemes.schemes)
      cout << s;
    cout << "checking primary quadratures on the triangle" << endl;
      for (auto &s: triquad_schemes.schemes)
      cout << s;
     cout << "checking secondary quadratures on the triangle" << endl;
   for (auto &s: triquad_schemes.secondary_schemes)
      cout << s;
    return 0;
  }

  // 

  if (N < 0)
    p.die(die_N_msg);
  zernike zm;
  smart_input is(filename);
  if (!is)
    p.die(cannot_open_msg + is.name + " (" + strerror(errno) + ")");

  // Output header

  cout << "# Produced by zm (" << p.version_text << ") from file: " << is.name << endl;
  cout << "# Date: " << now() << endl;

  // Identify type of input file
  
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
    zm = zernike(N, zm2);
    zm.output = zm_output::real;
  }
  // it is an OFF file
  else if (filetype == "OFF" || filetype == "off") {
    if ((!p("a") && N > N_exact))
      p.die(die_N_msg);
    mesh m;
    string err = read_object(is, m, p("v"));
    if (!err.empty())
      p.die(err);

    double rad = m.radius();
    cout << "# Mesh: " << m.points.size() << " vertices, "
         << m.triangles.size() << " facets, "
         << "radius: " << rad << endl;

    if (rad > 1.001)
      p.warn(radius_warning);

    // compute moments
    if (p("a")) {
      if (approx_err <= 0)
        p.die(die_approx_msg);
      const double facet_error = approx_err / m.triangles.size();
      if (facet_error < 1e-13) {
        ostringstream out;
        out << scientific << facet_error;
        p.warn(approx_warning + out.str());
      }
      zm = mesh_approx_integrate(m, N, approx_err, triquad_schemes, gauss_schemes, p("v"));
      cout << "# approximation error estimate: " << zm.error << endl;
    }
    else
      zm = mesh_exact_integrate(m, N, triquad_schemes, gauss_schemes, p("v"));
  }
  else
    p.die(die_unknown_format + is.name);
  
  // Select normalization and apply output options -c -p --chop

  zm.normalize(make_norm(p("raw"), p("dual"), p("n")));
  zm.output = make_output(p("c"), p("p"));
  if (p("chop"))
    zm.chop(zm.error / 10);

  // option -d read a secondary zm file

  zernike zm2;
  if (p("d")) {
    string err = read_file(zm_filename, zm2, p("v"));
    if (!err.empty())
      p.die(err);
    zm2.normalize(zm_norm::ortho);
    cout << "# Substracted data from file " << zm_filename << endl;
  }

// command -m output moments
  if (p("m")) {
    if (p("d"))
      zm = zm - zm2;
    cout << zm;
  }
  // command -i output rotational invariants
  else if (p("i")) {
    rotational_invariants ri(zm);
    if (p("d")) {
      rotational_invariants ri2(zm2);
      ri = ri - ri2;
    }
    cout << ri;
  }
  // command -s output signature
  else if (p("s")) {
    signature_invariants si(zm);
    if (p("d")) {
      signature_invariants si2(zm2);
      si = si - si2;
    }
    cout << si;
  }

  if (p("v"))
    cerr << "zm used " << (int) (timer.seconds() * 100) / 100. << " seconds to run.\n"; 
}
