/** \file zm.cpp
  A standalone program to compute zernike moments from an OFF file.
  It can also computes rotational invariants.
*/
#include "arg_parse.hpp"
#include "iotools.hpp"
#include "zernike_def.hpp"
#include "moments.hpp"

using namespace std;
using namespace argparse;

const scheme_selector schemes;
const int Nmax = (ZER_MAX_N < schemes.max_order()) ? ZER_MAX_N : schemes.max_order();
const string nmax = to_string(Nmax);

string sh =
  "Computes Zernike moments or invariants derived from them.\n"
  "Input should be in OFF format.";
string eh = "Currently works up to N = " + nmax
            + " for the computation of the moments.";
string v = "version 0.1";
string t_help = "runs internal sanity checks and exits";
string m_help = "Computes Zernike moments";
string i_help = "Computes Zernike rotational invariants";
string s_help = "Computes signature invariants";
string n_help = "Normalizes the moments such that order 0 gives 1";
string c_help = "Chop to 0 Zenike moments less than 1e-14";
string z_help = "Reads Zernike moments in ZM format instead of computing them";
string d_help = "Reads Zernike moments in ZM format and substract them from the computed moments";

string FILE_help = "Reads FILE in OFF format (default is standard input)";
string N_help = "The maximum order of Zernike moments computed";
string die_N_msg ="N must be positive and no more than " + nmax + " to compute moments.";
string radius_warning =
  "Warning: shape radius is larger than one. Risks of imprecisions.";
string f_help = "Numerical precision in fixed notation";
string e_help = "Numerical precision in scientific notation";

parser p(sh, eh, v);
int N = 0;
int digit = 6;
s_vec coord = {0, 0, 0};

int main (int argc, char *argv[])
{
  string filename = "-";
  string zm_filename;
  p.prog_name = "zm";

  p.flag("t", "tests", t_help);
  p.flag("m", "moments", m_help);
  p.flag("i", "invariants", i_help);
  p.flag("s", "signatures", s_help);
  p.flag("n", "normalize", n_help);
  p.flag("c", "chop", c_help);
  p.flag("z", "zm", z_help);
  p.option("d", "diff", "ZMFILE", zm_filename, d_help);
  p.option("f", "", "DIGITS", digit, f_help);
  p.option("e", "", "DIGITS", digit, e_help);
  p.arg("N", N, N_help);
  p.opt_arg("FILE", filename, FILE_help);

  p.selection({"t", "m", "i", "s"});

  p.run(argc, argv);

  if (digit < 0)
    digit = 6;
  if (p("f"))
    cout << fixed << setprecision(digit);
  else if (p("e"))
    cout << scientific << setprecision(digit);

  if (p("t")) {
    for (auto &s: schemes.schemes)
      cout << s;
    return 0;
  }

  cout << "# Produced by zm (" << v << ") from file: " << filename << endl;
  cout << "# Date: " << now() << endl;
  
  zernike zm;
  if (p("z")) {
    if (N < 0)
      p.die(die_N_msg);
    zernike zm2;
    string err = read_file(filename, zm2);
    if (!err.empty())
      p.die(err);
    zm = zernike(N, zm2);
  }
  else {
    if (N < 0 || N > Nmax)
      p.die(die_N_msg);
    mesh m;
    string err = read_file(filename, m);
    if (!err.empty())
      p.die(err);

    double rad = m.radius();
    cout << "# Mesh: " << m.points.size() << " vertices, "
         << m.triangles.size() << " facets, "
         << "radius: " << rad << endl;

    if (rad > 1.001)
      p.warn(radius_warning);

    // compute moments
    const scheme &s = schemes.get_scheme(N);
    zernike_m_int zmi(N);
    mesh_simple_integrate(m, s, zmi);
    zm = zmi;
  }
  
  zm.normalize(zm_norm::ortho);
  if (p("n"))
    zm.rescale(1);
  if (p("c"))
    zm.chop(1e-14);

  zernike zm2;
  if (p("d")) {
    string err = read_file(zm_filename, zm2);
    if (!err.empty())
      p.die(err);
    zm2.normalize(zm_norm::ortho);
    cout << "# Substracted data from file " << zm_filename << endl;
  }

  if (p("m")) {
    if (p("d"))
      zm = zm - zm2;
    cout << zm;
  }
  else if (p("i")) {
    rotational_invariants ri(zm);
    if (p("d")) {
      rotational_invariants ri2(zm2);
      ri = ri - ri2;
    }
    cout << ri;
  }
  else if (p("s")) {
    signature_invariants si(zm);
    if (p("d")) {
      signature_invariants si2(zm2);
      si = si - si2;
    }
    cout << si;
  }
}
