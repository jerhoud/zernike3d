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
string r_help = "Does not regularized the mesh";
string a_help = "Approximate density evaluation, give n, low, high in quote (e.g. \"30 0.2 0.8\")\n"
                "Does only compute up to n if intermediate result is outside the given thresholds";
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

class cascading_eval
{
public:
  std::function<double(const vec &)> fast, slow;
  const double t_low, t_high;

  cascading_eval(std::function<double(const vec &)> f1, double tl, double hl, std::function<double(const vec &)> f2):
  fast(f1), slow(f2), t_low(tl), t_high(hl)
  {}

  double operator() (const vec &v)
  {
    const double result = fast(v);
    if (result < t_low || result > t_high)
      return result;
    else
      return slow(v);
  }
};

class approx_data
{
public:
  int n;
  double tl, hl;
};

istream &operator >>(istream &is, approx_data &ad)
{
  approx_data a;
  is >> a.n >> a.tl >> a.hl;
  if (is)
    ad = a;
  return is;
}

int main (int argc, char *argv[])
{
  elapsed timer;
  approx_data ad;
  string filename = "-";
  p.prog_name = "rzm";

  p.flag("v", "verbose", v_help);
  p.option("a", "approx", "n low high", ad, a_help);
  p.option("f", "", "DIGITS", digit, f_help);
  p.option("e", "", "DIGITS", digit, e_help);
  p.flag("r", "raw", r_help);
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

  zernike_eval f(zernike(N, zm));
  mesh m;
  if (p("a")) {
    zernike_eval fast(zernike(ad.n, zm));
    cascading_eval g(fast, ad.tl, ad.hl, f);
    m = marching_tetrahedra({-1, 1, res}, {-1, 1, res}, {-1, 1, res}, g, thresh, !p("r"), p("v"));
  }
  else
    m = marching_tetrahedra({-1, 1, res}, {-1, 1, res}, {-1, 1, res}, f, thresh, !p("r"), p("v"));
  
  cout << m;
  
  if (p("v"))
    cerr << "rzm used " << (int) (timer.seconds() * 100) / 100. << " seconds to run.\n"; 
}