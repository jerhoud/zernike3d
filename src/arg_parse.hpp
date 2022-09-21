/** \file arg_parse.hpp
  A simple command line argument parser in one header without dependencies.
  \author J. Houdayer
*/

#include <cstring>
#include <vector>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <ctime>

#define strgf(x) #x
/** expands a macro inside quotes. */
#define stringify(x) strgf(x)

#ifndef VERSION
/** VERSION is written when option --version is used. */
#define VERSION undefined version
#endif

namespace argparse {

using namespace std;

const string color_red = "\033[0;31m";
const string color_blue = "\033[0;34m";
const string color_reset = "\033[0m";

const string usage_msg = "Usage: ";
const string example_msg = "Examples:";
const string options_msg = "Options";
const string arguments_msg = "Arguments";
const string std_help_msg = "shows this help message and exits";
const string std_version_msg = "shows version information and exits";
const string missing_opt_name_msg = "Missing option name for argument ";
const string unknown_opt_msg = "Unknown option ";
const string unwanted_opt_arg_msg = "Unwanted argument for option ";
const string unwanted_arg_msg = "Unwanted argument ";
const string invalid_arg_msg = "Invalid value for argument ";
const string ambiguous_opt_msg = "Ambiguous option name ";
const string missing_opt_arg_msg = "Missing argument for option ";
const string invalid_opt_arg_msg = "Invalid argument for option ";
const string missing_arg_msg = "Missing argument ";
const string more_help_msg = "Use option -h for help";
const string exclusion_msg = "The following options are incompatible: ";
const string selection_msg = "You must use one of the following options: ";
const string bad_test_msg = "Bug: testing unkown option: ";
const string bad_opt_arg = "Bug: optional argument before non-optional argument: ";

/** returns a string which represents the current date and time. */
string now()
{
  time_t rawtime;
  struct tm *timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  string s = asctime(timeinfo);

  if (s.size() != 0 && s.back() == '\n')
    s.pop_back();

 return s;
}

/** describes the different command line elements. */
enum token_type { tok_option, tok_long_option, tok_option_arg, tok_arg };

/** a class to hold the type and name of a command element. */
class token
{
public:
  const token_type tp;
  const string s;
};

/** parses command line elements into options and argument.*/
void parse_tokens(int argc, const char *const argv[], vector<token> &toks)
{
  bool opt = true;
  for (int i = 0 ; i < argc ; ++i) {
    const string arg = argv[i];
    if (arg == "--") // stop parsing options
      opt = false;
    else if (!opt || arg == "" || arg == "-" || arg[0] != '-') // a regular argument
      toks.push_back({tok_arg, arg});
    else {
      const char c = arg[1];
      if (c == '.' || (c >= '0' && c <= '9')) // not an option: a negative number
        toks.push_back({tok_arg, arg});
      else {                  // an option
        const bool long_opt = c == '-';
        const int s = (long_opt) ? 2 : 1;
        int p;
        if (long_opt)
          p = arg.find("=", s);
        else
          p = arg.find_first_not_of(
                "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ", s);
        if (p == s) // a missing option name (like -! or --=7)
          toks.push_back({tok_option, ""});
        else {
          const string head = arg.substr(s, p - s);
          if (long_opt)    // a long option
            toks.push_back({tok_long_option, head});
          else             // a bunch of short options
            for (auto c: head)
              toks.push_back({tok_option, string(1, c)});
        }
        if (p != (int) string::npos) // what is left is an argument
          toks.push_back({tok_option_arg, arg.substr(p + (arg[p]=='='))});
      }
    }
  }
}

/** pretty prints a paragraph. with a left margin of size w.*/
void show_paragraph(ostream &os, int w, const string &h)
{
  int p = 0, pp;
  while ((pp = h.find("\n", p)) != (int) string::npos) {
    os << h.substr(p, pp - p) << endl;
    os << string(w, ' ');
    p = pp + 1;
  }
  os << h.substr(p, pp - p) << endl;
}

/* a class to store the use of an option in a recorder */
class opt_recorder
{
public:
  string name;
  int pos;

  opt_recorder(const string &n): name(n), pos(-1) {}
  opt_recorder(const string &n, int p): name(n), pos(p) {}
};

/** a base class representing an option.*/
class option_desc_base
{
public:
  const string short_name;
  const string long_name;
  const bool want_arg;
  const string help;
  string short_usage, short_print, long_print;
  bool used;

  option_desc_base(const string &s, const string &l, const string &arg_name, const string &h):
  short_name(s), long_name(l), want_arg(!arg_name.empty()), help(h), used(false)
  {
    if (!short_name.empty()) {
      short_usage = "-" + short_name;
      if (want_arg)
        short_usage += " " + arg_name;
    }
    short_print = short_usage;
    if (!long_name.empty()) {
      long_print = "--" + long_name;
      if (want_arg)
        long_print += " " + arg_name;
      if (!short_name.empty())
        short_print += " ";
    }
  }

  /** called when the option is found alone on the command line. */
  virtual void process() { used = true; }
  /** called when the option is found with an arg on the command line. */
  virtual bool process(const string &) { used = true; return true; }

  /** help output of the option.*/ 
  void show(ostream &os, int w_short, int w_long)
  {
    os << " " << left << setw(w_short) << short_print
       << setw(w_long) << long_print << "  ";
    show_paragraph(os, w_short + w_long + 3, help);
  }
};

/** reads a string into a variable of type T, returns true if successfull.*/ 
template <typename T>
bool parse_arg(const string &arg, T &var)
{
  istringstream iss(arg);
  T v;
  iss >> v;
  if (iss.fail() || !iss.eof())
    return false;
  var = v;
  return true;
}

template <>
bool parse_arg<string> (const string &arg, string &var)
{ var = arg; return true; }

/** tries reading a T variable and appends it to the vector.*/
template <typename T>
bool parse_list(const string &arg, vector<T> &var)
{
    T v;
    if (!parse_arg(arg, v))
      return false;
    var.push_back(v);
    return true;
}

/** a yes / no option that record its usage in a recorder */
class rec_flag_desc: public option_desc_base
{
public:
  std::vector<opt_recorder> &rec;

  rec_flag_desc(const string &s, const string &l, std::vector<opt_recorder> &r, const string &h):
  option_desc_base(s, l, "", h), rec(r) {}

  void process()
  {
    used = true;
    rec.push_back(opt_recorder((short_name.empty()) ? long_name : short_name));
  }
};

/** an option that sets a variable with its argument. */
template <typename T>
class set_option_desc: public option_desc_base
{
public:
  T &var;

  set_option_desc(const string &s, const string &l, const string &an, T &v, const string &h):
  option_desc_base(s, l, an, h), var(v) {}

  bool process(const string &arg)
  {
    used = true;
    return parse_arg(arg, var);
  }
};

/** an option that appends its argument to a vector each time it is met on the command line. */
template <typename T>
class list_option_desc: public option_desc_base
{
public:
  vector<T> &var;

  list_option_desc(const string &s, const string &l, const string &an, vector<T> &v, const string &h):
  option_desc_base(s, l, an, h), var(v) {}

  bool process(const string &arg)
  { used = true; return parse_list(arg, var); }
};

/** an option that appends and records its argument to a vector each time it is met on the command line. */
template <typename T>
class rec_list_option_desc: public option_desc_base
{
public:
  vector<T> &var;
  vector<opt_recorder> &rec;

  rec_list_option_desc(const string &s, const string &l, const string &an, vector<T> &v, vector<opt_recorder> &r, const string &h):
  option_desc_base(s, l, an, h), var(v), rec(r) {}

  bool process(const string &arg)
  {
    used = true;
    rec.push_back(opt_recorder((short_name.empty()) ? long_name : short_name, var.size()));
    return parse_list(arg, var);
  }
};

/** base class for non option arguments.*/
class arg_desc_base
{
public:
  const bool need;
  const bool rest;
  const string arg_name;
  const string help;
  bool used;

  arg_desc_base(const string &a, bool n, bool r, const string &h):
  need(n), rest(r), arg_name(a), help(h) , used(false) {}

  /** called when the argument is met on the command line.*/
  virtual bool process(const string &arg) = 0;

  /** help output for the argument.*/
  void show(ostream &os, int w)
  {
    os << " " << left << setw(w) << arg_name << "   ";
    show_paragraph(os, w + 5, help);
  }
};

/** a positional argument that sets a variable.*/
template<typename T>
class pos_arg_desc: public arg_desc_base
{
public:
  T &var;

  pos_arg_desc(const string &a, T &v, bool n, const string &h):
  arg_desc_base(a, n, false, h), var(v) {}

  bool process(const string &arg)
  { used = true; return parse_arg(arg, var); }
};

/** an argument that collects all the reminding argument on the command line.*/
template<typename T>
class rest_arg_desc: public arg_desc_base
{
public:
  vector<T> &var;

  rest_arg_desc(const string &a, vector<T> &v, const string &h):
  arg_desc_base(a, false, true, h), var(v) {}

  bool process(const string &arg)
  { used = true; return parse_list(arg, var); }
};

/** The main class to parse the command line.*/ 
class parser
{
public:
  const string start_help, end_help, example_help, version_text;
  string prog_name, missing_arg;
  vector<token> toks;
  vector<option_desc_base *> opts;
  vector<arg_desc_base *> args;
  vector<vector<string> > exclusions;
  vector<vector<string> > selections;

  /** The constructor
    
    @param sh The help string written before the arg / opt descriptions
    @param eh The help string written after the arg / opt descriptions
    @param v The version string (for opt --version) default is macro VERSION
   */
  parser(const string &sh, const string &eh, const string &ex="", const string &v=stringify(VERSION)):
  start_help(sh), end_help(eh), example_help(ex), version_text(v)
  {
    flag("h", "help", std_help_msg);
    flag("", "version", std_version_msg);
  }

  void run(int argc, const char *const argv[]);
  void parse(int argc, const char *const argv[]);
  void analysis();
  void standard() const;
  void check() const;

  bool operator()(const string &s) const;

  /** returns true if any of the given options was met on the command line.*/
  bool any(const vector<string> &s) const
  { return any_of(s.begin(), s.end(), *this); }

  /** returns true if all of the given options were met on the command line.*/
  bool all(const vector<string> &s) const
  { return all_of(s.begin(), s.end(), *this); }

  /** returns true if none of the given options was met on the command line.*/
  bool none(const vector<string> &s) const
  { return none_of(s.begin(), s.end(), *this); }

  /** counts how many of the given options were met on the command line.*/
  int count(const vector<string> &s) const
  { return count_if(s.begin(), s.end(), *this); }

  void show_usage(ostream &os) const;
  void show_help(ostream &os) const;

  /** prints version informations.*/
  void show_version(ostream &os) const
  { os << prog_name << " " << version_text << endl; }

  /** adds a yes / no option to the parser.*/
  void flag(const string &s, const string &l, const string &h)
  { opts.push_back(new option_desc_base(s, l, "", h)); }

  /** adds a yes / no recorded option to the parser.*/
  void rec_flag(const string &s, const string &l, vector<opt_recorder> &r, const string &h)
  { opts.push_back(new rec_flag_desc(s, l, r, h)); }

  /** adds an option with an argument to the parser.*/
  template<typename T>
  void option(const string &s, const string &l, const string &a, T &v, const string &h)
  { opts.push_back(new set_option_desc<T>(s, l, a, v, h)); }

  /** adds a repeatable option with an argument to the parser.*/
  template<typename T>
  void list_option(const string &s, const string &l, const string &a, vector<T> &v, const string &h)
  { opts.push_back(new list_option_desc<T>(s, l, a, v, h)); }

  /** adds a repeatable recorded option with an argument to the parser.*/
  template<typename T>
  void rec_list_option(const string &s, const string &l, const string &a, vector<T> &v, vector<opt_recorder> &r, const string &h)
  { opts.push_back(new rec_list_option_desc<T>(s, l, a, v, r, h)); }

  /** adds a positional argument to the parser.*/
  template<typename T>
  void arg(const string &a, T &v, const string &h)
  {
    if (!args.empty() && !args.back()->need)
      die(bad_opt_arg + a);
    args.push_back(new pos_arg_desc<T>(a, v, true, h));
  }

  /** adds an optional positional argument to the parser, should be at the end.*/
  template<typename T>
  void opt_arg(const string &a, T &v, const string &h)
  { args.push_back(new pos_arg_desc<T>(a, v, false, h)); }

  /** adds a list of arguments to the parser, should be at the end.*/
  template<typename T>
  void rest_arg(const string &a, std::vector<T> &v, const string &h)
  { args.push_back(new rest_arg_desc<T>(a, v, h)); }

  /** sets mutally exclusive options or arguments.*/
  void exclusion(const vector<string> &s)
  { exclusions.push_back(s); }

  /** sets a list of options or arguments from which one and only one must be used.*/
  void selection(const vector<string> &s)
  { selections.push_back(s); }

  /** prints a warning.*/
  void warn(const string &w) const
  { cerr << prog_name << ": " << color_blue << w << color_reset << endl; }

  [[ noreturn ]] void die(const string &err) const;
  [[ noreturn ]] void die(const string &err, const vector<string> &l) const;
};

/** Prints the usage message triggered by an error.*/
void parser::show_usage(ostream &os) const
{
  string flag_str, opt_str, arg_str;
  for (auto &opt: opts)
    if (!opt->short_usage.empty()) {
      if (opt->want_arg) {
        if (!opt_str.empty())
          opt_str += "|";
        opt_str += opt->short_usage;
      }
      else
        flag_str += opt->short_name;
    }
  if (!flag_str.empty()) {
    flag_str = "-" + flag_str;
    if (opt_str.empty())
      opt_str = flag_str;
    else
      opt_str = flag_str + "|" + opt_str;
  }
  for (auto &arg: args) {
    if (!arg_str.empty())
      arg_str += " ";
    if (arg->rest)
      arg_str += arg->arg_name + "...";
    else if (arg->need)
      arg_str += arg->arg_name;
    else
      arg_str += "[" + arg->arg_name + "]";
  }
  if (opt_str.length()>30)
    opt_str = options_msg;
  if (arg_str.length()>30)
    arg_str = arguments_msg;
  os << usage_msg << prog_name;
  if (!opt_str.empty())
    os << " [" << opt_str << "]";
  if (!arg_str.empty())
    os << " " << arg_str;
  os << endl << endl;
  if (example_help != "") {
    os << example_msg << endl << " ";
    show_paragraph(os, 1, example_help);
    os << endl;
  }
}

/** prints the help message triggered by --help.*/
void parser::show_help(ostream &os) const
{
  int w_short = 0, w_long = 0, w_arg = 0;
  for (auto &opt: opts) {
    w_short = max(w_short, (int)opt->short_print.length());
    w_long = max(w_long, (int)opt->long_print.length());
  }
  if (w_short)
    w_short++;
  if (w_long)
    w_long++;
  for (auto &arg: args)
    w_arg = max(w_arg, (int)arg->arg_name.length());

  os << start_help << endl << endl;
  show_usage(os);
  if (!args.empty()) {
    os << arguments_msg << ":" << endl;
    for (auto &arg: args)
      arg->show(os, w_arg);
    os << endl;
  }
  if (!opts.empty()) {
    os << options_msg << ":" << endl;
    for (auto &opt: opts)
      opt->show(os, w_short, w_long);
    os << endl;
  }
  if (!end_help.empty())
    os << end_help << endl << endl;
}

/** parses the command line and executes actions linked with options.*/
void parser::run(int argc, const char *const argv[])
{
  parse(argc, argv);
  analysis();
  standard();
  check();
}

/** parses the command line.*/
void parser::parse(int argc, const char *const argv[])
{
  if (prog_name.empty())
    prog_name = argv[0];
  parse_tokens(argc - 1, argv + 1, toks);
}

/** verifies syntax and execute option actions.*/
void parser::analysis()
{
  int mi = toks.size(), j = 0, mj = args.size();
  for (int i = 0 ; i < mi ; i++) {
    token &tk = toks[i];
    if (tk.tp == tok_option_arg)
      die(unwanted_opt_arg_msg + toks[i-1].s + ": " + tk.s);
    if (tk.tp == tok_arg) {
      if (j >= mj)
        die(unwanted_arg_msg + tk.s);
      arg_desc_base *arg = args[j];
      if (!arg->process(tk.s))
        die(invalid_arg_msg + arg->arg_name + ": " + tk.s);
      if (!arg->rest)
        j++;
      continue;
    }
    option_desc_base *opt = 0;
    if (tk.tp == tok_option) {
      if (tk.s.empty())
        die(missing_opt_name_msg + toks[i + 1].s);
      for (auto &o: opts)
        if (o->short_name == tk.s) {
          opt = o;
          break;
        }
    }
    else {
      for (auto &o: opts)
        if (o->long_name.find(tk.s) == 0) {
          if (opt)
            die(ambiguous_opt_msg + tk.s);
          else
            opt = o;
        }
    }
    if (!opt)
      die(unknown_opt_msg + tk.s);
    if (!opt->want_arg)
      opt->process();
    else {
      i++;
      if (i >= mi || toks[i].tp == tok_option || toks[i].tp == tok_long_option)
        die(missing_opt_arg_msg + tk.s);
      if (!opt->process(toks[i].s))
        die(invalid_opt_arg_msg + tk.s + ": " + toks[i].s);
    }
  }
  if (j != mj && args[j]->need)
    missing_arg = args[j]->arg_name;
}

/** checks and executes option --help and --version.*/
void parser::standard() const
{
  if (operator()("help")) {
    show_help(cout);
    exit(0);
  }
  if (operator()("version")) {
    show_version(cout);
    exit(0);
  }
}

/** checks for missing args and constraints.*/
void parser::check() const
{
  if (!missing_arg.empty())
    die(missing_arg_msg + missing_arg);
  for (auto &excl: exclusions)
    if (count(excl) > 1)
      die(exclusion_msg, excl);
  for(auto &sel: selections)
    if (count(sel) != 1)
      die(selection_msg, sel);
}

/** returns true if the options was met on the command line.*/
bool parser::operator()(const string &s) const
{
  for (auto &opt: opts)
    if (s == opt->short_name || s == opt->long_name)
      return opt->used;
  for (auto &arg: args)
    if (s == arg->arg_name)
      return arg->used;
  die(bad_test_msg + s);
}

/** prints an error message and exits.*/
void parser::die(const string &err) const
{
  cerr << prog_name << ": " << color_red << err << color_reset << endl;
  show_usage(cerr);
  if (!more_help_msg.empty())
    cerr << more_help_msg << endl;
  exit(-1);
}

/** prints an error message and exits.*/
void parser::die(const string &err, const vector<string> &l) const
{
  string errl = err;
  if (!l.empty()) {
    for (auto i = l.begin() ; i != l.end() - 1 ; i++)
      errl += *i + ", ";
    errl += l.back();
  }
  die(errl);
}

}
