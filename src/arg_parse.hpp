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

namespace argparse {

using namespace std;

string color_red = "\033[0;31m";
string color_blue = "\033[0;34m";
string color_reset = "\033[0m";

string usage_msg = "Usage: ";
string options_msg = "Options";
string arguments_msg = "Arguments";
string std_help_msg = "shows this help message and exits";
string std_version_msg = "shows version information and exits";
string missing_opt_name_msg = "Missing option name for argument ";
string unknown_opt_msg = "Unknown option ";
string unwanted_opt_arg_msg = "Unwanted argument for option ";
string unwanted_arg_msg = "Unwanted argument ";
string invalid_arg_msg = "Invalid value for argument ";
string ambiguous_opt_msg = "Ambiguous option name ";
string missing_opt_arg_msg = "Missing argument for option ";
string invalid_opt_arg_msg = "Invalid argument for option ";
string missing_arg_msg = "Missing argument ";
string more_help_msg = "Use option -h for help";
string exclusion_msg = "The following options are incompatible: ";
string selection_msg = "You must use one of the following options: ";
string bad_test_msg = "Bug: testing unkown option: ";
string bad_opt_arg = "Bug: optional argument before non-optional argument: ";

string now()
{
  time_t rawtime;
  struct tm *timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  string s = asctime(timeinfo);

  if (s.size()!=0 && s.back()=='\n')
    s.pop_back();

 return s;
}

enum token_type { tok_option, tok_long_option, tok_option_arg, tok_arg };

class token {
public:
  token_type tp;
  string s;
};

void add_token(const string &arg, vector<token> &toks, bool long_opt)
{
  int s, p;
  if (long_opt) {
    s = 2;
    p = arg.find("=", s);
  }
  else {
    s = 1;
    p = arg.find_first_not_of(
      "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ", s);
  }
  if (p == s)
    toks.push_back({tok_option, ""});
  else {
    string head = arg.substr(s, p - s);
    if (long_opt)
      toks.push_back({tok_long_option, head});
    else
      for (auto c: head)
        toks.push_back({tok_option, string(1, c)});
  }
  if (p != (int) string::npos)
    toks.push_back({tok_option_arg, arg.substr(p + (arg[p]=='='))});
}

void parse_tokens(int argc, char *argv[], vector<token> &toks)
{
  bool opt = true;
  for (int i = 0 ; i < argc ; ++i) {
    string arg = argv[i];
    if (arg == "--") // stop parsing options
      opt = false;
    else if (!opt || arg == "" || arg == "-" || arg[0] != '-') // a regular argument
      toks.push_back({tok_arg, arg});
    else {
      char c = arg[1];
      if (c == '.' || (c >= '0' && c <= '9')) // not an option: a negative number
        toks.push_back({tok_arg, arg});
      else if (c == '-') // a long option starting with --
        add_token(arg, toks, true);
      else // a short option starting with -
        add_token(arg, toks, false);
    }
  }
}

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

class option_desc_base {
public:
  string short_name;
  string long_name;
  bool want_arg;
  string help;
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
        short_print += ",";
    }
  }

  virtual void process() { used = true; }
  virtual bool process(const string &arg) { used = true; return true; }

  void show(ostream &os, int w_short, int w_long)
  {
    os << "  " << left << setw(w_short) << short_print
       << setw(w_long) << long_print << "  ";
    show_paragraph(os, w_short + w_long + 4, help);
  }
};

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

template <typename T>
bool parse_list(const string &arg, vector<T> &var)
{
    T v;
    if (!parse_arg(arg, v))
      return false;
    var.push_back(v);
    return true;
}

template <typename T>
class set_option_desc: public option_desc_base {
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

template <typename T>
class list_option_desc: public option_desc_base {
public:
  vector<T> &var;

  list_option_desc(const string &s, const string &l, const string &an, vector<T> &v, const string &h):
  option_desc_base(s, l, an, h), var(v) {}

  bool process(const string &arg)
  { used = true; return parse_list(arg, var); }
};

class arg_desc_base {
public:
  bool need;
  bool rest;
  string arg_name;
  string help;
  bool used;

  arg_desc_base(const string &a, bool n, bool r, const string &h):
  need(n), rest(r), arg_name(a), help(h) , used(false) {}

  virtual bool process(const string &arg) = 0;
  void show(ostream &os, int w)
  {
    os << "  " << left << setw(w) << arg_name << "   ";
    show_paragraph(os, w + 5, help);
  }
};

template<typename T>
class pos_arg_desc: public arg_desc_base {
public:
  T &var;

  pos_arg_desc(const string &a, T &v, bool n, const string &h):
  arg_desc_base(a, n, false, h), var(v) {}

  bool process(const string &arg)
  { used = true; return parse_arg(arg, var); }
};

template<typename T>
class rest_arg_desc: public arg_desc_base {
public:
  vector<T> &var;

  rest_arg_desc(const string &a, vector<T> &v, const string &h):
  arg_desc_base(a, false, true, h), var(v) {}

  bool process(const string &arg)
  { used = true; return parse_list(arg, var); }
};

class parser {
public:
  string start_help, end_help, version_text, prog_name;
  string missing_arg;
  vector<token> toks;
  vector<option_desc_base *> opts;
  vector<arg_desc_base *> args;
  vector<vector<string> > exclusions;
  vector<vector<string> > selections;

  parser(const string &sh, const string &eh, const string &v, bool std_help_version = true):
  start_help(sh), end_help(eh), version_text(v)
  {
    if (std_help_version) {
      flag("h", "help", std_help_msg);
      flag("", "version", std_version_msg);
    }
  }

  void run(int argc, char *argv[]);
  void parse(int argc, char *argv[]);
  void analysis();
  void standard();
  void check();

  bool operator()(const string &s);

  bool any(const vector<string> &s)
  { return any_of(s.begin(), s.end(), *this); }

  bool all(const vector<string> &s)
  { return all_of(s.begin(), s.end(), *this); }

  bool none(const vector<string> &s)
  { return none_of(s.begin(), s.end(), *this); }

  int count(const vector<string> &s)
  { return count_if(s.begin(), s.end(), *this); }

  void show_usage(ostream &os);
  void show_help(ostream &os);
  void show_version(ostream &os)
  { os << prog_name << ": " << version_text << endl; }

  void flag(const string &s, const string &l, const string &h)
  { opts.push_back(new option_desc_base(s, l, "", h)); }

  template<typename T>
  void option(const string &s, const string &l, const string &a, T &v, const string &h)
  { opts.push_back(new set_option_desc<T>(s, l, a, v, h)); }

  template<typename T>
  void list_option(const string &s, const string &l, const string &a, vector<T> &v, const string &h)
  { opts.push_back(new list_option_desc<T>(s, l, a, v, h)); }

  template<typename T>
  void arg(const string &a, T &v, const string &h)
  {
    if (!args.empty() && !args.back()->need)
      die(bad_opt_arg + a);
    args.push_back(new pos_arg_desc<T>(a, v, true, h));
  }

  template<typename T>
  void opt_arg(const string &a, T &v, const string &h)
  { args.push_back(new pos_arg_desc<T>(a, v, false, h)); }

  template<typename T>
  void rest_arg(const string &a, std::vector<T> &v, const string &h)
  { args.push_back(new rest_arg_desc<T>(a, v, h)); }

  void exclusion(const vector<string> &s)
  { exclusions.push_back(s); }

  void selection(const vector<string> &s)
  { selections.push_back(s); }

  void warn(const string &w)
  { cerr << prog_name << ": " << color_blue << w << color_reset << endl; }

  [[ noreturn ]] void die(const string &err);
  [[ noreturn ]] void die(const string &err, const vector<string> &l);
};


void parser::show_usage(ostream &os)
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
  os << endl;
}

void parser::show_help(ostream &os)
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

  show_usage(os);
  os << endl << start_help << endl << endl;
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

void parser::run(int argc, char *argv[])
{
  parse(argc, argv);
  analysis();
  standard();
  check();
}

void parser::parse(int argc, char *argv[])
{
  if (prog_name.empty())
    prog_name = argv[0];
  parse_tokens(argc - 1, argv + 1, toks);
}

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
        if (o->long_name.find(tk.s)==0) {
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

void parser::standard()
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

void parser::check()
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

bool parser::operator()(const string &s)
{
  for (auto &opt: opts)
    if (s == opt->short_name || s == opt->long_name)
      return opt->used;
  for (auto &arg: args)
    if (s == arg->arg_name)
      return arg->used;
  die(bad_test_msg + s);
}

void parser::die(const string &err)
{
  cerr << prog_name << ": " << color_red << err << color_reset << endl;
  show_usage(cerr);
  if (!more_help_msg.empty())
    cerr << more_help_msg << endl;
  exit(-1);
}

void parser::die(const string &err, const vector<string> &l)
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
