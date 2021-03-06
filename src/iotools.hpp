/** \file iotools.hpp
  Small tools to help reading and writing files.
  \author J. Houdayer
*/

#ifndef IOTOOLS_HPP
#define IOTOOLS_HPP

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <chrono>


extern const std::string cannot_open_msg;
extern const std::string invalid_file_msg;
extern const std::string bad_file_msg;
extern const std::string unexpect_eof_msg;

std::istream &failed(std::istream &is);

class elapsed
{
public:
  const std::chrono::time_point<std::chrono::system_clock> start_time;

  elapsed();
  double seconds() const;
};

class progression
{
public:
  elapsed timer;
  size_t size, step;
  int old_percent, old_rest;
  bool silent;

  progression(size_t sz, bool show);
  ~progression();
  void progress(const std::string &s ="");
};

class smart_input
{
public:
  smart_input(const std::string &name);
  smart_input(std::istream &is, const std::string &name);
  ~smart_input();
  smart_input(const smart_input &) = delete;
  smart_input &operator=(const smart_input &) = delete;

  bool bad() const
  { return input->bad(); }

  bool fail() const
  {return input->fail(); }

  bool eof() const
  { return input->eof(); }

  void clear()
  { input->clear(); }

  explicit operator bool() const;
  smart_input &next_line(std::istringstream &iss);
  smart_input &peek_line(std::istringstream &iss);
  smart_input &failed();

  size_t line_count;
  std::string name;
  std::istream *input;

private:
  bool resend;
  std::string line;
  std::ifstream *file;
};

template <typename T>
smart_input &operator>>(smart_input &is, T &x)
{
  *is.input >> x;
  return is;
}

template<typename T>
std::string read_object(smart_input &is, T &x,  bool verbose = false)
{
  if (!is)
    return cannot_open_msg + is.name + " (" + strerror(errno) + ")";
  if (verbose)
    std::cerr << "Reading file " << is.name << "...";
  is >> x;
  if (verbose)
    std::cerr << "Done" << std::endl;
  if (is.bad())
    return bad_file_msg + is.name + " (" + strerror(errno) + ")";
  if (is.eof())
    return unexpect_eof_msg + is.name;
  if (is.fail())
    return invalid_file_msg + is.name + " at line " + std::to_string(is.line_count);
  return "";
}

template <typename T>
std::string read_file(const std::string &filename, T &x,  bool verbose = false)
{
  smart_input is(filename);
  return read_object(is, x, verbose);
}

#endif