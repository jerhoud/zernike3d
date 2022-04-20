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


extern int line_count;
extern const std::string cannot_open_msg;
extern const std::string invalid_file_msg;
extern const std::string bad_file_msg;
extern const std::string unexpect_eof_msg;

std::istream &next_line(std::istream &is, std::istringstream &iss);
std::istream &failed(std::istream &is);

/** reads a stream into an object.
  It reads the object using operator >>.
  It reports errors with line count if possible.
  @param name is the name of the stream for error reporting
  @param is the stream to read
  @param x is the variable to read the stream in
  @return "" if no error happens and a helpfull message otherwise
*/
template<typename T>
std::string read_stream(const std::string &name, std::istream &is, T &x, bool verbose = false)
{
  if (verbose)
    std::cerr << "Reading file " << name << "...";
  line_count = 0;
  if (!is)
    return cannot_open_msg + name + " (" + strerror(errno) + ")";
  is >> x;
  if (is.bad())
    return bad_file_msg + name + " (" + strerror(errno) + ")";
  if (is.eof())
    return unexpect_eof_msg + name;
  if (is.fail())
    return invalid_file_msg + name + " at line " + std::to_string(line_count);
  if (verbose)
    std::cerr << " Done\n"; 
  return "";
}

/** reads a file into an object of type T.
  It reads the object using operator >>.
  It reports errors with line count if possible.
  @param filename is the name of the file to read ("-" for standard input)
  @param x is the variable to read the stream in
  @return "" if no error happens and a helpfull message otherwise
*/

template<typename T>
std::string read_file(const std::string &filename, T &x, bool verbose = false)
{
  if (filename=="-")
    return read_stream("standard input", std::cin, x, verbose);
  else {
    std::ifstream file(filename);
    std::string err = read_stream(filename, file, x, verbose);
    file.close();
    return err;
  }
}

class progression
{
public:
  time_t start_time, current_time;
  size_t size, step;
  bool silent;

  progression(size_t sz, bool show);
  ~progression();
  void progress(const std::string &s ="");
};

#endif