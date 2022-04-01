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

template<typename T>
std::string read_stream(const std::string &name, std::istream &is, T &x)
{
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
  return "";
}

template<typename T>
std::string read_file(const std::string &filename, T &x)
{
  if (filename=="-")
    return read_stream("standard input", std::cin, x);
  else {
    std::ifstream file(filename);
    std::string err = read_stream(filename, file, x);
    file.close();
    return err;
  }
}

#endif