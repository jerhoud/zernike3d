#include "iotools.hpp"

int line_count = 0;

const std::string cannot_open_msg = "Cannot open file ";
const std::string invalid_file_msg = "Cannot read file ";
const std::string bad_file_msg = "Something bad happens when reading file ";
const std::string unexpect_eof_msg = "Unexpected end of file ";

/** Reads a line from a stream.
  Ignores initial white spaces, blank lines and lines starting with '#'.
  @param is The stream to read.
  @param iss The stream where to put the line.
*/
std::istream &next_line(std::istream &is, std::istringstream &iss)
{
  std::string s;
  while (getline(is, s)) {
    line_count++;
    s.push_back(' ');
    iss.clear();
    iss.str(s);
    iss >> std::ws;
    if (!iss.eof() && iss.peek()!='#')
      break;
  }
  return is;
}

/** Marks the stream as failed and returns it. */
std::istream &failed(std::istream &is)
{
  is.setstate(std::ios_base::failbit);
  return is;
}
