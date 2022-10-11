/** \file iotools.cpp
  Implementation of iotools.hpp
  \author J. Houdayer
*/

#include "iotools.hpp"

const std::string cannot_open_msg = "Cannot open file ";
const std::string invalid_file_msg = "Cannot read file ";
const std::string bad_file_msg = "Something bad happens when reading file ";
const std::string unexpect_eof_msg = "Unexpected end of file ";

/** Marks a stream as failed and returns it. */
std::istream &failed(std::istream &is)
{
  is.setstate(std::ios_base::failbit);
  return is;
}

elapsed::elapsed():
start_time(std::chrono::system_clock::now())
{}

/** returns the time in seconds since the creation of the object. */ 
double elapsed::seconds() const
{
  std::chrono::duration<double> diff_time =  std::chrono::system_clock::now() - start_time;
  return diff_time.count();
}

/** Creates a progression object.
 @param sz The number of steps in the progression.
 @param show Whether to show the progression.
 */
progression::progression(size_t sz, bool show):
timer(), size(sz), step(0),
old_sec(-1), silent(!show)
{
  if (!silent)
    std::cerr << "Starting 0/" << sz;
}

progression::~progression()
{
  if (!silent) {
    const double sec = timer.seconds();

    std::cerr << "\r";
    std::cerr << "Finished " << size << " steps in ";
    if (sec >= 1)
      std::cerr << ((int) (100 * sec)) / 100. << " s";
    else
      std::cerr << (int) (1000 * sec) << " ms";
    std::cerr << "\033[K" << std::endl;
  }
}

/** Advances the progression by one step.*/
void progression::progress(const std::string &s)
{
  ++step;
  if (!silent) {
    const double sec = timer.seconds();
    const int rest = (int) ((sec / step) * (size - step) + 0.5);
    const int percent = (int)(100 * double(step) / size + 0.5);
    if (sec > old_sec + 0.1) {
      std::cerr << "\r";
      std::cerr << step << "/" << size << ": "
                << percent
                << "% (" << rest << " s)";
      std::cerr << s << "\033[K";
      old_sec = sec;
    }
  }
}

/** Creates a smart_input from a filename.
 Uses cin if filename is set to "-".
 The created file is properly closed at destruction.
*/
smart_input::smart_input(const std::string &n):
line_count(0), name(n), input(NULL), resend(false), line(""), file(NULL)
{
  if (name == "-") {
    name = "standard input";
    input = &std::cin;
  }
  else {
    file = new std::ifstream(name);
    input = file;
  }
}

/** Make a smart_input from an istream with the given name (for error messages).*/
smart_input::smart_input(std::istream &is, const std::string &n):
line_count(0), name(n), input(&is), resend(false), line(""), file(NULL)
{}

smart_input::~smart_input()
{
  if (file != NULL) {
    file->close();
    delete file;
  }
}

/** To check whether the input is Ok.*/
smart_input::operator bool() const
{
  return bool(*input);
}

/** Returns the next line of the file in iss (properly reset).
 Ignores empty lines and lines starting with '#'. 
 Always check the status of the smart_input before using it.
*/
smart_input &smart_input::next_line(std::istringstream &iss)
{
  if (resend) {
    resend = false;
    iss.clear();
    iss.str(line);
    iss >> std::ws;
  }
  else
    while (getline(*input, line)) {
      line_count++;
      line.push_back(' ');
      iss.clear();
      iss.str(line);
      iss >> std::ws;
      if (!iss.eof() && iss.peek()!='#')
        break;
    }
  return *this;
}

/** Look at the next line of the file without actually removing it. */
smart_input &smart_input::peek_line(std::istringstream &iss)
{
  next_line(iss);
  resend = true;
  return *this;
}

/** marks the smart_input as failed. */
smart_input &smart_input::failed()
{
  input->setstate(std::ios_base::failbit);
  return *this;
}

/** Creates a smart_output from a filename.
 Uses cout if filename is set to "-".
 The created file is properly closed at destruction.
*/
smart_output::smart_output(const std::string &n):
name(n), output(NULL), file(NULL)
{
  if (name == "-") {
    name = "standard output";
    output = &std::cout;
  }
  else {
    file = new std::ofstream(name);
    output = file;
  }
}

smart_output::~smart_output()
{
  if (file != NULL) {
    file->close();
    delete file;
  }
}

/** To check whether the output is Ok.*/
smart_output::operator bool() const
{
  return bool(*output);
}

