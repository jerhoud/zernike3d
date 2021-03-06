#include "iotools.hpp"

const std::string cannot_open_msg = "Cannot open file ";
const std::string invalid_file_msg = "Cannot read file ";
const std::string bad_file_msg = "Something bad happens when reading file ";
const std::string unexpect_eof_msg = "Unexpected end of file ";

/** Marks the stream as failed and returns it. */
std::istream &failed(std::istream &is)
{
  is.setstate(std::ios_base::failbit);
  return is;
}

elapsed::elapsed():
start_time(std::chrono::system_clock::now())
{}

double elapsed::seconds() const
{
  std::chrono::duration<double> diff_time =  std::chrono::system_clock::now() - start_time;
  return diff_time.count();
}

progression::progression(size_t sz, bool show):
timer(), size(sz), step(0),
old_percent(-1), old_rest(-1), silent(!show)
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

void progression::progress(const std::string &s)
{
  ++step;
  if (!silent) {
    const double sec = timer.seconds();
    const int rest = (int) ((sec / step) * (size - step) + 0.5);
    const int percent = (int)(100 * double(step) / size + 0.5);
    if ((rest != old_rest) || (percent != old_percent)) {
      std::cerr << "\r";
      std::cerr << step << "/" << size << ": "
                << percent
                << "% (" << rest << " s)";
      std::cerr << s << "\033[K";
      old_rest = rest;
      old_percent = percent;
    }
  }
}

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

smart_input::operator bool() const
{
  return bool(*input);
}

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

smart_input &smart_input::peek_line(std::istringstream &iss)
{
  next_line(iss);
  resend = true;
  return *this;
}

smart_input &smart_input::failed()
{
  input->setstate(std::ios_base::failbit);
  return *this;
}
