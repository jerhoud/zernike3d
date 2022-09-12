/** \file parallel.hpp
  Classes to work with arrays in parallel, with progression bars.
  \author J. Houdayer
*/

#include <vector>
#include <functional>
#include <thread>
#include <mutex>
#include "iotools.hpp"

template<typename T>
void parallel_eval(int nt, std::vector<T> &v, std::function<T(size_t)> f, bool verbose = false)
{
  progression prog(v.size(), verbose);
  if (nt <= 1) // not parallel
    for (size_t i = 0 ; i < v.size() ; i++) {
      v[i] = f(i);
      prog.progress();
    }
  else { // parallel
    std::vector<std::thread> threads(nt);
    size_t idx = 0;
    std::mutex mtx;
    for (auto &t: threads)
      t = std::thread([&]{
        mtx.lock();
        size_t i = idx++;
        mtx.unlock();
        while (i < v.size()) {
          v[i] = f(i);
          mtx.lock();
          i = idx++;
          prog.progress();
          mtx.unlock();
        }
      });
    for (auto &t: threads)
      t.join();
  }
}

template<typename T, typename C>
const C parallel_collect(int nt, const std::vector<T> &v, const C &collector, bool verbose = false)
{
  progression prog(v.size(), verbose);
  if (nt <= 1) { // not parallel
    C c(collector);
    for (auto &x: v) {
      std::string s = c.collect(x);
      prog.progress(s);
    }
    return c;
  }
  else { // parallel
    std::vector<std::thread> threads(nt);
    std::vector<C> collectors(nt, collector);
    size_t idx = 0;
    std::mutex mtx;
    for (int it = 0 ; it < nt ; it++) {
      threads[it] = std::thread([&]{
        mtx.lock();
        size_t i = idx++;
        mtx.unlock();
        while (i < v.size()) {
          std::string s = collectors[it].collect(v[i]);
          mtx.lock();
          i = idx++;
          prog.progress(s);
          mtx.unlock();
        }
      });
    }
    for (auto &t: threads)
      t.join();
    C total(collector);
    for (auto &c: collectors)
      total.collect(c);
    return total;
  }
}