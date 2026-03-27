#include <iomanip>
#include <sstream>
#include <iostream>
#include "myomp.h"
#include "progressbar.h"

#ifdef USING_R
#include <Rcpp.h>
#endif

Progress::Progress(std::size_t total, std::string prefix, double interval)
  : done_(0),
    total_(total),
    prefix_(std::move(prefix)),
    last_percent_(0),
    last_width_(0),
    start_(clock::now()),
    last_print_(start_),
    min_interval_(interval)
{}

Progress::~Progress()
{
  finalize();
}

// ultra-cheap for workers
void Progress::tick() noexcept
{
  done_.fetch_add(1, std::memory_order_relaxed);
  update();
}

bool Progress::check_interrupt()
{
  if (omp_get_thread_num() != 0)
    return false; // # nocov

  const std::size_t done = done_.load(std::memory_order_relaxed);
  if (done % 10000 != 0) return false;

  #ifdef USING_R
  try
  {
    Rcpp::checkUserInterrupt();
  }
  catch (Rcpp::internal::InterruptedException e)
  {
    return true;
  }
  #endif

  return false;
}

void Progress::update()
{
  if (omp_get_thread_num() != 0)
    return;

  const auto now = clock::now();
  const double elapsed = std::chrono::duration<double>(now - start_).count();

  if (elapsed < min_interval_)
    return;

  const std::size_t done = done_.load(std::memory_order_relaxed);
  if (done == 0)
    return;

  const std::size_t percent = (done * 100) / total_;
  if (percent == last_percent_)
    return;

  last_percent_ = percent;
  last_print_ = now;

  // ETA
  const double rate = done / elapsed;
  const double remaining = (percent < 100) ? (total_ - done) / rate : 0;

  std::ostringstream os;
  os << prefix_ << ": "
     << percent << "% | "
     << std::fixed << std::setprecision(1)
     << remaining << "s ETA ("
     << omp_get_num_threads() << " threads)\r";

  const std::string msg = os.str();
  last_width_ = msg.size();

  #ifdef USING_R
    Rcpp::Rcout << msg;
    Rcpp::Rcout.flush();
  #else
    std::cout << msg;
    std::cout.flush();
  #endif
}

void Progress::finalize()
{
  if (omp_get_thread_num() != 0) return;
  if (finalized_) return;
  finalized_ = true;

  #ifdef USING_R
    Rcpp::Rcout << '\r' << std::string(last_width_, ' ') << '\r';
    Rcpp::Rcout.flush();
  #else
    std::cout << '\r' << std::string(last_width_, ' ') << '\r';
    std::cout.flush();
  #endif
}
