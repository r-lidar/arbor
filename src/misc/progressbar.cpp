#include <Rcpp.h>
#include <iomanip>
#include <sstream>
#include "myomp.h"
#include "progressbar.h"

Progress::Progress(std::size_t total, std::string prefix)
  : done_(0),
    total_(total),
    prefix_(std::move(prefix)),
    last_percent_(0),
    start_(clock::now()),
    last_print_(start_),
    min_interval_(0.5)
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
  if(omp_get_thread_num() != 0)
    return false; // # nocov

  const std::size_t done = done_.load(std::memory_order_relaxed);
  if(done % 10000 != 0) return false;

  try
  {
    Rcpp::checkUserInterrupt();
  }
  catch(Rcpp::internal::InterruptedException e)
  {
    return true;
  }

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
  os << prefix_
     << std::setw(3) << percent << "% | "
     << std::fixed << std::setprecision(1)
     << remaining << "s ETA ("
     << omp_get_num_threads() << " threads)\r";

  Rcpp::Rcout << os.str();
  Rcpp::Rcout.flush();
}

void Progress::finalize()
{
  if (omp_get_thread_num() != 0)
    return;

  Rcpp::Rcout << "\n";
  Rcpp::Rcout.flush();
}

