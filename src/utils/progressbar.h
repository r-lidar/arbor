#ifndef Progress_H
#define Progress_H

#include <atomic>
#include <string>
#include <chrono>
#include "services.h"

class Progress : public ProgressBar
{
public:
  Progress(std::size_t total, std::string prefix = "", double interval = 0.5);
  ~Progress() override;

  // Called by ALL threads
  void tick() noexcept override;

  // Called by master thread only
  void update() override;
  void finalize() override;
  bool check_interrupt() override;

private:
  using clock = std::chrono::steady_clock;

  std::atomic<std::size_t> done_;
  const std::size_t total_;

  std::string prefix_ = "";
  std::size_t last_percent_ = 0;
  std::size_t last_width_ = 0;
  bool finalized_ = false;

  clock::time_point start_;
  clock::time_point last_print_;

  double min_interval_; // seconds
};

#endif
