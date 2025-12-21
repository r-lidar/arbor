#include <atomic>
#include <string>
#include <chrono>

class Progress
{
public:
  Progress(std::size_t total, std::string prefix = "", double interval = 0.5);
  ~Progress();

  // Called by ALL threads
  void tick() noexcept;

  // Called by master thread only
  void update();
  void finalize();
  bool check_interrupt();

private:
  using clock = std::chrono::steady_clock;

  std::atomic<std::size_t> done_;
  const std::size_t total_;

  std::string prefix_;
  std::size_t last_percent_;

  clock::time_point start_;
  clock::time_point last_print_;

  double min_interval_; // seconds
};

