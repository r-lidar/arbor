#ifndef SERVICES_H
#define SERVICES_H

#include <functional>
#include <memory>
#include <string>
#include <cstddef>

using Logger = std::function<void(const std::string&)>;

class ProgressBar
{
public:
  virtual ~ProgressBar() = default;

  // Called by ALL threads
  virtual void tick() noexcept = 0;

  // Called by master thread only
  virtual void update() = 0;
  virtual void finalize() = 0;
  virtual bool check_interrupt() = 0;
};

using ProgressBarFactory = std::function<std::unique_ptr<ProgressBar>(std::size_t total, const std::string& prefix, double interval)>;

class ServiceLocator
{
public:
  static void register_logger(Logger logger);
  static void register_progress_bar(ProgressBarFactory factory);

  static const Logger& logger();
  static std::unique_ptr<ProgressBar> make_progress(std::size_t total, const std::string& prefix = "", double interval = 0.5);
};

#endif
