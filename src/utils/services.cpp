#include "services.h"
#include "progressbar.h"

namespace {

Logger& logger_store()
{
  static Logger l = [](const std::string&) {};
  return l;
}

ProgressBarFactory& factory_store()
{
  static ProgressBarFactory f = [](std::size_t total, const std::string& prefix, double interval) -> std::unique_ptr<ProgressBar> {
    return std::make_unique<Progress>(total, prefix, interval);
  };
  return f;
}

} // namespace

void ServiceLocator::register_logger(Logger logger)
{
  logger_store() = std::move(logger);
}

void ServiceLocator::register_progress_bar(ProgressBarFactory factory)
{
  factory_store() = std::move(factory);
}

const Logger& ServiceLocator::logger()
{
  return logger_store();
}

std::unique_ptr<ProgressBar> ServiceLocator::make_progress(std::size_t total, const std::string& prefix, double interval)
{
  return factory_store()(total, prefix, interval);
}
