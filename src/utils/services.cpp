/**
 * @file services.cpp
 * Project: Arbor
 * 
 * Copyright (C) 2026 Jean-Romain Roussel (r-lidar) <info @ r-lidar.com>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

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
