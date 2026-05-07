/**
 * @file services.h
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
