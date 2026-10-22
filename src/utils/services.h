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

#include <cstddef>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

// This file is the interface between the arbor core and the host application
// (R package, Qt GUI, command line tool). The core reports what it is doing
// through a Reporter registered by the host. The core itself is instrumented
// with arbor::Task, arbor::Counter and arbor::logging (see progress.h).
//
// Threading: every Reporter callback is invoked from the thread that runs the
// computation (never from OpenMP worker threads). A GUI running the computation
// in a worker thread must marshal the data to its UI thread (e.g. queued signals).

namespace arbor
{

enum class LogLevel { Debug, Info, Warning };

// One level of the task hierarchy
struct ProgressFrame
{
  std::string label;        // Name of the task (e.g. "Semantic segmentation", "SOR")
  std::string step;         // Label of the current step. Empty for counters or before the first step
  std::size_t current = 0;  // Current step (1-based) for tasks, items processed for counters
  std::size_t total = 0;    // Number of steps for tasks, number of items for counters
  double fraction = 0;      // Progress of this level in [0,1], including its sub-tasks
};

struct ProgressInfo
{
  double fraction = 0;              // Overall progress of the root task in [0,1]. Never decreases.
  double elapsed = 0;               // Seconds since the root task started
  double eta = -1;                  // Estimated remaining seconds. Negative if unknown
  std::vector<ProgressFrame> stack; // Active tasks, root first, innermost last
};

class Reporter
{
public:
  virtual ~Reporter() = default;

  // A root task started, progressed or ended. progress() is throttled by the core
  // but is always called when a step starts.
  virtual void begin(const ProgressInfo&) {}
  virtual void progress(const ProgressInfo& info) = 0;
  virtual void end(const ProgressInfo&, bool /*success*/) {}

  // Leveled message
  virtual void message(LogLevel level, const std::string& msg) = 0;

  // Polled regularly while a task is running. Returning true aborts the
  // computation, which then throws arbor::Cancelled.
  virtual bool cancel_requested() { return false; }
};

// Thrown by the core when Reporter::cancel_requested() returned true
class Cancelled : public std::runtime_error
{
public:
  Cancelled() : std::runtime_error("Operation cancelled by user") {}
};

// Reports nothing
class NullReporter : public Reporter
{
public:
  void progress(const ProgressInfo&) override {}
  void message(LogLevel, const std::string&) override {}
};

// Text mode reporter. Displays a single live line updated in place, e.g.
// "Semantic segmentation:  34% ETA 12s   | Accumulating passages"
// Messages are printed on their own timestamped lines above the live line.
class ConsoleReporter : public Reporter
{
public:
  explicit ConsoleReporter(std::ostream& os, LogLevel min_level = LogLevel::Info, bool color = true, std::size_t width = 80);

  void progress(const ProgressInfo& info) override;
  void end(const ProgressInfo& info, bool success) override;
  void message(LogLevel level, const std::string& msg) override;

  // Formats a progress state as a one-line string of at most `width` characters (0 = unlimited):
  // the root task, the percentage, the ETA and the current step of the root task
  static std::string format(const ProgressInfo& info, std::size_t width = 0);

protected:
  void clear_line();

  std::ostream& os;
  LogLevel min_level;
  bool color;
  std::size_t width;
  std::string line;  // Live line currently displayed
};

}

class ServiceLocator
{
public:
  static void register_reporter(std::shared_ptr<arbor::Reporter> reporter);
  static arbor::Reporter& reporter();
};

#endif
