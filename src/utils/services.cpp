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

#include <chrono>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <sstream>

#ifndef USING_R
#include <iostream>
#endif

namespace {

std::shared_ptr<arbor::Reporter>& reporter_store()
{
  // R forbids writing to stdout. The R package registers its own reporter.
  #ifdef USING_R
  static std::shared_ptr<arbor::Reporter> r = std::make_shared<arbor::NullReporter>();
  #else
  static std::shared_ptr<arbor::Reporter> r = std::make_shared<arbor::ConsoleReporter>(std::cout);
  #endif
  return r;
}

std::string format_duration(double s)
{
  long t = std::lround(s);
  std::ostringstream os;
  if (t < 60)
    os << t << "s";
  else if (t < 3600)
    os << t / 60 << "m" << std::setfill('0') << std::setw(2) << t % 60 << "s";
  else
    os << t / 3600 << "h" << std::setfill('0') << std::setw(2) << (t % 3600) / 60 << "m";
  return os.str();
}

} // namespace

void ServiceLocator::register_reporter(std::shared_ptr<arbor::Reporter> reporter)
{
  reporter_store() = reporter ? std::move(reporter) : std::make_shared<arbor::NullReporter>();
}

arbor::Reporter& ServiceLocator::reporter()
{
  return *reporter_store();
}

namespace arbor {

ConsoleReporter::ConsoleReporter(std::ostream& os, LogLevel min_level, bool color, std::size_t width)
  : os(os), min_level(min_level), color(color), width(width)
{}

std::string ConsoleReporter::format(const ProgressInfo& info, std::size_t width)
{
  if (info.stack.empty()) return "";

  // Fixed-width progress first so the percentage does not move, e.g.
  // "Semantic segmentation: 34% ETA 12s   | Accumulating passages"
  // Only the root task and its step are displayed: nested steps change too fast to be read.
  const ProgressFrame& root = info.stack[0];

  std::ostringstream os;
  os << root.label << ": "
     << std::setw(3) << static_cast<int>(info.fraction * 100) << "%"
     << " ETA " << std::left << std::setw(6) << (info.eta >= 0 ? format_duration(info.eta) : "--") << std::right;

  if (!root.step.empty()) os << "| " << root.step;

  std::string line = os.str();
  if (width > 3 && line.size() > width) line = line.substr(0, width - 3) + "...";
  return line;
}

void ConsoleReporter::progress(const ProgressInfo& info)
{
  // Keep one character free: writing in the last column wraps on some terminals
  std::string s = format(info, width > 0 ? width - 1 : 0);
  if (s == line) return;

  os << '\r' << s;
  if (s.size() < line.size()) os << std::string(line.size() - s.size(), ' ');
  os.flush();

  line = std::move(s);
}

void ConsoleReporter::end(const ProgressInfo&, bool)
{
  clear_line();
}

void ConsoleReporter::message(LogLevel level, const std::string& msg)
{
  if (level < min_level) return;

  std::string live = line;
  clear_line();

  std::time_t t = std::time(nullptr);
  std::tm* tm = std::localtime(&t);

  os << "[" << std::setfill('0')
     << std::setw(2) << tm->tm_hour << ":"
     << std::setw(2) << tm->tm_min  << ":"
     << std::setw(2) << tm->tm_sec  << "] ";

  if (level == LogLevel::Warning)
  {
    if (color) os << "\033[33m";
    os << "Warning: " << msg;
    if (color) os << "\033[0m";
  }
  else
  {
    os << msg;
  }

  os << '\n';

  // Redraw the live line below the message
  if (!live.empty())
  {
    os << live;
    line = std::move(live);
  }

  os.flush();
}

void ConsoleReporter::clear_line()
{
  if (line.empty()) return;
  os << '\r' << std::string(line.size(), ' ') << '\r';
  os.flush();
  line.clear();
}

}
