/**
 * @file progressbar.h
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
