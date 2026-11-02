/**
 * @file progress.h
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

#ifndef PROGRESS_H
#define PROGRESS_H

#include <atomic>
#include <memory>
#include <string>
#include <vector>
#include "services.h"

// Instrumentation of the core. Tasks and counters are RAII scopes that nest:
// a task or counter opened while a step of a task is running takes the share of
// the global progress allotted to this step. The first scope opened is the root
// and its progress is what the host receives in ProgressInfo::fraction.
//
//   Task task("Semantic segmentation", {1, 60, 8, 8, 10, 8, 1});
//   task.step("Partitioning");
//   ...
//   task.step("Accumulating passages");
//   accumulate_passages(...);     // opens its own Task, mapped into this step
//
//   Counter counter("SOR", n);
//   #pragma omp parallel for
//   for (size_t i = 0 ; i < n ; i++)
//   {
//     if (counter.cancelled()) continue;
//     counter.tick_at(i);
//     ...
//   }
//   counter.throw_if_cancelled();
//
// Rules
// - Tasks and counters must be created outside OpenMP parallel regions. Those created
//   inside a parallel region or inside a counter are inert (e.g. qsm() called per tree
//   in qsf()): they report nothing and never throw.
// - Task::step() throws arbor::Cancelled when the host requested a cancellation.
// - Counter::tick(), Counter::tick_at() and Counter::cancelled() can be called from any thread.
// - In hot loops use tick_at(i), which does something only every Counter::stride iterations.
//   tick() costs an atomic increment and is meant for coarse loops (e.g. one iteration per tree).
// - Messages from OpenMP worker threads are discarded.

namespace arbor
{

namespace detail { struct Frame; }

// A unit of work made of sequential labelled steps
class Task
{
public:
  // nsteps steps of equal weight
  explicit Task(std::string label, std::size_t nsteps = 1);
  // One step per weight. Weights are relative and do not need to sum to 1
  Task(std::string label, std::vector<double> weights);
  ~Task();

  Task(const Task&) = delete;
  Task& operator=(const Task&) = delete;

  // Enters the next step. Throws Cancelled if a cancellation was requested
  void step(const std::string& label);

  // Sets the step (0-based) entered by the next call to step(), e.g. to replay
  // steps on a retry. The global progress reported never decreases.
  void rewind(std::size_t next_step = 0);

private:
  std::unique_ptr<detail::Frame> frame; // nullptr if inert
};

// A loop of known size
class Counter
{
public:
  Counter(std::string label, std::size_t total);
  ~Counter();

  Counter(const Counter&) = delete;
  Counter& operator=(const Counter&) = delete;

  // Number of iterations between two updates in hot loops. Power of 2 for a cheap test.
  static constexpr std::size_t stride = std::size_t(1) << 16;

  // Hot loops. Call at every iteration with the iteration index i. Any thread.
  // Advances the counter by `stride` once every `stride` iterations and does nothing otherwise.
  void tick_at(std::size_t i) noexcept { if ((i & (stride - 1)) == 0) tick(stride); }

  // Coarse loops. Advances the counter by n. Any thread.
  void tick(std::size_t n = 1) noexcept;

  // Any thread. Cheap: reads a flag refreshed by tick()
  bool cancelled() const noexcept { return cancel.load(std::memory_order_relaxed); }

  // Outside parallel regions. Throws Cancelled if a cancellation was requested
  void throw_if_cancelled() const;

private:
  std::unique_ptr<detail::Frame> frame; // nullptr if inert
  std::atomic<bool> cancel{false};
};

namespace logging
{
  void debug(const std::string& msg);
  void info(const std::string& msg);
  void warning(const std::string& msg);

  // Discards all messages while alive
  class Mute
  {
  public:
    Mute();
    ~Mute();
    Mute(const Mute&) = delete;
    Mute& operator=(const Mute&) = delete;
  };
}

}

#endif
