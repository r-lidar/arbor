/**
 * @file progress.cpp
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

#include <algorithm>
#include <atomic>
#include <chrono>
#include <exception>
#include "myomp.h"
#include "progress.h"

namespace arbor {

namespace detail {

using clock = std::chrono::steady_clock;

// Minimal delay between two progress reports and two cancellation polls (seconds)
constexpr double interval = 0.1;

struct Root
{
  clock::time_point start;
  clock::time_point last_emit;
  double last_fraction = 0;
  int uncaught = 0;
  std::atomic<bool> cancelled{false};
};

struct Frame
{
  std::string label;
  bool counter = false;

  // Share of the global progress allotted to this frame
  double lo = 0;
  double hi = 1;

  // Tasks
  std::vector<double> cum; // Cumulative normalized weights. size = nsteps + 1
  long index = -1;
  std::string step;

  // Counters
  std::size_t total = 0;
  std::atomic<std::size_t> done{0};

  Root* root = nullptr;
  std::unique_ptr<Root> owned_root; // Only for the root frame

  std::size_t nsteps() const { return cum.size() - 1; }
  std::size_t current_step() const { return static_cast<std::size_t>(std::clamp<long>(index, 0, static_cast<long>(nsteps()) - 1)); }

  double local() const
  {
    if (counter) return total ? std::min(1.0, static_cast<double>(done.load(std::memory_order_relaxed)) / total) : 0.0;
    if (index < 0) return 0.0;
    return cum[current_step()];
  }

  double position() const { return lo + (hi - lo) * local(); }
};

// One task hierarchy per host thread. OpenMP thread 0 of a parallel region is the
// thread that opened the region, so it sees the stack of the computation.
thread_local std::vector<Frame*> stack;
thread_local int muted = 0;

static ProgressInfo snapshot(Root& root)
{
  ProgressInfo info;

  const double global = stack.empty() ? root.last_fraction : stack.back()->position();
  root.last_fraction = std::max(root.last_fraction, global);

  info.fraction = root.last_fraction;
  info.elapsed  = std::chrono::duration<double>(clock::now() - root.start).count();
  info.eta      = (info.fraction > 0.01 && info.elapsed > 1) ? info.elapsed * (1 - info.fraction) / info.fraction : -1;

  info.stack.reserve(stack.size());
  for (const Frame* f : stack)
  {
    ProgressFrame pf;
    pf.label = f->label;
    pf.step  = f->step;
    if (f->counter)
    {
      pf.current = std::min(f->done.load(std::memory_order_relaxed), f->total);
      pf.total   = f->total;
    }
    else
    {
      pf.current = (f->index < 0) ? 0 : f->current_step() + 1;
      pf.total   = f->nsteps();
    }
    pf.fraction = (f->hi > f->lo) ? std::clamp((global - f->lo) / (f->hi - f->lo), 0.0, 1.0) : 0.0;
    info.stack.push_back(std::move(pf));
  }

  return info;
}

static void emit(Root& root, bool force)
{
  const auto now = clock::now();
  if (!force && std::chrono::duration<double>(now - root.last_emit).count() < interval) return;
  root.last_emit = now;

  Reporter& reporter = ServiceLocator::reporter();
  if (reporter.cancel_requested()) root.cancelled = true;
  reporter.progress(snapshot(root));
}

// Returns nullptr if the frame must be inert
static std::unique_ptr<Frame> make_frame(std::string&& label, bool counter)
{
  if (omp_in_parallel()) return nullptr;
  if (!stack.empty() && stack.back()->counter) return nullptr;

  auto f = std::make_unique<Frame>();
  f->label   = std::move(label);
  f->counter = counter;

  if (stack.empty())
  {
    f->owned_root = std::make_unique<Root>();
    f->owned_root->start = clock::now();
    f->owned_root->last_emit = f->owned_root->start;
    f->owned_root->uncaught = std::uncaught_exceptions();
    f->root = f->owned_root.get();
  }
  else
  {
    const Frame& p = *stack.back();
    const double w = p.hi - p.lo;
    const std::size_t i = p.current_step();
    f->root = p.root;
    f->lo = p.lo + w * p.cum[i];
    f->hi = p.lo + w * p.cum[i + 1];
  }

  return f;
}

static void push(Frame* f)
{
  stack.push_back(f);

  try
  {
    if (f->owned_root)
      ServiceLocator::reporter().begin(snapshot(*f->root));
    if (f->counter)
      emit(*f->root, true);
  }
  catch (...) {}
}

static void pop(Frame* f) noexcept
{
  auto it = std::find(stack.begin(), stack.end(), f);
  if (it == stack.end()) return;
  stack.erase(it, stack.end());

  if (!f->owned_root) return;

  try
  {
    Root& root = *f->root;
    const bool success = std::uncaught_exceptions() <= root.uncaught;
    ProgressInfo info = snapshot(root);
    if (success)
    {
      info.fraction = 1;
      info.eta = 0;
    }
    info.stack.push_back({f->label, f->step, 0, 0, info.fraction});
    ServiceLocator::reporter().end(info, success);
  }
  catch (...) {}
}

} // namespace detail

// ---- Task -------------------------------------------------------------

Task::Task(std::string label, std::size_t nsteps)
  : Task(std::move(label), std::vector<double>(std::max<std::size_t>(nsteps, 1), 1.0))
{}

Task::Task(std::string label, std::vector<double> weights)
{
  frame = detail::make_frame(std::move(label), false);
  if (!frame) return;

  if (weights.empty()) weights.push_back(1);

  double sum = 0;
  for (double w : weights) sum += std::max(w, 0.0);

  frame->cum.resize(weights.size() + 1, 0.0);
  for (std::size_t i = 0; i < weights.size(); ++i)
  {
    const double w = (sum > 0) ? std::max(weights[i], 0.0) / sum : 1.0 / weights.size();
    frame->cum[i + 1] = frame->cum[i] + w;
  }
  frame->cum.back() = 1;

  detail::push(frame.get());
}

Task::~Task()
{
  if (frame) detail::pop(frame.get());
}

void Task::step(const std::string& label)
{
  if (!frame) return;

  frame->index++;
  frame->step = label;

  detail::emit(*frame->root, true);
  if (frame->root->cancelled) throw Cancelled();
}

void Task::rewind(std::size_t next_step)
{
  if (!frame) return;
  frame->index = static_cast<long>(next_step) - 1;
}

// ---- Counter ----------------------------------------------------------

Counter::Counter(std::string label, std::size_t total)
{
  frame = detail::make_frame(std::move(label), true);
  if (!frame) return;

  frame->total = total;
  detail::push(frame.get());
}

Counter::~Counter()
{
  if (frame) detail::pop(frame.get());
}

void Counter::tick(std::size_t n) noexcept
{
  if (!frame) return;

  frame->done.fetch_add(n, std::memory_order_relaxed);

  if (omp_get_thread_num() == 0)
  {
    try
    {
      detail::emit(*frame->root, false);
    }
    catch (...) {}
  }

  if (frame->root->cancelled.load(std::memory_order_relaxed))
    cancel.store(true, std::memory_order_relaxed);
}

void Counter::throw_if_cancelled() const
{
  if (cancelled()) throw Cancelled();
}

// ---- Logging ----------------------------------------------------------

namespace logging {

static void emit(LogLevel level, const std::string& msg)
{
  if (detail::muted > 0 || omp_get_thread_num() != 0) return;
  ServiceLocator::reporter().message(level, msg);
}

void debug(const std::string& msg)   { emit(LogLevel::Debug, msg); }
void info(const std::string& msg)    { emit(LogLevel::Info, msg); }
void warning(const std::string& msg) { emit(LogLevel::Warning, msg); }

Mute::Mute()  { ++detail::muted; }
Mute::~Mute() { --detail::muted; }

}

}
