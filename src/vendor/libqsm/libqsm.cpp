/**
 * @file libqsm.cpp
 * Project: Arbor
 *
 * MIT License
 *
 * Copyright (C) 2026 Jean-Romain Roussel (r-lidar) <info @ r-lidar.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "libqsm.h"

#include <fstream>
#include <sstream>
#include <iomanip>    // std::setprecision, std::fixed
#include <ctime>      // std::time, std::gmtime, std::strftime
#include <cstring>    // memset
#include <algorithm>
#include <cassert>
#include <functional>


namespace libqsm {

QSMnode read_node(std::ifstream& in, const QSMheader& hdr)
{
  float fx  = 0.0f, fy = 0.0f, fz = 0.0f;

  in.read(reinterpret_cast<char*>(&fx),  sizeof(float));
  in.read(reinterpret_cast<char*>(&fy),  sizeof(float));
  in.read(reinterpret_cast<char*>(&fz),  sizeof(float));

  QSMnode n;
  n.x  = static_cast<double>(fx) + hdr.x_offset;
  n.y  = static_cast<double>(fy) + hdr.y_offset;
  n.z  = static_cast<double>(fz) + hdr.z_offset;
  return n;
}

void write_node(std::ofstream& out, const QSMnode& n, const QSMheader& hdr)
{
  const float   fx = static_cast<float>(n.x - hdr.x_offset);
  const float   fy = static_cast<float>(n.y - hdr.y_offset);
  const float   fz = static_cast<float>(n.z - hdr.z_offset);

  out.write(reinterpret_cast<const char*>(&fx), sizeof(float));
  out.write(reinterpret_cast<const char*>(&fy), sizeof(float));
  out.write(reinterpret_cast<const char*>(&fz), sizeof(float));
}

namespace v1_0
{
  QSMedge read_edge(std::ifstream& in, const QSMheader& /*hdr*/)
  {
    uint32_t src = 0, tgt = 0;
    float   r = 0.0f;
    uint8_t q = 0;

    in.read(reinterpret_cast<char*>(&src),  sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&tgt),  sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&r),    sizeof(float));
    in.read(reinterpret_cast<char*>(&q),    sizeof(uint8_t));

    QSMedge e;
    e.source           = src;
    e.target           = tgt;
    e.radius           = r;
    return e;
  }

  void write_edge(std::ofstream& out, const QSMedge& e, const QSMheader& /*hdr*/)
  {
    const uint32_t src = e.source;
    const uint32_t tgt = e.target;

    out.write(reinterpret_cast<const char*>(&src),              sizeof(int32_t));
    out.write(reinterpret_cast<const char*>(&tgt),              sizeof(int32_t));
    out.write(reinterpret_cast<const char*>(&e.quality),        sizeof(uint8_t));
  }
} // namespace v1_1

namespace v1_1
{
  QSMedge read_edge(std::ifstream& in, const QSMheader& /*hdr*/)
  {
    uint32_t src = 0, tgt = 0, axis = 0;
    float   r   = 0.0f;
    uint8_t bo  = 0, q = 0;

    in.read(reinterpret_cast<char*>(&src),  sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&tgt),  sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&r),    sizeof(float));
    in.read(reinterpret_cast<char*>(&q),    sizeof(uint8_t));
    in.read(reinterpret_cast<char*>(&axis), sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&bo),   sizeof(uint8_t));

    QSMedge e;
    e.source           = src;
    e.target           = tgt;
    e.radius           = r;
    e.axis_id          = axis;
    e.branch_order     = bo;
    e.quality          = q;
    return e;
  }

  void write_edge(std::ofstream& out, const QSMedge& e, const QSMheader& /*hdr*/)
  {
    const uint32_t src = e.source;
    const uint32_t tgt = e.target;

    out.write(reinterpret_cast<const char*>(&src),              sizeof(int32_t));
    out.write(reinterpret_cast<const char*>(&tgt),              sizeof(int32_t));
    out.write(reinterpret_cast<const char*>(&e.quality),        sizeof(uint8_t));
    out.write(reinterpret_cast<const char*>(&e.axis_id),        sizeof(int32_t));
    out.write(reinterpret_cast<const char*>(&e.branch_order),   sizeof(uint8_t));
  }
} // namespace v1_1


namespace v1_2
{
  QSMedge read_edge(std::ifstream& in, const QSMheader& /*hdr*/)
  {
    uint32_t src = 0, tgt = 0, axis = 0;
    float   r   = 0.0f, sl = 0.0f, dtr = 0.0f;
    uint8_t bo  = 0, q = 0;

    in.read(reinterpret_cast<char*>(&src),  sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&tgt),  sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&r),    sizeof(float));
    in.read(reinterpret_cast<char*>(&sl),   sizeof(float));
    in.read(reinterpret_cast<char*>(&dtr),  sizeof(float));
    in.read(reinterpret_cast<char*>(&axis), sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&bo),   sizeof(uint8_t));
    in.read(reinterpret_cast<char*>(&q),    sizeof(uint8_t));

    QSMedge e;
    e.source           = src;
    e.target           = tgt;
    e.radius           = r;
    e.subtree_length   = sl;
    e.distance_to_root = dtr;
    e.axis_id          = axis;
    e.branch_order     = bo;
    e.quality          = q;
    // Fields added in future versions are left at their sentinel defaults.
    return e;
  }

  void write_edge(std::ofstream& out, const QSMedge& e, const QSMheader& /*hdr*/)
  {
    const int32_t src = e.source;
    const int32_t tgt = e.target;

    out.write(reinterpret_cast<const char*>(&src),              sizeof(int32_t));
    out.write(reinterpret_cast<const char*>(&tgt),              sizeof(int32_t));
    out.write(reinterpret_cast<const char*>(&e.radius),         sizeof(float));
    out.write(reinterpret_cast<const char*>(&e.subtree_length), sizeof(float));
    out.write(reinterpret_cast<const char*>(&e.distance_to_root), sizeof(float));
    out.write(reinterpret_cast<const char*>(&e.axis_id),        sizeof(int32_t));
    out.write(reinterpret_cast<const char*>(&e.branch_order),   sizeof(uint8_t));
    out.write(reinterpret_cast<const char*>(&e.quality),        sizeof(uint8_t));
  }
} // namespace v1_2

const std::vector<QSMversionSpec>& version_registry()
{
  static const std::vector<QSMversionSpec> reg = {
    {1, 0, 12, 13, { read_node,  v1_0::read_edge, write_node, v1_0::write_edge }},
    {1, 1, 12, 18, { read_node,  v1_1::read_edge, write_node, v1_1::write_edge }},
    {1, 2, 12, 26, { read_node,  v1_2::read_edge, write_node, v1_2::write_edge }}
  };
  return reg;
}

const QSMversionSpec& resolve_version(uint8_t major, uint8_t minor,  std::string* warning_out)
{
  const auto& reg = version_registry();

  for (const auto& spec : reg)
  {
    if (spec.major == major && spec.minor == minor)
      return spec;
  }

  // Known major, unknown minor -> use highest known minor, emit warning
  const QSMversionSpec* best = nullptr;
  for (const auto& spec : reg)
  {
    if (spec.major == major)
    {
      if (!best || spec.minor > best->minor)
        best = &spec;
    }
  }

  if (best)
  {
    std::string warn =
      "libqsm: file VERSION " + std::to_string(major) + "." + std::to_string(minor) +
      " is ahead of the highest known minor version " +
      std::to_string(best->major) + "." + std::to_string(best->minor) +
      ". Extra fields will be skipped. Consider upgrading libqsm.";

    if (warning_out)
      *warning_out = std::move(warn);

    return *best;
  }

  // Unknown major - cannot interpret safely
  throw std::runtime_error(
      "libqsm: file VERSION " + std::to_string(major) + "." + std::to_string(minor) +
        " uses an unknown major version. This reader supports major version(s): 1.");
}


// ===========================================================================
// Internal helpers (anonymous namespace)
// ===========================================================================

namespace {

static std::string utc_timestamp()
{
  char buf[32] = "unknown";
  const std::time_t t   = std::time(nullptr);
  std::tm* const    gmt = std::gmtime(&t);
  if (gmt) std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", gmt);
  return buf;
}

static QSMheader parse_header(std::ifstream& in)
{
  QSMheader   hdr;
  bool        data_seen = false;
  std::string line;

  while (std::getline(in, line))
  {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty() || line[0] == '#') continue;

    const auto        sp  = line.find(' ');
    const std::string key = (sp == std::string::npos) ? line          : line.substr(0, sp);
    const std::string val = (sp == std::string::npos) ? std::string{} : line.substr(sp + 1);

    if (key == "DATA")
    {
      if (val != "binary")
        throw std::runtime_error("libqsm: only 'DATA binary' is supported, got: DATA " + val);
      data_seen = true;
      break;
    }
    else if (key == "VERSION")
    {
      const std::size_t dot = val.find('.');
      if (dot == std::string::npos)
        throw std::runtime_error("libqsm: invalid VERSION format: " + val);
      try {
        hdr.version_major = static_cast<uint8_t>(std::stoi(val.substr(0, dot)));
        hdr.version_minor = static_cast<uint8_t>(std::stoi(val.substr(dot + 1)));
      } catch (const std::exception&) {
        throw std::runtime_error("libqsm: invalid VERSION number: " + val);
      }
    }
    else if (key == "NODES")     { hdr.n_nodes   = std::stoull(val); }
    else if (key == "EDGES")     { hdr.n_edges   = std::stoull(val); }
    else if (key == "NODE_SIZE") { hdr.node_size = static_cast<uint32_t>(std::stoul(val)); }
    else if (key == "EDGE_SIZE") { hdr.edge_size = static_cast<uint32_t>(std::stoul(val)); }
    else if (key == "X_OFFSET")  { hdr.x_offset  = std::stod(val); }
    else if (key == "Y_OFFSET")  { hdr.y_offset  = std::stod(val); }
    else if (key == "Z_OFFSET")  { hdr.z_offset  = std::stod(val); }
    else if (key == "SOFTWARE")  { hdr.software  = val; }
    else if (key == "CREATED")   { hdr.created   = val; }
    else if (key == "CRS")       { hdr.crs       = val; }
    else if (key == "MESSAGE")
    {
      std::string msg;
      msg.reserve(val.size());
      for (std::size_t i = 0; i < val.size(); ++i)
      {
        if (val[i] == '\\' && i + 1 < val.size() && val[i + 1] == 'n')
        { msg += '\n'; ++i; }
        else
        { msg += val[i]; }
      }
      hdr.messages.push_back(std::move(msg));
    }
    // Unknown keys are silently ignored (forward-compatibility).
  }

  if (!data_seen)
    throw std::runtime_error("libqsm: 'DATA binary' sentinel not found in header.");

  if (hdr.node_size == 0)
    throw std::runtime_error("libqsm: NODE_SIZE not specified in header.");

  if (hdr.edge_size == 0)
    throw std::runtime_error("libqsm: EDGE_SIZE not specified in header.");

  // Validate against the version-specific minimums from the registry.
  // resolve_version() will warn (not throw) on an unknown minor version.
  const QSMversionSpec& spec = resolve_version(hdr.version_major, hdr.version_minor);

  if (hdr.node_size < spec.min_node_size)
    throw std::runtime_error(
        "libqsm: NODE_SIZE " + std::to_string(hdr.node_size) +
          " is below the minimum of " + std::to_string(spec.min_node_size) +
          " bytes for version " +
          std::to_string(hdr.version_major) + "." + std::to_string(hdr.version_minor) + ".");

  if (hdr.edge_size < spec.min_edge_size)
    throw std::runtime_error(
        "libqsm: EDGE_SIZE " + std::to_string(hdr.edge_size) +
          " is below the minimum of " + std::to_string(spec.min_edge_size) +
          " bytes for version " +
          std::to_string(hdr.version_major) + "." + std::to_string(hdr.version_minor) + ".");

  return hdr;
}


static void read_binary(
      std::ifstream& in,
      const QSMheader& hdr,
      const std::function<void(const QSMnode&)>& on_node,
      const std::function<void(const QSMedge&)>& on_edge)
{
  const QSMversionSpec& spec = resolve_version(hdr.version_major, hdr.version_minor);

  const uint32_t node_extra = hdr.node_size - spec.min_node_size;
  const uint32_t edge_extra = hdr.edge_size - spec.min_edge_size;

  for (uint64_t i = 0; i < hdr.n_nodes; ++i)
  {
    const QSMnode n = spec.dispatch.read_node(in, hdr);
    if (node_extra > 0) in.ignore(node_extra);

    if (!in)
      throw std::runtime_error("libqsm: unexpected EOF reading node record " + std::to_string(i) + " of " + std::to_string(hdr.n_nodes) + ".");

    if (on_node) on_node(n);
  }

  for (uint64_t i = 0; i < hdr.n_edges; ++i)
  {
    const QSMedge e = spec.dispatch.read_edge(in, hdr);
    if (edge_extra > 0) in.ignore(edge_extra);

    if (!in)
      throw std::runtime_error("libqsm: unexpected EOF reading edge record " + std::to_string(i) + " of " + std::to_string(hdr.n_edges) + ".");

    if (on_edge) on_edge(e);
  }
}

} // anonymous namespace


// ===========================================================================
// QSMreader
// ===========================================================================

QSMreader::QSMreader(const std::string& filename) : filename_(filename)
{
  parse();
}

void QSMreader::parse()
{
  std::ifstream in(filename_, std::ios::binary);
  if (!in.is_open())
    throw std::runtime_error("libqsm: cannot open file: " + filename_);

  header      = parse_header(in);
  data_start_ = in.tellg();

  nodes_.reserve(static_cast<std::size_t>(header.n_nodes));
  edges_.reserve(static_cast<std::size_t>(header.n_edges));

  read_binary(in, header, [this](const QSMnode& n) { nodes_.push_back(n); }, [this](const QSMedge& e) { edges_.push_back(e); });

  if (nodes_.size() != header.n_nodes)
    throw std::runtime_error("libqsm: node count mismatch: header declared " + std::to_string(header.n_nodes) + ", read " + std::to_string(nodes_.size()) + ".");

  if (edges_.size() != header.n_edges)
    throw std::runtime_error("libqsm: edge count mismatch: header declared " + std::to_string(header.n_edges) + ", read " + std::to_string(edges_.size()) + ".");
}


// ===========================================================================
// QSMwriter
// ===========================================================================

QSMwriter::QSMwriter(const std::string& filename) : filename_(filename)
{}

void QSMwriter::write()
{
  // Sync record counts from the vectors - they are the source of truth.
  header.n_nodes = nodes_.size();
  header.n_edges = edges_.size();

  const QSMversionSpec& spec = resolve_version(header.version_major, header.version_minor);

  std::ofstream out(filename_, std::ios::binary);
  if (!out.is_open())
    throw std::runtime_error("libqsm: cannot open file for writing: " + filename_);

  write_header(out, spec);
  write_nodes (out, spec);
  write_edges (out, spec);

  if (!out)
    throw std::runtime_error("libqsm: I/O error while writing: " + filename_);
}

void QSMwriter::write_header(std::ofstream& out, const QSMversionSpec& spec) const
{
  const std::string created = header.created.empty() ? utc_timestamp() : header.created;

  out << "# QSM Binary Format\n";
  out << "VERSION "   << std::to_string(header.version_major) << "." << std::to_string(header.version_minor) << "\n";
  if (!header.software.empty())
    out << "SOFTWARE " << header.software << "\n";
  out << "CREATED "   << created          << "\n";
  out << "NODES "     << header.n_nodes   << "\n";
  out << "EDGES "     << header.n_edges   << "\n";
  // Use the version-registry minimum sizes, not magic constants.
  out << "NODE_SIZE " << spec.min_node_size << "\n";
  out << "EDGE_SIZE " << spec.min_edge_size << "\n";
  out << std::fixed << std::setprecision(4);
  out << "X_OFFSET "  << header.x_offset  << "\n";
  out << "Y_OFFSET "  << header.y_offset  << "\n";
  out << "Z_OFFSET "  << header.z_offset  << "\n";
  if (!header.crs.empty())
    out << "CRS "     << header.crs       << "\n";
  for (const std::string& msg : header.messages)
  {
    out << "MESSAGE ";
    for (const char c : msg)
    {
      if (c == '\n') out << "\\n";
      else           out << c;
    }
    out << "\n";
  }
  out << "DATA binary\n";
}

void QSMwriter::write_nodes(std::ofstream& out, const QSMversionSpec& spec) const
{
  for (const QSMnode& n : nodes_)
    spec.dispatch.write_node(out, n, header);
}

void QSMwriter::write_edges(std::ofstream& out, const QSMversionSpec& spec) const
{
  for (const QSMedge& e : edges_)
    spec.dispatch.write_edge(out, e, header);
}

} // namespace libqsm
