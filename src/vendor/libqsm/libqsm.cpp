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

// ===========================================================================
// ExtraField helpers
// ===========================================================================

uint32_t extra_field_type_size(ExtraFieldType t) noexcept
{
  switch (t) {
  case ExtraFieldType::I8:  return 1;
  case ExtraFieldType::U8:  return 1;
  case ExtraFieldType::I16: return 2;
  case ExtraFieldType::U16: return 2;
  case ExtraFieldType::I32: return 4;
  case ExtraFieldType::U32: return 4;
  case ExtraFieldType::F32: return 4;
  case ExtraFieldType::I64: return 8;
  case ExtraFieldType::U64: return 8;
  case ExtraFieldType::F64: return 8;
  }
  return 0; // unreachable
}

const char* extra_field_type_name(ExtraFieldType t) noexcept
{
  switch (t) {
  case ExtraFieldType::I8:  return "I8";
  case ExtraFieldType::U8:  return "U8";
  case ExtraFieldType::I16: return "I16";
  case ExtraFieldType::U16: return "U16";
  case ExtraFieldType::I32: return "I32";
  case ExtraFieldType::U32: return "U32";
  case ExtraFieldType::F32: return "F32";
  case ExtraFieldType::I64: return "I64";
  case ExtraFieldType::U64: return "U64";
  case ExtraFieldType::F64: return "F64";
  }
  return "??"; // unreachable
}

ExtraField parse_extra_field(const std::string& value)
{
  const auto sp = value.find(' ');
  if (sp == std::string::npos || sp == 0 || sp + 1 >= value.size())
    throw std::runtime_error("libqsm: malformed EXTRABYTES value: '" + value + "'. Expected '<TYPE> <name>'.");

  const std::string token = value.substr(0, sp);
  const std::string name  = value.substr(sp + 1);

  if (name.empty() || name.find(' ') != std::string::npos)
    throw std::runtime_error("libqsm: EXTRABYTES field name must be a single non-empty ASCII token, got: '" + name + "'.");

  ExtraFieldType type;
  if      (token == "I8" ) type = ExtraFieldType::I8;
  else if (token == "U8" ) type = ExtraFieldType::U8;
  else if (token == "I16") type = ExtraFieldType::I16;
  else if (token == "U16") type = ExtraFieldType::U16;
  else if (token == "I32") type = ExtraFieldType::I32;
  else if (token == "U32") type = ExtraFieldType::U32;
  else if (token == "F32") type = ExtraFieldType::F32;
  else if (token == "I64") type = ExtraFieldType::I64;
  else if (token == "U64") type = ExtraFieldType::U64;
  else if (token == "F64") type = ExtraFieldType::F64;
  else
    throw std::runtime_error("libqsm: unrecognised EXTRABYTES type token '" + token + "'. Cannot determine field size; aborting to avoid record misalignment.");

  return ExtraField{ type, name };
}

uint32_t QSMheader::computed_node_extra_bytes() const noexcept
{
  uint32_t total = 0;
  for (const auto& f : node_extra_fields) total += extra_field_type_size(f.type);
  return total;
}

uint32_t QSMheader::computed_edge_extra_bytes() const noexcept
{
  uint32_t total = 0;
  for (const auto& f : edge_extra_fields) total += extra_field_type_size(f.type);
  return total;
}

QSMnode read_node(std::ifstream& in, const QSMheader& hdr)
{
  float fx = 0.0f, fy = 0.0f, fz = 0.0f;

  in.read(reinterpret_cast<char*>(&fx), sizeof(float));
  in.read(reinterpret_cast<char*>(&fy), sizeof(float));
  in.read(reinterpret_cast<char*>(&fz), sizeof(float));

  QSMnode n;
  n.x = static_cast<double>(fx) + hdr.x_offset;
  n.y = static_cast<double>(fy) + hdr.y_offset;
  n.z = static_cast<double>(fz) + hdr.z_offset;
  return n;
}

void write_node(std::ofstream& out, const QSMnode& n, const QSMheader& hdr)
{
  const float fx = static_cast<float>(n.x - hdr.x_offset);
  const float fy = static_cast<float>(n.y - hdr.y_offset);
  const float fz = static_cast<float>(n.z - hdr.z_offset);

  out.write(reinterpret_cast<const char*>(&fx), sizeof(float));
  out.write(reinterpret_cast<const char*>(&fy), sizeof(float));
  out.write(reinterpret_cast<const char*>(&fz), sizeof(float));
}

// ===========================================================================
// Edge readers / writers, one namespace per format level
//
// Binary layout per spec (all formats are cumulative):
//
//  Format 0  (13 bytes):
//    offset  0  uint32_t  source_id
//    offset  4  uint32_t  target_id
//    offset  8  float     radius
//    offset 12  uint8_t   quality_level
//
//  Format 1  (+5 bytes = 18 total):
//    offset 13  uint32_t  axis_id
//    offset 17  uint8_t   branch_order
//
//  Format 2  (+8 bytes = 26 total):
//    offset 18  float     subtree_length
//    offset 22  float     distance_to_root
// ===========================================================================

namespace fmt0
{
  QSMedge read_edge(std::ifstream& in, const QSMheader& /*hdr*/)
  {
    uint32_t src = 0, tgt = 0;
    float    r   = 0.0f;
    uint8_t  q   = 0;

    in.read(reinterpret_cast<char*>(&src), sizeof(uint32_t));
    in.read(reinterpret_cast<char*>(&tgt), sizeof(uint32_t));
    in.read(reinterpret_cast<char*>(&r),   sizeof(float));
    in.read(reinterpret_cast<char*>(&q),   sizeof(uint8_t));

    QSMedge e;
    e.source  = src;
    e.target  = tgt;
    e.radius  = r;
    e.quality = q;
    // format 1+ fields left at sentinel defaults
    return e;
  }

  void write_edge(std::ofstream& out, const QSMedge& e, const QSMheader& /*hdr*/)
  {
    out.write(reinterpret_cast<const char*>(&e.source),  sizeof(uint32_t));
    out.write(reinterpret_cast<const char*>(&e.target),  sizeof(uint32_t));
    out.write(reinterpret_cast<const char*>(&e.radius),  sizeof(float));
    out.write(reinterpret_cast<const char*>(&e.quality), sizeof(uint8_t));
  }
} // namespace fmt0

namespace fmt1
{
  QSMedge read_edge(std::ifstream& in, const QSMheader& /*hdr*/)
  {
    uint32_t src = 0, tgt = 0, axis = 0;
    float    r   = 0.0f;
    uint8_t  q   = 0, bo = 0;

    in.read(reinterpret_cast<char*>(&src),  sizeof(uint32_t));
    in.read(reinterpret_cast<char*>(&tgt),  sizeof(uint32_t));
    in.read(reinterpret_cast<char*>(&r),    sizeof(float));
    in.read(reinterpret_cast<char*>(&q),    sizeof(uint8_t));
    in.read(reinterpret_cast<char*>(&axis), sizeof(uint32_t));
    in.read(reinterpret_cast<char*>(&bo),   sizeof(uint8_t));

    QSMedge e;
    e.source       = src;
    e.target       = tgt;
    e.radius       = r;
    e.quality      = q;
    e.axis_id      = axis;
    e.branch_order = bo;
    // format 2+ fields left at sentinel defaults
    return e;
  }

  void write_edge(std::ofstream& out, const QSMedge& e, const QSMheader& /*hdr*/)
  {
    out.write(reinterpret_cast<const char*>(&e.source),       sizeof(uint32_t));
    out.write(reinterpret_cast<const char*>(&e.target),       sizeof(uint32_t));
    out.write(reinterpret_cast<const char*>(&e.radius),       sizeof(float));
    out.write(reinterpret_cast<const char*>(&e.quality),      sizeof(uint8_t));
    out.write(reinterpret_cast<const char*>(&e.axis_id),      sizeof(uint32_t));
    out.write(reinterpret_cast<const char*>(&e.branch_order), sizeof(uint8_t));
  }
} // namespace fmt1

namespace fmt2
{
  QSMedge read_edge(std::ifstream& in, const QSMheader& /*hdr*/)
  {
    uint32_t src = 0, tgt = 0, axis = 0;
    float    r   = 0.0f, sl = 0.0f, dtr = 0.0f;
    uint8_t  q   = 0, bo = 0;

    in.read(reinterpret_cast<char*>(&src),  sizeof(uint32_t));
    in.read(reinterpret_cast<char*>(&tgt),  sizeof(uint32_t));
    in.read(reinterpret_cast<char*>(&r),    sizeof(float));
    in.read(reinterpret_cast<char*>(&q),    sizeof(uint8_t));
    in.read(reinterpret_cast<char*>(&axis), sizeof(uint32_t));
    in.read(reinterpret_cast<char*>(&bo),   sizeof(uint8_t));
    in.read(reinterpret_cast<char*>(&sl),   sizeof(float));
    in.read(reinterpret_cast<char*>(&dtr),  sizeof(float));

    QSMedge e;
    e.source           = src;
    e.target           = tgt;
    e.radius           = r;
    e.quality          = q;
    e.axis_id          = axis;
    e.branch_order     = bo;
    e.subtree_length   = sl;
    e.distance_to_root = dtr;
    // Fields added in future formats are left at their sentinel defaults.
    return e;
  }

  void write_edge(std::ofstream& out, const QSMedge& e, const QSMheader& /*hdr*/)
  {
    out.write(reinterpret_cast<const char*>(&e.source),           sizeof(uint32_t));
    out.write(reinterpret_cast<const char*>(&e.target),           sizeof(uint32_t));
    out.write(reinterpret_cast<const char*>(&e.radius),           sizeof(float));
    out.write(reinterpret_cast<const char*>(&e.quality),          sizeof(uint8_t));
    out.write(reinterpret_cast<const char*>(&e.axis_id),          sizeof(uint32_t));
    out.write(reinterpret_cast<const char*>(&e.branch_order),     sizeof(uint8_t));
    out.write(reinterpret_cast<const char*>(&e.subtree_length),   sizeof(float));
    out.write(reinterpret_cast<const char*>(&e.distance_to_root), sizeof(float));
  }
} // namespace fmt2

// ===========================================================================
// Format registry
// ===========================================================================

const std::vector<QSMformatSpec>& format_registry()
{
  static const std::vector<QSMformatSpec> reg = {
    { 0, 12, 13, { read_node, fmt0::read_edge, write_node, fmt0::write_edge } },
    { 1, 12, 18, { read_node, fmt1::read_edge, write_node, fmt1::write_edge } },
    { 2, 12, 26, { read_node, fmt2::read_edge, write_node, fmt2::write_edge } },
  };
  return reg;
}

const QSMformatSpec& resolve_format(uint8_t format, std::string* warning_out)
{
  const auto& reg = format_registry();

  for (const auto& spec : reg)
  {
    if (spec.format == format)
      return spec;
  }

  // Unknown format  fall back to the highest known one and emit a warning.
  const QSMformatSpec& best = reg.back();

  if (warning_out)
  {
    *warning_out =
      "libqsm: FORMAT " + std::to_string(format) +
      " is ahead of the highest known format " + std::to_string(best.format) +
      ". Extra fields will be skipped. Consider upgrading libqsm.";
  }

  return best;
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
    bool        data_seen      = false;
    bool        signature_seen = false;
    bool        format_seen    = false;
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
      else if (key == "SIGNATURE")
      {
        if (val != "QSMF")
          throw std::runtime_error("libqsm: invalid SIGNATURE '" + val + "', expected 'QSMF'.");

        signature_seen = true;
      }
      else if (key == "VERSION")
      {
        const std::size_t dot = val.find('.');
        if (dot == std::string::npos)
          throw std::runtime_error("libqsm: invalid VERSION format: " + val);

        try
        {
          hdr.version_major = static_cast<uint8_t>(std::stoi(val.substr(0, dot)));
          hdr.version_minor = static_cast<uint8_t>(std::stoi(val.substr(dot + 1)));
        }
        catch (const std::exception&)
        {
          throw std::runtime_error("libqsm: invalid VERSION number: " + val);
        }

        if (hdr.version_major != 1)
          throw std::runtime_error("libqsm: unsupported major VERSION " + std::to_string(hdr.version_major) +". This reader supports major version 1 only.");
      }
      else if (key == "FORMAT")
      {
        try {
          hdr.format = static_cast<uint8_t>(std::stoul(val));
        } catch (const std::exception&) {
          throw std::runtime_error("libqsm: invalid FORMAT value: " + val);
        }
        format_seen = true;
      }
      else if (key == "NODES")           { hdr.n_nodes  = std::stoull(val); }
      else if (key == "EDGES")           { hdr.n_edges  = std::stoull(val); }
      else if (key == "X_OFFSET")        { hdr.x_offset = std::stod(val); }
      else if (key == "Y_OFFSET")        { hdr.y_offset = std::stod(val); }
      else if (key == "Z_OFFSET")        { hdr.z_offset = std::stod(val); }
      else if (key == "XMIN")            { hdr.xmin     = std::stod(val); }
      else if (key == "YMIN")            { hdr.ymin     = std::stod(val); }
      else if (key == "ZMIN")            { hdr.zmin     = std::stod(val); }
      else if (key == "XMAX")            { hdr.xmax     = std::stod(val); }
      else if (key == "YMAX")            { hdr.ymax     = std::stod(val); }
      else if (key == "ZMAX")            { hdr.zmax     = std::stod(val); }
      else if (key == "TREEID")          { hdr.treeid   = std::stoul(val); }
      else if (key == "TREENAME")        { hdr.treename = val; }
      else if (key == "SOFTWARE")        { hdr.software = val; }
      else if (key == "CREATED")         { hdr.created  = val; }
      else if (key == "CRS")             { hdr.crs      = val; }
      else if (key == "NODE_EXTRABYTES") { hdr.node_extra_fields.push_back(parse_extra_field(val)); }
      else if (key == "EDGE_EXTRABYTES") { hdr.edge_extra_fields.push_back(parse_extra_field(val)); }
      // NODE_SIZE / EDGE_SIZE from older files are silently ignored; sizes are now derived.
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
      else
      {
        // Unknown keys are stored in a map.
        hdr.extra_keys[key] = val;
      }
    }

    if (!data_seen)
      throw std::runtime_error("libqsm: 'DATA binary' sentinel not found in header.");

    if (!signature_seen)
      throw std::runtime_error("libqsm: mandatory SIGNATURE key not found in header.");

    if (!format_seen)
      throw std::runtime_error("libqsm: mandatory FORMAT key not found in header.");

    // Resolve format and validate that the declared extra bytes do not produce a
    // record smaller than the format's own mandatory minimum.
    // resolve_format() will warn (not throw) for an unrecognised format value.
    const QSMformatSpec& spec = resolve_format(hdr.format);

    const uint32_t eff_node = spec.min_node_size + hdr.computed_node_extra_bytes();
    const uint32_t eff_edge = spec.min_edge_size + hdr.computed_edge_extra_bytes();

    // Hard lower bounds from the spec (these can only fail if the format registry
    // itself is wrong, but we check defensively).
    if (eff_node < 12)
      throw std::runtime_error("libqsm: effective node record size " + std::to_string(eff_node) + " is below the absolute minimum of 12 bytes. File is non-conforming.");

    if (eff_edge < 13)
      throw std::runtime_error("libqsm: effective edge record size " + std::to_string(eff_edge) + " is below the absolute minimum of 13 bytes. File is non-conforming.");

    return hdr;
  }


  static void read_binary(
      std::ifstream& in,
      const QSMheader& hdr,
      const std::function<void(const QSMnode&)>& on_node,
      const std::function<void(const QSMedge&)>& on_edge)
  {
    const QSMformatSpec& spec = resolve_format(hdr.format);

    const uint32_t node_extra = hdr.computed_node_extra_bytes();
    const uint32_t edge_extra = hdr.computed_edge_extra_bytes();

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

  read_binary(in, header,
              [this](const QSMnode& n) { nodes_.push_back(n); },
              [this](const QSMedge& e) { edges_.push_back(e); });

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

  // Compute bounding box in global coordinates (mandatory; must reflect actual data).
  if (!nodes_.empty())
  {
    header.xmin = header.xmax = nodes_[0].x;
    header.ymin = header.ymax = nodes_[0].y;
    header.zmin = header.zmax = nodes_[0].z;
    for (const auto& n : nodes_)
    {
      header.xmin = std::min(header.xmin, n.x);
      header.ymin = std::min(header.ymin, n.y);
      header.zmin = std::min(header.zmin, n.z);
      header.xmax = std::max(header.xmax, n.x);
      header.ymax = std::max(header.ymax, n.y);
      header.zmax = std::max(header.zmax, n.z);
    }
  }

  const QSMformatSpec& spec = resolve_format(header.format);

  std::ofstream out(filename_, std::ios::binary);
  if (!out.is_open())
    throw std::runtime_error("libqsm: cannot open file for writing: " + filename_);

  write_header(out, spec);
  write_nodes (out, spec);
  write_edges (out, spec);

  if (!out)
    throw std::runtime_error("libqsm: I/O error while writing: " + filename_);
}

void QSMwriter::write_header(std::ofstream& out, const QSMformatSpec& /*spec*/) const
{
  const std::string created = header.created.empty() ? utc_timestamp() : header.created;

  out << "# QSM Binary Format\n";
  out << "SIGNATURE QSMF\n";
  out << "VERSION "   << std::to_string(header.version_major) << "." << std::to_string(header.version_minor) << "\n";
  out << "FORMAT "    << std::to_string(header.format)        << "\n";
  if (!header.software.empty())
    out << "SOFTWARE " << header.software << "\n";
  out << "CREATED "   << created          << "\n";
  out << "NODES "     << header.n_nodes   << "\n";
  out << "EDGES "     << header.n_edges   << "\n";
  out << std::fixed << std::setprecision(4);
  out << "X_OFFSET "  << header.x_offset  << "\n";
  out << "Y_OFFSET "  << header.y_offset  << "\n";
  out << "Z_OFFSET "  << header.z_offset  << "\n";
  out << "XMIN "      << header.xmin      << "\n";
  out << "YMIN "      << header.ymin      << "\n";
  out << "ZMIN "      << header.zmin      << "\n";
  out << "XMAX "      << header.xmax      << "\n";
  out << "YMAX "      << header.ymax      << "\n";
  out << "ZMAX "      << header.zmax      << "\n";
  if (!header.crs.empty())
    out << "CRS "     << header.crs       << "\n";
  for (const ExtraField& f : header.node_extra_fields)
    out << "NODE_EXTRABYTES " << extra_field_type_name(f.type) << " " << f.name << "\n";
  for (const ExtraField& f : header.edge_extra_fields)
    out << "EDGE_EXTRABYTES " << extra_field_type_name(f.type) << " " << f.name << "\n";
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
  for (const auto& kv : header.extra_keys)
    out << kv.first << " " << kv.second << "\n";
  out << "DATA binary\n";
}

void QSMwriter::write_nodes(std::ofstream& out, const QSMformatSpec& spec) const
{
  for (const QSMnode& n : nodes_)
    spec.dispatch.write_node(out, n, header);
}

void QSMwriter::write_edges(std::ofstream& out, const QSMformatSpec& spec) const
{
  for (const QSMedge& e : edges_)
    spec.dispatch.write_edge(out, e, header);
}

} // namespace libqsm
