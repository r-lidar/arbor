/**
 * @file libqsf.cpp
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

#include "libqsf.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <ctime>

namespace libqsf {

const char* kind_name(QSFKind k) noexcept
{
  switch (k)
  {
    case QSFKind::VIRTUAL:  return "VIRTUAL";
    case QSFKind::EMBEDDED: return "EMBEDDED";
    default:                return "UNKNOWN";
  }
}

QSFKind parse_kind(const std::string& s)
{
  if (s == "VIRTUAL")  return QSFKind::VIRTUAL;
  if (s == "EMBEDDED") return QSFKind::EMBEDDED;
  throw std::runtime_error("libqsf: unknown KIND '" + s + "'. Expected 'VIRTUAL' or 'EMBEDDED'.");
}

namespace {

std::string utc_timestamp()
{
  char buf[32] = "unknown";
  const std::time_t t   = std::time(nullptr);
  std::tm* const    gmt = std::gmtime(&t);
  if (gmt) std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", gmt);
  return buf;
}

// Splits a whitespace-separated line into 'ntokens' tokens, with the last
// token containing the remainder of the line (so names may contain spaces).
std::vector<std::string> split_line(const std::string& line, std::size_t ntokens)
{
  std::vector<std::string> tokens;
  std::size_t pos = 0;

  while (tokens.size() + 1 < ntokens)
  {
    while (pos < line.size() && line[pos] == ' ') ++pos;
    std::size_t next = line.find(' ', pos);
    if (next == std::string::npos)
    {
      tokens.push_back(line.substr(pos));
      pos = line.size();
      break;
    }
    tokens.push_back(line.substr(pos, next - pos));
    pos = next;
  }

  while (pos < line.size() && line[pos] == ' ') ++pos;
  tokens.push_back(line.substr(pos));

  return tokens;
}

} // anonymous namespace

// ===========================================================================
// QSFwriter
// ===========================================================================

QSFwriter::QSFwriter(const std::string& filename) : filename_(filename) {}

void QSFwriter::write() const
{
  std::ofstream out(filename_);
  if (!out.is_open())
    throw std::runtime_error("libqsf: cannot open file for writing: " + filename_);

  const std::string created = header.created.empty() ? utc_timestamp() : header.created;

  out << "# QSF Manifest Format\n";
  out << "SIGNATURE QSFM\n";
  out << "VERSION "  << std::to_string(header.version_major) << "." << std::to_string(header.version_minor) << "\n";
  out << "KIND "     << kind_name(header.kind) << "\n";
  if (!header.software.empty())
    out << "SOFTWARE " << header.software << "\n";
  out << "CREATED "  << created << "\n";
  if (!header.crs.empty())
    out << "CRS "      << header.crs << "\n";
  out << "COUNT "    << entries_.size() << "\n";

  out << std::fixed << std::setprecision(4);
  for (const QSFEntry& e : entries_)
  {
    out << "ENTRY " << e.path << " " << e.id << " " << (e.has_bbox ? 1 : 0) << " "
        << e.xmin << " " << e.ymin << " " << e.zmin << " "
        << e.xmax << " " << e.ymax << " " << e.zmax << " "
        << (e.name.empty() ? "-" : e.name) << "\n";
  }
  out << "END\n";

  if (!out)
    throw std::runtime_error("libqsf: I/O error while writing: " + filename_);
}

// ===========================================================================
// QSFreader
// ===========================================================================

QSFreader::QSFreader(const std::string& filename) : filename_(filename)
{
  parse();
}

void QSFreader::parse()
{
  std::ifstream in(filename_, std::ios::binary);
  if (!in.is_open())
    throw std::runtime_error("libqsf: cannot open file for reading: " + filename_);

  bool signature_seen = false;
  bool kind_seen       = false;
  bool end_seen        = false;
  std::string line;

  while (std::getline(in, line))
  {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty() || line[0] == '#') continue;

    const auto        sp  = line.find(' ');
    const std::string key = (sp == std::string::npos) ? line          : line.substr(0, sp);
    const std::string val = (sp == std::string::npos) ? std::string{} : line.substr(sp + 1);

    if (key == "END")
    {
      end_seen = true;
      break;
    }
    else if (key == "SIGNATURE")
    {
      if (val != "QSFM")
        throw std::runtime_error("libqsf: invalid SIGNATURE '" + val + "', expected 'QSFM'.");
      signature_seen = true;
    }
    else if (key == "VERSION")
    {
      const std::size_t dot = val.find('.');
      if (dot == std::string::npos)
        throw std::runtime_error("libqsf: invalid VERSION format: " + val);

      try
      {
        header.version_major = static_cast<uint8_t>(std::stoi(val.substr(0, dot)));
        header.version_minor = static_cast<uint8_t>(std::stoi(val.substr(dot + 1)));
      }
      catch (const std::exception&)
      {
        throw std::runtime_error("libqsf: invalid VERSION number: " + val);
      }

      if (header.version_major != 1)
        throw std::runtime_error("libqsf: unsupported major VERSION " + std::to_string(header.version_major) + ". This reader supports major version 1 only.");
    }
    else if (key == "KIND")
    {
      header.kind = parse_kind(val);
      kind_seen = true;
    }
    else if (key == "SOFTWARE") { header.software = val; }
    else if (key == "CREATED")  { header.created  = val; }
    else if (key == "CRS")      { header.crs      = val; }
    else if (key == "COUNT")    { /* informational only; entries are counted as parsed */ }
    else if (key == "ENTRY")
    {
      // ENTRY <path> <id> <has_bbox> <xmin> <ymin> <zmin> <xmax> <ymax> <zmax> <name...>
      // 9 fixed fields, then the remainder of the line is the (possibly
      // space-containing) name, hence a 10-token split.
      std::vector<std::string> tok = split_line(val, 10);
      if (tok.size() < 10)
        throw std::runtime_error("libqsf: malformed ENTRY line: " + line);

      QSFEntry e;
      e.path     = tok[0];
      e.id       = std::stoi(tok[1]);
      e.has_bbox = std::stoi(tok[2]) != 0;
      e.xmin     = std::stod(tok[3]);
      e.ymin     = std::stod(tok[4]);
      e.zmin     = std::stod(tok[5]);
      e.xmax     = std::stod(tok[6]);
      e.ymax     = std::stod(tok[7]);
      e.zmax     = std::stod(tok[8]);
      e.name     = (tok[9] == "-") ? std::string{} : tok[9];

      entries_.push_back(std::move(e));
    }
    else
    {
      // Unknown keys are ignored for forward compatibility.
    }
  }

  if (!signature_seen)
    throw std::runtime_error("libqsf: mandatory SIGNATURE key not found in header.");

  if (!kind_seen)
    throw std::runtime_error("libqsf: mandatory KIND key not found in header.");

  if (!end_seen)
    throw std::runtime_error("libqsf: 'END' sentinel not found.");

  // KIND dispatch: only VIRTUAL is currently readable. EMBEDDED is reserved
  // for a future self-contained variant and must fail loudly rather than be
  // silently misinterpreted as VIRTUAL.
  if (header.kind == QSFKind::VIRTUAL)
  {
    // Nothing more to do: entries already parsed as relative file references.
  }
  else if (header.kind == QSFKind::EMBEDDED)
  {
    throw std::runtime_error("libqsf: KIND EMBEDDED is not yet supported by this reader.");
  }
  else
  {
    throw std::runtime_error("libqsf: unsupported KIND value.");
  }
}

} // namespace libqsf
