/**
 * @file PointCloudAttr.h
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

#ifndef POINTCLOUDATTR_H
#define POINTCLOUDATTR_H

#include <cstdint>
#include <string>

// Optional point attributes, as a bitmask. XYZ is always present and has no flag.
// Used to choose what a point cloud allocates (constructor, allocate(), release()), what subset() copies and what
// merge() keeps. Both PointCloudDefault and PointCloudDataFrame use the same flags.
enum class Attr : uint32_t
{
  None            = 0,
  RGB             = 1u << 0,
  Intensity       = 1u << 1,
  GPStime         = 1u << 2,
  ReturnNumber    = 1u << 3,
  NumberOfReturns = 1u << 4,
  Classification  = 1u << 5,
  UserData        = 1u << 6,
  HAG             = 1u << 7,
  PWood           = 1u << 8,
  Foliage         = 1u << 9,
  Passage         = 1u << 10,
  TreeID          = 1u << 11,

  // Presets
  LAS   = RGB | Intensity | GPStime | ReturnNumber | NumberOfReturns | Classification | UserData,
  Arbor = Classification | UserData | HAG | PWood | Foliage | Passage | TreeID,
  All   = LAS | Arbor
};

constexpr Attr operator|(Attr a, Attr b) { return static_cast<Attr>(static_cast<uint32_t>(a) | static_cast<uint32_t>(b)); }
constexpr Attr operator&(Attr a, Attr b) { return static_cast<Attr>(static_cast<uint32_t>(a) & static_cast<uint32_t>(b)); }
constexpr Attr operator^(Attr a, Attr b) { return static_cast<Attr>(static_cast<uint32_t>(a) ^ static_cast<uint32_t>(b)); }
constexpr Attr operator~(Attr a) { return static_cast<Attr>(~static_cast<uint32_t>(a)); }
constexpr Attr& operator|=(Attr& a, Attr b) { return a = a | b; }
constexpr Attr& operator&=(Attr& a, Attr b) { return a = a & b; }
constexpr bool any(Attr a) { return a != Attr::None; }

// Human-readable list of the flags of a, e.g. "treeID, hag". Used in error messages.
inline std::string describe(Attr a)
{
  static constexpr struct { Attr flag; const char* name; } names[] = {
    {Attr::RGB, "RGB"}, {Attr::Intensity, "Intensity"}, {Attr::GPStime, "gpstime"},
    {Attr::ReturnNumber, "ReturnNumber"}, {Attr::NumberOfReturns, "NumberOfReturns"},
    {Attr::Classification, "Classification"}, {Attr::UserData, "UserData"}, {Attr::HAG, "hag"},
    {Attr::PWood, "pwood"}, {Attr::Foliage, "foliage"}, {Attr::Passage, "passage"}, {Attr::TreeID, "treeID"}
  };

  std::string out;
  for (const auto& n : names)
  {
    if (!any(a & n.flag)) continue;
    if (!out.empty()) out += ", ";
    out += n.name;
  }
  return out.empty() ? "none" : out;
}

#endif
