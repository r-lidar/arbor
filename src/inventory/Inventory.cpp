/**
 * @file Inventory.cpp
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

#include "Inventory.h"

#include <cmath>
#include <stdexcept>

namespace arbor::inventory {

static const char* RESERVED_NAMES[] = {"id", "x", "y", "z", "radius", "pom", "species"};

static const char* to_string(PropertyType type) noexcept
{
  switch (type)
  {
  case PropertyType::Double: return "double";
  case PropertyType::Int:    return "int";
  case PropertyType::String: return "string";
  }
  return "unknown";
}

static PropertyValue default_of(PropertyType type)
{
  switch (type)
  {
  case PropertyType::Double: return 0.0;
  case PropertyType::Int:    return 0;
  case PropertyType::String: return std::string();
  }
  return 0.0;
}

static void check_id(int id)
{
  if (id <= 0)
    throw std::invalid_argument("Inventory: invalid tree id " + std::to_string(id) + ". IDs must be > 0");
}

static void check_position(int id, double x, double y, double z)
{
  if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
    throw std::invalid_argument("Inventory: tree " + std::to_string(id) + " has a non finite position");
}

static void check_positive(int id, const char* what, double value)
{
  if (!std::isfinite(value) || value <= 0)
    throw std::invalid_argument("Inventory: tree " + std::to_string(id) + " has an invalid " + what + " (" + std::to_string(value) + "). It must be > 0");
}

// ---- Access ----

size_t Inventory::index_of(int id) const
{
  auto it = index_.find(id);
  if (it == index_.end())
    throw std::out_of_range("Inventory: no tree with id " + std::to_string(id));
  return it->second;
}

const InventoryTree& Inventory::at(int id) const
{
  return trees_[index_of(id)];
}

const InventoryTree* Inventory::find(int id) const noexcept
{
  auto it = index_.find(id);
  return it == index_.end() ? nullptr : &trees_[it->second];
}

// ---- Trees ----

const InventoryTree& Inventory::add_tree(int id, double x, double y, double z, double radius, double pom, const std::string& species)
{
  check_id(id);
  check_position(id, x, y, z);
  check_positive(id, "radius", radius);
  check_positive(id, "pom", pom);

  if (contains(id))
    throw std::invalid_argument("Inventory: a tree with id " + std::to_string(id) + " already exists");

  InventoryTree tree{id, x, y, z, radius, pom, species, {}};
  tree.properties.reserve(schema_.size());
  for (const PropertyDef& def : schema_)
    tree.properties.push_back(default_of(def.type));

  index_[id] = trees_.size();
  trees_.push_back(std::move(tree));
  return trees_.back();
}

bool Inventory::remove_tree(int id)
{
  auto it = index_.find(id);
  if (it == index_.end()) return false;

  size_t i = it->second;
  trees_.erase(trees_.begin() + i);
  index_.erase(it);
  for (auto& [key, pos] : index_)
    if (pos > i) --pos;
  return true;
}

void Inventory::set_position(int id, double x, double y, double z)
{
  check_position(id, x, y, z);
  InventoryTree& t = trees_[index_of(id)];
  t.x = x;
  t.y = y;
  t.z = z;
}

void Inventory::set_radius(int id, double radius)
{
  check_positive(id, "radius", radius);
  trees_[index_of(id)].radius = radius;
}

void Inventory::set_pom(int id, double pom)
{
  check_positive(id, "pom", pom);
  trees_[index_of(id)].pom = pom;
}

void Inventory::set_species(int id, const std::string& species)
{
  trees_[index_of(id)].species = species;
}

// ---- Properties ----

bool Inventory::has_property(const std::string& name) const noexcept
{
  for (const PropertyDef& def : schema_)
    if (def.name == name) return true;
  return false;
}

size_t Inventory::property_index(const std::string& name) const
{
  for (size_t i = 0; i < schema_.size(); ++i)
    if (schema_[i].name == name) return i;
  throw std::out_of_range("Inventory: no property named '" + name + "'");
}

void Inventory::add_property(const std::string& name, PropertyType type)
{
  if (name.empty())
    throw std::invalid_argument("Inventory: property name cannot be empty");
  for (const char* r : RESERVED_NAMES)
    if (name == r) throw std::invalid_argument("Inventory: '" + name + "' is a reserved name");
  if (has_property(name))
    throw std::invalid_argument("Inventory: property '" + name + "' already exists");

  schema_.push_back({name, type});
  for (InventoryTree& t : trees_)
    t.properties.push_back(default_of(type));
}

void Inventory::set_property(int id, const std::string& name, const PropertyValue& value)
{
  size_t i = property_index(name);
  PropertyType type = schema_[i].type;
  PropertyType vtype = static_cast<PropertyType>(value.index());
  PropertyValue& dst = trees_[index_of(id)].properties[i];

  if (vtype == type)
    dst = value;
  else if (type == PropertyType::Double && vtype == PropertyType::Int)
    dst = static_cast<double>(std::get<int>(value));
  else
    throw std::invalid_argument(std::string("Inventory: property '") + name + "' is of type " + to_string(type) + ", got " + to_string(vtype));
}

const PropertyValue& Inventory::get_property(int id, const std::string& name) const
{
  return at(id).properties[property_index(name)];
}

} // namespace arbor::inventory
