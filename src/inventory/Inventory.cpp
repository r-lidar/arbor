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
#include <limits>

namespace arbor::inventory {

static const char* RESERVED_NAMES[] = {"id", "x", "y", "z", "radius", "pom", "species"};

const char* to_string(PropertyType type) noexcept
{
  switch (type)
  {
  case PropertyType::Double: return "double";
  case PropertyType::Int:    return "int";
  case PropertyType::String: return "string";
  }
  return "unknown";
}

static PropertyType type_of(const PropertyValue& value) noexcept
{
  return static_cast<PropertyType>(value.index());
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

// ---- Size / access ----

void Inventory::clear()
{
  trees_.clear();
  index_.clear();
}

const InventoryTree& Inventory::at_index(size_t i) const
{
  if (i >= trees_.size())
    throw std::out_of_range("Inventory: index " + std::to_string(i) + " out of range (size " + std::to_string(trees_.size()) + ")");
  return trees_[i];
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

size_t Inventory::index_of(int id) const
{
  auto it = index_.find(id);
  if (it == index_.end())
    throw std::out_of_range("Inventory: no tree with id " + std::to_string(id));
  return it->second;
}

InventoryTree& Inventory::mutable_tree(int id)
{
  return trees_[index_of(id)];
}

void Inventory::rebuild_index()
{
  index_.clear();
  index_.reserve(trees_.size());
  for (size_t i = 0; i < trees_.size(); ++i)
    index_[trees_[i].id] = i;
}

// ---- Tree CRUD ----

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

void Inventory::check_tree(const InventoryTree& tree)
{
  check_id(tree.id);
  check_position(tree.id, tree.x, tree.y, tree.z);
  check_positive(tree.id, "radius", tree.radius);
  check_positive(tree.id, "pom", tree.pom);
}

const InventoryTree& Inventory::add_tree(int id, double x, double y, double z, double radius, double pom, const std::string& species)
{
  InventoryTree tree;
  tree.id = id;
  tree.x = x;
  tree.y = y;
  tree.z = z;
  tree.radius = radius;
  tree.pom = pom;
  tree.species = species;
  return add_tree(std::move(tree));
}

const InventoryTree& Inventory::add_tree(InventoryTree tree)
{
  check_tree(tree);

  if (contains(tree.id))
    throw std::invalid_argument("Inventory: a tree with id " + std::to_string(tree.id) + " already exists");

  if (tree.properties.empty())
  {
    tree.properties.reserve(schema_.size());
    for (const PropertyDef& def : schema_)
      tree.properties.push_back(def.default_value);
  }
  else
  {
    if (tree.properties.size() != schema_.size())
      throw std::invalid_argument("Inventory: tree " + std::to_string(tree.id) + " has " + std::to_string(tree.properties.size()) + " properties, schema has " + std::to_string(schema_.size()));

    for (size_t i = 0; i < schema_.size(); ++i)
      tree.properties[i] = coerce(tree.properties[i], schema_[i].type, schema_[i].name);
  }

  index_[tree.id] = trees_.size();
  trees_.push_back(std::move(tree));
  return trees_.back();
}

bool Inventory::remove_tree(int id)
{
  auto it = index_.find(id);
  if (it == index_.end()) return false;
  trees_.erase(trees_.begin() + it->second);
  rebuild_index();
  return true;
}

int Inventory::next_id() const noexcept
{
  if (trees_.empty()) return 1;
  int max_id = trees_.front().id;
  for (const InventoryTree& t : trees_) max_id = std::max(max_id, t.id);
  return max_id + 1;
}

void Inventory::set_id(int id, int new_id)
{
  if (id == new_id) { index_of(id); return; }
  check_id(new_id);
  if (contains(new_id))
    throw std::invalid_argument("Inventory: a tree with id " + std::to_string(new_id) + " already exists");

  size_t i = index_of(id);
  trees_[i].id = new_id;
  index_.erase(id);
  index_[new_id] = i;
}

void Inventory::set_position(int id, double x, double y, double z)
{
  check_position(id, x, y, z);
  InventoryTree& t = mutable_tree(id);
  t.x = x;
  t.y = y;
  t.z = z;
}

void Inventory::set_radius(int id, double radius)
{
  check_positive(id, "radius", radius);
  mutable_tree(id).radius = radius;
}

void Inventory::set_pom(int id, double pom)
{
  check_positive(id, "pom", pom);
  mutable_tree(id).pom = pom;
}

void Inventory::set_species(int id, const std::string& species)
{
  mutable_tree(id).species = species;
}

// ---- Schema ----

bool Inventory::is_reserved_name(const std::string& name) noexcept
{
  for (const char* r : RESERVED_NAMES)
    if (name == r) return true;
  return false;
}

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

PropertyType Inventory::property_type(const std::string& name) const
{
  return schema_[property_index(name)].type;
}

void Inventory::check_new_property_name(const std::string& name) const
{
  if (name.empty())
    throw std::invalid_argument("Inventory: property name cannot be empty");
  if (is_reserved_name(name))
    throw std::invalid_argument("Inventory: '" + name + "' is a reserved name");
  if (has_property(name))
    throw std::invalid_argument("Inventory: property '" + name + "' already exists");
}

void Inventory::add_property(const std::string& name, PropertyType type)
{
  add_property(name, default_of(type));
}

void Inventory::add_property(const std::string& name, const PropertyValue& default_value)
{
  check_new_property_name(name);
  schema_.push_back({name, type_of(default_value), default_value});
  for (InventoryTree& t : trees_)
    t.properties.push_back(default_value);
}

void Inventory::remove_property(const std::string& name)
{
  size_t i = property_index(name);
  schema_.erase(schema_.begin() + i);
  for (InventoryTree& t : trees_)
    t.properties.erase(t.properties.begin() + i);
}

void Inventory::rename_property(const std::string& name, const std::string& new_name)
{
  size_t i = property_index(name);
  if (name == new_name) return;
  check_new_property_name(new_name);
  schema_[i].name = new_name;
}

// ---- Property values ----

PropertyValue Inventory::coerce(const PropertyValue& value, PropertyType type, const std::string& name) const
{
  PropertyType vtype = type_of(value);
  if (vtype == type) return value;
  if (type == PropertyType::Double && vtype == PropertyType::Int)
    return static_cast<double>(std::get<int>(value));

  throw std::invalid_argument(std::string("Inventory: property '") + name + "' is of type " + to_string(type) + ", got " + to_string(vtype));
}

void Inventory::set_property(int id, const std::string& name, const PropertyValue& value)
{
  size_t i = property_index(name);
  PropertyValue v = coerce(value, schema_[i].type, name);
  mutable_tree(id).properties[i] = std::move(v);
}

void Inventory::fill_property(const std::string& name, const PropertyValue& value)
{
  size_t i = property_index(name);
  PropertyValue v = coerce(value, schema_[i].type, name);
  for (InventoryTree& t : trees_)
    t.properties[i] = v;
}

const PropertyValue& Inventory::get_property(int id, const std::string& name) const
{
  size_t i = property_index(name);
  return at(id).properties[i];
}

// ---- Filter / sort ----

void Inventory::sort_by(const std::string& field, bool ascending)
{
  auto order = [ascending](const auto& a, const auto& b) { return ascending ? a < b : b < a; };

  if (field == "id")           sort([&](const InventoryTree& a, const InventoryTree& b) { return order(a.id, b.id); });
  else if (field == "x")       sort([&](const InventoryTree& a, const InventoryTree& b) { return order(a.x, b.x); });
  else if (field == "y")       sort([&](const InventoryTree& a, const InventoryTree& b) { return order(a.y, b.y); });
  else if (field == "z")       sort([&](const InventoryTree& a, const InventoryTree& b) { return order(a.z, b.z); });
  else if (field == "radius")  sort([&](const InventoryTree& a, const InventoryTree& b) { return order(a.radius, b.radius); });
  else if (field == "pom")     sort([&](const InventoryTree& a, const InventoryTree& b) { return order(a.pom, b.pom); });
  else if (field == "species") sort([&](const InventoryTree& a, const InventoryTree& b) { return order(a.species, b.species); });
  else
  {
    size_t i = property_index(field);
    sort([&](const InventoryTree& a, const InventoryTree& b) { return order(a.properties[i], b.properties[i]); });
  }
}

// ---- Spatial queries ----

std::optional<int> Inventory::nearest(double x, double y) const
{
  return nearest(x, y, std::numeric_limits<double>::infinity());
}

std::optional<int> Inventory::nearest(double x, double y, double max_distance) const
{
  std::optional<int> best;
  double best_d2 = max_distance * max_distance;
  for (const InventoryTree& t : trees_)
  {
    double dx = t.x - x;
    double dy = t.y - y;
    double d2 = dx*dx + dy*dy;
    if (d2 <= best_d2)
    {
      if (best && d2 == best_d2) continue; // keep the first tree on ties
      best_d2 = d2;
      best = t.id;
    }
  }
  return best;
}

std::vector<int> Inventory::within(double x, double y, double distance) const
{
  std::vector<std::pair<double, int>> hits;
  double r2 = distance * distance;
  for (const InventoryTree& t : trees_)
  {
    double dx = t.x - x;
    double dy = t.y - y;
    double d2 = dx*dx + dy*dy;
    if (d2 <= r2) hits.emplace_back(d2, t.id);
  }

  std::stable_sort(hits.begin(), hits.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

  std::vector<int> ids;
  ids.reserve(hits.size());
  for (const auto& h : hits) ids.push_back(h.second);
  return ids;
}

} // namespace arbor::inventory
