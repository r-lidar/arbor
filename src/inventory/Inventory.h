/**
 * @file Inventory.h
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

#ifndef INVENTORY_H
#define INVENTORY_H

#include <cstddef>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

namespace arbor::inventory {

// Extra fields of an inventory (columns other than the built-in ones), e.g. read from a file.
enum class PropertyType { Double, Int, String };
using PropertyValue = std::variant<double, int, std::string>;
struct PropertyDef { std::string name; PropertyType type; };

// A single tree. 'properties' is aligned with the schema of the owning Inventory:
// properties[i] holds the value of Inventory::schema()[i] and has its declared type.
struct InventoryTree
{
  int id = 0;
  double x = 0;
  double y = 0;
  double z = 0;
  double radius = 0;
  double pom = 0; // height of the point of measurement above ground
  std::string species;
  std::vector<PropertyValue> properties;
};

// A collection of uniquely identified trees sharing a typed property schema.
// Trees are exposed read-only; every edit goes through the Inventory so that
// ID uniqueness and property types are always consistent.
class Inventory
{
public:
  using const_iterator = std::vector<InventoryTree>::const_iterator;

  // ---- Access ----
  size_t size() const noexcept { return trees_.size(); }
  bool empty() const noexcept { return trees_.empty(); }
  void reserve(size_t n) { trees_.reserve(n); }

  const_iterator begin() const noexcept { return trees_.begin(); }
  const_iterator end() const noexcept { return trees_.end(); }
  const std::vector<InventoryTree>& trees() const noexcept { return trees_; }

  const InventoryTree& operator[](size_t i) const { return trees_[i]; }
  const InventoryTree& at(int id) const;            // throws if absent
  const InventoryTree* find(int id) const noexcept; // nullptr if absent
  bool contains(int id) const noexcept { return index_.count(id) > 0; }

  // ---- Trees ----
  // id must be > 0, coordinates finite, radius and pom finite and > 0. Throws otherwise.
  // Properties are initialized to the default value of their type.
  const InventoryTree& add_tree(int id, double x, double y, double z, double radius, double pom, const std::string& species = "");
  bool remove_tree(int id);

  void set_position(int id, double x, double y, double z);
  void set_radius(int id, double radius);
  void set_pom(int id, double pom);
  void set_species(int id, const std::string& species);

  // ---- Properties ----
  const std::vector<PropertyDef>& schema() const noexcept { return schema_; }
  bool has_property(const std::string& name) const noexcept;
  void add_property(const std::string& name, PropertyType type);

  // int values are promoted to double for Double properties; any other type mismatch throws.
  void set_property(int id, const std::string& name, const PropertyValue& value);
  const PropertyValue& get_property(int id, const std::string& name) const;

private:
  std::vector<InventoryTree> trees_;
  std::vector<PropertyDef> schema_;
  std::unordered_map<int, size_t> index_; // id -> position in trees_

  size_t index_of(int id) const; // throws if absent
  size_t property_index(const std::string& name) const; // throws if absent
};

} // namespace arbor::inventory

#endif // INVENTORY_H
