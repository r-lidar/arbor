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

#include <algorithm>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

namespace arbor::inventory {

enum class PropertyType { Double, Int, String };

// Order of alternatives matches PropertyType so that variant::index() == PropertyType
using PropertyValue = std::variant<double, int, std::string>;

struct PropertyDef
{
  std::string name;
  PropertyType type;
  PropertyValue default_value;
};

// A single tree. 'properties' is aligned with the schema of the owning Inventory:
// properties[i] holds the value of Inventory::schema()[i] and has its declared type.
struct InventoryTree
{
  int id = 0;
  double x = 0;
  double y = 0;
  double z = 0;
  double radius = 0;
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

  Inventory() = default;

  // ---- Size / access ----
  size_t size() const noexcept { return trees_.size(); }
  bool empty() const noexcept { return trees_.empty(); }
  void reserve(size_t n) { trees_.reserve(n); }
  void clear();

  const_iterator begin() const noexcept { return trees_.begin(); }
  const_iterator end() const noexcept { return trees_.end(); }
  const std::vector<InventoryTree>& trees() const noexcept { return trees_; }

  const InventoryTree& operator[](size_t i) const { return trees_[i]; }
  const InventoryTree& at_index(size_t i) const;
  const InventoryTree& at(int id) const;            // throws if absent
  const InventoryTree* find(int id) const noexcept; // nullptr if absent
  bool contains(int id) const noexcept { return index_.count(id) > 0; }
  size_t index_of(int id) const;                    // throws if absent

  // ---- Tree CRUD ----
  const InventoryTree& add_tree(int id, double x, double y, double z, double radius = 0, const std::string& species = "");
  const InventoryTree& add_tree(InventoryTree tree);
  bool remove_tree(int id);
  int next_id() const noexcept; // max(id) + 1, or 1 if empty

  void set_id(int id, int new_id);
  void set_position(int id, double x, double y, double z);
  void set_radius(int id, double radius);
  void set_species(int id, const std::string& species);

  // ---- Schema ----
  const std::vector<PropertyDef>& schema() const noexcept { return schema_; }
  bool has_property(const std::string& name) const noexcept;
  size_t property_index(const std::string& name) const; // throws if absent
  PropertyType property_type(const std::string& name) const;

  void add_property(const std::string& name, PropertyType type);
  void add_property(const std::string& name, const PropertyValue& default_value);
  void remove_property(const std::string& name);
  void rename_property(const std::string& name, const std::string& new_name);

  // ---- Property values ----
  // int values are promoted to double for Double properties; any other type mismatch throws.
  void set_property(int id, const std::string& name, const PropertyValue& value);
  // Sets the same value for all trees
  void fill_property(const std::string& name, const PropertyValue& value);
  const PropertyValue& get_property(int id, const std::string& name) const;

  template<typename T>
  const T& get(int id, const std::string& name) const
  {
    const PropertyValue& v = get_property(id, name);
    if (const T* p = std::get_if<T>(&v)) return *p;
    throw std::runtime_error("Inventory: property '" + name + "' is not of the requested type");
  }

  // ---- Filter / sort ----
  // Pred: bool(const InventoryTree&)
  template<typename Pred>
  Inventory filter(Pred pred) const
  {
    Inventory out;
    out.schema_ = schema_;
    for (const InventoryTree& t : trees_)
      if (pred(t)) out.trees_.push_back(t);
    out.rebuild_index();
    return out;
  }

  // Pred: bool(const InventoryTree&). Returns the number of trees removed.
  template<typename Pred>
  size_t remove_if(Pred pred)
  {
    size_t n = trees_.size();
    trees_.erase(std::remove_if(trees_.begin(), trees_.end(), pred), trees_.end());
    rebuild_index();
    return n - trees_.size();
  }

  // Compare: bool(const InventoryTree&, const InventoryTree&)
  template<typename Compare>
  void sort(Compare comp)
  {
    std::stable_sort(trees_.begin(), trees_.end(), comp);
    rebuild_index();
  }

  // Sorts by a built-in field (id, x, y, z, radius, species) or an extra property
  void sort_by(const std::string& field, bool ascending = true);

  // ---- Spatial queries (2D, horizontal distance, full scan) ----
  std::optional<int> nearest(double x, double y) const;
  std::optional<int> nearest(double x, double y, double max_distance) const;
  std::vector<int> within(double x, double y, double distance) const;

  static bool is_reserved_name(const std::string& name) noexcept;

private:
  std::vector<InventoryTree> trees_;
  std::vector<PropertyDef> schema_;
  std::unordered_map<int, size_t> index_; // id -> position in trees_

  InventoryTree& mutable_tree(int id);
  void rebuild_index();
  void check_new_property_name(const std::string& name) const;
  PropertyValue coerce(const PropertyValue& value, PropertyType type, const std::string& name) const;
};

const char* to_string(PropertyType type) noexcept;

} // namespace arbor::inventory

#endif // INVENTORY_H
