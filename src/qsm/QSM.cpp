#include "QSM.h"

#include <stdexcept>

QSMcylinder& QSM::get_cylinder_by_id(int cyl_id)
{
  auto it = cylinders_.find(cyl_id);
  if (it == cylinders_.end())
  {
    throw std::runtime_error("Internal error: No cylinder found with parent_ID.");
  }
  return it->second;
}

// Returns the ordered main axis as pointers (root -> tip)
std::vector<const QSMcylinder*> QSM::main_axis() const
{
  std::vector<const QSMcylinder*> axis;
  axis.reserve(cylinders_.size());

  for (const auto &kv : cylinders_)
  {
    const QSMcylinder* c = &kv.second;
    if (c->axis_ID == 1)
    {
      axis.push_back(c);
    }
  }

  std::sort(axis.begin(), axis.end(),  [](const QSMcylinder* a, const QSMcylinder* b)
  {
    return a->cyl_ID < b->cyl_ID;
  });

  return axis;
}

void QSM::add_cylinder(const QSMcylinder& c)
{
  cylinders_[c.cyl_ID] = c;

  // Ensure parent exists in the child map
  children_map_[c.parent_ID].push_back(c.cyl_ID);

  // Ensure the child also has a children entry even if empty
  if (!children_map_.count(c.cyl_ID))
    children_map_[c.cyl_ID] = {};
}

