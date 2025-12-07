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
