#include "QSM.h"

void QSM::add_cylinder(const QSMcylinder& c)
{
  cylinders_[c.cyl_ID] = c;

  // Ensure parent exists in the child map
  children_map_[c.parent_ID].push_back(c.cyl_ID);

  // Ensure the child also has a children entry even if empty
  if (!children_map_.count(c.cyl_ID))
    children_map_[c.cyl_ID] = {};
}
