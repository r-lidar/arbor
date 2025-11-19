#include "QSM.h"

void QSM::build_from_cylinders(const std::vector<QSMcylinder>& input)
{
  cylinders_.clear();
  children_map_.clear();

  // Insert cylinders and build children links
  for (const auto& c : input)
  {
    cylinders_[c.cyl_ID] = c;

    // Ensure parent exists in the child map
    children_map_[c.parent_ID].push_back(c.cyl_ID);

    // Ensure the child also has a children entry even if empty
    if (!children_map_.count(c.cyl_ID))
      children_map_[c.cyl_ID] = {};
  }
}
