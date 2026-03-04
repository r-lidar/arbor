#include "QSMbuilder.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

void QSMbuilder::prolongate(double d, double L)
{
  if (d <= 0.0) return;

  // Get ordered main axis (trunk)
  std::vector<const QSMcylinder*> axis = qsm.main_axis();
  const size_t n = axis.size();

  if (n < 2) return;

  // 2. Compute cyl lengths and cumulative lengths
  std::vector<double> lens(n);
  std::vector<double> cum(n);

  double total = 0.0;
  for (size_t i = 0; i < n; ++i)
  {
    lens[i] = axis[i]->length();
    total += lens[i];
    cum[i] = total;
  }
  if (total == 0.0) return;

  // 3. take only few first cylinders
  // first 10% of axis or 30 cm
  double cutoff = 0.1 * total;
  size_t k = 0;
  while (k < n && cum[k] <= cutoff) k++;
  if (k == 0) k = 1; // minimal 1 cyl
  if (k < n)
  {
    // Ensure sum < 0.3 : add next cylinder
    double s = 0.0;
    for (size_t i = 0; i < k; i++) s += lens[i];
    if (s < 0.3) k++;
  }
  if (k > n) k = n;

  const QSMcylinder* root = axis[0];
  const QSMcylinder* last  = axis[k-1];

  if (root->parent_ID != 0)
    throw std::runtime_error("Invalid QSM, the selected root parent ID is not 0");

  // Orientation estimated
  double dx = last->endX - root->startX;
  double dy = last->endY - root->startY;
  double dz = last->endZ - root->startZ;
  double N  = std::sqrt(dx*dx + dy*dy + dz*dz);
  if (N <= 0.0) return;
  double ox = dx / N;
  double oy = dy / N;
  double oz = dz / N;

  // Adjust prolongation distance by angle
  double d_adj;
  if (std::abs(oz) < 1e-9)
    d_adj = d;
  else
    d_adj = d / oz;

  // Start/end points of prolongation is the root
  double endX = root->startX;
  double endY = root->startY;
  double endZ = root->startZ;

  // Create new subdivided segments
  int nseg = std::max(1, int(std::ceil(d_adj / L)));
  double actual_L = d_adj / nseg;

  double root_subtree = root->subtree_length;
  if (root_subtree == SUBTREE_LENGTH_UNSET)
    throw std::runtime_error("Invalid QSM, the root has not subtree length");

  // Find safe negative cyl_ID
  int next_id = -1;
  int prev_id = 0;

  for (int i = 1; i <= nseg; i++)
  {
    double f1 = double(i - 1) / nseg;
    double f2 = double(i)     / nseg;

    double x1 = endX - ox * d_adj * f1;
    double y1 = endY - oy * d_adj * f1;
    double z1 = endZ - oz * d_adj * f1;

    double x2 = endX - ox * d_adj * f2;
    double y2 = endY - oy * d_adj * f2;
    double z2 = endZ - oz * d_adj * f2;

    QSMcylinder c;
    c.startX = x1;
    c.startY = y1;
    c.startZ = z1;
    c.endX   = x2;
    c.endY   = y2;
    c.endZ   = z2;
    c.cyl_ID     = next_id;
    c.parent_ID  = (i == 1 ? root->cyl_ID : prev_id);
    c.axis_ID    = 1;
    c.branch_order = 1;
    c.subtree_length = root_subtree + d_adj - actual_L * (nseg - i + 1);

    qsm.add_cylinder(c);

    prev_id = next_id;
    next_id--;
  }
}

void QSMbuilder::estimate_prolongation(const PointCloud& tree)
{
  prolongation_distance = 0.0;

  // Check if the attribute exists and QSM is not empty
  if (!tree.has_hag() || qsm.size() == 0) { return; }

  // Compute the minimum startZ in the QSM
  double min_start_z = std::numeric_limits<double>::max();
  for (const auto& pair : qsm)
  {
    const QSMcylinder& cyl = pair.second;
    if (cyl.startZ < min_start_z) {
      min_start_z = cyl.startZ;
    }
  }

  // Find max(hag) for points where tree.Z <= min_start_z
  double max_hag = -std::numeric_limits<double>::max();
  bool point_found = false;

  for (size_t i = 0; i < tree.size(); ++i)
  {
    if (tree.get_z(i) <= min_start_z)
    {
      double current_hag = tree.get_hag(i);
      if (current_hag > max_hag)
      {
        max_hag = current_hag;
        point_found = true;
      }
    }
  }

  if (point_found)
  {
    prolongation_distance = max_hag;
  }
}
