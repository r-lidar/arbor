#include <vector>
#include <memory>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <queue>
#include <map>

#include "SeedDetector.h"
#include "ransac.h"
#include "Grid3D.h"

namespace
{

// Helper function to validate if a detected circle is a valid tree
bool is_valid_circle(double radius, double angle_range, double pinlier, double pinside)
{
  if (radius > 2.0) return false;
  if (radius < 0.02) return false;
  if (radius < 0.05) return (angle_range > 180 && pinlier > 30);
  if (pinside > 20) return false;
  if (radius < 0.10) return (angle_range > 130 && pinlier > 60);
  return (angle_range > 140 && pinlier > 40);
}

// Structure to represent a cluster of points
struct Cluster
{
  int cluster_id;
  std::vector<size_t> indices; // indices into original point cloud
  Cluster() : cluster_id(0) {}
  Cluster(int id) : cluster_id(id) {}
};

std::unique_ptr<SeedDetectorGeometries::Circle> fit_circle_to_cluster(const Cluster& cluster, const PointCloud& point_cloud, int num_iterations = 400,  double inlier_threshold = 0.02)
{
  int id = cluster.cluster_id;

  if (cluster.indices.size() < 20) return nullptr;

  RansacCircle ransac(num_iterations, inlier_threshold);
  for (size_t idx : cluster.indices)
  {
    double x = point_cloud.get_x(idx);
    double y = point_cloud.get_y(idx);
    double z = point_cloud.get_z(idx);
    ransac.add_point(x, y, z);
  }
  ransac.find_circle();

  double radius = ransac.get_radius();
  double covered_arc_degree = ransac.get_arc_coverage();
  double percentage_inlier = ransac.get_inlier_percentage() * 100.0;
  double percentage_inside = ransac.get_inside_percentage() * 100.0;

  bool valid = is_valid_circle(radius, covered_arc_degree, percentage_inlier, percentage_inside);

  if (valid)
  {
    auto center = ransac.get_center();
    return std::make_unique<SeedDetectorGeometries::Circle>(center[0], center[1],  center[2], radius, id);
  }

  return nullptr;
}

void assign_cluster_ids(std::vector<SeedDetectorGeometries::Circle>& circles)
{
  int n = circles.size();
  if (n == 0) return;

  // 1. Build an adjacency list
  // We check if distance between centers < sum of radii
  // Optimization: Use squared distance to avoid expensive sqrt() calls
  std::vector<std::vector<int>> adj(n);
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double dx = circles[i].X - circles[j].X;
      double dy = circles[i].Y - circles[j].Y;
      double distSq = dx * dx + dy * dy;
      double radiiSum = circles[i].R + circles[j].R;

      if (distSq < (radiiSum * radiiSum)) {
        adj[i].push_back(j);
        adj[j].push_back(i);
      }
    }
  }

  // 2. BFS to find connected components
  std::vector<bool> visited(n, false);
  int current_gid = 0;

  for (int i = 0; i < n; ++i)
  {
    if (!visited[i]) {
      current_gid++;
      std::queue<int> q;
      q.push(i);
      visited[i] = true;

      while (!q.empty()) {
        int u = q.front();
        q.pop();
        circles[u].id = current_gid; // Assign group ID

        for (int neighbor : adj[u]) {
          if (!visited[neighbor]) {
            visited[neighbor] = true;
            q.push(neighbor);
          }
        }
      }
    }
  }
}

}

// Find group of circles in the wood slices
// -----------------------------------------------
// Main function to detect tree circles from wood points
std::vector<SeedDetectorGeometries::Circle> SeedDetector::detect_tree_circles(const PointCloud& wood, double resolution, int connectivity, int num_ransac_iterations, double inlier_threshold, size_t min_cluster_size)
{
  size_t n = wood.point_count();

  // Create Grid3D and compute connected components
  Grid3D grid(wood, resolution);
  std::vector<int> cluster_ids = grid.connected_components(connectivity);

  // Group points by cluster ID (skip cluster 0)
  std::map<int, Cluster> clusters_map;
  for (size_t i = 0; i < n; ++i)
  {
    int cluster_id = cluster_ids[i];
    if (cluster_id == 0) continue;
    if (clusters_map.find(cluster_id) == clusters_map.end())
    {
      clusters_map[cluster_id] = Cluster(cluster_id);
    }
    clusters_map[cluster_id].indices.push_back(i);
  }

  // Filter clusters by minimum size and convert to vector
  std::vector<Cluster> clusters;
  for (auto& pair : clusters_map)
  {
    if (pair.second.indices.size() >= min_cluster_size)
    {
      clusters.push_back(std::move(pair.second));
    }
  }

  // Process each cluster
  size_t total_clusters = clusters.size();
  std::vector<SeedDetectorGeometries::Circle> circles;
  for (size_t i = 0; i < total_clusters; ++i)
  {
    // Progress indicator
    /*if ((i + 1) % 50 == 0)
    {
      std::cout << "\r  Processed " << (i + 1) << " / " << total_clusters << std::flush;
    }*/

    // Try to fit a circle to this cluster
    auto result = fit_circle_to_cluster(clusters[i], wood, num_ransac_iterations, inlier_threshold);

    // If valid circle found, add to results
    if (result != nullptr)
    {
      circles.push_back(*result);
    }
  }

  assign_cluster_ids(circles);

  return circles;
}
