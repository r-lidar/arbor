/**
 * @file MemoryUtils.h
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

#ifndef MEMORY_UTILS_H
#define MEMORY_UTILS_H

#include "Graph.h"
#include <string>
#include <sstream>
#include <iomanip>

namespace arbor::segment {

class MemoryUtils
{
public:
  // Calculate memory size of AdjacencyList in bytes
  static size_t calculate_adjacency_list_size(const Graph::AdjacencyList& adj_list)
  {
    size_t total = sizeof(adj_list);  // vector overhead
    
    for (const auto& neighbors : adj_list)
    {
      total += sizeof(neighbors);  // vector overhead per node
      total += neighbors.capacity() * sizeof(Graph::Node);  // actual node storage
    }
    
    return total;
  }

  // Calculate memory size of GraphCache in bytes
  static size_t calculate_graph_cache_size(const Graph::GraphCache& cache)
  {
    const auto& [distances, predecessors] = cache;
    
    size_t total = 0;
    total += sizeof(distances) + distances.capacity() * sizeof(Graph::Cost);
    
    // PredecessorMap is an unordered_map
    total += sizeof(predecessors);
    total += predecessors.size() * (sizeof(Graph::NodeId) * 2 + sizeof(void*));  // approximate overhead
    total += predecessors.bucket_count() * sizeof(void*);  // bucket array
    
    return total;
  }

  // Calculate total memory size of a Graph in bytes
  static size_t calculate_graph_size(const Graph* graph)
  {
    if (!graph) return 0;
    
    size_t total = sizeof(Graph);  // object overhead
    total += calculate_adjacency_list_size(graph->adjacency_list);
    
    return total;
  }

  // Format bytes into human-readable string (KB, MB, GB)
  static std::string format_bytes(size_t bytes)
  {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2);
    
    if (bytes < 1024)
    {
      oss << bytes << " B";
    }
    else if (bytes < 1024 * 1024)
    {
      oss << (bytes / 1024.0) << " KB";
    }
    else if (bytes < 1024 * 1024 * 1024)
    {
      oss << (bytes / (1024.0 * 1024.0)) << " MB";
    }
    else
    {
      oss << (bytes / (1024.0 * 1024.0 * 1024.0)) << " GB";
    }
    
    return oss.str();
  }

  // Print memory statistics for a Graph
  static std::string print_graph_memory(const Graph* graph)
  {
    if (!graph)
      return "Graph is null";
    
    size_t adj_size = calculate_adjacency_list_size(graph->adjacency_list);
    size_t total_size = calculate_graph_size(graph);
    
    std::ostringstream oss;
    oss << "Graph Memory Usage:\n";
    oss << "  Adjacency List: " << format_bytes(adj_size) << "\n";
    oss << "  Total Graph: " << format_bytes(total_size) << "\n";
    oss << "  Node Count: " << graph->adjacency_list.size();
    
    return oss.str();
  }

  // Print memory statistics for a GraphCache
  static std::string print_cache_memory(const Graph::GraphCache& cache)
  {
    size_t cache_size = calculate_graph_cache_size(cache);
    
    std::ostringstream oss;
    oss << "GraphCache Memory Usage: " << format_bytes(cache_size);
    
    return oss.str();
  }
};

}

#endif // MEMORY_UTILS_H
