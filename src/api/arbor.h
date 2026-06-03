/**
 * @file arbor.h
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

#ifndef API_H
#define API_H

#include "PointCloud.h"
#include "QSM.h"
#include "QSF.h"
#include "services.h"

#include <vector>

#define ARBOR_VERSION_MAJOR 0
#define ARBOR_VERSION_MINOR 13
#define ARBOR_VERSION_PATCH 0

namespace arbor
{
  namespace settings
  {
    struct GlobalParameters
    {
      double cut_above_ground = 0.25;

      template<typename V> void visit(V& v)
      {
        v("cut_above_ground", cut_above_ground);
      }
    };

    struct WoodlikelihoodParameters
    {
      int k = 80;

      template<typename V> void visit(V& v)
      {
        v("k", k);
      }
    };

    struct GraphParameters
    {
      bool   downward   = false;
      int    k          = 10;
      int    k_seed     = 100;
      double decimation = 0.05;
      double space_res  = 0.2;
      double max_gap    = 1.0;
      double power      = 3.0;
      double wood2wood  = 0.1;
      double leaf2leaf  = 20.0;
      double wood2leaf  = 1000.0;
      std::vector<float> angle_penalty;

      GraphParameters() : angle_penalty(181)
      {
        for (int x = 0; x <= 180; ++x)
        {
          float y = std::exp(0.046051f * x);
          angle_penalty[x] = (x > 100) ? 100.0f : y;
        }
      }

      template<typename V> void visit_pathfinder(V& v)
      {
        v("k_neighborhood_connectivity", k);
        v("k_seed_connectivity",         k_seed);
        v("decimation",                  decimation);
        v("space_res",                   space_res);
        v("max_gap",                     max_gap);
        v("distance_power",              power);
        v("downward",                    downward);
      }

      template<typename V> void visit_instance(V& v)
      {
        v("wood2leaf_factor", wood2leaf);
        v("leaf2leaf_factor", leaf2leaf);
        v("wood2wood_factor", wood2wood);
      }
    };

    struct SemanticParameters
    {
      int    min_passage                  = 3;
      double high_pwood_threshold         = 0.9;
      double medium_pwood_threshold       = 0.75;
      double connected_components_res     = 0.05;
      int    connected_components_min     = 2000;
      int    wood_assignation_k           = 50;
      double wood_assignation_dist        = 0.05;
      int    wood_extra_reasignation_k    = 10;
      double wood_extra_reasignation_dist = 0.03;
      int    medium_pwood_sor_k           = 50;
      double medium_pwood_sor_m           = 0.05;
      double ground_res                   = 0.2;

      template<typename V> void visit(V& v)
      {
        v("min_passage",                  min_passage);
        v("high_pwood_threshold",         high_pwood_threshold);
        v("medium_pwood_threshold",       medium_pwood_threshold);
        v("connected_components_res",     connected_components_res);
        v("connected_components_min",     connected_components_min);
        v("wood_assignation_k",           wood_assignation_k);
        v("wood_assignation_dist",        wood_assignation_dist);
        v("wood_extra_reasignation_k",    wood_extra_reasignation_k);
        v("wood_extra_reasignation_dist", wood_extra_reasignation_dist);
        v("medium_pwood_sor_k",           medium_pwood_sor_k);
        v("medium_pwood_sor_m",           medium_pwood_sor_m);
        v("ground_res",                   ground_res);
      }
    };

    struct SeedParameters
    {
      std::vector<double> slice_at       = {0.7, 0.9};
      double              slice_thickness = 0.05;
      int                 min_passage    = 15;
      double              safe_zone      = 0.2;

      template<typename V> void visit(V& v)
      {
        v("slice_thickness", slice_thickness);
        v("min_passage",     min_passage);
        v("safe_zone",       safe_zone);
        // slice_at handled separately (vector needs Rcpp::wrap)
      }
    };

    struct QsmParameters
    {
      float       skeleton_node_distance   = 0.1;
      float       dbscan_eps_distance      = 0.1;
      float       max_d                    = 0.1;
      float       apex_radius              = 0.0025;
      int         smooth_steps             = 15;
      float       min_measurable_dbh       = 0.06;
      float       min_measurable_radius    = 0.03;
      bool        broken_detection_enabled = true;
      float       allometry_scale          = 1.0f;
      std::string allometry_name           = "Griese2025";

      template<typename V> void visit(V& v)
      {
        v("skeleton_node_distance",   skeleton_node_distance);
        v("dbscan_eps_distance",      dbscan_eps_distance);
        v("max_d",                    max_d);
        v("apex_radius",              apex_radius);
        v("smooth_steps",             smooth_steps);
        v("min_measurable_dbh",       min_measurable_dbh);
        v("min_measurable_radius",    min_measurable_radius);
        v("broken_detection_enabled", broken_detection_enabled);
        v("allometry_name",           allometry_name);
        v("allometry_scale",          allometry_scale);
      }
    };

    struct ArborParameters
    {
      GlobalParameters         global;
      WoodlikelihoodParameters woodlikelihood;
      GraphParameters          pathfinder;
      SemanticParameters       semantic;
      SeedParameters           seeds;
      QsmParameters            qsm;
    };
  }
  namespace segment
  {
    void segment_ground  (PointCloud& core,                          const settings::ArborParameters& params);
    void segment_instance(PointCloud& core, const PointCloud& seeds, const settings::ArborParameters& params);
    void segment_semantic(PointCloud& core, const PointCloud& dtm,   const settings::ArborParameters& params);

    std::vector<float> dist2root(const PointCloud& core, const PointCloud& dtm, const settings::GraphParameters& params);
  }

  namespace seeds
  {
    PointCloud find_seeds(const PointCloud&,  const settings::ArborParameters& params);
  }

  namespace qsm
  {
    QSM qsm(const PointCloud&,  const settings::ArborParameters& params);
    QSF qsf(const PointCloud&, double min_height, const settings::ArborParameters& params);
  }

  namespace dtm
  {
    PointCloud dtm(const PointCloud&);
  }


  namespace utils
  {
    std::vector<bool> homogeneization(const PointCloud& pc, double res, bool hybrid = true);
    std::vector<bool> sor(const PointCloud& pc, unsigned int k, double m);
    std::vector<float> anisotropy(const PointCloud& pc, int k);
    PointCloud smooth3d(const PointCloud& pc, double radius, int ncores);
  }
}

#endif
