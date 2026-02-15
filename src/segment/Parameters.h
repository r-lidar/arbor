#pragma once

#include <vector>

struct GraphBuilderParams
{
  bool downward = false;
  int k = 10;
  int k_seed = 100;
  double decimation = 0.05;
  double space_res = 0.2;
  double max_gap = 0.2;
  double power = 3.0;
  double wood2wood = 0.1;
  double leaf2leaf = 20;
  double wood2leaf = 100;
  std::vector<float> angle_penalty;
  GraphBuilderParams(): angle_penalty(181)
  {
    for (int x = 0; x <= 180; ++x)
    {
      float y = std::exp(0.046051f * x);
      angle_penalty[x] = (x > 100) ? 100.0f : y;
    }
  }
};

struct SemanticParams
{
  int   min_passage = 3;
  double high_pwood_threshold   = 0.9;
  double medium_pwood_threshold = 0.75;
  double connected_components_res = 0.05;
  int    connected_components_min = 2000;
  double wood_assignation_dist = 0.05f;
  int    wood_extra_reasignation_k = 10;
  double wood_extra_reasignation_dist = 0.03;
  int    medium_pwood_sor_k = 50;
  double medium_pwood_sor_m = 0.05;
  double ground_res = 0.2;
};

