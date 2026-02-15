#pragma once

#include <vector>

struct GraphBuilderParams
{
  bool downward = false;
  int k = 10;
  int k_seed = 100;
  float decimation = 0.05f;
  float space_res = 0.2;
  float max_gap = 0.2f;
  float power = 3.0f;
  float wood2wood = 0.1;
  float leaf2leaf = 20;
  float wood2leaf = 100;
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

  float high_pwood_threshold   = 0.9f;
  float medium_pwood_threshold = 0.75f;

  float connected_components_res = 0.05f;
  int   connected_components_min = 2000;

  float wood_assignation_dist = 0.05f;

  int   wood_extra_reasignation_k = 10;
  float wood_extra_reasignation_dist = 0.03f;

  int   medium_pwood_sor_k = 50;
  float medium_pwood_sor_m = 0.05f;

  float ground_res = 0.2f;
};

