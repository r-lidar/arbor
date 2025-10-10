#' @export
default_parameters = list(
  decimation = list(
    barycentric_predecimation_resolution = 0.05,
    cut_above_ground = 0.25
  ),
  anistotropy = list(
    k = 50
  ),
  path_finder = list(
    k_neighborhood_connectivity = 10,
    max_gap = 1,
    min_passage = 2,
    angle_penalty_function = function(x){ y = exp(log(100)/100*x); y[x>100]=100; y }
  ),
  semantic = list(
    z_scale = 0.8,
    high_anisotropy_threshold = 0.9,
    medium_anisotropy_thresold = 0.75,
    connected_components_res = 0.05,
    connected_components_min = 2000,
    wood_assignation_k = 50,
    wood_assignation_dist = 0.05,
    wood_extra_reasignation_k = 10,
    wood_extra_reasignation_dist = 0.03,
    medium_anisotropy_sor_k = 50,
    medium_anisotropy_sor_m = 0.05
  ),
  seed = list (
    slice_at = c(1,2,3),
    slice_thickness = 0.03,
    sor_k = 10,
    sor_m = 0.5,
    connected_components_res = 0.06,
    connected_components_min = 1,
    safe_zone = 0.2
  ),
  instance = list(
    z_scale = 0.8,
    wood2wood_factor = 0.1,
    leaf2leaf_factor = 20,
    wood2leaf_factor = 1000
  ),
  post_prod = list(
    isolated_wood = list(
      z_scale = 0.1,
      connected_components_res = 0.05,
      connected_components_min = 200
    )
  )
)
