default_parameters_v1 = list(
  anistotropy = list(
    k = 80
  ),
  path_finder = list(
    k_neighborhood_connectivity = 10,
    k_seed_connectivity = 100,
    decimation = 0.05,
    space_res = 0.2,
    max_gap = 0.2,
    distance_power = 3,
    angle_penalty = function(x) { y = exp(0.046051*x); ifelse(x > 100, 100, y) }
  ),
  semantic = list(
    min_passage = 3,
    high_pwood_threshold = 0.9,
    medium_pwood_thresold = 0.75,
    connected_components_res = 0.05,
    connected_components_min = 2000,
    wood_assignation_k = 50,
    wood_assignation_dist = 0.05,
    wood_extra_reasignation_k = 10,
    wood_extra_reasignation_dist = 0.03,
    medium_pwood_sor_k = 50,
    medium_pwood_sor_m = 0.05,
    ground_res = 0.2
  ),
  seed = list (
    slice_at = c(0.7,0.9),
    slice_thickness = 0.05,
    sor_k = 10,
    sor_m = 0.5,
    min_passage = 15,
    safe_zone = 0.2
  ),
  instance = list(
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

#' Parameters
#' @rdname parameters
#' @name parameters
#' @export
default_arbor_parameters = default_parameters_v1
