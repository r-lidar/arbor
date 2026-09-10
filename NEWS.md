# 1.1.0

### NEW FEATURES

- `segment_instance()` now detects and fixes over-segmentation, i.e. single
  trees split across two or more treeIDs. This occurs occasionally for large
  trees and especially often for buttress trees. The module is enabled by
  default and can be turned off with:
  ```r
    arbor_parameters_default$instance$oversegmentation_solver_enabled = FALSE
  ```
  If disabled at segmentation time, it can still be run afterwards — e.g. in
  post-processing, on an already segmented point cloud — with:
  ```r
    las <- resolve_oversegmentation(las)
  ```
  
### BUG FIXES

- Fix #15: Prevent QSM generation failure on non-constructible skeletons by defaulting to a zero-volume placeholder edge.


# 1.0.0

Public release
