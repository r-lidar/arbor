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

- `qsf_write()` can now merge an entire forest into a single combined
  `.stl` file, joining the existing single-file support for `.obj` and
  `.ply`. Documentation for `qsf_write()` now explains, format by format, how 
  well individual trees stay distinguishable when several QSMs are combined
  into one file.

- New QSF file format. When writing QSMs to `.qsm` files, `qsf_write()` now also creates a `.qsf` file that indexes all the `.qsm` files. A QSF can be reloaded with `qsf_read()` by reading the `.qsf` index file.
  ```r
  qsf_write(qsf, "forest.qsf")
  qsf <- qsf_read("forest.qsf")
  ```

### ENHANCES

- `find_seeds()` no longer throw the error "No circle detected in wood slices", 
  No shape detection is no longer an issue.
  
### BUG FIXES

- Fix #15: Prevent QSM generation failure on non-constructible skeletons by 
  defaulting to a zero-volume placeholder edge.


# 1.0.0

Public release
