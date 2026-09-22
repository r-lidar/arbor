# v1.1.0

### NEW FEATURES

1. `segment_instance()` now detects and fixes over-segmentation, i.e. single
  trees split across two or more tree IDs. This occurs occasionally for large
  trees and especially often for buttress trees. The module is enabled by
  default and can be turned off with:
    ```r
    arbor_parameters_default$instance$oversegmentation_solver_enabled = FALSE
    ```
    If disabled at segmentation time, it can still be run afterwards, e.g. in
  post-processing, on an already segmented point cloud, with:
    ```r
    las <- resolve_oversegmentation(las)
    ```

2. `qsf_write()` can now merge an entire forest into a single combined
  mesh file. Documentation for `qsf_write()` now explains, format by format, how 
  well individual trees stay distinguishable when several QSMs are combined
  into one mesh file.
    ```r
    qsf_write(qsf, "forest.obj")
    ```

3. New QSF file format. When writing QSMs to `.qsm` files, `qsf_write()` now also 
  creates a `.qsf` file that indexes all the `.qsm` files. A QSF can be reloaded with 
  `qsf_read()` by reading the `.qsf` index file.
    ```r
    qsf_write(qsf, "forest.qsf")
    qsf <- qsf_read("forest.qsf")
    ```

4. `qsf_log()` has been redesigned. It is easier to read, and easier to manipulated.

5. New function `qsf_select()` that allows to filter QSMs in the QSF object by message,
  by message type, by DBH, by height.
    ```r
    # Keep only trees flagged with code W2
    qsf_select(qsf, code = c("W0", "W3"))
    # Exclude trees flagged with code W2 or W3, while keeping trees with DBH > 20 cm
    qsf_select(qsf, code = c("W2", "W3"), invert = TRUE, dbh = c(0.2, Inf))
    ```

### ENHANCES

1. `find_seeds()` no longer throw the error "No circle detected in wood slices", 
  No shape detection is no longer an issue.
  
### BUG FIXES

1. Fix #15: Prevent QSM generation failure on non-constructible skeletons by 
  defaulting to a zero-volume placeholder edge.
  
2. Fix: QSM messages were not written in `.qsm` files. They are now recorded 
  properly.


# v1.0.0

Public release
