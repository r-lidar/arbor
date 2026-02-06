# arbor 0.8.0

## New features
- Introduced the `qsm` class with dedicated `print()` and `plot()` methods.
- Added the `qsf()` function, replacing `qsm_batch()` and introducing the Quantitative Structure Forest (QSF) concept and class.
- Added `qsf_write()` and a suite of `qsf_*()` utility functions.
- Added experimental function `segment_semantic_from_qsm()`.
- Added experimental function `extract_tree_context()`.
- Added experimental function `editor_seeds()`.

## Enhancements
- Significantly improved rendering performance of `plot()` / `plot_qsm()`: QSMs now render instantaneously instead of taking several seconds.
- Extended `plot()` to support rendering QSFs (entire forests).

## Bug fixes
- Fixed `colorize_trees()` when `darken_foliage = FALSE`.
- Fixed `qsm()` for trees taller than 36 m

# arbor 0.7.2

- Fix QSM. There is no longer any irrelevant and biologically impossible branches.
- Polynomial fitting applied only on missing diameters and outliers. Raw measurements preserved.

# arbor 0.7.1

- Default breast height is 1.30 m instead of 1.37 m
- `qsm_stats()` can pass arguments to `qsm_dbh()`.

# arbor 0.7.0

- `compute_anisotopy()` is twice faster
- `connected_component()` (internal function) is 5 times faster which provides a boost to several functions.
- SOR filter is 7 times faster and brings a significant boost to `segment_foliage()` 

# arbor 0.6.0

- QSMs in C++ are faster
- new function `qsm_dbh()`
- new function `qsm_stats()`

# arbor 0.5.0

- Improved QSM

# arbor 0.3.0

- Added QSM in arbor
