# arbor 0.9.1

- `colorize_trees()` in C++
- NEW: `qsf_treemap()`

# arbor 0.9.0

- Change: a lot of C++. 
- New: better QSM meshing
- New: better QSM smoothing
- Enhancement: drastic speed up of `qsf()`
- Change: no need to compute a DTM
- Change: `segment_semantic()` no longer requires a DTM
- Change: `cut_above_ground` is a parameters no longer a variable
- Change: `remove_small_trees()` no longer clips. It assigns NAs. Up to the user to clip.

# arbor 0.8.2

- Enhancement: `find_seeds()` full C++
- Fix: bug `filter_range()`

# arbor 0.8.1

- Enhancement: there are no longer two versions of ransac fitting in R and C++. All computations are done in C++
- Enhancement: lot of C++ change
- Enhancement: `segment_semantic()` 20% faster
- Fix: bug with small trees in `qsm()`
- Fix: bug when rendering `qsf` with `plot()`
- Fix: bugs in arbor studio

# arbor 0.8.0

## New features

- Introduced the `qsm` class with dedicated `print()` and `plot()` methods.
- Added the `qsf()` function, replacing `qsm_batch()` and introducing the Quantitative Structure Forest (QSF) concept and class.
- Added `qsf_write()` and a suite of `qsf_*()` utility functions.
- Added experimental function `qsf_segment_semantic()`.
- Added experimental function `extract_tree_context()`.
- Added experimental function `arbor_studio_seeds()`.
- Added experimental function `arbor_studio_instance()`.

## Enhancements

- Significantly improved rendering performance of `plot()` / `plot_qsm()`: QSMs now render instantaneously instead of taking several seconds.
- Extended `plot()` to support rendering QSFs (entire forests).

## Bug fixes

- Fixed `colorize_trees()` when `darken_foliage = FALSE`.
- Fixed `qsm()` for trees taller than 36 m

## Changes

- Change: `wood_likelihood()` replaces `compute_anisotropy()`
- Change: `plot_likelihood()` replaces `plot_anisotropy()`
- Change: `segment_semantic()` replaces `segment_foliage()`
- Change: `segment_instance()` replaces `segment_vegetation()`

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
