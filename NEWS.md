# arbor 0.7.3

- New: QSMs now have a class `qsm` and dedicated methods `print()` and `plot()`. `plot_qsm()` has been removed.
- Enhance: `plot()` renders QSMs instantaneously. It no longer takes several seconds.
- Fix `colorize_trees()` with `darken_foliage = FALSE`
- Fix `qsm()` for trees > 36 m

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
