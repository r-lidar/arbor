# arbor 0.8.0

- New: QSMs now have a class `qsm` and dedicated methods `print()` and `plot()`.
- New: function `qsf` replaces `qsm_batch()` and introduces the concept of Quantitative Structure Forest
- New: functions `qsf_write()` and other `qsf_*()` tools.
- Enhance: `plot()` or `plot_qsm()` render QSMs instantaneously. It no longer takes several seconds. `plot()` can render QSF to render the entire forest.
- Fix `colorize_trees()` with `darken_foliage = FALSE`
- Fix `qsm()` for trees > 36 m
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
