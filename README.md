# Arbor

### From raw forest scans to thousands of QSMs on your laptop, in minutes.

[![License: GPL-3](https://img.shields.io/badge/License-GPL3-green.svg)](LICENSE)

[📖 Book & Tutorial](https://r-lidar.github.io/arborBook/) · [🌐 r-lidar.com](https://www.r-lidar.com/)

---

Arbor is an R package with a fast C++ core designed for processing Mobile Laser Scanning (MLS) point clouds in forestry. It takes you through the entire workflow: from initial ground classification to final QSM reconstruction, including wood/foliage semantic segmentation and individual tree instance segmentation.

It handles real-world, noisy data, leaf-on scans, automatically without requiring parameter tuning or manual intervention.


[![Arbor demo](https://raw.githubusercontent.com/r-lidar/arbor/main/man/figures/yt-play-readme-500px.png)](https://youtu.be/GIbaMfSo6dc)

## Why Arbor?

- **Built to scale:** Process 900 m² in about a minute, a quarter-hectare in 8 minutes. Compute roughly 700 QSMs per minute. All on a standard laptop.
- **Consistent results out of the box:** Run your datasets through a clear, repeatable pipeline using default settings that just work.
- **Well documented:** Explore step-by-step guides, theoretical context, and hands-on examples in the [Arbor Book](https://r-lidar.github.io/arborBook/).
- **Grounded in real data:** Validated on an extensive collection of real-world forest scans.

## Installation

```r
install.packages('arbor', repos = c('https://r-lidar.r-universe.dev', 'https://cloud.r-project.org'))
```

## Learn more

The **[Arbor Book](https://r-lidar.github.io/arborBook/)** is a complete, illustrated guide to the full pipeline with example datasets you can download and run right now.

Minimal Reproducible Example:

```r
las <- readTLS("forest.laz", select = "xyz", filter = "-keep_random_fraction 0.25")
las <- hybrid_homogeneization(las)
las <- segment_ground(las)
las <- wood_likelihood(las)
las <- segment_semantic(las)
see <- find_seeds(las)
las <- segment_instance(las, see)
qsf <- qsf(las)
plot(qsf, pal = "chocolate4")
```
