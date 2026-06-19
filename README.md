# Arbor

### From raw forest scan to thousands of QSMs, on a laptop, in minutes.

[![License: GPL-3](https://img.shields.io/badge/License-GPL3-green.svg)](LICENSE)

[📖 Book & Tutorial](https://r-lidar.github.io/arborBook/) · [🌐 r-lidar.com](https://www.r-lidar.com/)

---

Arbor is a **production-grade** C++ library with an R API for processing Mobile and Terrestrial Laser Scanning point clouds of forests. It covers the full pipeline — ground classification, wood/foliage segmentation, individual tree detection, and QSM reconstruction — automatically, without parameter tuning, on real-world noisy data and at an unmatched speed.

No supervision. No tweaking. Works leaf-on.

## What makes it different?

**It scales.** 900 m² in 1 minute. 5 000 m² in 8 minutes. 10 QSMs per second. On a laptop.

**It is easy.** Every dataset goes through the same pipeline with the same default parameters and comes out the other end with accurate 3D models of the forest.

**It is documented.** Arbor is deeply documented in a [dedicated book](https://r-lidar.github.io/arborBook/).

**It is validated.** On an exceptionally large dataset of real data.

## Installation

```r
install.packages("arbor")
```

Or if not yet on CRAN:

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
qsf <- qsf(las)          # one QSM per tree, in parallel
plot(qsf, pal = "chocolate4")
```

## Acknowledgements

Developed by [r-lidar inc.](https://www.r-lidar.com/) with support from Sherbrooke University, Natural Resources Canada, IRD Montpellier, the Québec Ministry of Natural Resources and Forests, and Ontario Forestry Futures Trust.
