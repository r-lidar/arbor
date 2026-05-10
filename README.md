# Arbor

High-performance tree segmentation and structural modeling for TLS/MLS point clouds.

**arbor** is a production-grade C++ library for large-scale processing of terrestrial and mobile laser scanning data of forests. It introduces a new generation of algorithms for tree segmentation and modeling that significantly advance the state of the art in forest point-cloud analysis.

Designed for both **accuracy and scalability**, arbor can process very large datasets covering thousands of square meters within few minutes, compute dozens of Quantitative Structure Models (QSMs) per second, and robustly segment forests ranging from simple stands to extremely complex tropical ecosystems.

Arbor maintains high segmentation accuracy in dense understory and low vegetation, enabling reliable reconstruction of complete forest structure from mature trees to saplings.

📖 See the [Arbor Book](https://r-lidar.github.io/arborBook/) for tutorials and examples.

---

# Features

- 🌲 **Advanced semantic segmentation:** classifies wood, foliage, and ground with high accuracy even in dense forests and complex understory environments.
- 🌳 **Instance segmentation: automatically** detects and separates individual trees in large point clouds.
- 🪵 **Quantitative Structure Models (QSM):** reconstructs detailed tree skeletons and branching architecture at high throughput (dozens of QSMs per second).
- 🌐 **Quantitative Structure Forest (QSF):** extends QSM reconstruction to entire forest scenes, enabling structural analysis at ecosystem scale.
-  ⚡ **Designed for massive datasets:** optimized algorithms allow arbor to process very large TLS/MLS point clouds efficiently, making it suitable for production pipelines and large forest inventories.
- 📦 **Interoperable outputs:** support results to CSV, OBJ, PLY, and STL for analysis and 3D visualization.
- 🖥 **Fast 3D visualization:** built-in tools to quickly render and inspect forest scenes.

---

# Installation

```r
# Install from GitHub
devtools::install_github("r-lidar-lab/arbor")
```
