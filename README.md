# 🌍 geoRDDprep

[![PyPI version](https://badge.fury.io/py/geoRDDprep.svg)](https://badge.fury.io/py/geoRDDprep)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)

**Streamline your Geographical Regression Discontinuity Design (GeoRDD) workflow.**

`geoRDDprep` is a powerful Python toolkit designed to take the pain out of spatial data preparation. Whether you are an economist, political scientist, or data analyst, this package helps you assign points to boundaries, clean up messy polygons, and implement rigorous spatial algorithms with ease.

## 🚀 Why geoRDDprep?

*   **⚡️ Fast & Efficient**: Optimized spatial joins and geometric operations using `geopandas` and `shapely`.
*   **📐 Turner Algorithm Ready**: Out-of-the-box implementation of the orthogonal distance criteria from *Turner et al. (2014)*.
*   **🧹 Data Cleaning**: Automatically remove sliver polygons and tiny, noisy line segments that mess up your analysis.
*   **🛠️ Easy Integration**: Works seamlessly with your existing `pandas` and `geopandas` workflows.

## 📦 Installation

Install directly from PyPI:

```bash
pip install geoRDDprep
```

## 🛠️ Usage Examples

### 1. Assign Addresses to Districts
Map millions of points to their administrative regions in seconds.

```python
import geopandas as gpd
from geoRDDprep import points_in_polygon

# Load your data
points = gpd.read_file("addresses.geojson")
districts = gpd.read_file("school_districts.geojson")

# 🪄 Magic happens here
result = points_in_polygon(points, districts, suffix_name="_district")

print(result.head())
```

### 2. The Turner Algorithm (2014)
Assign points to boundaries only if they are within a strict orthogonal distance—crucial for valid RDD analysis.

```python
from geoRDDprep import poly_to_line, drop_tiny_lines, turner

# 1. Convert polygons to boundary lines
lines = poly_to_line(districts)

# 2. Clean up noise (remove lines < 500m)
clean_lines = drop_tiny_lines(lines, method='length', meters=500)

# 3. Match points to boundaries (within 15m)
matched_data = turner(points, clean_lines, orth_distance=15)

# Check which points passed the test
print(matched_data['turner_pass'].value_counts())
```

### 3. Clean Messy Polygons
Got "slivers" or gaps in your map? Fix them automatically.

```python
from geoRDDprep import remove_sliver

# Merge slivers into their largest neighbors
clean_polygons = remove_sliver(messy_polygons, boundary_clip)
```

## 🤝 Contributing

We love contributions!
1.  Fork the repo.
2.  Create a branch (`git checkout -b feature/AmazingFeature`).
3.  Commit your changes (`git commit -m 'Add AmazingFeature'`).
4.  Push to the branch (`git push origin feature/AmazingFeature`).
5.  Open a Pull Request.

## 📄 License

Distributed under the MIT License. See `LICENSE` for more information.

---

*Built with ❤️ for the spatial data community.*
