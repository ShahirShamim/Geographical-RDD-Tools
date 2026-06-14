# 🌍 geoRDDprep

[![PyPI version](https://badge.fury.io/py/geoRDDprep.svg)](https://badge.fury.io/py/geoRDDprep)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)

**Streamline your Geographical Regression Discontinuity Design (GeoRDD) workflow.**

`geoRDDprep` is a high-performance Python toolkit designed to simplify spatial data preparation for boundary-analysis. Whether you are an economist, political scientist, or data analyst, this package helps you assign points to districts, clean up messy polygons, and implement rigorous spatial algorithms (such as the Turner orthogonal distance criteria) with ease.

---

## 🚀 Key Features

*   **⚡️ Vectorized & Fast**: Rewritten to use fully vectorized operations via `geopandas` and `shapely`, delivering up to **18x speedups** compared to standard row-by-row iteration.
*   **📐 Bug-Free Turner Algorithm**: Out-of-the-box implementation of the orthogonal distance criteria from *Turner et al. (2014)*, with robust math handling both horizontal and vertical boundary projections.
*   **🧹 Smart Sliver Cleaning**: Merge sliver polygons and boundary gaps using Voronoi diagrams with **automatic ID mapping** and **dynamic padding** (making it compatible with both degree and metric coordinate systems).
*   **🔄 Automatic CRS Alignment**: Automatically detects Coordinate Reference System (CRS) mismatches and re-projects inputs dynamically (issuing a warning instead of crashing). Fully supports naive geometries (where `crs = None`).
*   **🛠️ Easy Integration**: Integrates seamlessly with your existing `pandas` and `geopandas` pipelines.

---

## 📦 Installation

Install directly from PyPI:

```bash
pip install geoRDDprep
```

---

## 🛠️ API Reference

### 1. `points_in_polygon`
Assigns polygon characteristics to points that fall within them.
```python
def points_in_polygon(
    points_gdf: gpd.GeoDataFrame, 
    polygons_gdf: gpd.GeoDataFrame, 
    suffix_name: str
) -> gpd.GeoDataFrame
```
*   **`points_gdf`**: Point geometries.
*   **`polygons_gdf`**: Polygon geometries with attributes to join.
*   **`suffix_name`**: Suffix to append to joined polygon columns. Overlapping column names are cleanly renamed with a single underscore (e.g., `id_district` instead of `id__district`).

### 2. `poly_to_line`
Converts polygon boundaries into LineStrings.
```python
def poly_to_line(polygon_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame
```
*   **`polygon_gdf`**: Polygons or MultiPolygons to convert. Returns exploded boundary `LineString` elements, preserving all original attributes.

### 3. `turner`
Verifies if points satisfy the Turner et al. (2014) orthogonal distance criteria relative to boundaries.
```python
def turner(
    points_gdf: gpd.GeoDataFrame, 
    boundaries_gdf: gpd.GeoDataFrame, 
    *,
    orth_distance: float = 15.0,
    reduced: bool = True,
    unit_crs: int = 3857
) -> gpd.GeoDataFrame
```
*   **`points_gdf`**: Point geometries.
*   **`boundaries_gdf`**: LineString boundary geometries.
*   **`orth_distance`** *(keyword-only)*: Orthogonal distance threshold (in meters). Default `15`.
*   **`reduced`** *(keyword-only)*: If `True`, returns only original columns plus the `turner_pass` boolean result. Default `True`.
*   **`unit_crs`** *(keyword-only)*: EPSG code for metric distance calculation. Default `3857` (Web Mercator).

### 4. `drop_tiny_lines`
Filters out small boundary LineStrings to reduce map noise.
```python
def drop_tiny_lines(
    boundaries_gdf: gpd.GeoDataFrame, 
    method: str = 'percentile', 
    *,
    percentile: float = 0.01,
    num_dev: float = 2.0,
    meters: float = 500.0,
    reduced: bool = True,
    unit_crs: int = 3857
) -> gpd.GeoDataFrame
```
*   **`method`**: Threshold method (`'percentile'`, `'number_of_std'`, or `'length'`).
*   **`percentile`** *(keyword-only)*: Quantile threshold (0-1). Default `0.01`.
*   **`num_dev`** *(keyword-only)*: Number of standard deviations below the mean. Default `2.0`.
*   **`meters`** *(keyword-only)*: Length cutoff in meters. Default `500.0`.

### 5. `remove_sliver`
Cleans sliver polygons and gaps by assigning them to their nearest neighbor using a Voronoi diagram.
```python
def remove_sliver(
    polygons_gdf: gpd.GeoDataFrame, 
    boundary_gdf: gpd.GeoDataFrame,
    *,
    id_col: Optional[str] = None
) -> gpd.GeoDataFrame
```
*   **`polygons_gdf`**: Input polygons to clean.
*   **`boundary_gdf`**: Bounding geometry to clip the output.
*   **`id_col`** *(keyword-only)*: Name of the unique identifier column. If `None`, automatically searches for `'id'`, checks the index name, or uses default indexing.

### 6. `remove_overlaps`
Removes overlapping line segments from `df1` that are present in `df2`.
```python
def remove_overlaps(df1: gpd.GeoDataFrame, df2: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
```
*   **`df1`**: The GeoDataFrame containing LineStrings to clean.
*   **`df2`**: The GeoDataFrame containing geometries to subtract.

---

## 🛠️ Usage Examples

### 1. Assign Addresses to Districts
```python
import geopandas as gpd
from geoRDDprep import points_in_polygon

points = gpd.read_file("addresses.geojson")
districts = gpd.read_file("school_districts.geojson")

# Merges district characteristics into matching points
result = points_in_polygon(points, districts, suffix_name="_district")
print(result.head())
```

### 2. The Turner Algorithm (2014)
Check if points are within 15 meters orthogonal distance of school boundaries and not close to vertices or endpoints.
```python
from geoRDDprep import poly_to_line, drop_tiny_lines, turner

# 1. Convert school district polygons to boundary lines
lines = poly_to_line(districts)

# 2. Remove tiny boundary segments (less than 500m) to reduce noise
clean_lines = drop_tiny_lines(lines, method='length', meters=500)

# 3. Match points to boundaries (within 15m)
matched_data = turner(points, clean_lines, orth_distance=15)

# Check which points passed the Turner check
print(matched_data['turner_pass'].value_counts())
```

### 3. Clean Slivers and Gaps
```python
from geoRDDprep import remove_sliver

# Merge gaps/slivers into neighbor polygons using a custom identifier column
clean_polygons = remove_sliver(messy_polygons, boundary_clip, id_col="district_code")
```

---

## 🧪 Running Tests

To verify package modifications, you can run the test suite using `pytest`.

1. Install test dependencies:
   ```bash
   pip install pytest
   ```
2. Run tests:
   ```bash
   pytest tests/
   ```

---

## 🤝 Contributing

We welcome contributions!
1. Fork the repository.
2. Create a feature branch (`git checkout -b feature/AmazingFeature`).
3. Commit your changes (`git commit -m 'Add AmazingFeature'`).
4. Push to the branch (`git push origin feature/AmazingFeature`).
5. Open a Pull Request.

---

## 📄 License

Distributed under the MIT License. See `LICENSE` for more information.
