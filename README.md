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

### 7. `calculate_signed_distance`
Calculates the absolute and signed distance from points to a boundary. Points inside a treatment polygon receive a positive distance; points outside receive a negative distance.
```python
def calculate_signed_distance(
    points_gdf: gpd.GeoDataFrame,
    boundary_gdf: gpd.GeoDataFrame,
    treatment_gdf: gpd.GeoDataFrame,
    *,
    distance_col: str = 'distance',
    signed_distance_col: str = 'signed_distance',
    treatment_col: str = 'is_treated',
    unit_crs: int = 3857
) -> gpd.GeoDataFrame
```
*   **`points_gdf`**: Point geometries.
*   **`boundary_gdf`**: Boundary LineString geometries.
*   **`treatment_gdf`**: Treatment polygon geometry defining the treated area.
*   **`distance_col`** *(keyword-only)*: Column name for absolute distance. Default `'distance'`.
*   **`signed_distance_col`** *(keyword-only)*: Column name for signed distance. Default `'signed_distance'`.
*   **`treatment_col`** *(keyword-only)*: Column name for the treatment indicator. Default `'is_treated'`.
*   **`unit_crs`** *(keyword-only)*: EPSG code for metric distance calculation. Default `3857`.

### 8. `extract_shared_boundaries`
Extracts shared border lines (interfaces) between adjacent polygons in a GeoDataFrame.
```python
def extract_shared_boundaries(
    gdf: gpd.GeoDataFrame,
    *,
    id_col: Optional[str] = None
) -> gpd.GeoDataFrame
```
*   **`gdf`**: GeoDataFrame containing Polygon/MultiPolygon geometries.
*   **`id_col`** *(keyword-only)*: Name of the ID column. If `None`, uses index.

### 9. `shift_boundary_placebo`
Shifts/translates a boundary by a specified offset in meters. Useful for generating placebo borders.
```python
def shift_boundary_placebo(
    boundary_gdf: gpd.GeoDataFrame,
    xoff: float = 0.0,
    yoff: float = 0.0,
    *,
    unit_crs: int = 3857
) -> gpd.GeoDataFrame
```
*   **`boundary_gdf`**: Boundary geometries.
*   **`xoff`**: Translation offset in the X direction (meters). Default `0.0`.
*   **`yoff`**: Translation offset in the Y direction (meters). Default `0.0`.
*   **`unit_crs`** *(keyword-only)*: EPSG code for metric calculation. Default `3857`.

### 10. `filter_by_boundary_distance`
Filters points to only those within a specified bandwidth (distance threshold) of the boundary.
```python
def filter_by_boundary_distance(
    points_gdf: gpd.GeoDataFrame,
    boundary_gdf: gpd.GeoDataFrame,
    max_distance: float,
    *,
    unit_crs: int = 3857
) -> gpd.GeoDataFrame
```
*   **`points_gdf`**: Point geometries to filter.
*   **`boundary_gdf`**: Boundary geometries.
*   **`max_distance`**: Maximum distance threshold (meters).
*   **`unit_crs`** *(keyword-only)*: EPSG code for metric calculation. Default `3857`.

### 11. `assign_nearest_boundary`
Assigns each point to its nearest boundary feature, attaching the boundary's identifier and the metric distance. Pair it with `segment_boundary` to build **boundary-segment fixed effects**.
```python
def assign_nearest_boundary(
    points_gdf: gpd.GeoDataFrame,
    boundary_gdf: gpd.GeoDataFrame,
    *,
    id_col: Optional[str] = None,
    boundary_id_col: str = 'boundary_id',
    distance_col: str = 'boundary_distance',
    unit_crs: int = 3857
) -> gpd.GeoDataFrame
```
*   **`points_gdf`**: Point geometries.
*   **`boundary_gdf`**: Boundary geometries.
*   **`id_col`** *(keyword-only)*: Column in `boundary_gdf` holding the identifier to attach. If `None`, the boundary index is used.
*   **`boundary_id_col`** *(keyword-only)*: Output column for the boundary identifier. Default `'boundary_id'`.
*   **`distance_col`** *(keyword-only)*: Output column for the distance to the nearest boundary (meters). Default `'boundary_distance'`.
*   **`unit_crs`** *(keyword-only)*: EPSG code for metric calculation. Default `3857`.

### 12. `snap_points_to_boundary`
Projects each point onto its nearest boundary, returning the snapped (projected) point geometry. Useful for the "boundary point" RD approach and for visualization.
```python
def snap_points_to_boundary(
    points_gdf: gpd.GeoDataFrame,
    boundary_gdf: gpd.GeoDataFrame,
    *,
    snapped_col: str = 'snapped_geometry',
    distance_col: Optional[str] = None,
    unit_crs: int = 3857
) -> gpd.GeoDataFrame
```
*   **`points_gdf`**: Point geometries.
*   **`boundary_gdf`**: Boundary geometries.
*   **`snapped_col`** *(keyword-only)*: Output column for the projected `Point` geometries. Default `'snapped_geometry'`.
*   **`distance_col`** *(keyword-only)*: If provided, also stores the distance from each point to its snapped location (meters).
*   **`unit_crs`** *(keyword-only)*: EPSG code for metric calculation. Default `3857`.

### 13. `segment_boundary`
Splits boundary `LineString`s into consecutive equal-length segments (at most `segment_length` meters each), assigning a unique `segment_id`. The standard prep step for boundary-segment fixed effects.
```python
def segment_boundary(
    boundary_gdf: gpd.GeoDataFrame,
    segment_length: float,
    *,
    segment_id_col: str = 'segment_id',
    unit_crs: int = 3857
) -> gpd.GeoDataFrame
```
*   **`boundary_gdf`**: `LineString` boundary geometries.
*   **`segment_length`**: Target maximum segment length (meters).
*   **`segment_id_col`** *(keyword-only)*: Output column for the unique segment id. Default `'segment_id'`.
*   **`unit_crs`** *(keyword-only)*: EPSG code used for length-based splitting. Default `3857`. Original (non-geometry) attributes of each source line are carried over to its segments.

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

### 4. RDD Distances, Bandwidth Filtering & Placebo Tests
```python
import geopandas as gpd
from geoRDDprep import (
    calculate_signed_distance, 
    filter_by_boundary_distance, 
    shift_boundary_placebo
)

# 1. Calculate absolute and signed distances to a school boundary
# Points inside the treatment area polygon receive a positive distance; points outside are negative.
prepared_points = calculate_signed_distance(
    points_gdf=points,
    boundary_gdf=clean_lines,
    treatment_gdf=treatment_area,
    distance_col='dist_to_border',
    signed_distance_col='running_var',
    treatment_col='is_treated'
)

# 2. Select a subset of points within a 1000m RDD bandwidth for local estimation
rdd_sample = filter_by_boundary_distance(
    points_gdf=prepared_points,
    boundary_gdf=clean_lines,
    max_distance=1000.0
)

# 3. Create a placebo boundary shifted by 500 meters north
placebo_boundary = shift_boundary_placebo(clean_lines, xoff=0.0, yoff=500.0)
```

### 5. Boundary-Segment Fixed Effects
Split the border into local segments and tag each observation with the segment it sits closest to (Keele & Titiunik, 2015). You can also snap points onto the border for the "boundary point" approach.
```python
from geoRDDprep import segment_boundary, assign_nearest_boundary, snap_points_to_boundary

# 1. Cut the boundary into ~1km segments, each with a unique segment_id
segments = segment_boundary(clean_lines, segment_length=1000.0)

# 2. Assign each point to its nearest segment (for segment fixed effects)
points_fe = assign_nearest_boundary(points, segments, id_col='segment_id')
print(points_fe[['segment_id', 'boundary_distance']].head())

# 3. (Optional) Project points onto the border for boundary-point analysis
snapped = snap_points_to_boundary(points, clean_lines, distance_col='dist_to_border')
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
