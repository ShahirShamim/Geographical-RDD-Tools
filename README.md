# geoRDDprep

**geoRDDprep** is a Python package designed to streamline the data preparation process for **Geographical Regression Discontinuity Design (GeoRDD)**. It provides efficient tools for spatial joins, polygon-to-line conversions, and implementing the Turner et al. (2014) algorithm for assigning points to boundaries.

## Features

*   **`points_in_polygon`**: Efficiently assign points to polygons (e.g., addresses to school districts).
*   **`turner`**: Assign points to LineStrings based on orthogonal distance criteria (Turner et al., 2014).
*   **`poly_to_line`**: Convert Polygon geometries to LineStrings for boundary analysis.
*   **`drop_tiny_lines`**: Filter out small, noisy line segments to improve analysis quality.
*   **`remove_sliver`**: Clean up sliver polygons using Voronoi diagrams.
*   **`remove_overlaps`**: Remove overlapping segments between line datasets.

## Installation

You can install the package directly from the source:

```bash
pip install .
```

Or, if you are developing:

```bash
pip install -e .
```

## Usage

### 1. Assign Points to Polygons

Map addresses or other points to their respective administrative regions.

```python
import geopandas as gpd
from geoRDDprep import points_in_polygon

# Load your data
points = gpd.read_file("addresses.geojson")
districts = gpd.read_file("districts.geojson")

# Assign points to districts
# The resulting GeoDataFrame will have columns from 'districts' suffixed with '_district'
result = points_in_polygon(points, districts, suffix_name="_district")
```

### 2. Prepare Boundaries (Polygons to Lines)

Convert polygon boundaries into lines for distance analysis.

```python
from geoRDDprep import poly_to_line, drop_tiny_lines

# Convert polygons to lines
lines = poly_to_line(districts)

# Clean up noise by dropping very short lines (e.g., < 500 meters)
clean_lines = drop_tiny_lines(lines, method='length', meters=500)
```

### 3. Turner Algorithm

Assign points to boundaries based on distance and orthogonality.

```python
from geoRDDprep import turner

# Match points to the nearest boundary within 15 meters
# 'turner_pass' column will be True if the point satisfies the criteria
matched_data = turner(points, clean_lines, orth_distance=15)
```

## Requirements

*   `geopandas`
*   `shapely`
*   `numpy`
*   `pandas`
*   `scipy`

## License

[MIT License](LICENSE)
