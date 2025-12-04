

---

## points_in_polygon Function

This function performs a spatial join operation to find points that fall within a polygon.

### Parameters:

- `pts`: GeoDataFrame
  - Input GeoDataFrame containing Point geometries.
  
- `cha`: GeoDataFrame
  - Input GeoDataFrame containing Polygon geometries.
  
- `cha_name`: str
  - A suffix to distinguish columns from the `cha` GeoDataFrame after the spatial join.

### Returns:

- `joined`: GeoDataFrame
  - GeoDataFrame containing the result of the spatial join operation.

### Operation:

- The function performs a spatial join between the `pts` and `cha` GeoDataFrames.
- It checks whether each point falls within any polygon in the `cha` GeoDataFrame.
- The spatial join is performed using the 'within' predicate, meaning only points falling completely within polygons are considered.
- The result is a GeoDataFrame (`joined`) containing the input points along with additional columns from the `cha` GeoDataFrame.

### Example Usage:

```python
result = points_in_polygon(points_data, polygons_data, "polygon_data")
```

This example performs a spatial join to find points within polygons in `polygons_data` GeoDataFrame and adds the suffix "_polygon_data" to the columns from the `polygons_data` GeoDataFrame in the result.

---



## turner Function

This function performs a geometric operation to determine whether points lie inside LineString geometries in a GeoDataFrame.

### Parameters

- `gp_pts`: GeoDataFrame
  - Input GeoDataFrame containing Point geometries.
  
- `gp_bds`: GeoDataFrame
  - Input GeoDataFrame containing LineString geometries.
  
- `**kwargs`: Additional keyword arguments
  - `orth_distance`: int, optional (default: 15)
    - Distance in meters used to calculate orthogonal points from the LineString.
  - `reduced`: bool, optional (default: True)
    - If set to `True`, the output GeoDataFrame will contain reduced columns.

### Returns

- `gdf`: GeoDataFrame
  - GeoDataFrame containing the result of the geometric operation.

### Operation

1. Filtering input GeoDataFrames to only contain Point (`'Point'`) and LineString (`'LineString'`) geometries.
2. Finding the nearest LineString geometry to each Point geometry and calculating the distance between them.
3. Calculating the shortest line between each Point and its corresponding LineString geometry.
4. Finding orthogonal points from the LineString, at a specified distance from the nearest Point.
5. Determining if both orthogonal points lie inside the LineString geometry.
6. Optionally reducing the columns in the output GeoDataFrame.
7. Converting the output GeoDataFrame back to the original coordinate reference system.

### Notes

- The function relies on GeoPandas (`gpd`) library for geospatial operations.
- It assumes the input GeoDataFrames (`gp_pts` and `gp_bds`) are in the same coordinate reference system.
- Orthogonal points are calculated based on a fixed distance (`orth_distance`) from the LineString.
- The function returns a GeoDataFrame (`gdf`) with the original Point geometries and additional columns indicating the result of the operation.

### Example Usage

```python
result = turner(points_data, lines_data, orth_distance=20, reduced=False)
```

This example performs the 'turner' operation on `points_data` GeoDataFrame with `lines_data` GeoDataFrame using an orthogonal distance of 20 meters and keeps all columns in the output GeoDataFrame.

---



## drop_tiny_lines Function

This function is designed to filter out small LineString geometries from a GeoDataFrame based on specified criteria.

### Parameters

- `bds`: GeoDataFrame
  - Input GeoDataFrame containing LineString geometries.
  
- `method`: str, optional (default: 'percentile')
  - Method to determine the threshold for filtering out small LineString geometries. Available options are:
    - 'percentile': Filter out geometries below a certain percentile of length.
    - 'number_of_std': Filter out geometries below a certain number of standard deviations from the mean length.
    - 'length': Filter out geometries below a specified length threshold in meters.
    
- `**kwargs`: Additional keyword arguments
  - Parameters specific to each method:
    - For 'percentile' method:
      - `percentile`: float, optional (default: 0.01)
        - Percentile threshold for filtering geometries based on length.
    - For 'number_of_std' method:
      - `num_dev`: float, optional (default: 2)
        - Number of standard deviations from the mean length to use as the threshold.
    - For 'length' method:
      - `meters`: int, optional (default: 500)
        - Length threshold in meters.
        
### Returns

- `df_n`: GeoDataFrame
  - Filtered GeoDataFrame containing LineString geometries that meet the specified criteria.

### Notes

- The function first converts the input GeoDataFrame to the WGS 84 coordinate reference system (EPSG:4326) to ensure consistent length calculations.
- It then calculates the length of each LineString geometry and filters out geometries based on the specified method and threshold.
- If `reduced` is set to `True`, the function removes the 'length' column from the output GeoDataFrame.
- Finally, the filtered GeoDataFrame is converted back to its original coordinate reference system and returned.

### Example Usage

```python
filtered_data = drop_tiny_lines(input_data, method='percentile', percentile=0.05)
```

This example filters out LineString geometries below the 5th percentile of length from the `input_data` GeoDataFrame.

---

## remove_sliver Function

This function removes small, sliver-like polygons from a GeoDataFrame by creating a Voronoi diagram based on the centroids of the input polygons and then performing geometric operations to clip and refine the boundaries.

### Parameters:

- `polygons`: GeoDataFrame
  - Input GeoDataFrame containing Polygon geometries.
  
- `bound`: GeoDataFrame
  - GeoDataFrame containing a boundary for clipping the polygons.
  
### Returns:

- `final_output`: GeoDataFrame
  - GeoDataFrame containing the refined polygons with slivers removed.

### Operation:

1. Converts the `polygons` and `bound` GeoDataFrames to a common CRS (EPSG:4326) for processing.
2. Computes centroids for each polygon and creates a convex hull and Voronoi diagram based on these centroids.
3. Performs geometric operations to create a bounding box or envelope and then clips the Voronoi polygons with this geometry.
4. Finds the nearest polygon ID for each centroid and updates the Voronoi polygons accordingly.
5. Computes the difference between the clipped Voronoi polygons and the original polygons.
6. Merges the original polygons with the refined polygons and performs intersection with the clipping boundary.
7. Returns the final GeoDataFrame with the refined polygons.

### Example Usage:

```python
refined_polygons = remove_sliver(polygons_data, boundary_data)
```

This example removes sliver-like polygons from the `polygons_data` GeoDataFrame based on the `boundary_data` GeoDataFrame.

---

## poly_to_line Function

This function converts Polygon and MultiPolygon geometries in a GeoDataFrame into LineString geometries.

### Parameters:

- `GDF_polygon`: GeoDataFrame
  - Input GeoDataFrame containing Polygon and MultiPolygon geometries.

### Returns:

- `output`: GeoDataFrame
  - GeoDataFrame containing the converted LineString geometries.

### Operation:

1. Converts the input GeoDataFrame to a common CRS (EPSG:4326) for consistent processing.
2. Defines a helper function to convert Polygon or MultiPolygon geometries to LineString geometries.
3. Applies the conversion function to each geometry in the GeoDataFrame.
4. Explodes any MultiLineString geometries into individual LineStrings and resets the index.
5. Converts the GeoDataFrame back to the original CRS.

### Example Usage:

```python
lines_data = poly_to_line(polygon_data)
```

This example converts Polygon and MultiPolygon geometries in the `polygon_data` GeoDataFrame to LineString geometries.

---

## remove_overlaps Function

This function removes overlapping segments between two GeoDataFrames containing LineString geometries.

### Parameters:

- `df1`: GeoDataFrame
  - Input GeoDataFrame containing LineString geometries.
  
- `df2`: GeoDataFrame
  - Input GeoDataFrame containing LineString geometries to check for overlaps.

### Returns:

- `gdf1_cleaned`: GeoDataFrame
  - GeoDataFrame containing LineString geometries from `df1` with overlapping segments removed.

### Operation:

1. Iterates through each `LineString` in `df1` and compares it with each `LineString` in `df2`.
2. Subtracts overlapping segments from each `LineString` in `df1`.
3. Ensures that the resulting geometries are all `LineString` objects.
4. Creates a new GeoDataFrame with the cleaned `LineString` geometries.

### Example Usage:

```python
cleaned_lines = remove_overlaps(lines_data1, lines_data2)
```

This example removes overlapping segments from `lines_data1` GeoDataFrame using `lines_data2` GeoDataFrame.


# `remove_overlaps` Function Documentation

This function resolves overlapping polygons within a GeoDataFrame by creating non-overlapping geometries.

---

## Parameters:

- `gdf`: GeoDataFrame  
  Input GeoDataFrame containing Polygon geometries.

- `id_field`: str  
  The name of the column in the GeoDataFrame that contains unique identifiers for the polygons.

---

## Returns:

- `result`: GeoDataFrame  
  A GeoDataFrame with non-overlapping polygons.  
  Each polygon retains the original identifier from the `id_field`.

---

## Operation:

1. Identifies overlapping areas among the polygons in the input GeoDataFrame.
2. Resolves these overlaps by splitting the polygons into non-overlapping regions.
3. Each resulting region is assigned to the polygon with the highest priority based on its `id_field`.

---

## Example Usage:

```python
# Removing overlaps in a GeoDataFrame of administrative boundaries
cleaned_gdf = remove_overlaps(admin_boundaries, id_field="region_id")

