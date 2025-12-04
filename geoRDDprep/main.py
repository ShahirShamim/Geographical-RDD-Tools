import geopandas as gpd
from shapely.geometry import LineString, Point, Polygon, MultiPolygon
from shapely.ops import nearest_points, linemerge
from scipy.spatial import Voronoi, ConvexHull
import numpy as np
import pandas as pd
from typing import Optional, Union, List, Tuple

def points_in_polygon(
    points_gdf: gpd.GeoDataFrame, 
    polygons_gdf: gpd.GeoDataFrame, 
    suffix_name: str
) -> gpd.GeoDataFrame:
    """
    Assigns characteristics to points based on the polygon they fall within.
    Performs a spatial join between points and polygons.

    Args:
        points_gdf (gpd.GeoDataFrame): GeoDataFrame containing the points.
        polygons_gdf (gpd.GeoDataFrame): GeoDataFrame containing the polygons with characteristics.
        suffix_name (str): Suffix to append to columns from the polygons_gdf in the result.

    Returns:
        gpd.GeoDataFrame: A GeoDataFrame with points and their assigned polygon characteristics.
    """
    # Perform spatial join with op='within'
    joined = gpd.sjoin(
        points_gdf, 
        polygons_gdf, 
        how="left", 
        lsuffix='', 
        rsuffix=suffix_name, 
        predicate="within"
    )
    return joined

def poly_to_line(polygon_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Converts Polygons and MultiPolygons in a GeoDataFrame to LineStrings.
    Useful for preparing polygon boundaries for the 'turner' function.

    Args:
        polygon_gdf (gpd.GeoDataFrame): Input GeoDataFrame containing Polygon or MultiPolygon geometries.

    Returns:
        gpd.GeoDataFrame: GeoDataFrame with geometries converted to LineStrings.
    """
    og_crs = polygon_gdf.crs
    # Work in WGS84 for consistency if needed, though mostly geometric operation
    polygon_gdf = polygon_gdf.to_crs(4326)

    def polygon_to_linestring(geometry: Union[Polygon, MultiPolygon]) -> Union[LineString, List[LineString], None]:
        if isinstance(geometry, Polygon):
            return LineString(geometry.exterior)
        elif isinstance(geometry, MultiPolygon):
            # Convert each polygon in the MultiPolygon to a LineString
            return [LineString(p.exterior) for p in geometry]
        else:
            return None

    # Apply the conversion
    polygon_gdf['geometry'] = polygon_gdf['geometry'].apply(polygon_to_linestring)
    
    # Explode the list of LineStrings (if any) into separate rows
    output = polygon_gdf.explode(index_parts=False).reset_index(drop=True)
    
    if og_crs:
        output = output.to_crs(og_crs)
        
    return output

def turner(
    points_gdf: gpd.GeoDataFrame, 
    boundaries_gdf: gpd.GeoDataFrame, 
    **kwargs
) -> gpd.GeoDataFrame:
    """
    Matches points to LineString boundaries based on the criteria provided in 
    Landuse Regulation and Welfare, Turner et al (2014).
    
    Checks if points are within a certain orthogonal distance from a line segment.

    Args:
        points_gdf (gpd.GeoDataFrame): GeoDataFrame containing Point geometries.
        boundaries_gdf (gpd.GeoDataFrame): GeoDataFrame containing LineString geometries.
        **kwargs:
            orth_distance (float): Orthogonal distance threshold in meters. Default is 15.
            reduced (bool): If True, returns only original columns plus result. Default is True.
            unit_crs (int): EPSG code for metric distance calculation. Default is 3857 (Web Mercator).

    Returns:
        gpd.GeoDataFrame: GeoDataFrame with points and a 'turner_pass' boolean column.
    """
    # Default parameters
    orth_distance = kwargs.get('orth_distance', 15)
    reduced = kwargs.get('reduced', True)
    unit_crs = kwargs.get('unit_crs', 3857)

    # Filter for valid geometries
    pts = points_gdf[points_gdf['geometry'].geom_type == 'Point'].copy()
    bds = boundaries_gdf[boundaries_gdf['geometry'].geom_type == 'LineString'].copy()
    
    orig_crs = points_gdf.crs
    
    # Convert to metric CRS for accurate distance calculation
    pts = pts.to_crs(unit_crs)
    bds = bds.to_crs(unit_crs)
    
    # Find nearest linestring and calculate distance
    # sjoin_nearest finds the nearest boundary for each point
    df_n = gpd.sjoin_nearest(pts, bds).merge(bds, left_on="index_right", right_index=True)
    
    # Calculate exact distance to the nearest boundary
    df_n["distance"] = df_n.apply(lambda r: r["geometry_x"].distance(r["geometry_y"]), axis=1)

    # Find the shortest line between point and nearest linestring
    p = gpd.GeoSeries(df_n["geometry_x"])
    l = gpd.GeoSeries(df_n["geometry_y"])
    shortest_line = p.shortest_line(l, align=True)
    df_n["shortest_line"] = shortest_line

    def orthogonal_points(line: LineString, distance_in_meters: float) -> Optional[Tuple[Point, Point]]:
        """Calculates two points orthogonal to the line segment at the given distance."""
        dist = distance_in_meters
        # The second point of the shortest line is on the boundary
        # The first point is the original point
        # We want to check orthogonality relative to the boundary? 
        # Actually, looking at original logic:
        # line is the shortest line connecting point to boundary.
        # line.coords[1] is the point ON the boundary.
        # line.coords[0] is the original point.
        
        if len(line.coords) < 2:
            return None
            
        point_on_boundary = Point(line.coords[1][0], line.coords[1][1])
        
        # Calculate the slope of the shortest line
        dx = line.coords[1][0] - line.coords[0][0]
        dy = line.coords[1][1] - line.coords[0][1]
        
        if dx == 0:  # If the line is vertical, swap dx and dy for perpendicular calculation logic below
             dx, dy = dy, dx
        
        if dx == 0: # Still 0 means points are identical
            return None
            
        slope = dy / dx

        # Calculate the angle of the line
        angle = np.arctan(slope)

        # Calculate the coordinates of the orthogonal points
        # These points are 'dist' away from the point_on_boundary, perpendicular to the shortest line?
        # Wait, the original code adds pi/2 to the angle of the shortest line.
        # This creates points perpendicular to the shortest line, centered at point_on_boundary.
        # This effectively checks if we are "within the width" of the segment if we consider the segment as a thick line?
        
        x1 = point_on_boundary.x + dist * np.cos(angle + np.pi / 2)
        y1 = point_on_boundary.y + dist * np.sin(angle + np.pi / 2)
        x2 = point_on_boundary.x + dist * np.cos(angle - np.pi / 2)
        y2 = point_on_boundary.y + dist * np.sin(angle - np.pi / 2)

        return Point(x1, y1), Point(x2, y2)

    df_n['orthogonal_points'] = df_n.apply(
        lambda row: orthogonal_points(row['shortest_line'], orth_distance), axis=1)

    # Expand the tuple of points into columns
    df_n[['orthogonal_point1', 'orthogonal_point2']] = df_n['orthogonal_points'].apply(pd.Series)

    def both_points_inside(geo: LineString, op1: Point, op2: Point) -> bool:
        """Checks if both orthogonal points lie on the LineString geometry."""
        if not isinstance(geo, LineString) or op1 is None or op2 is None:
            return False
        # Check if distance is effectively zero (point is on the line)
        if geo.distance(op1) < 1e-6 and geo.distance(op2) < 1e-6:
            return True
        else:
            return False

    df_n['turner_pass'] = df_n.apply(
        lambda row: both_points_inside(
            row['geometry_y'], 
            row['orthogonal_point1'],
            row['orthogonal_point2']
        ), axis=1
    )

    if reduced:
        # Keep only original geometry and the result
        # Assuming geometry_x is the original point geometry
        columns_to_drop = df_n.loc[:, 'geometry_y':'orthogonal_point2'].columns
        df_n = df_n.drop(columns_to_drop, axis=1)

    gdf = gpd.GeoDataFrame(df_n, geometry='geometry_x', crs=unit_crs)
    
    if orig_crs:
        gdf = gdf.to_crs(orig_crs)
        
    return gdf

def drop_tiny_lines(
    boundaries_gdf: gpd.GeoDataFrame, 
    method: str = 'percentile', 
    **kwargs
) -> gpd.GeoDataFrame:
    """
    Filters out small LineString geometries to reduce noise.

    Args:
        boundaries_gdf (gpd.GeoDataFrame): GeoDataFrame containing LineStrings.
        method (str): Method to determine threshold ('percentile', 'number_of_std', 'length').
        **kwargs:
            percentile (float): Percentile threshold (0-1). Default 0.01.
            num_dev (float): Number of standard deviations below mean. Default 2.
            meters (float): Length threshold in meters. Default 500.
            reduced (bool): If True, removes the temporary 'length' column. Default True.
            unit_crs (int): EPSG code for metric calculation. Default 3857.

    Returns:
        gpd.GeoDataFrame: Filtered GeoDataFrame.
    """
    reduced = kwargs.get('reduced', True)
    unit_crs = kwargs.get('unit_crs', 3857)
    
    orig_crs = boundaries_gdf.crs
    
    # Filter for LineStrings
    bds = boundaries_gdf[boundaries_gdf['geometry'].geom_type == 'LineString'].copy()
    
    # Convert to metric CRS for length calculation
    df = bds.to_crs(unit_crs)
    df['length'] = df.length

    if method == 'percentile':
        percentile = kwargs.get('percentile', 0.01)
        cut_off = df.length.quantile(percentile)
        df = df[df['length'] >= cut_off]

    elif method == 'number_of_std':
        num_dev = kwargs.get('num_dev', 2)
        cut_off = df.length.mean() - num_dev * df.length.std()
        df = df[df['length'] >= cut_off]

    elif method == 'length':
        meters = kwargs.get('meters', 500)
        # Note: In the original code, this was meters / 100000. 
        # Assuming the user meant meters directly if using a metric CRS like 3857.
        # If the original intent was degrees, 500/100000 = 0.005 degrees approx 500m.
        # Since we converted to unit_crs (metric), we should use meters directly.
        cut_off = meters
        df = df[df['length'] >= cut_off]

    if reduced:
        df = df.drop(columns=['length'])

    if orig_crs:
        df_n = df.to_crs(orig_crs)
    else:
        df_n = df
        
    return df_n

def remove_sliver(
    polygons_gdf: gpd.GeoDataFrame, 
    boundary_gdf: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Removes sliver polygons by merging them into neighbors using Voronoi diagrams.

    Args:
        polygons_gdf (gpd.GeoDataFrame): Input polygons to clean.
        boundary_gdf (gpd.GeoDataFrame): Boundary to clip the result.

    Returns:
        gpd.GeoDataFrame: Cleaned polygons.
    """
    # Ensure we are working in WGS84 for consistency in this algorithm
    area = polygons_gdf.to_crs(4326)
    clip = boundary_gdf.to_crs(4326)

    def find_nearest_id(point: Point, centers: gpd.GeoDataFrame) -> int:
        nearest_center = centers.geometry.distance(point).idxmin()
        return centers.loc[nearest_center, 'id']

    # 1. Create Voronoi regions from centroids
    centroids = area.geometry.centroid
    center = gpd.GeoDataFrame(area[['id']].copy(), geometry=centroids)
    
    points = np.array(list(center.geometry.apply(lambda geom: (geom.x, geom.y))))
    # Remove duplicates
    points = np.unique(points, axis=0)
    
    if len(points) < 3:
        # ConvexHull needs at least 3 points
        return polygons_gdf

    hull = ConvexHull(points)
    hull_points = points[hull.vertices]
    
    # Add dummy points to bound the Voronoi diagram
    padding = 0.1 
    min_x, min_y = hull_points.min(axis=0)
    max_x, max_y = hull_points.max(axis=0)
    min_x -= padding; max_x += padding
    min_y -= padding; max_y += padding
    
    dummy_points = np.array([
        [min_x, min_y], [min_x, max_y], [max_x, min_y], [max_x, max_y],
        [(min_x + max_x) / 2, min_y], [(min_x + max_x) / 2, max_y],
        [min_x, (min_y + max_y) / 2], [max_x, (min_y + max_y) / 2]
    ])
    all_points = np.vstack([points, dummy_points])
    
    vor = Voronoi(all_points)
    
    voronoi_polygons = []
    for point_index, region_index in enumerate(vor.point_region):
        region = vor.regions[region_index]
        if -1 not in region and len(region) > 0:
            polygon = Polygon([vor.vertices[i] for i in region])
            if polygon.is_valid:
                # Match the centroid back to its original id
                # Only consider original points, not dummy ones
                if point_index < len(center):
                    voronoi_polygons.append({'geometry': polygon, 'id': center.iloc[point_index]['id']})
                    
    voronoi_gdf = gpd.GeoDataFrame(voronoi_polygons, crs=area.crs)
    
    # 2. Clip Voronoi regions to the external boundary of the original area
    merged_geometry = area.unary_union
    external_boundary = merged_geometry.convex_hull
    external_boundary_gdf = gpd.GeoDataFrame([{'geometry': external_boundary}], crs=area.crs)
    
    # Intersect with external boundary
    clipped_voronoi = gpd.overlay(voronoi_gdf, external_boundary_gdf, how='intersection')

    # Re-assign IDs based on nearest original centroid (to ensure correctness)
    clipped_voronoi['centroid'] = clipped_voronoi.geometry.centroid
    clipped_voronoi['id'] = clipped_voronoi['centroid'].apply(lambda x: find_nearest_id(x, center))
    clipped_voronoi.drop(columns='centroid', inplace=True)

    # 3. Calculate difference between new Voronoi regions and original polygons
    # This identifies the "slivers" or expanded areas
    difference_polygons = []
    
    # This iteration might be slow for large datasets
    for index, row in clipped_voronoi.iterrows():
        poly = row['geometry']
        poly_id = row['id']
        # Subtract all original polygons from this Voronoi region
        # Optimization: Subtract only the original polygon with the SAME ID?
        # The logic seems to be: fill gaps. 
        # So we take the Voronoi region (which covers the gap) and subtract existing polygons.
        # The remainder is the gap/sliver, which we assign to this ID.
        
        # Using unary_union of area might be faster than iterating
        poly = poly.difference(merged_geometry)
        
        if poly.is_valid and not poly.is_empty:
            difference_polygons.append({'geometry': poly, 'id': poly_id})
            
    difference_gdf = gpd.GeoDataFrame(difference_polygons, crs=area.crs)
    
    # 4. Merge the difference (slivers) back into the original polygons
    merged_gdf = area.merge(difference_gdf, on='id', suffixes=('_area', '_difference'), how='left')
    
    combined_polygons = []
    for index, row in merged_gdf.iterrows():
        geom_area = row['geometry_area']
        geom_diff = row['geometry_difference'] if not pd.isna(row['geometry_difference']) else None
        
        if geom_diff:
            combined_geometry = geom_area.union(geom_diff)
        else:
            combined_geometry = geom_area
            
        combined_polygons.append({'geometry': combined_geometry, 'id': row['id']})
        
    combined_gdf = gpd.GeoDataFrame(combined_polygons, crs=area.crs)
    
    # 5. Final clip with the provided boundary
    final_output = gpd.overlay(combined_gdf, clip, how='intersection')
    
    # Cleanup columns
    if 'id_2' in final_output.columns:
        final_output = final_output.drop('id_2', axis=1)
        
    return final_output

def remove_overlaps(df1: gpd.GeoDataFrame, df2: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Removes overlapping segments from df1 that are present in df2.
    
    Args:
        df1 (gpd.GeoDataFrame): The GeoDataFrame to clean.
        df2 (gpd.GeoDataFrame): The GeoDataFrame containing geometries to subtract.

    Returns:
        gpd.GeoDataFrame: df1 with overlapping segments removed.
    """
    result_geometries = []
    
    # Spatial index could speed this up, but sticking to logic for now
    for geom1 in df1.geometry:
        non_overlapping_geom = geom1

        # Check against all geometries in df2
        # Optimization: Use spatial index to only check intersecting geometries
        possible_matches_index = list(df2.sindex.intersection(geom1.bounds))
        possible_matches = df2.iloc[possible_matches_index]

        for geom2 in possible_matches.geometry:
            if non_overlapping_geom.intersects(geom2):
                non_overlapping_geom = non_overlapping_geom.difference(geom2)
                if non_overlapping_geom.is_empty:
                    break

        if not non_overlapping_geom.is_empty:
            if isinstance(non_overlapping_geom, LineString):
                result_geometries.append(non_overlapping_geom)
            elif hasattr(non_overlapping_geom, 'geoms'): # MultiLineString or GeometryCollection
                result_geometries.extend(non_overlapping_geom.geoms)
            else:
                 result_geometries.append(non_overlapping_geom)

    exploded_gdf = gpd.GeoDataFrame(geometry=result_geometries, crs=df1.crs)
    return exploded_gdf
