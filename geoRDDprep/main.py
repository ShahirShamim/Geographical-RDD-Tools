import warnings
import geopandas as gpd
from shapely.geometry import LineString, Point, Polygon, MultiPolygon
from shapely.ops import nearest_points, linemerge, substring
import shapely
from scipy.spatial import Voronoi, ConvexHull
import numpy as np
import pandas as pd
from typing import Optional, Union, List, Tuple

def _align_crs(gdf1: gpd.GeoDataFrame, gdf2: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Helper to check and align CRS between two GeoDataFrames."""
    if gdf1.crs is None or gdf2.crs is None:
        return gdf2
    if gdf1.crs != gdf2.crs:
        warnings.warn(
            f"CRS mismatch: reprojecting second GeoDataFrame from {gdf2.crs} to {gdf1.crs}.",
            UserWarning
        )
        return gdf2.to_crs(gdf1.crs)
    return gdf2

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
    polygons_gdf = _align_crs(points_gdf, polygons_gdf)
    
    # Perform spatial join with predicate='within'
    joined = gpd.sjoin(
        points_gdf, 
        polygons_gdf, 
        how="left", 
        lsuffix='', 
        rsuffix=suffix_name, 
        predicate="within"
    )
    
    # Clean up double underscores in columns if suffix starts with an underscore
    if suffix_name.startswith('_'):
        double_suffix = f"_{suffix_name}"
        rename_dict = {}
        for col in joined.columns:
            if col.endswith(double_suffix):
                cleaned_name = col[:-len(double_suffix)] + suffix_name
                rename_dict[col] = cleaned_name
        if rename_dict:
            joined = joined.rename(columns=rename_dict)
            
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
    output = polygon_gdf.copy()
    output['geometry'] = output.boundary
    output = output.explode(index_parts=False).reset_index(drop=True)
    return output

def turner(
    points_gdf: gpd.GeoDataFrame, 
    boundaries_gdf: gpd.GeoDataFrame, 
    *,
    orth_distance: float = 15.0,
    reduced: bool = True,
    unit_crs: int = 3857
) -> gpd.GeoDataFrame:
    """
    Matches points to LineString boundaries based on the criteria provided in 
    Landuse Regulation and Welfare, Turner et al (2014).
    
    Checks if points are within a certain orthogonal distance from a line segment.

    Args:
        points_gdf (gpd.GeoDataFrame): GeoDataFrame containing Point geometries.
        boundaries_gdf (gpd.GeoDataFrame): GeoDataFrame containing LineString geometries.
        orth_distance (float): Orthogonal distance threshold in meters. Default is 15.
        reduced (bool): If True, returns only original columns plus result. Default is True.
        unit_crs (int): EPSG code for metric distance calculation. Default is 3857 (Web Mercator).

    Returns:
        gpd.GeoDataFrame: GeoDataFrame with points and a 'turner_pass' boolean column.
    """
    # Align CRS
    boundaries_gdf = _align_crs(points_gdf, boundaries_gdf)

    # Filter for valid geometries
    pts = points_gdf[points_gdf['geometry'].geom_type == 'Point'].copy()
    bds = boundaries_gdf[boundaries_gdf['geometry'].geom_type == 'LineString'].copy()
    
    orig_crs = points_gdf.crs
    
    # Convert to metric CRS only if CRS is defined
    if orig_crs:
        pts = pts.to_crs(unit_crs)
        bds = bds.to_crs(unit_crs)
    
    # Spatial join + merge ONLY the geometry column of boundaries to avoid name clashes
    df_n = gpd.sjoin_nearest(pts, bds).merge(bds[['geometry']], left_on="index_right", right_index=True)
    
    # Extract coordinates in a vectorized manner
    shortest_line = gpd.GeoSeries(df_n["geometry_x"]).shortest_line(gpd.GeoSeries(df_n["geometry_y"]), align=True)
    coords = shapely.get_coordinates(shortest_line)
    P_coords = coords[0::2]
    B_coords = coords[1::2]
    
    dx = B_coords[:, 0] - P_coords[:, 0]
    dy = B_coords[:, 1] - P_coords[:, 1]
    L = np.hypot(dx, dy)
    
    # Vectorized orthogonal offset calculations (avoiding trig functions & division-by-zero)
    valid = L > 1e-9
    op1_x = np.full_like(L, np.nan)
    op1_y = np.full_like(L, np.nan)
    op2_x = np.full_like(L, np.nan)
    op2_y = np.full_like(L, np.nan)
    
    # Perpendicular unit vector (-dy/L, dx/L) applied at boundary points B
    op1_x[valid] = B_coords[valid, 0] - orth_distance * (dy[valid] / L[valid])
    op1_y[valid] = B_coords[valid, 1] + orth_distance * (dx[valid] / L[valid])
    op2_x[valid] = B_coords[valid, 0] + orth_distance * (dy[valid] / L[valid])
    op2_y[valid] = B_coords[valid, 1] - orth_distance * (dx[valid] / L[valid])
    
    op1 = shapely.points(op1_x, op1_y)
    op2 = shapely.points(op2_x, op2_y)
    
    # Element-wise distance check
    boundary_geoms = df_n["geometry_y"].values
    pass1 = np.zeros(len(df_n), dtype=bool)
    pass2 = np.zeros(len(df_n), dtype=bool)
    pass1[valid] = shapely.distance(boundary_geoms[valid], op1[valid]) < 1e-6
    pass2[valid] = shapely.distance(boundary_geoms[valid], op2[valid]) < 1e-6
    
    df_n['turner_pass'] = pass1 & pass2
    
    if reduced:
        df_n = df_n.drop(columns=['geometry_y'])
        
    gdf = gpd.GeoDataFrame(df_n, geometry='geometry_x', crs=unit_crs if orig_crs else None)
    
    if orig_crs:
        gdf = gdf.to_crs(orig_crs)
        
    return gdf

def drop_tiny_lines(
    boundaries_gdf: gpd.GeoDataFrame, 
    method: str = 'percentile', 
    *,
    percentile: float = 0.01,
    num_dev: float = 2.0,
    meters: float = 500.0,
    reduced: bool = True,
    unit_crs: int = 3857
) -> gpd.GeoDataFrame:
    """
    Filters out small LineString geometries to reduce noise.

    Args:
        boundaries_gdf (gpd.GeoDataFrame): GeoDataFrame containing LineStrings.
        method (str): Method to determine threshold ('percentile', 'number_of_std', 'length').
        percentile (float): Percentile threshold (0-1). Default 0.01.
        num_dev (float): Number of standard deviations below mean. Default 2.
        meters (float): Length threshold in meters. Default 500.
        reduced (bool): If True, removes the temporary 'length' column. Default True.
        unit_crs (int): EPSG code for metric calculation. Default 3857.

    Returns:
        gpd.GeoDataFrame: Filtered GeoDataFrame.
    """
    orig_crs = boundaries_gdf.crs
    
    # Filter for LineStrings
    bds = boundaries_gdf[boundaries_gdf['geometry'].geom_type == 'LineString'].copy()
    if bds.empty:
        return bds
        
    # Calculate lengths in the metric CRS
    if orig_crs is None or orig_crs.is_projected:
        lengths = bds.geometry.length
    else:
        lengths = bds.to_crs(unit_crs).geometry.length

    if method == 'percentile':
        cut_off = lengths.quantile(percentile)
        mask = lengths >= cut_off

    elif method == 'number_of_std':
        cut_off = lengths.mean() - num_dev * lengths.std()
        mask = lengths >= cut_off

    elif method == 'length':
        cut_off = meters
        mask = lengths >= cut_off
    else:
        mask = np.ones(len(bds), dtype=bool)

    # Filter original boundaries directly to avoid re-projection overhead
    filtered_bds = bds[mask].copy()
    
    if not reduced:
        filtered_bds['length'] = lengths[mask]
        
    return filtered_bds

def remove_sliver(
    polygons_gdf: gpd.GeoDataFrame, 
    boundary_gdf: gpd.GeoDataFrame,
    *,
    id_col: Optional[str] = None
) -> gpd.GeoDataFrame:
    """
    Removes sliver polygons by merging them into neighbors using Voronoi diagrams.

    Args:
        polygons_gdf (gpd.GeoDataFrame): Input polygons to clean.
        boundary_gdf (gpd.GeoDataFrame): Boundary to clip the result.
        id_col (str, optional): Name of the column containing unique identifiers.
                                If None, defaults to 'id' or the index name.

    Returns:
        gpd.GeoDataFrame: Cleaned polygons.
    """
    # Align CRS
    boundary_gdf = _align_crs(polygons_gdf, boundary_gdf)

    if id_col is None:
        if 'id' in polygons_gdf.columns:
            id_col = 'id'
        else:
            id_col = polygons_gdf.index.name or 'index'

    # If id_col is not in columns, reset index to make it a column
    polygons_temp = polygons_gdf.copy()
    if id_col not in polygons_temp.columns:
        polygons_temp = polygons_temp.reset_index()
        if polygons_gdf.index.name and polygons_gdf.index.name in polygons_temp.columns:
            id_col = polygons_gdf.index.name
        else:
            id_col = 'index'

    # Rename id_col to 'id' internally for consistent processing
    polygons_temp = polygons_temp.rename(columns={id_col: 'id'})

    orig_crs = polygons_temp.crs
    if orig_crs is None:
        area = polygons_temp
        clip = boundary_gdf.copy()
    elif orig_crs.is_projected:
        area = polygons_temp
        clip = boundary_gdf.copy()
    else:
        area = polygons_temp.to_crs(3857)
        clip = boundary_gdf.to_crs(3857)

    # De-duplicate centroids preserving original order
    centroids = area.geometry.centroid
    coords = np.column_stack([centroids.x, centroids.y])
    _, unique_indices = np.unique(coords, axis=0, return_index=True)
    unique_indices = np.sort(unique_indices)
    
    area_unique = area.iloc[unique_indices].copy()
    num_unique = len(area_unique)
    
    if num_unique < 3:
        return polygons_gdf.copy()

    points = coords[unique_indices]
    
    hull = ConvexHull(points)
    hull_points = points[hull.vertices]
    
    min_x, min_y = hull_points.min(axis=0)
    max_x, max_y = hull_points.max(axis=0)
    padding = max(max_x - min_x, max_y - min_y, 1.0) * 10
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
    # vor.point_region aligns with all_points
    for point_index in range(num_unique):
        region_index = vor.point_region[point_index]
        region = vor.regions[region_index]
        if -1 not in region and len(region) > 0:
            polygon = Polygon([vor.vertices[i] for i in region])
            if polygon.is_valid:
                # Direct ID assignment from the aligned unique area row
                voronoi_polygons.append({
                    'geometry': polygon, 
                    'id': area_unique.iloc[point_index]['id']
                })
                    
    voronoi_gdf = gpd.GeoDataFrame(voronoi_polygons, crs=area.crs)
    
    merged_geometry = area.union_all()
    external_boundary = merged_geometry.convex_hull
    external_boundary_gdf = gpd.GeoDataFrame([{'geometry': external_boundary}], crs=area.crs)
    
    # Clip to external boundary
    clipped_voronoi = gpd.overlay(voronoi_gdf, external_boundary_gdf, how='intersection')

    # Vectorized gap/sliver calculation
    gaps = external_boundary.difference(merged_geometry)
    gaps_gdf = gpd.GeoDataFrame([{'geometry': gaps}], crs=area.crs)
    difference_gdf = gpd.overlay(clipped_voronoi, gaps_gdf, how='intersection')
    
    # Merge difference back to the original area
    difference_gdf = difference_gdf.rename(columns={'geometry': 'geometry_difference'})
    merged_gdf = area.merge(difference_gdf[['id', 'geometry_difference']], on='id', how='left')
    
    # Vectorized union
    geom_area = merged_gdf['geometry'].values
    geom_diff = merged_gdf['geometry_difference'].values
    has_diff = pd.notna(geom_diff)
    
    combined_geoms = np.where(has_diff, shapely.union(geom_area, geom_diff), geom_area)
    combined_gdf = gpd.GeoDataFrame({'id': merged_gdf['id']}, geometry=combined_geoms, crs=area.crs)
    
    final_output = gpd.overlay(combined_gdf, clip, how='intersection')
    if 'id_2' in final_output.columns:
        final_output = final_output.drop('id_2', axis=1)
        
    # Rename 'id' back to id_col
    final_output = final_output.rename(columns={'id': id_col})

    # Restore index if it was converted from index
    if id_col == 'index' and 'index' in final_output.columns:
        final_output = final_output.set_index('index')
        final_output.index.name = polygons_gdf.index.name
    elif id_col == polygons_gdf.index.name and id_col in final_output.columns:
        final_output = final_output.set_index(id_col)
        
    if orig_crs is not None and not orig_crs.is_projected:
        final_output = final_output.to_crs(orig_crs)
        
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
    # Align CRS
    df2 = _align_crs(df1, df2)

    if df1.empty or df2.empty:
        return df1.copy()
        
    union_df2 = df2.union_all()
    diff = df1.geometry.difference(union_df2)
    
    # Explode to convert MultiLineStrings to LineStrings and filter out empty geometries
    exploded = diff.explode(index_parts=False)
    exploded = exploded[~exploded.is_empty]
    
    return gpd.GeoDataFrame(geometry=exploded, crs=df1.crs).reset_index(drop=True)

def calculate_signed_distance(
    points_gdf: gpd.GeoDataFrame,
    boundary_gdf: gpd.GeoDataFrame,
    treatment_gdf: gpd.GeoDataFrame,
    *,
    distance_col: str = 'distance',
    signed_distance_col: str = 'signed_distance',
    treatment_col: str = 'is_treated',
    unit_crs: int = 3857
) -> gpd.GeoDataFrame:
    """
    Calculates the shortest distance and signed distance from points to a boundary.
    Points inside the treatment area (treatment_gdf) receive a positive distance,
    while points outside (control area) receive a negative distance.
    
    Args:
        points_gdf (gpd.GeoDataFrame): GeoDataFrame containing the Point geometries.
        boundary_gdf (gpd.GeoDataFrame): GeoDataFrame containing the boundary LineString geometries.
        treatment_gdf (gpd.GeoDataFrame): GeoDataFrame containing the treatment polygon(s).
        distance_col (str): Name of the column for absolute distance. Default is 'distance'.
        signed_distance_col (str): Name of the column for signed distance. Default is 'signed_distance'.
        treatment_col (str): Name of the column for treatment indicator. Default is 'is_treated'.
        unit_crs (int): EPSG code for metric distance calculation. Default is 3857.
        
    Returns:
        gpd.GeoDataFrame: A copy of points_gdf with distance, signed distance, and treatment indicator columns added.
    """
    if boundary_gdf.empty:
        raise ValueError("boundary_gdf cannot be empty.")
    if treatment_gdf.empty:
        raise ValueError("treatment_gdf cannot be empty.")

    # Align CRS
    boundary_gdf = _align_crs(points_gdf, boundary_gdf)
    treatment_gdf = _align_crs(points_gdf, treatment_gdf)

    orig_crs = points_gdf.crs
    pts = points_gdf.copy()
    bds = boundary_gdf.copy()

    if orig_crs:
        pts = pts.to_crs(unit_crs)
        bds = bds.to_crs(unit_crs)
        treatment_proj = treatment_gdf.to_crs(unit_crs)
    else:
        treatment_proj = treatment_gdf

    # Calculate distance using sjoin_nearest
    joined = gpd.sjoin_nearest(pts, bds, how='left', distance_col=distance_col)
    joined = joined[~joined.index.duplicated(keep='first')]

    # Calculate treatment status
    treatment_union = treatment_proj.union_all()
    is_treated = pts.geometry.within(treatment_union)

    # Compile results
    res = points_gdf.copy()
    res[treatment_col] = is_treated
    res[distance_col] = joined[distance_col]
    res[signed_distance_col] = np.where(res[treatment_col], res[distance_col], -res[distance_col])

    return res

def extract_shared_boundaries(
    gdf: gpd.GeoDataFrame,
    *,
    id_col: Optional[str] = None
) -> gpd.GeoDataFrame:
    """
    Extracts the shared boundaries (borders) between adjacent polygons in a GeoDataFrame.
    
    Args:
        gdf (gpd.GeoDataFrame): GeoDataFrame containing Polygon or MultiPolygon geometries.
        id_col (str, optional): Column name to identify polygons. If None, uses index.
        
    Returns:
        gpd.GeoDataFrame: A GeoDataFrame where each row is a shared boundary segment (LineString or MultiLineString)
                          between two adjacent polygons, with columns identifying the two neighbors.
    """
    if gdf.empty:
        return gpd.GeoDataFrame(columns=['left_id', 'right_id', 'geometry'], crs=gdf.crs)

    # Spatial join to find touching polygons
    adj = gpd.sjoin(gdf, gdf, predicate='touches')
    
    # Filter to unique pairs where left index < right index
    left_index = adj.index
    right_index = adj['index_right']
    
    mask = left_index < right_index
    adj_unique = adj[mask]

    if adj_unique.empty:
        return gpd.GeoDataFrame(columns=['left_id', 'right_id', 'geometry'], crs=gdf.crs)

    geom_left = gdf.loc[adj_unique.index, 'geometry'].values
    geom_right = gdf.loc[adj_unique['index_right'], 'geometry'].values
    
    # Vectorized intersection
    shared_geoms = shapely.intersection(geom_left, geom_right)
    
    # Keep only linear components (dimension 1: LineString, MultiLineString)
    is_linear = ~shapely.is_empty(shared_geoms) & (shapely.get_dimensions(shared_geoms) == 1)

    if id_col is not None:
        left_ids = gdf.loc[adj_unique.index, id_col].values
        right_ids = gdf.loc[adj_unique['index_right'], id_col].values
    else:
        left_ids = adj_unique.index.values
        right_ids = adj_unique['index_right'].values

    result_gdf = gpd.GeoDataFrame({
        'left_id': left_ids[is_linear],
        'right_id': right_ids[is_linear],
        'geometry': shared_geoms[is_linear]
    }, crs=gdf.crs)

    return result_gdf.reset_index(drop=True)

def shift_boundary_placebo(
    boundary_gdf: gpd.GeoDataFrame,
    xoff: float = 0.0,
    yoff: float = 0.0,
    *,
    unit_crs: int = 3857
) -> gpd.GeoDataFrame:
    """
    Shifts/translates a boundary by a specified offset in meters.
    Useful for creating placebo boundaries in GeoRDD analysis.
    
    Args:
        boundary_gdf (gpd.GeoDataFrame): GeoDataFrame containing the boundary geometries.
        xoff (float): Offset in the X direction (east-west) in meters.
        yoff (float): Offset in the Y direction (north-south) in meters.
        unit_crs (int): EPSG code for metric translation. Default is 3857.
        
    Returns:
        gpd.GeoDataFrame: A new GeoDataFrame with the translated boundary in the original CRS.
    """
    orig_crs = boundary_gdf.crs
    gdf_temp = boundary_gdf.copy()
    if orig_crs:
        gdf_temp = gdf_temp.to_crs(unit_crs)
        
    gdf_temp['geometry'] = gdf_temp['geometry'].translate(xoff=xoff, yoff=yoff)
    
    if orig_crs:
        gdf_temp = gdf_temp.to_crs(orig_crs)
        
    return gdf_temp

def filter_by_boundary_distance(
    points_gdf: gpd.GeoDataFrame,
    boundary_gdf: gpd.GeoDataFrame,
    max_distance: float,
    *,
    unit_crs: int = 3857
) -> gpd.GeoDataFrame:
    """
    Filters points to only those within a specified maximum distance (bandwidth) from the boundary.
    
    Args:
        points_gdf (gpd.GeoDataFrame): GeoDataFrame containing the Point geometries.
        boundary_gdf (gpd.GeoDataFrame): GeoDataFrame containing the boundary geometries.
        max_distance (float): Maximum distance threshold in meters.
        unit_crs (int): EPSG code for metric distance calculation. Default is 3857.
        
    Returns:
        gpd.GeoDataFrame: A filtered GeoDataFrame containing only the points within the distance threshold.
    """
    if boundary_gdf.empty:
        return points_gdf.copy().iloc[0:0]

    # Align CRS
    boundary_gdf = _align_crs(points_gdf, boundary_gdf)

    orig_crs = points_gdf.crs
    pts = points_gdf.copy()
    bds = boundary_gdf.copy()
    if orig_crs:
        pts = pts.to_crs(unit_crs)
        bds = bds.to_crs(unit_crs)
        
    joined = gpd.sjoin_nearest(pts, bds, how='left', distance_col='temp_dist')
    joined = joined[~joined.index.duplicated(keep='first')]

    mask = (joined['temp_dist'] <= max_distance) & (joined['temp_dist'].notna())

    return points_gdf[mask].copy()

def assign_nearest_boundary(
    points_gdf: gpd.GeoDataFrame,
    boundary_gdf: gpd.GeoDataFrame,
    *,
    id_col: Optional[str] = None,
    boundary_id_col: str = 'boundary_id',
    distance_col: str = 'boundary_distance',
    unit_crs: int = 3857
) -> gpd.GeoDataFrame:
    """
    Assigns each point to its nearest boundary feature, attaching the boundary's
    identifier and the metric distance to it.

    This is the building block for boundary-segment fixed effects in GeoRDD: when
    there are multiple boundary segments (e.g. the output of ``segment_boundary``),
    each observation can be tagged with the segment it is closest to.

    Args:
        points_gdf (gpd.GeoDataFrame): GeoDataFrame containing the Point geometries.
        boundary_gdf (gpd.GeoDataFrame): GeoDataFrame containing the boundary geometries.
        id_col (str, optional): Column in ``boundary_gdf`` holding the identifier to
                                attach. If None, the boundary's index is used.
        boundary_id_col (str): Name of the output column for the boundary identifier.
                               Default is 'boundary_id'.
        distance_col (str): Name of the output column for the distance to the nearest
                            boundary (in meters). Default is 'boundary_distance'.
        unit_crs (int): EPSG code for metric distance calculation. Default is 3857.

    Returns:
        gpd.GeoDataFrame: A copy of points_gdf with the boundary identifier and distance
                          columns added (in the original CRS).
    """
    if boundary_gdf.empty:
        raise ValueError("boundary_gdf cannot be empty.")

    boundary_gdf = _align_crs(points_gdf, boundary_gdf)

    orig_crs = points_gdf.crs
    pts = points_gdf.copy()
    bds = boundary_gdf.copy()
    if orig_crs:
        pts = pts.to_crs(unit_crs)
        bds = bds.to_crs(unit_crs)

    # Reset to a clean positional index so 'index_right' can be mapped reliably.
    bds_reset = bds.reset_index(drop=True)
    if id_col is not None:
        bid_values = bds_reset[id_col].values
    else:
        bid_values = boundary_gdf.index.values

    joined = gpd.sjoin_nearest(pts, bds_reset, how='left', distance_col=distance_col)
    joined = joined[~joined.index.duplicated(keep='first')]

    res = points_gdf.copy()
    res[distance_col] = joined[distance_col]

    idx_right = joined['index_right']
    mapped = pd.Series(np.nan, index=joined.index, dtype=object)
    ok = idx_right.notna()
    mapped.loc[ok] = bid_values[idx_right[ok].astype(int).values]
    res[boundary_id_col] = mapped

    return res

def snap_points_to_boundary(
    points_gdf: gpd.GeoDataFrame,
    boundary_gdf: gpd.GeoDataFrame,
    *,
    snapped_col: str = 'snapped_geometry',
    distance_col: Optional[str] = None,
    unit_crs: int = 3857
) -> gpd.GeoDataFrame:
    """
    Projects each point onto its nearest boundary, returning the snapped (projected)
    point geometry. Useful for the "boundary point" RD approach and for visualizing
    which part of the border each observation maps to.

    Args:
        points_gdf (gpd.GeoDataFrame): GeoDataFrame containing the Point geometries.
        boundary_gdf (gpd.GeoDataFrame): GeoDataFrame containing the boundary geometries.
        snapped_col (str): Name of the output column holding the projected Point
                           geometries. Default is 'snapped_geometry'.
        distance_col (str, optional): If provided, also store the distance from each
                                      point to its snapped location (in meters).
        unit_crs (int): EPSG code for metric distance calculation. Default is 3857.

    Returns:
        gpd.GeoDataFrame: A copy of points_gdf with the snapped geometry column added.
                          The active geometry remains the original points; the snapped
                          geometries are returned in the original CRS.
    """
    if boundary_gdf.empty:
        raise ValueError("boundary_gdf cannot be empty.")

    boundary_gdf = _align_crs(points_gdf, boundary_gdf)

    orig_crs = points_gdf.crs
    pts = points_gdf.copy()
    bds = boundary_gdf.copy()
    if orig_crs:
        pts = pts.to_crs(unit_crs)
        bds = bds.to_crs(unit_crs)

    bds_reset = bds.reset_index(drop=True)
    joined = gpd.sjoin_nearest(pts, bds_reset, how='left', distance_col='_snap_dist')
    joined = joined[~joined.index.duplicated(keep='first')]

    matched_geom = bds_reset.geometry.values[joined['index_right'].astype(int).values]
    point_geom = joined.geometry.values

    # The shortest line runs point -> boundary; its second coordinate is the
    # projected location on the boundary.
    lines = shapely.shortest_line(point_geom, matched_geom)
    coords = shapely.get_coordinates(lines)
    snapped_pts = shapely.points(coords[1::2])

    snapped_series = gpd.GeoSeries(
        snapped_pts, index=joined.index, crs=unit_crs if orig_crs else None
    )
    if orig_crs:
        snapped_series = snapped_series.to_crs(orig_crs)

    res = points_gdf.copy()
    res[snapped_col] = snapped_series
    if distance_col is not None:
        res[distance_col] = joined['_snap_dist']

    return res

def segment_boundary(
    boundary_gdf: gpd.GeoDataFrame,
    segment_length: float,
    *,
    segment_id_col: str = 'segment_id',
    unit_crs: int = 3857
) -> gpd.GeoDataFrame:
    """
    Splits boundary LineStrings into consecutive segments of approximately equal
    length (at most ``segment_length`` meters each). Each output segment receives a
    unique identifier, enabling boundary-segment fixed effects in GeoRDD analysis.

    Each line is divided into ``ceil(length / segment_length)`` equal pieces, so the
    pieces of a given line are all the same length and never exceed ``segment_length``.
    Original (non-geometry) attributes of each source line are carried over to its
    segments.

    Args:
        boundary_gdf (gpd.GeoDataFrame): GeoDataFrame containing LineString geometries.
        segment_length (float): Target maximum segment length in meters.
        segment_id_col (str): Name of the output column for the unique segment id.
                              Default is 'segment_id'.
        unit_crs (int): EPSG code used for length-based splitting when the input is
                        geographic. Default is 3857.

    Returns:
        gpd.GeoDataFrame: A GeoDataFrame of segment LineStrings (in the original CRS)
                          with a unique ``segment_id_col`` column.
    """
    if segment_length <= 0:
        raise ValueError("segment_length must be a positive number.")

    orig_crs = boundary_gdf.crs
    bds = boundary_gdf[boundary_gdf.geometry.geom_type == 'LineString'].copy()
    if bds.empty:
        return gpd.GeoDataFrame(columns=[segment_id_col, 'geometry'], crs=orig_crs)

    # Split using metric lengths; reproject only when the CRS is geographic.
    if orig_crs is not None and not orig_crs.is_projected:
        work = bds.to_crs(unit_crs)
    else:
        work = bds

    work_reset = work.reset_index(drop=True)
    attr_cols = [c for c in work_reset.columns if c != 'geometry']

    records = []
    for pos, geom in enumerate(work_reset.geometry.values):
        if geom is None or geom.is_empty:
            continue
        total = geom.length
        if total == 0:
            continue
        n = max(1, int(np.ceil(total / segment_length)))
        cuts = np.linspace(0.0, total, n + 1)
        base = {c: work_reset.iloc[pos][c] for c in attr_cols}
        for i in range(n):
            sub = substring(geom, cuts[i], cuts[i + 1])
            if sub.is_empty:
                continue
            rec = dict(base)
            rec['geometry'] = sub
            records.append(rec)

    seg_gdf = gpd.GeoDataFrame(records, crs=work.crs)
    seg_gdf.insert(0, segment_id_col, np.arange(len(seg_gdf)))

    if orig_crs is not None and not orig_crs.is_projected:
        seg_gdf = seg_gdf.to_crs(orig_crs)

    return seg_gdf
