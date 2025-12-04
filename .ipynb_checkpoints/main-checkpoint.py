import geopandas as gpd
from shapely.ops import nearest_points, linemerge
from shapely.geometry import LineString, Point, Polygon, MultiPolygon
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely import LineString, Point
from scipy.spatial import Voronoi, ConvexHull


# Assign characteristics to points
## first input is the shape file containing the points on the map
## second input is the shape file containing the polygons for the characteristics
## e,g assign_characteristics(addresses, school_districs)
def points_in_polygon(pts, cha, cha_name):
    # Perform spatial join with op='within'
    joined = gpd.sjoin(pts, cha, how="left", lsuffix='', rsuffix=cha_name, predicate="within")
    return joined



# Turner algorithm to match points to boundries
## Boundries must be linestrings
##

# Input orth distace based on pts.crs units of measurement
def turner(gp_pts, gp_bds, **kwargs):
	# Default parameters
	orth_distance = kwargs.get('orth_distance', 15)
	reduced = kwargs.get('reduced', True)

	# Only keeping Points and lineStrings
	pts = gp_pts[gp_pts['geometry'].geom_type == 'Point']
	bds = gp_bds[gp_bds['geometry'].geom_type == 'LineString']
	orig_crs = gp_pts.crs
	# Finding nearest linestring and distance
	gpd.options.display_precision = 10
	pts = pts.to_crs(4326)
	bds = bds.to_crs(4326)
	df_n = gpd.sjoin_nearest(pts, bds).merge(bds, left_on="index_right", right_index=True)
	df_n["distance"] = df_n.apply(lambda r: r["geometry_x"].distance(r["geometry_y"]), axis=1)

	# Shortest line between point and linestring
	p = gpd.GeoSeries(df_n["geometry_x"])
	l = gpd.GeoSeries(df_n["geometry_y"])
	shortest_line = p.shortest_line(l, align=True)
	df_n["shortest_line"] = shortest_line

	def orthogonal_points(line, distance_in_meters):
		dist = distance_in_meters / 100000
		point = Point(line.coords[1][0], line.coords[1][1])
		# Calculate the slope of the line
		dx = line.coords[1][0] - line.coords[0][0]
		dy = line.coords[1][1] - line.coords[0][1]
		if dx == 0:  # If the line is vertical, swap dx and dy
			dx, dy = dy, dx
		if dx == 0:  # Accounting for divide by 0
			return
		slope = dy / dx

		# Calculate the angle of the line
		angle = np.arctan(slope)

		# Calculate the coordinates of the orthogonal points
		x1 = point.x + dist * np.cos(angle + np.pi / 2)
		y1 = point.y + dist * np.sin(angle + np.pi / 2)
		x2 = point.x + dist * np.cos(angle - np.pi / 2)
		y2 = point.y + dist * np.sin(angle - np.pi / 2)

		return Point(x1, y1), Point(x2, y2)

	df_n['orthogonal_points'] = df_n.apply(
		lambda row: orthogonal_points(LineString(row['shortest_line']), orth_distance), axis=1)

	df_n[['orthogonal_point1', 'orthogonal_point2']] = df_n['orthogonal_points'].apply(pd.Series)

	def both_points_inside(geo, op1, op2):
		if geo.geom_type != 'LineString':
			return False
		if geo.distance(op1) < 0.000001 and geo.distance(op2) < 0.000001:
			return True
		else:
			return False

	df_n['turner_pass'] = df_n.apply(
		lambda row: both_points_inside(LineString(row['geometry_y']), Point(row['orthogonal_point1']),
									   Point(row['orthogonal_point2'])), axis=1)

	if reduced == True:
		df_n = df_n.drop(df_n.loc[:, 'geometry_y':'orthogonal_point2'].columns, axis=1)

	# df_n.rename(columns={'geometry_x': 'geometry'})
	gdf = gpd.GeoDataFrame(df_n, geometry='geometry_x', crs="EPSG:4326")
	gdf = gdf.to_crs(orig_crs)
	return gdf


def drop_tiny_lines(bds, method='percentile', **kwargs):
	reduced = kwargs.get('reduced', True)
	orig_crs = bds.crs
	bds = bds[bds['geometry'].geom_type == 'LineString']
	df = bds.to_crs(4326)
	df['length'] = df.length

	if method == 'percentile':
		percentile = kwargs.get('percentile', 0.01)
		cut_off = df.length.quantile(percentile)
		df = df[df['length'] >= cut_off]

	if method == 'number_of_std':
		num_dev = kwargs.get('num_dev', 2)
		cut_off = df.length.mean() - num_dev * df.length.std()
		df = df[df['length'] >= cut_off]

	if method == 'length':
		meters = kwargs.get('meters', 500)
		cut_off = meters / 100000
		df = df[df['length'] >= cut_off]

	if reduced == True:
		df = df.drop(df.loc[:, 'length'.columns], axis=1)

	df_n = df.to_crs(orig_crs)
	return df_n

# Rename key value coumn to id before use
def remove_sliver(polygons, bound):
	orig_crs = polygons.crs
	clip_crs = bound.crs
	area = polygons.to_crs(4326)
	clip = bound.to_crs(4326)

	def find_nearest_id(point, centers):
		nearest_center = centers.geometry.distance(point).idxmin()
		return centers.loc[nearest_center, 'id']

	centroids = area.geometry.centroid
	centroids = centroids.to_crs(area.crs)
	center = gpd.GeoDataFrame(area[['id']].copy(), geometry=centroids)
	center.set_crs(area.crs, inplace=True)
	points = np.array(list(center.geometry.apply(lambda geom: (geom.x, geom.y))))
	points = np.unique(points, axis=0)
	hull = ConvexHull(points)
	hull_points = points[hull.vertices]
	padding = 0.1  # Adjust the padding as necessary
	min_x, min_y = hull_points.min(axis=0)
	max_x, max_y = hull_points.max(axis=0)

	min_x -= padding
	max_x += padding
	min_y -= padding
	max_y += padding
	dummy_points = np.array([
		[min_x, min_y], [min_x, max_y], [max_x, min_y], [max_x, max_y],  # corners
		[(min_x + max_x) / 2, min_y], [(min_x + max_x) / 2, max_y],  # midpoints of top and bottom edges
		[min_x, (min_y + max_y) / 2], [max_x, (min_y + max_y) / 2]  # midpoints of left and right edges
	])
	all_points = np.vstack([points, dummy_points])
	vor = Voronoi(all_points)
	voronoi_polygons = []
	for point_index, region_index in enumerate(vor.point_region):
		region = vor.regions[region_index]
		if not -1 in region and len(region) > 0:
			polygon = Polygon([vor.vertices[i] for i in region])
			if polygon.is_valid:
				# Match the centroid back to its original id
				if point_index < len(center):
					voronoi_polygons.append({'geometry': polygon, 'id': center.iloc[point_index]['id']})
	voronoi_gdf = gpd.GeoDataFrame(voronoi_polygons, crs=area.crs)
	merged_geometry = area.unary_union
	external_boundary = merged_geometry.convex_hull
	external_boundary_gdf = gpd.GeoDataFrame([{'geometry': external_boundary}], crs=area.crs)
	multi_poly = area.geometry.unary_union
	rotated_bounding_box = multi_poly.minimum_rotated_rectangle
	rotated_envelope = gpd.GeoDataFrame([{'geometry': rotated_bounding_box}], crs=area.crs)

	# Pick one or the other

	clipped_voronoi = gpd.overlay(voronoi_gdf, rotated_envelope, how='intersection')
	clipped_voronoi = gpd.overlay(voronoi_gdf, external_boundary_gdf, how='intersection')

	clipped_voronoi['centroid'] = clipped_voronoi.geometry.centroid
	clipped_voronoi['id'] = clipped_voronoi['centroid'].apply(lambda x: find_nearest_id(x, center))
	clipped_voronoi.drop(columns='centroid', inplace=True)
	difference_polygons = []
	for index, row in clipped_voronoi.iterrows():
		poly = row['geometry']
		poly_id = row['id']
		for area_poly in area.geometry:
			poly = poly.difference(area_poly)
		if poly.is_valid and not poly.is_empty:
			difference_polygons.append({'geometry': poly, 'id': poly_id})
	difference_gdf = gpd.GeoDataFrame(difference_polygons, crs=area.crs)
	merged_gdf = area.merge(difference_gdf, on='id', suffixes=('_area', '_difference'))
	combined_polygons = []
	for index, row in merged_gdf.iterrows():
		# Combine the geometries from area and difference_gdf
		combined_geometry = row['geometry_area'].union(row['geometry_difference'])

		# Append the combined geometry and id to the results
		combined_polygons.append({'geometry': combined_geometry, 'id': row['id']})
	combined_gdf = gpd.GeoDataFrame(combined_polygons, crs=area.crs)
	final_output = gpd.overlay(combined_gdf, clip, how='intersection')
	final_output = final_output.drop('id_2', axis=1)
	return final_output


def poly_to_line(GDF_polygon):
    og_crs = GDF_polygon.crs
    GDF_polygon = GDF_polygon.to_crs(4326)
    # Function to convert Polygon or MultiPolygon to LineString
    def polygon_to_linestring(geometry):
        if isinstance(geometry, Polygon):
            return LineString(geometry.exterior)
        elif isinstance(geometry, MultiPolygon):
            # Convert each polygon in the MultiPolygon to a LineString and return a list of LineStrings
            return [LineString(p.exterior) for p in geometry]
        else:
            return None
    # Apply the conversion to each geometry in the GeoDataFrame
    GDF_polygon['geometry'] = GDF_polygon['geometry'].apply(polygon_to_linestring)
    
    # Explode the list of LineStrings (if any) into separate rows
    output = GDF_polygon.explode(index_parts=False).reset_index(drop=True)
    output = output.to_crs(og_crs)
    return output


def remove_overlaps(df1, df2):
    def remove_overlapping_segments(gdf1, gdf2):
        result_geometries = []
    
        for geom1 in gdf1.geometry:
            # Initialize the non-overlapping portion of the line
            non_overlapping_geom = geom1
    
            for geom2 in gdf2.geometry:
                if non_overlapping_geom.intersects(geom2):
                    # Subtract the overlapping portion
                    non_overlapping_geom = non_overlapping_geom.difference(geom2)
    
                    # If the result is empty, break the loop
                    if non_overlapping_geom.is_empty:
                        break
    
            if not non_overlapping_geom.is_empty:
                # Ensure the result is exploded into individual LineStrings
                if isinstance(non_overlapping_geom, LineString):
                    result_geometries.append(non_overlapping_geom)
                else:  # Handle MultiLineString case
                    result_geometries.extend(non_overlapping_geom.geoms)
    
        # Create a new GeoDataFrame with only LineStrings
        exploded_gdf = gpd.GeoDataFrame(geometry=result_geometries, crs=gdf1.crs)
    
        return exploded_gdf
    gdf1_cleaned = remove_overlapping_segments(df1, df2)
    return gdf1_cleaned

