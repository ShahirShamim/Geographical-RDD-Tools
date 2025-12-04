import geopandas as gpd
from shapely.geometry import Point, Polygon, LineString
import numpy as np
import os

def create_samples():
    print("Generating sample data...")
    
    # 1. Create Polygons (e.g., Districts)
    polygons = [
        Polygon([(0, 0), (0, 10), (10, 10), (10, 0)]),
        Polygon([(10, 0), (10, 10), (20, 10), (20, 0)]),
        Polygon([(0, 10), (0, 20), (10, 20), (10, 10)]),
        Polygon([(10, 10), (10, 20), (20, 20), (20, 10)])
    ]
    ids = [1, 2, 3, 4]
    names = ['District A', 'District B', 'District C', 'District D']
    
    gdf_poly = gpd.GeoDataFrame(
        {'id': ids, 'name': names, 'geometry': polygons}, 
        crs="EPSG:4326"
    )
    
    # 2. Create Points (e.g., Addresses)
    # Generate random points within the bounding box of polygons
    np.random.seed(42)
    x = np.random.uniform(0, 20, 50)
    y = np.random.uniform(0, 20, 50)
    points = [Point(xi, yi) for xi, yi in zip(x, y)]
    
    gdf_points = gpd.GeoDataFrame(
        {'id': range(len(points)), 'value': np.random.randint(1, 100, 50), 'geometry': points},
        crs="EPSG:4326"
    )
    
    # 3. Create Overlapping Lines (for remove_overlaps)
    line1 = LineString([(5, 5), (15, 5)])
    line2 = LineString([(10, 5), (20, 5)]) # Overlaps with line1 from (10,5) to (15,5)
    
    gdf_lines1 = gpd.GeoDataFrame({'id': [1], 'geometry': [line1]}, crs="EPSG:4326")
    gdf_lines2 = gpd.GeoDataFrame({'id': [2], 'geometry': [line2]}, crs="EPSG:4326")

    # Save to files
    output_dir = os.path.dirname(os.path.abspath(__file__))
    gdf_poly.to_file(os.path.join(output_dir, "polygons.geojson"), driver="GeoJSON")
    gdf_points.to_file(os.path.join(output_dir, "points.geojson"), driver="GeoJSON")
    gdf_lines1.to_file(os.path.join(output_dir, "lines1.geojson"), driver="GeoJSON")
    gdf_lines2.to_file(os.path.join(output_dir, "lines2.geojson"), driver="GeoJSON")
    
    print(f"Sample data saved to {output_dir}")

if __name__ == "__main__":
    create_samples()
