import geopandas as gpd
import os
import sys
from shapely.geometry import box

# Add parent directory to path to import geoRDDprep if running from source
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geoRDDprep import (
    points_in_polygon,
    poly_to_line,
    turner,
    drop_tiny_lines,
    remove_sliver,
    remove_overlaps,
    calculate_signed_distance,
    extract_shared_boundaries,
    shift_boundary_placebo,
    filter_by_boundary_distance,
    segment_boundary,
    assign_nearest_boundary,
    snap_points_to_boundary
)

def run_examples():
    print("Running geoRDDprep examples...\n")
    
    # Load Data
    data_dir = os.path.dirname(os.path.abspath(__file__))
    try:
        points = gpd.read_file(os.path.join(data_dir, "points.geojson"))
        polygons = gpd.read_file(os.path.join(data_dir, "polygons.geojson"))
        lines1 = gpd.read_file(os.path.join(data_dir, "lines1.geojson"))
        lines2 = gpd.read_file(os.path.join(data_dir, "lines2.geojson"))
    except Exception as e:
        print(f"Error loading data: {e}")
        print("Did you run generate_data.py first?")
        return

    # 1. points_in_polygon
    print("--- 1. points_in_polygon ---")
    joined = points_in_polygon(points, polygons, suffix_name="_district")
    print(f"Assigned {len(joined)} points to districts.")
    print(f"Columns: {joined.columns}")
    print(joined.head(3))
    print("\n")

    # 2. poly_to_line
    print("--- 2. poly_to_line ---")
    lines_from_poly = poly_to_line(polygons)
    print(f"Converted {len(polygons)} polygons to {len(lines_from_poly)} lines.")
    print("\n")

    # 3. drop_tiny_lines
    print("--- 3. drop_tiny_lines ---")
    # Our synthetic lines are long, so let's set a high threshold to see filtering
    # Note: These are in degrees (EPSG:4326), so length is small numbers.
    # The function converts to meters (EPSG:3857) internally.
    clean_lines = drop_tiny_lines(lines_from_poly, method='length', meters=1000) 
    print(f"Retained {len(clean_lines)} lines after filtering.")
    print("\n")

    # 4. turner
    print("--- 4. turner ---")
    # Match points to the lines we created
    # Using a large distance because our synthetic data is sparse
    matched = turner(points, lines_from_poly, orth_distance=500000) 
    print(f"Turner analysis complete. Columns: {matched.columns}")
    print(f"Points passing Turner criteria: {matched['turner_pass'].sum()}")
    print("\n")

    # 5. remove_overlaps
    print("--- 5. remove_overlaps ---")
    print(f"Line 1 Length: {lines1.length.iloc[0]:.2f}")
    print(f"Line 2 Length: {lines2.length.iloc[0]:.2f}")
    
    cleaned_lines1 = remove_overlaps(lines1, lines2)
    
    if not cleaned_lines1.empty:
        print(f"Line 1 Length after removing overlap with Line 2: {cleaned_lines1.length.iloc[0]:.2f}")
    else:
        print("Line 1 was completely removed (full overlap).")
    print("\n")

    # 6. segment_boundary
    print("--- 6. segment_boundary ---")
    # Split the boundary lines into roughly 10 segments total (adaptive to data scale).
    total_len_m = clean_lines.to_crs(3857).length.sum()
    seg_len = max(total_len_m / 10.0, 1.0)
    segments = segment_boundary(clean_lines, segment_length=seg_len)
    print(f"Split {len(clean_lines)} boundary line(s) into {len(segments)} segments "
          f"(target {seg_len:.0f} m each).")
    print(f"First segment ids: {segments['segment_id'].tolist()[:10]}")
    print("\n")

    # 7. assign_nearest_boundary (boundary-segment fixed effects)
    print("--- 7. assign_nearest_boundary ---")
    # id_col names the source column to read; the value lands in 'boundary_id'.
    pts_fe = assign_nearest_boundary(points, segments, id_col='segment_id')
    print(f"Assigned each of {len(pts_fe)} points to its nearest boundary segment.")
    print(pts_fe[['boundary_id', 'boundary_distance']].head(3))
    print("\n")

    # 8. snap_points_to_boundary
    print("--- 8. snap_points_to_boundary ---")
    snapped = snap_points_to_boundary(points, clean_lines, distance_col='dist_to_border')
    print(f"Projected {len(snapped)} points onto the nearest boundary.")
    print(snapped[['snapped_geometry', 'dist_to_border']].head(3))
    print("\n")

    # 9. remove_sliver
    print("--- 9. remove_sliver ---")
    # Clip to the polygon extent and merge any slivers/gaps into their neighbours.
    minx, miny, maxx, maxy = polygons.total_bounds
    clip = gpd.GeoDataFrame({'geometry': [box(minx, miny, maxx, maxy)]}, crs=polygons.crs)
    cleaned_polys = remove_sliver(polygons, clip, id_col='id')
    print(f"Cleaned {len(polygons)} polygons into {len(cleaned_polys)} (slivers/gaps merged).")
    print("\n")

    # 10. extract_shared_boundaries
    print("--- 10. extract_shared_boundaries ---")
    shared = extract_shared_boundaries(polygons, id_col='id')
    print(f"Found {len(shared)} shared border(s) between adjacent districts.")
    if not shared.empty:
        print(shared[['left_id', 'right_id']].head())
    print("\n")

    # 11. calculate_signed_distance
    print("--- 11. calculate_signed_distance ---")
    # Treat district id==1 as the treated area; points inside get a positive distance.
    treatment = polygons[polygons['id'] == 1]
    signed = calculate_signed_distance(points, clean_lines, treatment)
    print(f"Treated points: {int(signed['is_treated'].sum())} of {len(signed)}.")
    print(signed[['is_treated', 'distance', 'signed_distance']].head(3))
    print("\n")

    # 12. filter_by_boundary_distance
    print("--- 12. filter_by_boundary_distance ---")
    # Keep only points within a 500 km bandwidth of the boundary.
    nearby = filter_by_boundary_distance(points, clean_lines, max_distance=500000.0)
    print(f"{len(nearby)} of {len(points)} points fall within the 500 km bandwidth.")
    print("\n")

    # 13. shift_boundary_placebo
    print("--- 13. shift_boundary_placebo ---")
    # Shift the boundary 100 km north to build a placebo border.
    placebo = shift_boundary_placebo(clean_lines, xoff=0.0, yoff=100000.0)
    print(f"Created a placebo boundary ({len(placebo)} segment(s)) shifted 100 km north.")
    print("\n")

    print("Examples completed successfully!")

if __name__ == "__main__":
    run_examples()
