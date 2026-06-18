import pytest
import geopandas as gpd
import numpy as np
import shapely
from shapely.geometry import Point, LineString, Polygon, MultiPolygon

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
    filter_by_boundary_distance
)

def test_points_in_polygon():
    # Setup simple polygons (districts)
    polys = [
        Polygon([(0,0), (0,10), (10,10), (10,0)]),
        Polygon([(10,0), (10,10), (20,10), (20,0)])
    ]
    gdf_poly = gpd.GeoDataFrame({'id': [1, 2], 'name': ['District A', 'District B'], 'geometry': polys}, crs='EPSG:4326')
    
    # Points inside districts
    pts = [Point(5, 5), Point(15, 5), Point(25, 25)] # Last one is outside
    gdf_pts = gpd.GeoDataFrame({'id': [10, 20, 30], 'name': ['pt1', 'pt2', 'pt3'], 'geometry': pts}, crs='EPSG:4326')
    
    res = points_in_polygon(gdf_pts, gdf_poly, suffix_name='_district')
    
    # Column assertions (should be single-underscore after cleanup)
    assert 'name_district' in res.columns
    assert 'id_district' in res.columns
    
    # Value assertions
    assert res.loc[res['id_'] == 10, 'name_district'].values[0] == 'District A'
    assert res.loc[res['id_'] == 20, 'name_district'].values[0] == 'District B'
    assert np.isnan(res.loc[res['id_'] == 30, 'id_district'].values[0])


def test_poly_to_line():
    polys = [
        Polygon([(0,0), (0,10), (10,10), (10,0)]),
        MultiPolygon([
            Polygon([(10,0), (10,10), (20,10), (20,0)]),
            Polygon([(20,0), (20,10), (30,10), (30,0)])
        ])
    ]
    gdf_poly = gpd.GeoDataFrame({'id': [1, 2], 'geometry': polys}, crs='EPSG:4326')
    
    res = poly_to_line(gdf_poly)
    
    # Result should contain 3 linestrings (1 from polygon, 2 from multipolygon)
    assert len(res) == 3
    assert (res['geometry'].geom_type == 'LineString').all()
    assert res.crs == 'EPSG:4326'

    # Test with naive geometries (No CRS)
    gdf_naive = gpd.GeoDataFrame({'id': [1], 'geometry': [polys[0]]}, crs=None)
    res_naive = poly_to_line(gdf_naive)
    assert len(res_naive) == 1
    assert res_naive.crs is None


def test_turner_correctness_and_dx_zero():
    # Horizontal boundary line along y=0, from x=-10 to x=10
    boundary = LineString([(-10, 0), (10, 0)])
    gdf_bds = gpd.GeoDataFrame({'id': [1], 'geometry': [boundary]}, crs='EPSG:3857')
    
    # Point at (0, 5) - nearest point on boundary is (0, 0)
    # The shortest line from (0, 5) to (0, 0) is vertical (dx=0, dy=-5)
    # With orth_distance=5, the orthogonal points along boundary are (-5, 0) and (5, 0)
    # Since boundary extends from -10 to 10, these points lie inside the boundary segment,
    # so the point should pass the Turner check.
    gdf_pts_pass = gpd.GeoDataFrame({'geometry': [Point(0, 5)]}, crs='EPSG:3857')
    res_pass = turner(gdf_pts_pass, gdf_bds, orth_distance=5)
    assert res_pass['turner_pass'].values[0] == True
    
    # Point at (9, 5) - nearest point on boundary is (9, 0)
    # With orth_distance=5, orthogonal points are (4, 0) and (14, 0)
    # (14, 0) is outside the boundary (-10 to 10), so it should fail the Turner check.
    gdf_pts_fail = gpd.GeoDataFrame({'geometry': [Point(9, 5)]}, crs='EPSG:3857')
    res_fail = turner(gdf_pts_fail, gdf_bds, orth_distance=5)
    assert res_fail['turner_pass'].values[0] == False


def test_turner_naive_crs():
    # Ensure turner doesn't crash on naive CRS
    boundary = LineString([(-10, 0), (10, 0)])
    gdf_bds = gpd.GeoDataFrame({'id': [1], 'geometry': [boundary]}, crs=None)
    gdf_pts = gpd.GeoDataFrame({'geometry': [Point(0, 5)]}, crs=None)
    
    res = turner(gdf_pts, gdf_bds, orth_distance=5)
    assert 'turner_pass' in res.columns
    assert res.crs is None


def test_drop_tiny_lines():
    lines = [
        LineString([(0,0), (1,0)]),     # length 1
        LineString([(0,0), (10,0)]),    # length 10
        LineString([(0,0), (100,0)])    # length 100
    ]
    gdf_lines = gpd.GeoDataFrame({'id': [1, 2, 3], 'geometry': lines}, crs='EPSG:3857')
    
    # Test 'length' method
    res = drop_tiny_lines(gdf_lines, method='length', meters=50)
    assert len(res) == 1
    assert res['id'].values[0] == 3
    
    # Test 'percentile' method (quantile 0.5 cuts off length < 10)
    res_p = drop_tiny_lines(gdf_lines, method='percentile', percentile=0.5)
    assert len(res_p) == 2
    assert 1 not in res_p['id'].values


def test_remove_sliver_correctness():
    # Grid of 4 small polygons with gaps
    # Poly 1: (0, 0) to (2, 2) -> Centroid (1, 1)
    # Poly 2: (0, 3) to (2, 5) -> Centroid (1, 4)
    # Poly 3: (3, 0) to (5, 2) -> Centroid (4, 1)
    # Poly 4: (3, 3) to (5, 5) -> Centroid (4, 4)
    # Order them differently to verify ID preservation
    polys = [
        Polygon([(3, 3), (3, 5), (5, 5), (5, 3)]), # Centroid (4, 4) -> ID 4
        Polygon([(0, 3), (0, 5), (2, 5), (2, 3)]), # Centroid (1, 4) -> ID 2
        Polygon([(3, 0), (3, 2), (5, 2), (5, 0)]), # Centroid (4, 1) -> ID 3
        Polygon([(0, 0), (0, 2), (2, 2), (2, 0)])  # Centroid (1, 1) -> ID 1
    ]
    ids = [4, 2, 3, 1]
    
    gdf_poly = gpd.GeoDataFrame({'id': ids, 'geometry': polys}, crs='EPSG:4326')
    boundary = gpd.GeoDataFrame({'geometry': [Polygon([(-1, -1), (-1, 6), (6, 6), (6, -1)])]}, crs='EPSG:4326')
    
    res = remove_sliver(gdf_poly, boundary)
    
    # Ensure IDs are correct and geometries are merged without losing any polygons
    assert len(res) == 4
    for orig_id in ids:
        assert orig_id in res['id'].values
        
    # Check that they filled the gaps (total area should increase)
    assert res.to_crs(3857).area.sum() > gdf_poly.to_crs(3857).area.sum()


def test_remove_overlaps():
    # Line 1: (0,0) to (10,0)
    # Line 2: (5,0) to (15,0) (overlaps from 5 to 10)
    l1 = LineString([(0, 0), (10, 0)])
    l2 = LineString([(5, 0), (15, 0)])
    
    gdf1 = gpd.GeoDataFrame({'id': [1], 'geometry': [l1]}, crs='EPSG:3857')
    gdf2 = gpd.GeoDataFrame({'id': [2], 'geometry': [l2]}, crs='EPSG:3857')
    
    res = remove_overlaps(gdf1, gdf2)
    
    assert len(res) == 1
    assert res['geometry'].iloc[0].equals(LineString([(0, 0), (5, 0)]))


def test_crs_alignment_warning():
    boundary = LineString([(-10, 0), (10, 0)])
    gdf_bds = gpd.GeoDataFrame({'id': [1], 'geometry': [boundary]}, crs='EPSG:4326')
    gdf_pts = gpd.GeoDataFrame({'geometry': [Point(0, 5)]}, crs='EPSG:3857')
    
    # Should issue a UserWarning about CRS mismatch
    with pytest.warns(UserWarning, match="CRS mismatch: reprojecting second GeoDataFrame"):
        res = turner(gdf_pts, gdf_bds, orth_distance=5)
    
    # Should automatically align CRS to points_gdf CRS
    assert res.crs == 'EPSG:3857'


def test_remove_sliver_custom_id_col():
    polys = [
        Polygon([(3, 3), (3, 5), (5, 5), (5, 3)]),
        Polygon([(0, 3), (0, 5), (2, 5), (2, 3)]),
        Polygon([(3, 0), (3, 2), (5, 2), (5, 0)]),
        Polygon([(0, 0), (0, 2), (2, 2), (2, 0)])
    ]
    custom_ids = [40, 20, 30, 10]
    gdf_poly = gpd.GeoDataFrame({'district_id': custom_ids, 'geometry': polys}, crs='EPSG:4326')
    boundary = gpd.GeoDataFrame({'geometry': [Polygon([(-1, -1), (-1, 6), (6, 6), (6, -1)])]}, crs='EPSG:4326')
    
    res = remove_sliver(gdf_poly, boundary, id_col='district_id')
    assert 'district_id' in res.columns
    assert 'id' not in res.columns
    for cid in custom_ids:
        assert cid in res['district_id'].values


def test_calculate_signed_distance():
    # Treatment area
    treatment = gpd.GeoDataFrame(
        {'geometry': [Polygon([(0,0), (0,10), (10,10), (10,0)])]}, 
        crs='EPSG:3857'
    )
    # Boundary between treatment and control
    boundary = gpd.GeoDataFrame(
        {'geometry': [LineString([(10, 0), (10, 10)])]}, 
        crs='EPSG:3857'
    )
    # Points
    points = gpd.GeoDataFrame(
        {'id': [1, 2], 'geometry': [Point(5, 5), Point(15, 5)]},
        crs='EPSG:3857'
    )
    
    res = calculate_signed_distance(points, boundary, treatment)
    
    assert res.loc[res['id'] == 1, 'is_treated'].values[0] == True
    assert res.loc[res['id'] == 2, 'is_treated'].values[0] == False
    assert pytest.approx(res.loc[res['id'] == 1, 'distance'].values[0]) == 5.0
    assert pytest.approx(res.loc[res['id'] == 2, 'distance'].values[0]) == 5.0
    assert pytest.approx(res.loc[res['id'] == 1, 'signed_distance'].values[0]) == 5.0
    assert pytest.approx(res.loc[res['id'] == 2, 'signed_distance'].values[0]) == -5.0


def test_extract_shared_boundaries():
    polys = [
        Polygon([(0,0), (0,10), (10,10), (10,0)]),
        Polygon([(10,0), (10,10), (20,10), (20,0)])
    ]
    gdf = gpd.GeoDataFrame({'district_id': [101, 102], 'geometry': polys}, crs='EPSG:3857')
    
    shared = extract_shared_boundaries(gdf, id_col='district_id')
    
    assert len(shared) == 1
    row = shared.iloc[0]
    # left_id and right_id are determined by index ordering, so left_id should be 101, right_id 102
    assert row['left_id'] == 101
    assert row['right_id'] == 102
    
    # The geometry should be a LineString along x=10
    geom = row['geometry']
    assert isinstance(geom, LineString)
    assert geom.equals(LineString([(10, 0), (10, 10)])) or geom.equals(LineString([(10, 10), (10, 0)]))


def test_shift_boundary_placebo():
    boundary = gpd.GeoDataFrame(
        {'geometry': [LineString([(0, 0), (10, 0)])]}, 
        crs='EPSG:3857'
    )
    # Shift east by 5 meters, north by 10 meters
    shifted = shift_boundary_placebo(boundary, xoff=5.0, yoff=10.0)
    
    assert len(shifted) == 1
    geom = shifted['geometry'].iloc[0]
    assert isinstance(geom, LineString)
    assert geom.equals(LineString([(5, 10), (15, 10)]))


def test_filter_by_boundary_distance():
    boundary = gpd.GeoDataFrame(
        {'geometry': [LineString([(10, 0), (10, 10)])]}, 
        crs='EPSG:3857'
    )
    points = gpd.GeoDataFrame(
        {'id': [1, 2, 3], 'geometry': [Point(8, 5), Point(12, 5), Point(20, 5)]},
        crs='EPSG:3857'
    )
    
    # Max distance = 5. Point 3 is 10 units away, so it should be filtered out
    filtered = filter_by_boundary_distance(points, boundary, max_distance=5.0)
    
    assert len(filtered) == 2
    assert 1 in filtered['id'].values
    assert 2 in filtered['id'].values
    assert 3 not in filtered['id'].values
