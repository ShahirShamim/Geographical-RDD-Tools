import sys
import os

# Ensure we can import from the current directory
sys.path.append(os.getcwd())

try:
    import geoRDDprep
    print("Successfully imported geoRDDprep from root")
    
    functions = [
        'points_in_polygon',
        'poly_to_line',
        'turner',
        'drop_tiny_lines',
        'remove_sliver',
        'remove_overlaps'
    ]
    
    for func_name in functions:
        if hasattr(geoRDDprep, func_name):
            print(f"Function '{func_name}' is available.")
        else:
            print(f"ERROR: Function '{func_name}' is MISSING.")
            
except ImportError as e:
    print(f"Failed to import geoRDDprep: {e}")
