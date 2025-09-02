#!/usr/bin/env python3
"""
Calculate closest distances from earthquake points to slab contours from USGS Slab2.0 GRD files.

Usage:
    python eq_distances.py --grd sco_slab2_dep_02.23.18.grd --thk sco_slab2_thk_02.23.18.grd --eq earthquakes.csv --out distances.csv
"""

import argparse
import pandas as pd
import numpy as np
import xarray as xr
from scipy.spatial.distance import cdist
from scipy.interpolate import griddata, RegularGridInterpolator
import sys
import os


def load_grd_file(grd_path):
    """Load GRD file using xarray."""
    try:
        # Try loading with xarray (works for NetCDF format GRD files)
        ds = xr.open_dataset(grd_path)
        return ds
    except Exception as e:
        print(f"Error loading GRD file: {e}")
        print("Trying alternative loading methods...")
        
        try:
            # Alternative: try loading as a raster with rasterio if available
            import rasterio
            with rasterio.open(grd_path) as src:
                data = src.read(1)
                transform = src.transform
                # Convert to xarray format
                lon = np.linspace(src.bounds.left, src.bounds.right, src.width)
                lat = np.linspace(src.bounds.bottom, src.bounds.top, src.height)
                ds = xr.Dataset({
                    'z': (['lat', 'lon'], data)
                }, coords={'lat': lat, 'lon': lon})
                return ds
        except ImportError:
            print("rasterio not available. Please install with: pip install rasterio")
        except Exception as e2:
            print(f"Alternative loading failed: {e2}")
            
        sys.exit(1)


def haversine_distance(lon1, lat1, lon2, lat2):
    """
    Calculate haversine distance between points in kilometers.
    
    Parameters:
    lon1, lat1: longitude and latitude of first point(s)
    lon2, lat2: longitude and latitude of second point(s)
    
    Returns:
    distance in kilometers
    """
    # Convert to radians
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    
    # Haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    
    # Earth's radius in kilometers
    R = 6371.0
    
    return R * c


def calculate_distances_to_slab(eq_df, ds, depth_var, lon_coord, lat_coord, thickness_ds):
    """
    Calculate 3D distances from earthquake points to slab surface and center.
    Positive distance = above slab (shallower than slab)
    Negative distance = below slab (deeper than slab)
    
    Parameters:
    eq_df: DataFrame with earthquake data
    ds: xarray Dataset containing slab depth data
    depth_var, lon_coord, lat_coord: coordinate names
    thickness_ds: Dataset with slab thickness data
    
    Returns:
    surface_distances: array of signed distances to slab surface
    center_distances: array of signed distances to slab center
    """
    print(f"Calculating 3D distances from earthquakes to slab surface and center...")
    
    # Get slab data
    depth_data = ds[depth_var]
    lons = depth_data[lon_coord].values
    lats = depth_data[lat_coord].values
    slab_surface_depths = depth_data.values
    
    # Get thickness data
    print("Processing thickness data to calculate slab center...")
    thickness_data = None
    for var_name in thickness_ds.data_vars.keys():
        thickness_data = thickness_ds[var_name]
        print(f"Using thickness variable: {var_name}")
        break
    
    if thickness_data is None:
        print("Error: No thickness data found in thickness file")
        sys.exit(1)
    
    slab_thickness = thickness_data.values
    # Calculate slab center: surface_depth - thickness/2
    # Note: slab depths are negative, thickness is positive
    slab_center_depths = slab_surface_depths - (slab_thickness / 2)
    print("Calculating distances to both SURFACE and CENTER")
    
    # Convert earthquake coordinates if needed
    eq_lons = eq_df['lon360'].values.copy()
    eq_lats = eq_df['latitude'].values
    eq_depths = eq_df['depth_km'].values
    
    # Convert longitudes to match slab coordinate system
    if np.any(eq_lons > 180) and np.all(lons <= 180):
        print("Converting earthquake longitudes from 0-360 to -180-180 to match slab grid")
        eq_lons = np.where(eq_lons > 180, eq_lons - 360, eq_lons)
    elif np.any(eq_lons < 0) and np.all(lons >= 0):
        print("Converting earthquake longitudes from -180-180 to 0-360 to match slab grid")
        eq_lons = np.where(eq_lons < 0, eq_lons + 360, eq_lons)
    
    # Prepare coordinate arrays for interpolation
    if lats[0] > lats[-1]:  # If latitudes are decreasing, flip
        lats = lats[::-1]
        slab_surface_depths = slab_surface_depths[::-1, :]
        slab_center_depths = slab_center_depths[::-1, :]
    if lons[0] > lons[-1]:  # If longitudes are decreasing, flip  
        lons = lons[::-1]
        slab_surface_depths = slab_surface_depths[:, ::-1]
        slab_center_depths = slab_center_depths[:, ::-1]
    
    # Create interpolators
    surface_interpolator = RegularGridInterpolator(
        (lats, lons), 
        slab_surface_depths, 
        method='linear', 
        bounds_error=False, 
        fill_value=np.nan
    )
    
    center_interpolator = RegularGridInterpolator(
        (lats, lons), 
        slab_center_depths, 
        method='linear', 
        bounds_error=False, 
        fill_value=np.nan
    )
    
    surface_distances = []
    center_distances = []
    
    for i, (eq_lon, eq_lat, eq_depth) in enumerate(zip(eq_lons, eq_lats, eq_depths)):
        # Interpolate slab depths at earthquake location
        slab_surface_at_eq = surface_interpolator([eq_lat, eq_lon])[0]
        slab_center_at_eq = center_interpolator([eq_lat, eq_lon])[0]
        
        if np.isnan(slab_surface_at_eq):
            # If outside slab grid, find nearest slab point
            if i == 0:  # Only print this message once
                print("Some earthquakes outside slab grid, finding nearest slab points...")
            
            # Create meshgrid and find valid points
            lon_grid, lat_grid = np.meshgrid(lons, lats)
            valid_mask = ~np.isnan(slab_surface_depths)
            
            if np.any(valid_mask):
                valid_lons = lon_grid[valid_mask]
                valid_lats = lat_grid[valid_mask]
                valid_surface_depths = slab_surface_depths[valid_mask]
                valid_center_depths = slab_center_depths[valid_mask]
                
                # Find closest valid slab point
                distances_to_slab_points = haversine_distance(
                    eq_lon, eq_lat, valid_lons, valid_lats
                )
                closest_idx = np.argmin(distances_to_slab_points)
                
                horizontal_distance = distances_to_slab_points[closest_idx]
                slab_surface_at_closest = valid_surface_depths[closest_idx]
                slab_center_at_closest = valid_center_depths[closest_idx]
                
                # Calculate surface distance
                vertical_distance_surface = eq_depth - slab_surface_at_closest
                distance_3d_surface = np.sqrt(horizontal_distance**2 + vertical_distance_surface**2)
                
                if eq_depth > slab_surface_at_closest:
                    surface_distances.append(distance_3d_surface)  # Above surface
                else:
                    surface_distances.append(-distance_3d_surface)  # Below surface
                
                # Calculate center distance
                vertical_distance_center = eq_depth - slab_center_at_closest
                distance_3d_center = np.sqrt(horizontal_distance**2 + vertical_distance_center**2)
                
                if eq_depth > slab_center_at_closest:
                    center_distances.append(distance_3d_center)  # Above center
                else:
                    center_distances.append(-distance_3d_center)  # Below center
            else:
                surface_distances.append(np.nan)
                center_distances.append(np.nan)
                
        else:
            # Earthquake is within slab grid - calculate vertical distances
            
            # Surface distance
            vertical_distance_surface = eq_depth - slab_surface_at_eq
            if eq_depth > slab_surface_at_eq:
                surface_distances.append(abs(vertical_distance_surface))  # Above surface
            else:
                surface_distances.append(-abs(vertical_distance_surface))  # Below surface
            
            # Center distance
            if not np.isnan(slab_center_at_eq):
                vertical_distance_center = eq_depth - slab_center_at_eq
                if eq_depth > slab_center_at_eq:
                    center_distances.append(abs(vertical_distance_center))  # Above center
                else:
                    center_distances.append(-abs(vertical_distance_center))  # Below center
            else:
                center_distances.append(np.nan)
        
        if (i + 1) % 100 == 0 or (i + 1) == len(eq_lons):
            print(f"Processed {i+1}/{len(eq_lons)} earthquakes")
    
    return np.array(surface_distances), np.array(center_distances)


def main():
    parser = argparse.ArgumentParser(
        description='Calculate distances from earthquakes to slab surface and center'
    )
    parser.add_argument('--grd', required=True, help='Path to depth GRD file')
    parser.add_argument('--thk', required=True, help='Path to thickness GRD file')
    parser.add_argument('--eq', required=True, help='Path to earthquakes CSV file')
    parser.add_argument('--out', required=True, help='Output CSV file path')
    
    args = parser.parse_args()
    
    # Check if files exist
    if not os.path.exists(args.grd):
        print(f"Error: GRD file not found: {args.grd}")
        sys.exit(1)
    
    if not os.path.exists(args.thk):
        print(f"Error: Thickness file not found: {args.thk}")
        sys.exit(1)
    
    if not os.path.exists(args.eq):
        print(f"Error: Earthquake file not found: {args.eq}")
        sys.exit(1)
    
    # Load earthquake data
    print(f"Loading earthquake data from {args.eq}")
    try:
        eq_df = pd.read_csv(args.eq)
        required_cols = ['date_time', 'lon360', 'latitude', 'depth_km', 'magnitude']
        
        if not all(col in eq_df.columns for col in required_cols):
            print(f"Error: Earthquake file must contain columns: {required_cols}")
            print(f"Found columns: {list(eq_df.columns)}")
            sys.exit(1)
        
        print(f"Loaded {len(eq_df)} earthquakes")
        
    except Exception as e:
        print(f"Error loading earthquake file: {e}")
        sys.exit(1)
    
    # Load GRD files
    print(f"Loading depth GRD file from {args.grd}")
    ds = load_grd_file(args.grd)
    
    print(f"Loading thickness GRD file from {args.thk}")
    thickness_ds = load_grd_file(args.thk)
    
    # Get coordinate and variable information
    data_vars = list(ds.data_vars.keys())
    print(f"Available data variables: {data_vars}")
    
    if 'z' in data_vars:
        depth_var = 'z'
    elif len(data_vars) == 1:
        depth_var = data_vars[0]
    else:
        print(f"Multiple variables found: {data_vars}")
        depth_var = data_vars[0]
        print(f"Using variable: {depth_var}")
    
    depth_data = ds[depth_var]
    coords = list(depth_data.coords.keys())
    print(f"Available coordinates: {coords}")
    
    # Identify coordinates
    lat_coord = None
    lon_coord = None
    for coord in coords:
        coord_lower = coord.lower()
        if any(pattern in coord_lower for pattern in ['lat', 'y']):
            lat_coord = coord
        elif any(pattern in coord_lower for pattern in ['lon', 'x']):
            lon_coord = coord
    
    # If automatic detection fails, analyze coordinate values
    if lat_coord is None or lon_coord is None:
        print("Analyzing coordinate values for identification...")
        for coord in coords:
            values = depth_data[coord].values
            print(f"{coord}: min={np.min(values):.2f}, max={np.max(values):.2f}")
            
            if np.min(values) >= -90 and np.max(values) <= 90:
                lat_coord = coord
                print(f"Identified {coord} as latitude")
            elif (np.min(values) >= -180 and np.max(values) <= 180) or (np.min(values) >= 0 and np.max(values) <= 360):
                lon_coord = coord
                print(f"Identified {coord} as longitude")
    
    if lat_coord is None or lon_coord is None:
        print(f"Could not identify coordinates. Available: {coords}")
        sys.exit(1)
    
    print(f"Using coordinates: longitude={lon_coord}, latitude={lat_coord}, depth_variable={depth_var}")
    
    # Calculate 3D signed distances to slab
    print("Calculating signed distances to slab SURFACE and CENTER...")
    print("  Positive distance = earthquake above slab (shallower)")
    print("  Negative distance = earthquake below slab (deeper)")
    
    surface_distances, center_distances = calculate_distances_to_slab(
        eq_df, ds, depth_var, lon_coord, lat_coord, thickness_ds
    )
    
    # Create output dataframe
    output_df = eq_df.copy()
    output_df['distance_to_surface_km'] = surface_distances
    output_df['distance_to_center_km'] = center_distances
    
    # Save results
    print(f"Saving results to {args.out}")
    output_df.to_csv(args.out, index=False)
    
    # Print summary statistics
    valid_surface = surface_distances[~np.isnan(surface_distances)]
    valid_center = center_distances[~np.isnan(center_distances)]
    
    print("\nSummary Statistics:")
    print(f"Number of earthquakes processed: {len(output_df)}")
    print(f"Number with valid distances: {len(valid_surface)}")
    
    if len(valid_surface) > 0:
        print(f"\nSURFACE distances:")
        print(f"  Mean: {np.mean(valid_surface):.2f} km")
        print(f"  Median: {np.median(valid_surface):.2f} km")
        print(f"  Min: {np.min(valid_surface):.2f} km")
        print(f"  Max: {np.max(valid_surface):.2f} km")
        
        above_surface = valid_surface > 0
        below_surface = valid_surface < 0
        print(f"  Above surface (positive): {np.sum(above_surface)} ({100*np.sum(above_surface)/len(valid_surface):.1f}%)")
        print(f"  Below surface (negative): {np.sum(below_surface)} ({100*np.sum(below_surface)/len(valid_surface):.1f}%)")
    
    if len(valid_center) > 0:
        print(f"\nCENTER distances:")
        print(f"  Mean: {np.mean(valid_center):.2f} km")
        print(f"  Median: {np.median(valid_center):.2f} km")
        print(f"  Min: {np.min(valid_center):.2f} km")
        print(f"  Max: {np.max(valid_center):.2f} km")
        
        above_center = valid_center > 0
        below_center = valid_center < 0
        print(f"  Above center (positive): {np.sum(above_center)} ({100*np.sum(above_center)/len(valid_center):.1f}%)")
        print(f"  Below center (negative): {np.sum(below_center)} ({100*np.sum(below_center)/len(valid_center):.1f}%)")
    
    print(f"\nResults saved successfully to {args.out}")
    print("CSV contains: date_time, lon360, latitude, depth_km, magnitude, distance_to_surface_km, distance_to_center_km")


if __name__ == "__main__":
    main()