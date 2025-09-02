#!/usr/bin/env python3
"""
Slab 2.0 Cross-Section Extractor with Plotting and Earthquake Filtering

Extracts north-south or east-west cross-sections from Slab 2.0 grid files
and outputs top and bottom slab depths to CSV files, with optional plotting
and earthquake point filtering. Automatically crops the plot around earthquake
data when provided.

Usage:
python slab2_cross_section_extractor.py --input_dep van_slab2_dep_02.23.18.grd \
                                        --input_thk van_slab2_thk_02.23.18.grd \
                                        --orientation n_s \
                                        --center_lat 8 \
                                        --center_lon 126 \
                                        --output_plot van_cross_section.png \
                                        --dpi 300 \
                                        --title "Philippines Slab Cross-Section" \
                                        --earthquakes quakes.csv \
                                        --buffer_degrees 2
"""

import argparse
import numpy as np
import pandas as pd
try:
    import xarray as xr
except ImportError:
    print("Error: xarray is required. Install with: pip install xarray")
    exit(1)

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("Error: matplotlib is required for plotting. Install with: pip install matplotlib")
    exit(1)

def haversine_distance(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # Convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    r = 6371  # Radius of earth in kilometers
    return c * r

def filter_earthquakes(earthquakes_df, center_lon, center_lat, buffer_degrees):
    """
    Filter earthquakes within a buffer around the cross-section line
    """
    # Calculate distances from the cross-section line
    distances = haversine_distance(
        earthquakes_df['lon360'], 
        earthquakes_df['latitude'], 
        center_lon, 
        center_lat
    )
    
    # Convert buffer from degrees to kilometers (rough approximation)
    buffer_km = buffer_degrees * 111  # 1 degree is roughly 111 km
    
    # Filter earthquakes within the buffer
    filtered_quakes = earthquakes_df[distances <= buffer_km].copy()
    
    print(f"Total earthquakes: {len(earthquakes_df)}")
    print(f"Earthquakes within {buffer_degrees}° buffer: {len(filtered_quakes)}")
    
    return filtered_quakes

def inspect_grd_file(filename):
    """
    Inspect the structure of a GMT grid file
    """
    try:
        ds = xr.open_dataset(filename)
        print(f"\n=== Inspecting {filename} ===")
        print(f"Dimensions: {list(ds.dims.keys())}")
        print(f"Coordinates: {list(ds.coords.keys())}")
        print(f"Data variables: {list(ds.data_vars.keys())}")
        
        # Print coordinate ranges
        for coord in ds.coords:
            coord_data = ds.coords[coord]
            if coord_data.size > 1:
                print(f"{coord}: {coord_data.min().values:.3f} to {coord_data.max().values:.3f} ({coord_data.size} points)")
            else:
                print(f"{coord}: {coord_data.values}")
        
        # Print data variable info
        for var in ds.data_vars:
            var_data = ds[var]
            print(f"Data variable '{var}': shape {var_data.shape}, dims {var_data.dims}")
            valid_data = var_data.where(~np.isnan(var_data))
            if valid_data.size > 0:
                print(f"  Value range: {valid_data.min().values:.3f} to {valid_data.max().values:.3f}")
        
        return ds
    except Exception as e:
        print(f"Error inspecting {filename}: {e}")
        return None

def load_grd_file(filename):
    """
    Load a GMT grid file using xarray
    """
    try:
        ds = xr.open_dataset(filename)
        return ds
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return None

def get_coordinate_names(dataset):
    """
    Determine the coordinate names in the dataset
    """
    coords = list(dataset.coords.keys())
    dims = list(dataset.dims.keys())
    
    # Common coordinate name variations
    lat_names = ['lat', 'latitude', 'y', 'Y']
    lon_names = ['lon', 'longitude', 'x', 'X']
    
    lat_coord = None
    lon_coord = None
    
    # Find latitude coordinate
    for name in lat_names:
        if name in coords or name in dims:
            lat_coord = name
            break
    
    # Find longitude coordinate
    for name in lon_names:
        if name in coords or name in dims:
            lon_coord = name
            break
    
    # If not found, use the first two coordinates/dimensions
    if lat_coord is None or lon_coord is None:
        available = coords if coords else dims
        if len(available) >= 2:
            # Assume first is longitude, second is latitude (common GMT convention)
            lon_coord = available[0]
            lat_coord = available[1]
    
    return lat_coord, lon_coord

def extract_cross_section(dep_data, thk_data, orientation, center_lat, center_lon):
    """
    Extract cross-section data based on orientation
    
    Parameters:
    - dep_data: xarray dataset with depth data
    - thk_data: xarray dataset with thickness data  
    - orientation: 'n_s' for north-south, 'w_e' for east-west
    - center_lat: latitude of cross-section center
    - center_lon: longitude of cross-section center
    
    Returns:
    - DataFrame with lon360, latitude, top_depth, bottom_depth
    """
    
    # Get coordinate names
    lat_coord, lon_coord = get_coordinate_names(dep_data)
    print(f"Using coordinates: latitude='{lat_coord}', longitude='{lon_coord}'")
    
    # Get variable names
    dep_var = list(dep_data.data_vars)[0]
    thk_var = list(thk_data.data_vars)[0]
    print(f"Using variables: depth='{dep_var}', thickness='{thk_var}'")
    
    if orientation.lower() == 'n_s':
        # North-South cross-section: fix longitude, vary latitude
        # Find closest longitude index
        lon_idx = np.abs(dep_data[lon_coord] - center_lon).argmin()
        actual_lon = float(dep_data[lon_coord].isel({lon_coord: lon_idx}).values)
        print(f"Selected longitude: {actual_lon:.3f} (requested: {center_lon})")
        
        # Extract data along this longitude
        dep_slice = dep_data[dep_var].isel({lon_coord: lon_idx})
        thk_slice = thk_data[thk_var].isel({lon_coord: lon_idx})
        
        # Create output arrays
        lats = dep_slice[lat_coord].values
        lons = np.full_like(lats, actual_lon)
        x_coord = lats  # For plotting
        x_label = 'Latitude (°)'
        
    elif orientation.lower() == 'w_e':
        # East-West cross-section: fix latitude, vary longitude
        # Find closest latitude index
        lat_idx = np.abs(dep_data[lat_coord] - center_lat).argmin()
        actual_lat = float(dep_data[lat_coord].isel({lat_coord: lat_idx}).values)
        print(f"Selected latitude: {actual_lat:.3f} (requested: {center_lat})")
        
        # Extract data along this latitude
        dep_slice = dep_data[dep_var].isel({lat_coord: lat_idx})
        thk_slice = thk_data[thk_var].isel({lat_coord: lat_idx})
        
        # Create output arrays
        lons = dep_slice[lon_coord].values
        lats = np.full_like(lons, actual_lat)
        x_coord = lons  # For plotting
        x_label = 'Longitude (°)'
        
    else:
        raise ValueError("Orientation must be 'n_s' or 'w_e'")
    
    # Get depth and thickness values
    top_depths = dep_slice.values
    thicknesses = thk_slice.values
    
    # Calculate bottom depths (top - thickness)
    # Bottom should be deeper (more negative) than top
    bottom_depths = top_depths - thicknesses
    
    # Convert longitudes to 0-360 format if needed
    lons_360 = np.where(lons < 0, lons + 360, lons)
    
    # Create DataFrame
    df = pd.DataFrame({
        'lon360': lons_360,
        'latitude': lats,
        'top_depth_km': -np.abs(top_depths),  # Make negative (depth below surface)
        'bottom_depth_km': -np.abs(bottom_depths),  # Make negative (depth below surface)
        'x_coord': x_coord,  # For plotting
        'x_label': x_label   # For plotting
    })
    
    # Remove rows with NaN values
    df = df.dropna()
    
    # Verify that bottom is deeper than top
    if len(df) > 0:
        print(f"Top depth range: {df['top_depth_km'].min():.1f} to {df['top_depth_km'].max():.1f} km")
        print(f"Bottom depth range: {df['bottom_depth_km'].min():.1f} to {df['bottom_depth_km'].max():.1f} km")
        print(f"Thickness range: {(df['top_depth_km'] - df['bottom_depth_km']).min():.1f} to {(df['top_depth_km'] - df['bottom_depth_km']).max():.1f} km")
    
    return df

def calculate_plot_bounds(df, earthquakes_df=None, margin_percent=0.15):
    """
    Calculate plot bounds based on data, with optional cropping to earthquake extent
    
    Parameters:
    - df: DataFrame with slab data
    - earthquakes_df: DataFrame with earthquake data (optional)
    - margin_percent: Percentage margin to add around the data bounds
    
    Returns:
    - x_min, x_max, y_min, y_max
    """
    if earthquakes_df is not None and not earthquakes_df.empty:
        # Crop to earthquake extent with margin
        x_min_eq = earthquakes_df['x_coord'].min()
        x_max_eq = earthquakes_df['x_coord'].max()
        y_min_eq = earthquakes_df['depth_km'].min()
        y_max_eq = earthquakes_df['depth_km'].max()
        
        # Add margin to earthquake bounds
        x_range = x_max_eq - x_min_eq
        y_range = y_max_eq - y_min_eq
        
        x_margin = x_range * margin_percent if x_range > 0 else 1.0
        y_margin = y_range * margin_percent if y_range > 0 else 10.0
        
        x_min = x_min_eq - x_margin
        x_max = x_max_eq + x_margin
        y_min = y_min_eq - y_margin  # More negative (deeper)
        y_max = y_max_eq + y_margin  # Less negative (shallower)
        
        print(f"Cropping plot to earthquake extent:")
        print(f"  X-axis: {x_min:.2f} to {x_max:.2f} (margin: ±{x_margin:.2f})")
        print(f"  Y-axis: {y_min:.1f} to {y_max:.1f} km (margin: ±{y_margin:.1f} km)")
        
    else:
        # Use full slab data extent
        all_depths = np.concatenate([df['top_depth_km'], df['bottom_depth_km']])
        
        x_min = df['x_coord'].min()
        x_max = df['x_coord'].max()
        y_min = all_depths.min()  # Most negative (deepest)
        y_max = all_depths.max()  # Least negative (shallowest)
        
        # Add margin
        x_range = x_max - x_min
        y_range = y_max - y_min
        
        x_margin = x_range * margin_percent if x_range > 0 else 1.0
        y_margin = y_range * margin_percent if y_range > 0 else 10.0
        
        x_min -= x_margin
        x_max += x_margin
        y_min -= y_margin  # More negative (deeper)
        y_max += y_margin  # Less negative (shallower)
    
    return x_min, x_max, y_min, y_max

def create_plot(df, output_plot, earthquakes_df=None, dpi=300, title=None, crop_to_earthquakes=True):
    """
    Create a cross-section plot of the slab
    
    Parameters:
    - df: DataFrame with slab cross-section data
    - output_plot: Output filename for the plot
    - earthquakes_df: DataFrame with earthquake data (optional)
    - dpi: Resolution for the output image
    - title: Custom title for the plot
    - crop_to_earthquakes: Whether to crop the plot to earthquake extent
    """
    if len(df) == 0:
        print("No data to plot")
        return
    
    # Sort by x-coordinate for proper line plotting
    df_sorted = df.sort_values('x_coord')
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Calculate plot bounds
    if crop_to_earthquakes and earthquakes_df is not None and not earthquakes_df.empty:
        x_min, x_max, y_min, y_max = calculate_plot_bounds(df_sorted, earthquakes_df)
        
        # Filter slab data to plot bounds for cleaner visualization
        mask = ((df_sorted['x_coord'] >= x_min) & 
                (df_sorted['x_coord'] <= x_max))
        df_plot = df_sorted[mask]
    else:
        x_min, x_max, y_min, y_max = calculate_plot_bounds(df_sorted, None)
        df_plot = df_sorted
    
    # Plot top slab (line only, no points)
    if len(df_plot) > 0:
        ax.plot(df_plot['x_coord'], df_plot['top_depth_km'], 
                color='black', linewidth=1.5, label='Top of slab')
        
        # Plot bottom slab (line only, no points)
        ax.plot(df_plot['x_coord'], df_plot['bottom_depth_km'], 
                color='black', linewidth=1.5, label='Bottom of slab')
    
    # Plot earthquakes if provided
    if earthquakes_df is not None and not earthquakes_df.empty:
        # Sort earthquakes by depth to ensure smaller magnitude points are not hidden
        earthquakes_sorted = earthquakes_df.sort_values('magnitude')
        
        # Scatter plot of earthquakes
        scatter = ax.scatter(
            earthquakes_sorted['x_coord'], 
            earthquakes_sorted['depth_km'], 
            c=earthquakes_sorted['magnitude'], 
            cmap='viridis', 
            s=earthquakes_sorted['magnitude']*20,  # Size proportional to magnitude
            alpha=0.7,
            edgecolors='black',
            linewidth=0.5,
            label='Earthquakes'
        )
        
        # Add colorbar
        plt.colorbar(scatter, ax=ax, label='Magnitude')
    
    # Set plot bounds
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    
    # Set labels and formatting
    x_label = df['x_label'].iloc[0]  # Get the appropriate label
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel('Depth (km)', fontsize=12)
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Set title
    if title:
        ax.set_title(title, fontsize=14, fontweight='bold')
    else:
        if 'Latitude' in x_label:
            default_title = 'North-South Cross-Section'
        else:
            default_title = 'East-West Cross-Section'
        
        # Add cropping info to title if applicable
        if crop_to_earthquakes and earthquakes_df is not None and not earthquakes_df.empty:
            default_title += ' (Cropped to Earthquake Extent)'
        
        ax.set_title(default_title, fontsize=14, fontweight='bold')
    
    # Add legend
    ax.legend()
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(output_plot, dpi=dpi, bbox_inches='tight')
    print(f"Plot saved to: {output_plot} (DPI: {dpi})")
    
    # Close the figure to free memory
    plt.close()

def main():
    parser = argparse.ArgumentParser(
        description='Extract cross-sections from Slab 2.0 grid files with optional plotting and earthquake filtering',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('--input_dep', required=True,
                       help='Input depth grid file (e.g., van_slab2_dep_02.23.18.grd)')
    parser.add_argument('--input_thk', required=True,
                       help='Input thickness grid file (e.g., van_slab2_thk_02.23.18.grd)')
    parser.add_argument('--orientation', required=True, choices=['n_s', 'w_e'],
                       help='Cross-section orientation: n_s (north-south) or w_e (east-west)')
    parser.add_argument('--center_lat', type=float, required=True,
                       help='Center latitude for cross-section')
    parser.add_argument('--center_lon', type=float, required=True,
                       help='Center longitude for cross-section')
    parser.add_argument('--output_plot', 
                       help='Output PNG file for cross-section plot')
    parser.add_argument('--dpi', type=int, default=300,
                       help='DPI for output plot (default: 300)')
    parser.add_argument('--title', 
                       help='Custom title for the plot')
    parser.add_argument('--earthquakes',
                       help='CSV file with earthquake data')
    parser.add_argument('--buffer_degrees', type=float, default=2,
                       help='Buffer in degrees for filtering earthquakes (default: 2)')
    parser.add_argument('--no_crop', action='store_true',
                       help='Disable automatic cropping to earthquake extent')
    parser.add_argument('--crop_margin', type=float, default=0.15,
                       help='Margin percentage for cropping (default: 0.15 = 15%%)')
    parser.add_argument('--inspect', action='store_true',
                       help='Only inspect the grid files without extracting data')
    
    args = parser.parse_args()
    
    print(f"Loading depth data from: {args.input_dep}")
    dep_data = load_grd_file(args.input_dep)
    if dep_data is None:
        return
    
    print(f"Loading thickness data from: {args.input_thk}")
    thk_data = load_grd_file(args.input_thk)
    if thk_data is None:
        return
    
    # Inspect files if requested
    if args.inspect:
        inspect_grd_file(args.input_dep)
        inspect_grd_file(args.input_thk)
        return
    
    print(f"\nExtracting {args.orientation.upper()} cross-section at lat={args.center_lat}, lon={args.center_lon}")
    
    try:
        df = extract_cross_section(dep_data, thk_data, args.orientation, 
                                 args.center_lat, args.center_lon)
        
        print(f"\nExtracted {len(df)} valid data points")
        if len(df) > 0:
            #save_csv_files(df, args.output_top_csv, args.output_bottom_csv)
            
            # Process earthquakes if file provided
            earthquakes_df = None
            if args.earthquakes:
                # Read earthquake CSV
                earthquakes_df = pd.read_csv(args.earthquakes)
                
                # Adjust longitude if needed (subtract 180 if > 360)
                if earthquakes_df['lon360'].max() > 360:
                    earthquakes_df['lon360'] = earthquakes_df['lon360'] - 180
                
                # Filter earthquakes
                earthquakes_df = filter_earthquakes(
                    earthquakes_df, 
                    args.center_lon, 
                    args.center_lat, 
                    args.buffer_degrees
                )
                
                # Add x_coord for plotting based on cross-section orientation
                if args.orientation.lower() == 'n_s':
                    earthquakes_df['x_coord'] = earthquakes_df['latitude']
                else:
                    earthquakes_df['x_coord'] = earthquakes_df['lon360']
            
            # Create plot if requested
            if args.output_plot:
                crop_to_earthquakes = not args.no_crop and earthquakes_df is not None
                create_plot(df, args.output_plot, earthquakes_df, args.dpi, args.title, crop_to_earthquakes)
        else:
            print("No valid data points found in the cross-section")
        
    except Exception as e:
        print(f"Error during extraction: {e}")
        print("\nTry running with --inspect flag to see the file structure:")
        print(f"python {parser.prog} --input_dep {args.input_dep} --input_thk {args.input_thk} --inspect")
        return
    
    print("Cross-section extraction completed successfully!")

if __name__ == '__main__':
    main()