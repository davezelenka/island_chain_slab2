import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from io import StringIO

# Sample data (you can replace this with pd.read_csv('your_file.csv'))
sample_data = """Location,date_time,lon360,latitude,depth_km,magnitude,distance_to_surface_km,distance_to_center_km
Okinawa,2008-03-24T04:51:09.400Z,126.47,25.35,0.00,3.60,28.96,79.57
Okinawa,2008-04-02T16:22:26.070Z,126.56,24.30,-1.00,4.30,9.30,53.79
Okinawa,2010-02-28T02:51:41.590Z,128.48,25.95,-1.30,5.00,10.81,65.14
Okinawa,2008-09-19T10:01:43.870Z,128.17,28.30,-1.30,3.70,105.04,156.19
Okinawa,2008-07-13T07:45:39.080Z,127.74,24.67,-1.40,4.30,14.35,62.58
Okinawa,2008-07-13T07:45:39.080Z,127.74,24.67,30.00,4.30,14.35,62.58
Okinawa,2008-07-13T07:45:39.080Z,127.74,24.67,50.00,4.30,14.35,62.58
Mariana,2008-03-24T04:51:09.400Z,145.47,15.35,50.00,4.60,15.96,45.57
Mariana,2008-04-02T16:22:26.070Z,145.56,14.30,75.00,5.30,25.30,35.79
Mariana,2010-02-28T02:51:41.590Z,146.48,15.95,100.00,6.00,35.81,25.14
Tonga,2008-03-24T04:51:09.400Z,185.47,-15.35,200.00,5.60,45.96,15.57
Tonga,2008-04-02T16:22:26.070Z,185.56,-14.30,250.00,6.30,55.30,5.79
Tonga,2010-02-28T02:51:41.590Z,186.48,-15.95,300.00,7.00,65.81,-5.14"""

# Function to remove outliers using IQR method
def remove_outliers_iqr(data, multiplier=1.5):
    """
    Remove outliers using the Interquartile Range (IQR) method
    multiplier: typically 1.5 (standard) or 3.0 (more conservative)
    """
    Q1 = np.percentile(data, 25)
    Q3 = np.percentile(data, 75)
    IQR = Q3 - Q1
    lower_bound = Q1 - multiplier * IQR
    upper_bound = Q3 + multiplier * IQR
    return data[(data >= lower_bound) & (data <= upper_bound)]

# Function to remove outliers using Z-score method
def remove_outliers_zscore(data, threshold=3):
    """
    Remove outliers using Z-score method
    threshold: typically 2, 2.5, or 3 standard deviations
    """
    z_scores = np.abs((data - np.mean(data)) / np.std(data))
    return data[z_scores < threshold]

# Function to remove outliers using percentile method
def remove_outliers_percentile(data, lower_percentile=5, upper_percentile=95):
    """
    Remove outliers by trimming extreme percentiles
    Common ranges: 5-95%, 2.5-97.5%, 1-99%
    """
    lower_bound = np.percentile(data, lower_percentile)
    upper_bound = np.percentile(data, upper_percentile)
    return data[(data >= lower_bound) & (data <= upper_bound)]

# Load data
df = pd.read_csv('earthquake_data_full.csv')

# FILTER BY DEPTH: Remove earthquakes with depth < 25 km
MAX_DEPTH = -25  # Minimum depth in km
original_count = len(df)
df_filtered = df[df['depth_km'] <= MAX_DEPTH].copy()
filtered_count = len(df_filtered)
removed_shallow = original_count - filtered_count

print(f"Depth filtering applied: Removed {removed_shallow} earthquakes with depth < {MAX_DEPTH} km")
print(f"Original dataset: {original_count} earthquakes")
print(f"Filtered dataset: {filtered_count} earthquakes")

# Use filtered dataframe for the rest of the analysis
df = df_filtered

# CHOOSE YOUR OUTLIER REMOVAL METHOD HERE:
OUTLIER_METHOD = 'iqr'  # Options: 'iqr', 'zscore', 'percentile', 'none'
IQR_MULTIPLIER = 1.5    # 1.5 = standard, 3.0 = more conservative
ZSCORE_THRESHOLD = 3    # 2, 2.5, or 3 standard deviations
PERCENTILE_RANGE = (5, 95)  # (lower, upper) percentiles to keep

# Get unique locations (after filtering)
locations = df['Location'].unique()
num_locations = len(locations)

# Check if we have data for all locations after filtering
if num_locations == 0:
    print(f"Warning: No data remaining after depth filter (>= {MAX_DEPTH} km)")
    exit()

# Create figure and axis
fig, ax = plt.subplots(figsize=(15, 8))

# Prepare data for box plots
box_data_surface = []
box_data_center = []
box_positions_surface = []
box_positions_center = []
outlier_counts = {'surface': {}, 'center': {}}

for i, location in enumerate(locations):
    location_data = df[df['Location'] == location]
    
    # Skip locations with no data after filtering
    if len(location_data) == 0:
        continue
    
    # Extract distance data for this location
    surface_distances = location_data['distance_to_surface_km'].values
    center_distances = location_data['distance_to_center_km'].values
    
    # Store original counts
    original_surface_count = len(surface_distances)
    original_center_count = len(center_distances)
    
    # Apply outlier removal based on selected method
    if OUTLIER_METHOD == 'iqr':
        surface_distances_clean = remove_outliers_iqr(surface_distances, IQR_MULTIPLIER)
        center_distances_clean = remove_outliers_iqr(center_distances, IQR_MULTIPLIER)
    elif OUTLIER_METHOD == 'zscore':
        surface_distances_clean = remove_outliers_zscore(surface_distances, ZSCORE_THRESHOLD)
        center_distances_clean = remove_outliers_zscore(center_distances, ZSCORE_THRESHOLD)
    elif OUTLIER_METHOD == 'percentile':
        surface_distances_clean = remove_outliers_percentile(surface_distances, 
                                                           PERCENTILE_RANGE[0], PERCENTILE_RANGE[1])
        center_distances_clean = remove_outliers_percentile(center_distances, 
                                                          PERCENTILE_RANGE[0], PERCENTILE_RANGE[1])
    else:  # 'none' - no outlier removal
        surface_distances_clean = surface_distances
        center_distances_clean = center_distances
    
    # Store outlier counts for reporting
    outlier_counts['surface'][location] = original_surface_count - len(surface_distances_clean)
    outlier_counts['center'][location] = original_center_count - len(center_distances_clean)
    
    box_data_surface.append(surface_distances_clean)
    box_data_center.append(center_distances_clean)
    
    # Position boxes side by side for each location
    box_positions_surface.append(i * 2 + 0.8)  # Slightly left of center
    box_positions_center.append(i * 2 + 1.2)   # Slightly right of center

# Create box plots with outliers hidden (showfliers=False)
box1 = ax.boxplot(box_data_surface, positions=box_positions_surface, 
                  widths=0.3, patch_artist=True, showfliers=False,  # This hides outliers
                  boxprops=dict(facecolor='white', hatch='///', color='black'),
                  medianprops=dict(color='black', linewidth=2),
                  whiskerprops=dict(color='black'),
                  capprops=dict(color='black'))

box2 = ax.boxplot(box_data_center, positions=box_positions_center, 
                  widths=0.3, patch_artist=True, showfliers=False,  # This hides outliers
                  boxprops=dict(facecolor='black', hatch='\\\\\\', color='white'),
                  medianprops=dict(color='black', linewidth=2),
                  whiskerprops=dict(color='black'),
                  capprops=dict(color='black'))

# Add horizontal line at distance = 0
ax.axhline(y=0, color='black', linestyle='-', linewidth=1, alpha=0.7)

# Customize the plot
ax.set_xlabel('Location', fontsize=14, fontweight='bold')
ax.set_ylabel('Distance (km)', fontsize=14, fontweight='bold')

# Update title to reflect both depth filtering and outlier removal
title = f'Earthquake Distance Distributions by Location \nDistance to Slab Surface vs. Distance to Slab Center'
    
ax.set_title(title, fontsize=16, fontweight='bold', pad=20)

# Set x-axis labels
ax.set_xticks([i * 2 + 1 for i in range(len(locations))])
ax.set_xticklabels(locations, rotation=45, ha='right', fontsize=14)

# Add grid for better readability
ax.grid(True, alpha=0.3, axis='y')

# Create legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='white', edgecolor='black', hatch='///', label='Distance to Surface'),
                  Patch(facecolor='black', edgecolor='white', hatch='\\\\\\', label='Distance to Center')]
ax.legend(handles=legend_elements, loc='upper right', fontsize=12)

# Adjust layout to prevent label cutoff
plt.tight_layout()

# Add statistics including depth filtering and outlier information
total_outliers_surface = sum(outlier_counts['surface'].values())
total_outliers_center = sum(outlier_counts['center'].values())

stats_text = f"Earthquakes: {filtered_count}\nOutliers removed (Surface): {total_outliers_surface}\nOutliers removed (Center): {total_outliers_center}\nMethod: {OUTLIER_METHOD.upper()}"

ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# Show the plot
#plt.show()

# Optional: Save the plot
plt.savefig('charts/earthquake_distance_boxplots_filtered.png', dpi=600, bbox_inches='tight')

print("\nBox plot visualization complete!")
print(f"Depth filter: <= {MAX_DEPTH} km (removed {removed_shallow} shallow earthquakes)")
print(f"Outlier removal method: {OUTLIER_METHOD.upper()}")
if OUTLIER_METHOD != 'none':
    print(f"Total outliers removed - Surface: {total_outliers_surface}, Center: {total_outliers_center}")
    print("\nOutliers removed by location:")
    for location in locations:
        if location in outlier_counts['surface']:
            print(f"  {location}: Surface={outlier_counts['surface'][location]}, Center={outlier_counts['center'][location]}")

print("\nSummary statistics (after depth filtering and outlier removal):")
for location in locations:
    location_data = df[df['Location'] == location]
    if len(location_data) == 0:
        continue
        
    surface_data = location_data['distance_to_surface_km'].values
    center_data = location_data['distance_to_center_km'].values
    
    # Apply same outlier removal for statistics
    if OUTLIER_METHOD == 'iqr':
        surface_clean = remove_outliers_iqr(surface_data, IQR_MULTIPLIER)
        center_clean = remove_outliers_iqr(center_data, IQR_MULTIPLIER)
    elif OUTLIER_METHOD == 'zscore':
        surface_clean = remove_outliers_zscore(surface_data, ZSCORE_THRESHOLD)
        center_clean = remove_outliers_zscore(center_data, ZSCORE_THRESHOLD)
    elif OUTLIER_METHOD == 'percentile':
        surface_clean = remove_outliers_percentile(surface_data, PERCENTILE_RANGE[0], PERCENTILE_RANGE[1])
        center_clean = remove_outliers_percentile(center_data, PERCENTILE_RANGE[0], PERCENTILE_RANGE[1])
    else:
        surface_clean = surface_data
        center_clean = center_data
    
    print(f"\n{location}:")
    print(f"  Count after filtering: {len(location_data)} earthquakes")
    print(f"  Count after cleaning: Surface={len(surface_clean)}, Center={len(center_clean)}")
    print(f"  Depth range: {location_data['depth_km'].min():.1f} - {location_data['depth_km'].max():.1f} km")
    print(f"  Surface distance - Mean: {surface_clean.mean():.2f} km, Std: {surface_clean.std():.2f} km")
    print(f"  Center distance - Mean: {center_clean.mean():.2f} km, Std: {center_clean.std():.2f} km")