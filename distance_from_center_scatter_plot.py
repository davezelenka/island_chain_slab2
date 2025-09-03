import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

def create_scatter_plot(args):
    # Ensure output directory exists
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    # Read the CSV file
    df = pd.read_csv(args.csv)

    # Filter for specified location
    df_filtered = df[df['Location'] == args.location].copy()

    # Create the figure and axis
    plt.figure(figsize=(10, 6), dpi=args.dpi)

    # Create scatter plot
    plt.scatter(df_filtered['distance_to_center_km'], df_filtered['depth_km'], 
                color='black', alpha=0.7, s=20, edgecolors='none')

    # Customize the plot
    plt.title(args.title, fontsize=16, fontweight='bold')
    plt.xlabel('Distance from Slab Center (km)', fontsize=12)
    plt.ylabel('Depth (km)', fontsize=12)

    # Set x and y limits based on data with some padding
    x_min = df_filtered['distance_to_center_km'].min()
    x_max = df_filtered['distance_to_center_km'].max()
    y_min = df_filtered['depth_km'].min()
    y_max = df_filtered['depth_km'].max()

    # Add 5% padding to x and y limits
    x_padding = (x_max - x_min) * 0.05
    y_padding = (y_max - y_min) * 0.05

    plt.xlim(x_min - x_padding, x_max + x_padding)
    plt.ylim(y_min - y_padding, y_max + y_padding)

    # Add grid for readability
    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

    # Add dashed vertical line at x = 0
    plt.axvline(x=0, color='red', linestyle='--', linewidth=1, alpha=0.7)
    
    # Tight layout to use space efficiently
    plt.tight_layout()

    # Save the plot
    print(f"Saving plot to: {os.path.abspath(args.output)}")
    plt.savefig(args.output, dpi=args.dpi, bbox_inches='tight')
    plt.close()  # Close the figure to free memory

    # Print some basic statistics
    print(f"{args.location} Earthquake Data Summary:")
    print(f"Total earthquakes: {len(df_filtered)}")
    print(f"Depth range: {df_filtered['depth_km'].min():.2f} to {df_filtered['depth_km'].max():.2f} km")
    print(f"Distance from slab center range: {df_filtered['distance_to_center_km'].min():.2f} to {df_filtered['distance_to_center_km'].max():.2f} km")

    # Additional diagnostic information
    print("\nData Distribution:")
    print(df_filtered[['depth_km', 'distance_to_center_km']].describe())

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description='Create earthquake scatter plot')
    
    # Add arguments with default values
    parser.add_argument('--title', type=str, default='Earthquake Scatter Plot', 
                        help='Title of the plot')
    parser.add_argument('--location', type=str, required=True, 
                        help='Location to filter earthquakes')
    parser.add_argument('--csv', type=str, default='earthquake_data_full.csv', 
                        help='Path to input CSV file')
    parser.add_argument('--output', type=str, default='charts/earthquake_scatter.png', 
                        help='Output file path for the plot')
    parser.add_argument('--dpi', type=int, default=600, 
                        help='DPI (resolution) of the output image')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Create the scatter plot
    create_scatter_plot(args)

if __name__ == '__main__':
    main()