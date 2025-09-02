#!/usr/bin/env python3
"""
Publication-ready black and white histogram plotter for earthquake distance data.
Creates overlapping histograms with statistical markers and comprehensive legends.
Now includes Wilcoxon signed-rank tests, Cohen's d, and confidence intervals.
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
from pathlib import Path
from scipy import stats
from scipy.stats import wilcoxon
import warnings

def load_and_validate_data(csv_path,max_depth,min_depth):
    """Load CSV data and validate required columns."""
    try:
        df = pd.read_csv(csv_path)
        required_cols = ['distance_to_surface_km', 'distance_to_center_km', 'latitude', 'lon360', 'depth_km']
        
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
        
        # Remove any rows with NaN values in required columns
        df_clean = df[required_cols].dropna()
        df_clean = df_clean[df_clean['depth_km'] <= max_depth]
        df_clean = df_clean[df_clean['depth_km'] > min_depth]
        
        if len(df_clean) == 0:
            raise ValueError("No valid data rows found after removing NaN values")
        
        return df_clean
    
    except FileNotFoundError:
        raise FileNotFoundError(f"CSV file not found: {csv_path}")
    except Exception as e:
        raise Exception(f"Error loading data: {str(e)}")

def calculate_statistics(data):
    """Calculate median, mean, and standard deviation."""
    return {
        'median': np.median(data),
        'mean': np.mean(data),
        'std': np.std(data),
    }

def cohens_d(x, y):
    """Calculate Cohen's d effect size for two samples."""
    # Calculate the size of samples
    nx = len(x)
    ny = len(y)
    
    # Calculate the variance of the samples
    dof = nx + ny - 2
    pooled_std = np.sqrt(((nx-1)*np.var(x, ddof=1) + (ny-1)*np.var(y, ddof=1)) / dof)
    
    # Calculate Cohen's d
    d = (np.mean(x) - np.mean(y)) / pooled_std
    return d

def wilcoxon_test_analysis(surface_data, center_data):
    """
    Perform Wilcoxon signed-rank test comparing paired distances.
    Tests if earthquakes are systematically closer to interface vs center.
    """
    # Calculate differences: negative means closer to interface
    differences = surface_data - center_data
    
    # Remove zero differences (ties)
    non_zero_diffs = differences[differences != 0]
    
    if len(non_zero_diffs) < 5:
        return {
            'statistic': np.nan,
            'p_value': np.nan,
            'effect_size': np.nan,
            'n_pairs': len(differences),
            'n_zero_diffs': np.sum(differences == 0),
            'interpretation': 'Insufficient non-zero differences for test'
        }
    
    # Perform Wilcoxon signed-rank test
    # alternative='less' tests if surface distances are systematically smaller (closer to interface)
    try:
        statistic, p_value = wilcoxon(non_zero_diffs, alternative='less')
    except ValueError as e:
        return {
            'statistic': np.nan,
            'p_value': np.nan,
            'effect_size': np.nan,
            'n_pairs': len(differences),
            'n_zero_diffs': np.sum(differences == 0),
            'interpretation': f'Test failed: {str(e)}'
        }
    
    # Calculate effect size (r = Z / sqrt(N))
    # For large samples, we can approximate using normal distribution
    n = len(non_zero_diffs)
    z_score = stats.norm.ppf(1 - p_value) if p_value < 0.5 else -stats.norm.ppf(p_value)
    effect_size_r = abs(z_score) / np.sqrt(n)
    
    # Interpretation of effect size (Cohen's conventions for r)
    if effect_size_r < 0.1:
        effect_interpretation = 'negligible'
    elif effect_size_r < 0.3:
        effect_interpretation = 'small'
    elif effect_size_r < 0.5:
        effect_interpretation = 'medium'
    else:
        effect_interpretation = 'large'
    
    return {
        'statistic': statistic,
        'p_value': p_value,
        'effect_size': effect_size_r,
        'effect_interpretation': effect_interpretation,
        'n_pairs': len(differences),
        'n_zero_diffs': np.sum(differences == 0),
        'median_difference': np.median(differences),
        'mean_difference': np.mean(differences)
    }

def calculate_proportion_confidence_interval(successes, total, confidence_level=0.95):
    """
    Calculate confidence interval for a proportion using Wilson score interval.
    More accurate than normal approximation, especially for extreme proportions.
    """
    if total == 0:
        return 0.0, 0.0, 0.0
    
    p = successes / total
    alpha = 1 - confidence_level
    z = stats.norm.ppf(1 - alpha/2)
    
    # Wilson score interval
    denominator = 1 + z**2/total
    centre_adjusted_probability = (p + z**2/(2*total)) / denominator
    adjusted_standard_deviation = np.sqrt((p*(1-p) + z**2/(4*total)) / total) / denominator
    
    lower_bound = centre_adjusted_probability - z * adjusted_standard_deviation
    upper_bound = centre_adjusted_probability + z * adjusted_standard_deviation
    
    return p, max(0, lower_bound), min(1, upper_bound)

def comprehensive_statistical_analysis(df):
    """
    Perform comprehensive statistical analysis of earthquake distance data.
    """
    surface_data = df['distance_to_surface_km'].values
    center_data = df['distance_to_center_km'].values
 
    print("=" * 80)
    print("COMPREHENSIVE STATISTICAL ANALYSIS")
    print("=" * 80)

    # Basic descriptive statistics
    print("\n1. DESCRIPTIVE STATISTICS")
    print("-" * 40)
    surface_stats = calculate_statistics(surface_data)
    center_stats = calculate_statistics(center_data)
    
    print(f"Interface Distances:")
    print(f"  Median: {surface_stats['median']:.2f} km")
    print(f"  Mean: {surface_stats['mean']:.2f} km") 
    print(f"  Std Dev: {surface_stats['std']:.2f} km")
    
    print(f"\nCenter Distances:")
    print(f"  Median: {center_stats['median']:.2f} km")
    print(f"  Mean: {center_stats['mean']:.2f} km")
    print(f"  Std Dev: {center_stats['std']:.2f} km")
    
    # Proportion analysis with confidence intervals
    print("\n2. PROPORTION ANALYSIS")
    print("-" * 40)
    
    # Count earthquakes closer to interface (absolute distance comparison)
    closer_to_interface = np.abs(surface_data) < np.abs(center_data)
    n_closer = np.sum(closer_to_interface)
    total = len(surface_data)
    
    prop, ci_lower, ci_upper = calculate_proportion_confidence_interval(n_closer, total)
    
    print(f"Earthquakes closer to interface: {n_closer}/{total} ({prop:.1%})")
    print(f"95% Confidence Interval: [{ci_lower:.1%}, {ci_upper:.1%}]")
    
    # Alternative analysis: earthquakes above vs below center
    above_center = center_data > 0  # positive means above center
    n_above = np.sum(above_center)
    prop_above, ci_lower_above, ci_upper_above = calculate_proportion_confidence_interval(n_above, total)
    
    print(f"Earthquakes above slab center: {n_above}/{total} ({prop_above:.1%})")
    print(f"95% Confidence Interval: [{ci_lower_above:.1%}, {ci_upper_above:.1%}]")
    
    # Wilcoxon signed-rank test
    print("\n3. WILCOXON SIGNED-RANK TEST")
    print("-" * 40)
    wilcox_results = wilcoxon_test_analysis(surface_data, center_data)
    
    if not np.isnan(wilcox_results['p_value']):
        print(f"Test Statistic: {wilcox_results['statistic']:.0f}")
        print(f"P-value: {wilcox_results['p_value']:.2e}")
        print(f"Effect Size (r): {wilcox_results['effect_size']:.3f} ({wilcox_results['effect_interpretation']})")
        print(f"Median difference (Interface - Center): {wilcox_results['median_difference']:.2f} km")
        print(f"Mean difference (Interface - Center): {wilcox_results['mean_difference']:.2f} km")
        
        # Statistical significance interpretation
        if wilcox_results['p_value'] < 0.001:
            significance = "highly significant (p < 0.001)"
        elif wilcox_results['p_value'] < 0.01:
            significance = "very significant (p < 0.01)" 
        elif wilcox_results['p_value'] < 0.05:
            significance = "significant (p < 0.05)"
        else:
            significance = "not significant (p ≥ 0.05)"
        
        print(f"Result: {significance}")
        
        if wilcox_results['p_value'] < 0.05:
            print("CONCLUSION: Earthquakes are systematically closer to slab interface than center")
        else:
            print("CONCLUSION: No significant difference in distance to interface vs center")
    else:
        print(f"Test could not be performed: {wilcox_results['interpretation']}")
    
    print(f"Sample size: {wilcox_results['n_pairs']} paired observations")
    if wilcox_results['n_zero_diffs'] > 0:
        print(f"Tied differences removed: {wilcox_results['n_zero_diffs']}")
    
    # Cohen's d for unpaired comparison (supplementary)
    print("\n4. COHEN'S D EFFECT SIZE (supplementary)")
    print("-" * 40)
    cohens_d_value = cohens_d(surface_data, center_data)
    print(f"Cohen's d: {cohens_d_value:.3f}")
    
    if abs(cohens_d_value) < 0.2:
        d_interpretation = "negligible"
    elif abs(cohens_d_value) < 0.5:
        d_interpretation = "small"
    elif abs(cohens_d_value) < 0.8:
        d_interpretation = "medium"
    else:
        d_interpretation = "large"
    
    print(f"Interpretation: {d_interpretation} effect")
    
    # Summary for publication
    print("\n5. PUBLICATION SUMMARY")
    print("-" * 40)
    print(f"• {prop:.0%} of earthquakes occur closer to slab interface than center")
    print(f"• 95% CI: [{ci_lower:.0%}, {ci_upper:.0%}]")
    if not np.isnan(wilcox_results['p_value']):
        print(f"• Wilcoxon test: p = {wilcox_results['p_value']:.2e}, effect size r = {wilcox_results['effect_size']:.3f}")
        print(f"• Statistical significance: {significance}")
    
    return {
        'proportion_interface_proximal': prop,
        'proportion_ci_lower': ci_lower,
        'proportion_ci_upper': ci_upper,
        'wilcoxon_p_value': wilcox_results['p_value'],
        'wilcoxon_effect_size': wilcox_results['effect_size'],
        'cohens_d': cohens_d_value,
        'n_total': total,
        'n_closer_interface': n_closer
    }

def determine_bin_parameters(surface_data, center_data, n_bins=None):
    """Determine optimal bin parameters based on data range."""
    all_data = np.concatenate([surface_data, center_data])
    data_min, data_max = np.min(all_data), np.max(all_data)
    data_range = data_max - data_min
    
    if n_bins is None:
        # Use Freedman-Diaconis rule for bin width
        q75, q25 = np.percentile(all_data, [75, 25])
        iqr = q75 - q25
        bin_width = 2 * iqr / (len(all_data) ** (1/3))
        
        if bin_width == 0:  # Handle edge case
            bin_width = data_range / 50
        
        n_bins = max(10, min(100, int(data_range / bin_width)))
    
    # Extend range slightly for better visualization
    margin = data_range * 0.05
    bin_range = (data_min - margin, data_max + margin)
    
    return n_bins, bin_range

def create_histogram_plot(df, title, x_axis_label, n_bins=None):
    """Create the main histogram plot with all required elements."""
    
    # Extract data
    surface_data = df['distance_to_surface_km'].values
    center_data = df['distance_to_center_km'].values
    latitude_data = df['latitude'].values
    longitude_data = df['lon360'].values
    
    # Perform comprehensive statistical analysis
    stats_results = comprehensive_statistical_analysis(df)
    
    # Calculate statistics for plot labels
    surface_stats = calculate_statistics(surface_data)
    center_stats = calculate_statistics(center_data)
    
    # Calculate center coordinates (max-min)/2
    lat_range = np.max(latitude_data) - np.min(latitude_data)
    lon_range = np.max(longitude_data) - np.min(longitude_data)
    center_lat = np.min(latitude_data) + lat_range / 2
    center_lon = np.min(longitude_data) + lon_range / 2
    
    # Basic summary for console output
    print(f"\n--- PLOT SUMMARY ---")
    print(f"Center Latitude: {center_lat:.3f}°")
    print(f"Center Longitude: {center_lon:.3f}°")
    n = len(center_data)
    ngt0 = len(center_data[center_data > 0])
    print(f"N: {n}")
    print(f"N > Mid-slab: {ngt0}")
    print(f"Surface Median: {surface_stats['median']:.2f} km")
    print(f"Mid-slab Median: {center_stats['median']:.2f} km")
    
    # Determine optimal binning
    n_bins_final, bin_range = determine_bin_parameters(surface_data, center_data, n_bins)
    bins = np.linspace(bin_range[0], bin_range[1], n_bins_final + 1)
    
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Create overlapping histograms
    n_surface, _, patches_surface = ax.hist(
        surface_data, bins=bins, alpha=0.7, 
        color='white', edgecolor='black', linewidth=1,
        hatch='///', label='Slab Interface (Surface)'
    )
    
    n_center, _, patches_center = ax.hist(
        center_data, bins=bins, alpha=0.7,
        color='lightgray', edgecolor='black', linewidth=1,
        hatch='\\\\\\', label='Mid-slab (center)'
    )
    
    # Add vertical markers
    ax.axvline(x=0, color='black', linestyle='-', linewidth=2, label='x=0')
    ax.axvline(x=surface_stats['median'], color='black', linestyle='--', 
               linewidth=2, label=f"Interface Median ({surface_stats['median']:.1f})")
    ax.axvline(x=center_stats['median'], color='black', linestyle=':', 
               linewidth=2, label=f"Center Median ({center_stats['median']:.1f})")
    
    # Set labels and title (increased font sizes by 1.5x)
    ax.set_xlabel(x_axis_label, fontsize=21, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=21, fontweight='bold')
    ax.set_title(title, fontsize=24, fontweight='bold', pad=20)
    
    # Increase tick label font size
    ax.tick_params(axis='both', which='major', labelsize=15)
    
    # Add location info at top left - no box, no color
    location_text = f"Center Latitude: {center_lat:.2f}°\nCenter Longitude: {center_lon:.2f}°"
    ax.text(0.02, 0.98, location_text, transform=ax.transAxes,
            fontsize=16, ha='left', va='top')
    
    # Add n-value at top right - no box
    n_total = len(df)
    ax.text(0.98, 0.98, f'n = {n_total}', transform=ax.transAxes,
            fontsize=18, fontweight='bold', ha='right', va='top')
    
    # Add statistical summary box
    prop = stats_results['proportion_interface_proximal']
    ci_lower = stats_results['proportion_ci_lower'] 
    ci_upper = stats_results['proportion_ci_upper']
    p_val = stats_results['wilcoxon_p_value']
    
    # Format p-value for display
    if p_val < 0.001:
        p_str = "p < 0.001"
    elif p_val < 0.01:
        p_str = f"p < 0.01"
    elif p_val < 0.05:
        p_str = f"p < 0.05"
    else:
        p_str = f"p = {p_val:.3f}"
    
    stats_text = (
        f"Statistical Summary:\n"
        f"• {prop:.0%} interface-proximal\n"
        f"• 95% CI: [{ci_lower:.0%}-{ci_upper:.0%}]\n"
        f"• {p_str} (Wilcoxon)\n"
        f"• n = {stats_results['n_total']:,} pairs"
    )
    
    ax.text(0.02, 0.25, stats_text, transform=ax.transAxes,
            fontsize=12, ha='left', va='top',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white', 
                     edgecolor='black', alpha=0.9, linewidth=1))
    
    # Create legends with statistics
    # Left legend - Surface of Slab
    surface_legend_text = (
        f"Slab Interface (surface)\n"
        f"Median: {surface_stats['median']:.1f} km\n"
        f"Mean: {surface_stats['mean']:.1f} km\n"
        f"Std Dev: {surface_stats['std']:.1f} km"
    )
    
    ax.text(0.02, 0.85, surface_legend_text, transform=ax.transAxes,
            fontsize=16, ha='left', va='top',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', 
                     edgecolor='black', alpha=0.9))
    
    # Right legend - Center of Slab (moved up to mirror left side)
    center_legend_text = (
        f"Mid-slab (center)\n"
        f"Median: {center_stats['median']:.1f} km\n"
        f"Mean: {center_stats['mean']:.1f} km\n"
        f"Std Dev: {center_stats['std']:.1f} km"
    )
    
    ax.text(0.98, 0.85, center_legend_text, transform=ax.transAxes,
            fontsize=16, ha='right', va='top',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgray', 
                     edgecolor='black', alpha=0.9))
    
    # Improve layout
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)
    
    # Make sure all elements are visible
    plt.tight_layout()
    
    return fig, ax, stats_results

def main():
    """Main function to parse arguments and create the plot."""
    parser = argparse.ArgumentParser(
        description="Create publication-ready histogram of earthquake distance data with statistical analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python eq_plot_v2.py --csv distances_surface.csv --title "Philippines" --x-axis "Earthquake distances (km)"
    python eq_plot_v2.py --csv distances_surface.csv --title "Philippines" --x-axis "Earthquake distances (km)" --bins 30
        """
    )
    
    parser.add_argument('--csv', required=True, type=str,
                       help='Path to CSV file containing earthquake data')
    parser.add_argument('--title', required=True, type=str,
                       help='Title for the histogram plot')
    parser.add_argument('--x-axis', required=True, type=str,
                       help='Label for the x-axis')
    parser.add_argument('--bins', type=int, default=None,
                       help='Number of bins for histogram (optional, defaults to auto)')
    parser.add_argument('--output', type=str, default=None,
                       help='Output filename (optional, defaults to display)')
    parser.add_argument('--dpi', type=int, default=300,
                       help='Resolution for saved figure (default: 300)')
    parser.add_argument('--max', type=int, default=-70,
                       help='Max depth (default: -70)')
    parser.add_argument('--min', type=int, default=-900,
                       help='Min depth (default: -900)')
    parser.add_argument('--stats-only', action='store_true',
                       help='Only run statistical analysis, do not create plot')
    
    args = parser.parse_args()
    print(f"\nTitle: {args.title}")
    try:
        # Load and validate data
        df = load_and_validate_data(args.csv, args.max, args.min)
        print(f"Successfully loaded {len(df)} data points")
        
        if args.stats_only:
            # Just run statistical analysis
            comprehensive_statistical_analysis(df)
        else:
            # Create the plot with statistical analysis
            print("Creating histogram plot with statistical analysis...")
            fig, ax, stats_results = create_histogram_plot(df, args.title, args.x_axis, args.bins)
            
            # Save or display the plot
            if args.output:
                print(f"\nSaving plot to: {args.output}")
                plt.savefig(args.output, dpi=args.dpi, bbox_inches='tight', 
                           facecolor='white', edgecolor='none')
                print(f"Plot saved successfully with {args.dpi} DPI")
            else:
                print("\nDisplaying plot...")
                plt.show()
            
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()