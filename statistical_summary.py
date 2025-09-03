#!/usr/bin/env python3
"""
Create publication-ready statistical summary charts for multiple regions/depths.
Generates proportion bar charts with confidence intervals and significance markers.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path

def create_statistical_summary_chart(data_dict, chart_type='proportion'):
    """
    Create statistical summary chart from multiple region analyses.
    
    Parameters:
    -----------
    data_dict : dict
        Dictionary with region names as keys, statistical results as values
        Each value should have: proportion_interface_proximal, proportion_ci_lower, 
        proportion_ci_upper, wilcoxon_p_value, n_total
    chart_type : str
        'proportion' for percentage chart, 'forest' for effect size plot
    """
    
    # Extract data for plotting
    regions = [1,1]
    proportions = [0.5,0.5]
    ci_lowers = [0.5,0.5]
    ci_uppers = [0.5,0.5]
    p_values = [0.5,0.5]
    sample_sizes = [55,55]
    depth_categories = ["a","b"]
    
    for region_depth, stats in data_dict.items():
        regions.append(region_depth)
        proportions.append(stats['proportion_interface_proximal'] * 100)  # Convert to percentage
        ci_lowers.append(stats['proportion_ci_lower'] * 100)
        ci_uppers.append(stats['proportion_ci_upper'] * 100)
        p_values.append(stats['wilcoxon_p_value'])
        sample_sizes.append(stats['n_total'])
        
        # Determine depth category from region name
        if '25-70' in region_depth or 'shallow' in region_depth.lower():
            depth_categories.append('Shallow')
        else:
            depth_categories.append('Deep')
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(16, 10))
    
    if chart_type == 'proportion':
        return create_proportion_chart(ax, regions, proportions, ci_lowers, ci_uppers, 
                                     p_values, sample_sizes, depth_categories)
    else:
        # Future: implement forest plot
        raise NotImplementedError("Forest plot not yet implemented")

def create_proportion_chart(ax, regions, proportions, ci_lowers, ci_uppers, 
                          p_values, sample_sizes, depth_categories):
    """Create proportion bar chart with confidence intervals."""
    
    # Set up positions
    n_regions = len(regions)
    x_positions = np.arange(n_regions)
    
    # Color coding for depth
    colors = ['skyblue' if cat == 'Shallow' else 'lightcoral' for cat in depth_categories]
    edge_colors = ['darkblue' if cat == 'Shallow' else 'darkred' for cat in depth_categories]
    
    # Create bars
    bars = ax.bar(x_positions, proportions, color=colors, edgecolor=edge_colors, 
                  linewidth=1.5, alpha=0.8)
    
    # Add error bars (confidence intervals)
    error_lower = [prop - ci_low for prop, ci_low in zip(proportions, ci_lowers)]
    error_upper = [ci_up - prop for prop, ci_up in zip(proportions, ci_uppers)]
    
    ax.errorbar(x_positions, proportions, yerr=[error_lower, error_upper], 
                fmt='none', color='black', capsize=4, capthick=1.5, linewidth=1.5)
    
    # Add significance markers above bars
    max_height = max([prop + (ci_up - prop) for prop, ci_up in zip(proportions, ci_uppers)])
    for i, (p_val, prop, ci_up) in enumerate(zip(p_values, proportions, ci_uppers)):
        marker_height = ci_up + (max_height * 0.02)
        
        if p_val < 0.001:
            marker = '***'
        elif p_val < 0.01:
            marker = '**'
        elif p_val < 0.05:
            marker = '*'
        else:
            marker = 'ns'
        
        ax.text(i, marker_height, marker, ha='center', va='bottom', 
                fontsize=12, fontweight='bold')
    
    # Add horizontal reference line at 50% (random expectation)
    ax.axhline(y=50, color='gray', linestyle='--', linewidth=2, alpha=0.7, 
               label='Random expectation (50%)')
    
    # Customize the plot
    ax.set_xlabel('Subduction Zone - Depth Range', fontsize=14, fontweight='bold')
    ax.set_ylabel('Earthquakes Closer to Interface (%)', fontsize=14, fontweight='bold')
    ax.set_title('Interface-Proximal Earthquake Proportions Across Subduction Zones\n' +
                 'Error bars show 95% confidence intervals', 
                 fontsize=16, fontweight='bold', pad=20)
    
    # Set x-axis labels
    ax.set_xticks(x_positions)
    ax.set_xticklabels(regions, rotation=45, ha='right', fontsize=10)
    
    # Set y-axis limits
    ax.set_ylim(0, max_height * 1.1)
    
    # Add grid
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_axisbelow(True)
    
    # Create legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='skyblue', edgecolor='darkblue', label='Shallow (25-70 km)'),
        Patch(facecolor='lightcoral', edgecolor='darkred', label='Deep (>70 km)'),
        plt.Line2D([0], [0], color='gray', linestyle='--', label='Random expectation'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=12)
    
    # Add significance legend
    sig_text = ('Significance levels:\n'
                '*** p < 0.001\n'
                '**  p < 0.01\n'
                '*   p < 0.05\n'
                'ns  not significant')
    ax.text(0.02, 0.98, sig_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Add sample size information
    for i, (n, prop) in enumerate(zip(sample_sizes, proportions)):
        ax.text(i, prop/2, f'n={n:,}', ha='center', va='center', 
                fontsize=8, fontweight='bold', rotation=90)
    
    plt.tight_layout()
    return fig, ax

def load_multiple_analyses(csv_directory, analysis_function):
    """
    Load and analyze multiple CSV files from a directory.
    
    Parameters:
    -----------
    csv_directory : str or Path
        Directory containing CSV files
    analysis_function : callable
        Function that takes a CSV path and returns statistical results
    
    Returns:
    --------
    dict : Statistical results for each file
    """
    csv_dir = Path(csv_directory)
    results = {}
    
    for csv_file in csv_dir.glob('*.csv'):
        # Extract region/depth info from filename
        region_name = csv_file.stem.replace('_', ' ').title()
        
        try:
            # Run analysis (you'd need to import your analysis function)
            stats = analysis_function(csv_file)
            results[region_name] = stats
            print(f"Processed: {region_name}")
        except Exception as e:
            print(f"Error processing {csv_file}: {e}")
    
    return results

def example_usage():
    """Example of how to use this script with mock data."""
    
    # Mock data for demonstration (replace with your actual data)
    mock_data = {
        'Philippines Deep': {
            'proportion_interface_proximal': 0.755,
            'proportion_ci_lower': 0.734,
            'proportion_ci_upper': 0.775,
            'wilcoxon_p_value': 6.77e-275,
            'n_total': 1671
        },
        'Chile Shallow': {
            'proportion_interface_proximal': 1.00,
            'proportion_ci_lower': 0.986,
            'proportion_ci_upper': 1.00,
            'wilcoxon_p_value': 1e-100,
            'n_total': 267
        },
        'Vanuatu Deep': {
            'proportion_interface_proximal': 0.09,
            'proportion_ci_lower': 0.076,
            'proportion_ci_upper': 0.11,
            'wilcoxon_p_value': 1e-200,
            'n_total': 1516
        },
        'Japan Shallow': {
            'proportion_interface_proximal': 1.00,
            'proportion_ci_lower': 0.998,
            'proportion_ci_upper': 1.00,
            'wilcoxon_p_value': 0.0,
            'n_total': 2076
        },
        'Tonga Deep': {
            'proportion_interface_proximal': 0.57,
            'proportion_ci_lower': 0.54,
            'proportion_ci_upper': 0.60,
            'wilcoxon_p_value': 1e-50,
            'n_total': 1053
        }
    }
    
    fig, ax = create_statistical_summary_chart(mock_data, 'proportion')
    
    return fig, ax

def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(
        description="Create statistical summary charts for earthquake analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python statistical_summary.py --example
    python statistical_summary.py --data results.csv --output summary_chart.png
        """
    )
    
    parser.add_argument('--example', action='store_true',
                       help='Run example with mock data')
    parser.add_argument('--data', type=str,
                       help='CSV file with compiled statistical results')
    parser.add_argument('--output', type=str, default=None,
                       help='Output filename for chart')
    parser.add_argument('--dpi', type=int, default=300,
                       help='Resolution for saved figure')
    
    args = parser.parse_args()
    
    try:
        if args.example:
            print("Running example with mock data...")
            fig, ax = example_usage()
        else:
            raise NotImplementedError("Data loading not yet implemented. Use --example for now.")
        
        if args.output:
            print(f"Saving chart to: {args.output}")
            plt.savefig(args.output, dpi=args.dpi, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            print(f"Chart saved successfully")
        else:
            plt.show()
            
    except Exception as e:
        print(f"Error: {str(e)}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())