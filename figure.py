#!/usr/bin/env python3
"""
Plot disease simulation metrics vs sensitivity parameter (V) for different networks
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgb

# ============ CONFIGURATION ============
# Choose metric to plot: 'peak_prevalence', 'final_size', or 'min_edge_reduction'
METRIC = 'min_edge_reduction'  # Options: 'peak_prevalence', 'final_size', 'min_edge_reduction'

# Fixed uncertainty band width for each network and metric
# Format: {metric: {network: value}}
UNCERTAINTY_BANDS = {
    'final_size': {
        'Barabasi': 10/500,
        'Erdos-Renyi': 10/500,
        'Small World': 10/500
    },
    'peak_prevalence': {
        'Barabasi': 10/500,
        'Erdos-Renyi': 10/500,
        'Small World': 10/500
    },
    'min_edge_reduction': {
        'Barabasi': 0.01,
        'Erdos-Renyi': 0.01,
        'Small World': 0.01
    }
}
# =======================================

# Read the data
df = pd.read_csv('./Data/20251122115710_all_networks_bepidist_heatmap_data.csv')

# Filter for T=7
df = df[df['T'] == 7].copy()

# Compute the metric based on configuration
if METRIC == 'peak_prevalence':
    # Get the peak prevalence (max i_mean) for each combination of V and network
    def get_peak_mean(group):
        peak_idx = group['i_mean'].idxmax()
        return group.loc[peak_idx, 'i_mean']
    
    plot_data = df.groupby(['V', 'network']).apply(get_peak_mean).reset_index()
    plot_data.columns = ['sensitivity_parameter', 'network', 'mean_value']
    
    # Divide by 500 to get proportion
    plot_data['mean_value'] = plot_data['mean_value'] / 500
    
    # Get uncertainty bands for this metric
    uncertainty_bands = UNCERTAINTY_BANDS['peak_prevalence']
    
    y_label = 'Peak Prevalence'
    title = 'Peak Prevalence vs Sensitivity Parameter by Network'
    output_filename = 'peak_prevalence_by_sensitivity.png'
    plot_type = 'single'  # Single line per network
    
elif METRIC == 'final_size':
    # Filter for day 200 and get r_mean for each combination of V and network
    day_200 = df[df['day'] == 200].copy()
    plot_data = day_200.groupby(['V', 'network'])['r_mean'].mean().reset_index()
    plot_data.columns = ['sensitivity_parameter', 'network', 'mean_value']
    
    # Divide by 500 to get proportion
    plot_data['mean_value'] = plot_data['mean_value'] / 500
    
    # Get uncertainty bands for this metric
    uncertainty_bands = UNCERTAINTY_BANDS['final_size']
    
    y_label = 'Final Size'
    title = 'Final Size vs Sensitivity Parameter by Network'
    output_filename = 'final_size_by_sensitivity.png'
    plot_type = 'single'  # Single line per network
    
elif METRIC == 'min_edge_reduction':
    # Get minimum edgecount_mean and suscedgecount_mean for each V and network
    edge_data = df.groupby(['V', 'network']).agg({
        'edgecount_mean': 'min',
        'suscedgecount_mean': 'min'
    }).reset_index()
    
    plot_data = edge_data
    plot_data.columns = ['sensitivity_parameter', 'network', 'edgecount_min', 'suscedgecount_min']
    
    # Get uncertainty bands for this metric
    uncertainty_bands = UNCERTAINTY_BANDS['min_edge_reduction']
    
    y_label = 'Min Edge Reduction'
    title = 'Min Edge Reduction vs Sensitivity Parameter by Network'
    output_filename = 'min_edge_reduction_by_sensitivity.png'
    plot_type = 'dual'  # Two lines per network
    
else:
    raise ValueError(f"Invalid METRIC: {METRIC}. Choose 'peak_prevalence', 'final_size', or 'min_edge_reduction'")

# Get unique networks and sort them for consistent ordering
networks = sorted(plot_data['network'].unique())

# Create the plot
fig, ax = plt.subplots(figsize=(12, 7))

# Base color palette
base_colors = plt.cm.Set2(np.linspace(0, 1, len(networks)))

if plot_type == 'single':
    # Plot each network as a smooth line with shaded region
    for i, network in enumerate(networks):
        network_data = plot_data[plot_data['network'] == network].sort_values('sensitivity_parameter')
        
        x = network_data['sensitivity_parameter'].values
        y_mean = network_data['mean_value'].values
        
        # Get the fixed uncertainty band for this network
        uncertainty = uncertainty_bands.get(network, 0.01)  # Default to 0.01 if network not in config
        
        # Calculate upper and lower bounds (mean ± fixed uncertainty)
        y_upper = y_mean + uncertainty
        y_lower = y_mean - uncertainty
        
        # Plot the mean line (smooth, no markers)
        ax.plot(x, y_mean, linewidth=2.5, label=network, color=base_colors[i])
        
        # Fill the region between mean-uncertainty and mean+uncertainty
        ax.fill_between(x, y_lower, y_upper, alpha=0.3, color=base_colors[i])
        
elif plot_type == 'dual':
    # Plot two lines per network (edgecount and suscedgecount)
    for i, network in enumerate(networks):
        network_data = plot_data[plot_data['network'] == network].sort_values('sensitivity_parameter')
        
        x = network_data['sensitivity_parameter'].values
        y_edge = network_data['edgecount_min'].values
        y_susc = network_data['suscedgecount_min'].values
        
        # Get the fixed uncertainty band for this network
        uncertainty = uncertainty_bands.get(network, 0.01)
        
        # Create darker and lighter versions of the base color
        base_rgb = to_rgb(base_colors[i])
        # Darker version (multiply by 0.7)
        dark_color = tuple(c * 0.7 for c in base_rgb)
        # Lighter version (interpolate with white)
        light_color = tuple(c * 0.6 + 0.4 for c in base_rgb)
        
        # Plot edgecount_mean line (solid, darker)
        ax.plot(x, y_edge, linewidth=2.5, linestyle='-', 
                label=f'{network} (edge)', color=dark_color)
        y_upper_edge = y_edge + uncertainty
        y_lower_edge = y_edge - uncertainty
        ax.fill_between(x, y_lower_edge, y_upper_edge, alpha=0.3, color=dark_color)
        
        # Plot suscedgecount_mean line (dashed, lighter)
        ax.plot(x, y_susc, linewidth=2.5, linestyle='--', 
                label=f'{network} (susc)', color=light_color)
        y_upper_susc = y_susc + uncertainty
        y_lower_susc = y_susc - uncertainty
        ax.fill_between(x, y_lower_susc, y_upper_susc, alpha=0.3, color=light_color)

# Styling
ax.set_xlabel('Sensitivity Parameter (v)', fontsize=16, fontweight='bold')
ax.set_ylabel(y_label, fontsize=16, fontweight='bold')
ax.set_title(title, fontsize=16, fontweight='bold', pad=20)

# Make tick labels bold
ax.tick_params(axis='both', which='major', labelsize=16, width=2)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontweight('bold')

# Add grid for better readability
ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.8)

# Legend
ax.legend(fontsize=16, frameon=True, shadow=True, loc='best')

# Adjust layout to prevent label cutoff
plt.tight_layout()

# Save the figure
output_path = f'./Figures/{output_filename}'
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"Figure saved to: {output_path}")
print(f"Metric plotted: {METRIC}")
print(f"Filtered for T=7")
print(f"\nUncertainty bands used: {uncertainty_bands}")

# Also display some summary statistics
print("\nSummary Statistics:")
print("=" * 60)
if plot_type == 'single':
    for network in networks:
        network_data = plot_data[plot_data['network'] == network]
        uncertainty = uncertainty_bands.get(network, 0.01)
        print(f"\n{network}:")
        print(f"  Sensitivity range: {network_data['sensitivity_parameter'].min():.4f} - {network_data['sensitivity_parameter'].max():.4f}")
        print(f"  {y_label} (mean) range: {network_data['mean_value'].min():.4f} - {network_data['mean_value'].max():.4f}")
        print(f"  Fixed uncertainty band: ±{uncertainty:.4f}")
elif plot_type == 'dual':
    for network in networks:
        network_data = plot_data[plot_data['network'] == network]
        uncertainty = uncertainty_bands.get(network, 0.01)
        print(f"\n{network}:")
        print(f"  Sensitivity range: {network_data['sensitivity_parameter'].min():.4f} - {network_data['sensitivity_parameter'].max():.4f}")
        print(f"  Edgecount min range: {network_data['edgecount_min'].min():.4f} - {network_data['edgecount_min'].max():.4f}")
        print(f"  Suscedgecount min range: {network_data['suscedgecount_min'].min():.4f} - {network_data['suscedgecount_min'].max():.4f}")
        print(f"  Fixed uncertainty band: ±{uncertainty:.4f}")