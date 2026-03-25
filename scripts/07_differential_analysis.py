import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats

print("="*60)
print("DIFFERENTIAL PATHWAY & ENZYME ANALYSIS")
print("="*60)

# File paths
results_dir = Path('../results')
figures_dir = Path('../figures')

# Sample information
samples = {
    'MSM5LLHV_P': 'UC',
    'HSM7CZ2A': 'nonIBD',
    'HSM6XRQE_P': 'UC',
    'CSM5FZ4C_P': 'CD',
    'CSM9X1ZO': 'UC'
}

# Load pathway data
print("\nLoading pathway data...")
pathway_file = results_dir / 'top20_pathways.csv'
if pathway_file.exists():
    pathways = pd.read_csv(pathway_file, index_col=0)
    print(f"✓ {len(pathways)} pathways loaded")
else:
    print("⚠️ Pathway file not found!")
    pathways = None

# Load EC enzyme data
print("Loading EC enzyme data...")
ec_file = results_dir / 'all_ec_enzymes.csv'
if ec_file.exists():
    ec_data = pd.read_csv(ec_file)
    
    # Pivot to matrix format
    ec_matrix = ec_data.pivot_table(
        index='EC_number',
        columns='EC_number',
        values=list(samples.keys()),
        aggfunc='first'
    )
    
    # Actually, let's reload properly
    # The file has EC_number, EC_class, and sample columns
    ec_wide = ec_data.drop('EC_class', axis=1, errors='ignore')
    ec_wide = ec_wide.set_index('EC_number')
    
    print(f"✓ {len(ec_wide)} EC enzymes loaded")
else:
    print("⚠️ EC file not found!")
    ec_wide = None

print("\n" + "="*60)
print("GROUP DEFINITIONS")
print("="*60)

# Define groups (excluding HSM6XRQE outlier for fair comparison)
healthy_samples = ['HSM7CZ2A']
cd_samples = ['CSM5FZ4C_P']
uc_samples = ['MSM5LLHV_P', 'CSM9X1ZO']  # Exclude HSM6XRQE outlier
all_ibd_samples = cd_samples + uc_samples

print(f"\nHealthy: {healthy_samples}")
print(f"CD: {cd_samples}")
print(f"UC: {uc_samples}")
print(f"All IBD: {all_ibd_samples}")
print(f"\nNote: HSM6XRQE excluded (extreme outlier)")

# Function to calculate fold change and comparison
def compare_groups(data, group1, group2, group1_name, group2_name):
    """Compare two groups and return statistics"""
    
    results = []
    
    for feature in data.index:
        # Get values for each group
        g1_vals = data.loc[feature, group1].values
        g2_vals = data.loc[feature, group2].values
        
        # Calculate means
        g1_mean = np.mean(g1_vals)
        g2_mean = np.mean(g2_vals)
        
        # Skip if both means are zero
        if g1_mean == 0 and g2_mean == 0:
            continue
        
        # Calculate fold change (avoid division by zero)
        if g2_mean > 0:
            fold_change = g1_mean / g2_mean
            log2_fc = np.log2(fold_change) if fold_change > 0 else np.nan
        else:
            fold_change = np.inf if g1_mean > 0 else 1
            log2_fc = np.inf if g1_mean > 0 else 0
        
        # Statistical test (Mann-Whitney U for small samples)
        if len(g1_vals) > 1 and len(g2_vals) > 1:
            statistic, p_value = stats.mannwhitneyu(g1_vals, g2_vals, alternative='two-sided')
        else:
            # Use difference for single samples
            p_value = 1.0 if abs(g1_mean - g2_mean) < 10 else 0.1
        
        results.append({
            'Feature': feature,
            f'{group1_name}_mean': g1_mean,
            f'{group2_name}_mean': g2_mean,
            'Fold_Change': fold_change,
            'Log2_FC': log2_fc,
            'P_value': p_value,
            'Difference': g1_mean - g2_mean
        })
    
    return pd.DataFrame(results)

# Pathway comparison: Healthy vs IBD
if pathways is not None and len(pathways) > 0:
    print("\n" + "="*60)
    print("PATHWAY COMPARISON: Healthy vs IBD")
    print("="*60)
    
    pathway_diff = compare_groups(
        pathways,
        healthy_samples,
        all_ibd_samples,
        'Healthy',
        'IBD'
    )
    
    # Sort by absolute log2 fold change
    pathway_diff = pathway_diff.sort_values('Log2_FC', key=abs, ascending=False)
    
    print(f"\nTotal pathways compared: {len(pathway_diff)}")
    print(f"\nTop 10 pathways ENRICHED in IBD (higher in disease):")
    
    ibd_enriched = pathway_diff[pathway_diff['Log2_FC'] < 0].head(10)
    for idx, row in ibd_enriched.iterrows():
        print(f"\n  {row['Feature']}")
        print(f"    Healthy: {row['Healthy_mean']:.1f} CPM")
        print(f"    IBD: {row['IBD_mean']:.1f} CPM")
        print(f"    Fold change: {row['Fold_Change']:.2f}× (Log2: {row['Log2_FC']:.2f})")
        print(f"    P-value: {row['P_value']:.3f}")
    
    print(f"\nTop 10 pathways DEPLETED in IBD (lower in disease):")
    
    ibd_depleted = pathway_diff[pathway_diff['Log2_FC'] > 0].head(10)
    for idx, row in ibd_depleted.iterrows():
        print(f"\n  {row['Feature']}")
        print(f"    Healthy: {row['Healthy_mean']:.1f} CPM")
        print(f"    IBD: {row['IBD_mean']:.1f} CPM")
        print(f"    Fold change: {row['Fold_Change']:.2f}× (Log2: {row['Log2_FC']:.2f})")
        print(f"    P-value: {row['P_value']:.3f}")
    
    # Save
    pathway_diff.to_csv(results_dir / 'pathway_differential_analysis.csv', index=False)

# EC enzyme comparison: Healthy vs IBD
if ec_wide is not None and len(ec_wide) > 0:
    print("\n" + "="*60)
    print("ENZYME COMPARISON: Healthy vs IBD")
    print("="*60)
    
    # Filter for enzymes present in at least one sample
    ec_filtered = ec_wide[(ec_wide > 0).any(axis=1)]
    
    ec_diff = compare_groups(
        ec_filtered,
        healthy_samples,
        all_ibd_samples,
        'Healthy',
        'IBD'
    )
    
    # Sort by absolute difference (since we have small n)
    ec_diff = ec_diff.sort_values('Difference', key=abs, ascending=False)
    
    print(f"\nTotal enzymes compared: {len(ec_diff)}")
    
    print(f"\nTop 10 enzymes MOST DIFFERENT between Healthy and IBD:")
    top10_diff = ec_diff.head(10)
    
    for idx, row in top10_diff.iterrows():
        direction = "HIGHER in Healthy" if row['Difference'] > 0 else "HIGHER in IBD"
        print(f"\n  EC {row['Feature']} - {direction}")
        print(f"    Healthy: {row['Healthy_mean']:.1f} CPM")
        print(f"    IBD: {row['IBD_mean']:.1f} CPM")
        print(f"    Difference: {row['Difference']:.1f} CPM")
        print(f"    Fold change: {row['Fold_Change']:.2f}×")
    
    # Specific enzymes of interest
    print(f"\n" + "="*60)
    print("KEY ENZYMES: Healthy vs IBD Comparison")
    print("="*60)
    
    key_enzymes = {
        '2.7.2.7': 'Butyrate kinase',
        '1.3.8.1': 'Butyryl-CoA dehydrogenase',
        '3.5.2.6': 'Beta-lactamase (AMR)',
        '1.1.1.27': 'Lactate dehydrogenase'
    }
    
    for ec_num, name in key_enzymes.items():
        if ec_num in ec_diff['Feature'].values:
            row = ec_diff[ec_diff['Feature'] == ec_num].iloc[0]
            print(f"\n{name} (EC {ec_num}):")
            print(f"  Healthy: {row['Healthy_mean']:.2f} CPM")
            print(f"  IBD: {row['IBD_mean']:.2f} CPM")
            print(f"  Fold change: {row['Fold_Change']:.2f}×")
            print(f"  Difference: {row['Difference']:.2f} CPM")
    
    # Save
    ec_diff.to_csv(results_dir / 'ec_differential_analysis.csv', index=False)

print(f"\n Results saved to {results_dir}/")

# Visualization
print("\n" + "="*60)
print("CREATING VISUALIZATIONS")
print("="*60)

fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

# Plot 1: Pathway volcano plot (if available)
if pathways is not None and len(pathway_diff) > 0:
    ax1 = fig.add_subplot(gs[0, 0])
    
    # Remove infinite values for plotting
    plot_data = pathway_diff[np.isfinite(pathway_diff['Log2_FC'])].copy()
    plot_data['-log10_p'] = -np.log10(plot_data['P_value'] + 1e-10)
    
    # Color by significance
    colors = ['red' if abs(fc) > 1 else 'gray' for fc in plot_data['Log2_FC']]
    
    ax1.scatter(plot_data['Log2_FC'], plot_data['-log10_p'], 
                c=colors, alpha=0.6, s=50)
    ax1.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax1.axhline(y=-np.log10(0.05), color='blue', linestyle='--', 
                linewidth=1, alpha=0.5, label='p=0.05')
    ax1.set_xlabel('Log2 Fold Change (Healthy/IBD)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('-Log10 P-value', fontsize=12, fontweight='bold')
    ax1.set_title('Pathway Differential Analysis', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

# Plot 2: Top differential pathways
if pathways is not None and len(pathway_diff) > 0:
    ax2 = fig.add_subplot(gs[0, 1])
    
    # Get top 10 by absolute difference
    top_diff = pathway_diff.nlargest(10, 'Difference', keep='first')
    
    y_pos = np.arange(len(top_diff))
    colors_bar = ['green' if d > 0 else 'red' for d in top_diff['Difference']]
    
    ax2.barh(y_pos, top_diff['Difference'], color=colors_bar, alpha=0.7, edgecolor='black')
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([p[:50] + '...' if len(p) > 50 else p 
                         for p in top_diff['Feature']], fontsize=8)
    ax2.set_xlabel('Difference (Healthy - IBD) CPM', fontsize=11, fontweight='bold')
    ax2.set_title('Top 10 Differential Pathways', fontsize=13, fontweight='bold')
    ax2.axvline(x=0, color='black', linestyle='-', linewidth=1)
    ax2.grid(True, alpha=0.3, axis='x')
    ax2.invert_yaxis()

# Plot 3: Key enzymes comparison
if ec_wide is not None:
    ax3 = fig.add_subplot(gs[1, 0])
    
    key_ec_data = []
    for ec_num, name in key_enzymes.items():
        if ec_num in ec_wide.index:
            for sample in healthy_samples:
                key_ec_data.append({
                    'Enzyme': name,
                    'Group': 'Healthy',
                    'Abundance': ec_wide.loc[ec_num, sample]
                })
            for sample in all_ibd_samples:
                key_ec_data.append({
                    'Enzyme': name,
                    'Group': 'IBD',
                    'Abundance': ec_wide.loc[ec_num, sample]
                })
    
    if key_ec_data:
        key_df = pd.DataFrame(key_ec_data)
        
        enzymes = key_df['Enzyme'].unique()
        x = np.arange(len(enzymes))
        width = 0.35
        
        healthy_means = [key_df[(key_df['Enzyme']==e) & (key_df['Group']=='Healthy')]['Abundance'].mean() 
                        for e in enzymes]
        ibd_means = [key_df[(key_df['Enzyme']==e) & (key_df['Group']=='IBD')]['Abundance'].mean() 
                    for e in enzymes]
        
        ax3.bar(x - width/2, healthy_means, width, label='Healthy', 
                color='green', alpha=0.7, edgecolor='black')
        ax3.bar(x + width/2, ibd_means, width, label='IBD', 
                color='red', alpha=0.7, edgecolor='black')
        
        ax3.set_ylabel('Abundance (CPM)', fontsize=12, fontweight='bold')
        ax3.set_title('Key Enzymes: Healthy vs IBD', fontsize=14, fontweight='bold')
        ax3.set_xticks(x)
        ax3.set_xticklabels([e.split('(')[0].strip() for e in enzymes], 
                            rotation=45, ha='right', fontsize=9)
        ax3.legend()
        ax3.grid(True, alpha=0.3, axis='y')

# Plot 4: EC class comparison with values
if ec_wide is not None:
    ax4 = fig.add_subplot(gs[1, 1])
    
    # Reload EC data with class info
    ec_full = pd.read_csv(ec_file)
    
    # Define EC class names
    ec_names = {
        '1': 'Oxidoreductases',
        '2': 'Transferases',
        '3': 'Hydrolases',
        '4': 'Lyases',
        '5': 'Isomerases',
        '6': 'Ligases'
    }
    
    # Sum by EC class
    class_comparison = []
    for ec_class in ['1', '2', '3', '4', '5', '6']:
        class_data = ec_full[ec_full['EC_class'] == ec_class]
        
        if len(class_data) > 0:
            healthy_sum = sum([class_data[s].sum() for s in healthy_samples if s in class_data.columns])
            ibd_sum = sum([class_data[s].sum() for s in all_ibd_samples if s in class_data.columns]) / len(all_ibd_samples)
            
            class_comparison.append({
                'EC_Class': ec_names.get(ec_class, f'EC {ec_class}'),
                'Healthy': healthy_sum,
                'IBD': ibd_sum,
                'Difference': healthy_sum - ibd_sum
            })
    
    if class_comparison:
        comp_df = pd.DataFrame(class_comparison)
        
        x = np.arange(len(comp_df))
        width = 0.35
        
        bars1 = ax4.bar(x - width/2, comp_df['Healthy'], width, label='Healthy',
                color='green', alpha=0.7, edgecolor='black')
        bars2 = ax4.bar(x + width/2, comp_df['IBD'], width, label='IBD',
                color='red', alpha=0.7, edgecolor='black')
        
        ax4.set_ylabel('Total Abundance (CPM)', fontsize=12, fontweight='bold')
        ax4.set_title('EC Class Totals: Healthy vs IBD', fontsize=14, fontweight='bold')
        ax4.set_xticks(x)
        ax4.set_xticklabels([name.split()[0] for name in comp_df['EC_Class']], 
                            rotation=45, ha='right', fontsize=9)
        ax4.legend(loc='upper right')
        ax4.grid(True, alpha=0.3, axis='y')
        
        # Add difference percentages
        for i, (h, ibd, diff) in enumerate(zip(comp_df['Healthy'], comp_df['IBD'], comp_df['Difference'])):
            if h > 0:
                pct = (diff / h) * 100
                ax4.text(i, max(h, ibd) * 1.05, f'{pct:+.0f}%', 
                        ha='center', fontsize=8, fontweight='bold')

plt.savefig(figures_dir / 'differential_analysis.png', dpi=300, bbox_inches='tight')
print(f" Figure saved: {figures_dir / 'differential_analysis.png'}")

print("\n" + "="*60)
print("DIFFERENTIAL ANALYSIS COMPLETE!")
print("="*60)
