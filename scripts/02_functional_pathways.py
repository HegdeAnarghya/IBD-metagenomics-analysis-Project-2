import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import glob

print("="*60)
print("FUNCTIONAL PATHWAY ANALYSIS")
print("="*60)

# File paths
data_dir = Path('../data')
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

print(f"\nAnalyzing functional profiles for {len(samples)} samples")

# Read pathway abundance files
pathway_data = {}

for sample_id, diagnosis in samples.items():
    # Find the pathway file
    filepath = data_dir / f"{sample_id}_pathabundance_cpm.tsv"
    
    if filepath.exists():
        print(f"\nReading {sample_id} ({diagnosis})...")
        
        # Read pathway abundances
        df = pd.read_csv(filepath, sep='\t', comment='#')
        
        # First column is pathway name, second is abundance
        df.columns = ['Pathway', 'Abundance_CPM']
        
        # Remove stratified pathways (contain '|')
        df = df[~df['Pathway'].str.contains(r'\|', na=False)]
        
        # Remove unmapped/unintegrated
        df = df[~df['Pathway'].str.contains('UNMAPPED|UNINTEGRATED', na=False)]
        
        pathway_data[sample_id] = {
            'data': df,
            'diagnosis': diagnosis
        }
        
        print(f"  {len(df)} pathways detected")
        print(f"  Total abundance: {df['Abundance_CPM'].sum():.1f} CPM")
    else:
        print(f"  File not found: {filepath}")

# Create combined pathway matrix
print("\n" + "="*60)
print("CREATING PATHWAY ABUNDANCE MATRIX")
print("="*60)

all_pathways = pd.DataFrame()

for sample_id, info in pathway_data.items():
    temp = info['data'].copy()
    temp.columns = ['Pathway', sample_id]
    
    if all_pathways.empty:
        all_pathways = temp
    else:
        all_pathways = all_pathways.merge(temp, on='Pathway', how='outer')

all_pathways = all_pathways.fillna(0)
all_pathways = all_pathways.set_index('Pathway')

print(f"\nTotal unique pathways: {len(all_pathways)}")
print(f"Samples: {list(all_pathways.columns)}")

# Identify key pathways
print("\n" + "="*60)
print("KEY PATHWAYS OF INTEREST")
print("="*60)

# Butyrate synthesis pathways
butyrate_keywords = ['butyrate', 'butanoate', 'butyrat']
butyrate_pathways = all_pathways[
    all_pathways.index.str.contains('|'.join(butyrate_keywords), case=False, na=False)
]

print(f"\n1. BUTYRATE SYNTHESIS PATHWAYS ({len(butyrate_pathways)} found):")
if len(butyrate_pathways) > 0:
    for pathway in butyrate_pathways.index:
        print(f"\n  {pathway}")
        for sample in all_pathways.columns:
            diagnosis = samples[sample]
            abundance = all_pathways.loc[pathway, sample]
            print(f"    {sample} ({diagnosis}): {abundance:.2f} CPM")
else:
    print("  No butyrate pathways detected")

# Inflammation-related
inflam_keywords = ['lipopolysaccharide', 'LPS', 'endotoxin']
inflam_pathways = all_pathways[
    all_pathways.index.str.contains('|'.join(inflam_keywords), case=False, na=False)
]

print(f"\n2. INFLAMMATION-RELATED PATHWAYS ({len(inflam_pathways)} found):")
if len(inflam_pathways) > 0:
    for pathway in inflam_pathways.index[:5]:  # Show top 5
        print(f"  {pathway}")

# Short-chain fatty acid synthesis
scfa_keywords = ['propanoate', 'propionate', 'acetate']
scfa_pathways = all_pathways[
    all_pathways.index.str.contains('|'.join(scfa_keywords), case=False, na=False)
]

print(f"\n3. SHORT-CHAIN FATTY ACID PATHWAYS ({len(scfa_pathways)} found):")
if len(scfa_pathways) > 0:
    for pathway in scfa_pathways.index[:5]:
        print(f"  {pathway}")

# Get top 20 most abundant pathways
print("\n" + "="*60)
print("TOP 20 MOST ABUNDANT PATHWAYS")
print("="*60)

# Calculate mean abundance across samples
all_pathways['mean_abundance'] = all_pathways.mean(axis=1)
top20 = all_pathways.nlargest(20, 'mean_abundance')

for idx, pathway in enumerate(top20.index, 1):
    mean_abund = top20.loc[pathway, 'mean_abundance']
    print(f"\n{idx}. {pathway}")
    print(f"   Mean: {mean_abund:.2f} CPM")

# Save results
top20_for_save = top20.drop('mean_abundance', axis=1)
top20_for_save.to_csv(results_dir / 'top20_pathways.csv')

print(f"\n Results saved to {results_dir}/")

# Visualization
print("\n" + "="*60)
print("CREATING VISUALIZATIONS")
print("="*60)

# Plot top 15 pathways as heatmap
fig, ax = plt.subplots(figsize=(12, 10))

top15 = all_pathways.nlargest(15, 'mean_abundance').drop('mean_abundance', axis=1)

# Create heatmap
sns.heatmap(top15, cmap='YlOrRd', annot=True, fmt='.0f', 
            cbar_kws={'label': 'Abundance (CPM)'},
            linewidths=0.5, ax=ax)

ax.set_title('Top 15 Metabolic Pathways (Abundance in CPM)', 
             fontsize=14, fontweight='bold', pad=20)
ax.set_xlabel('Sample', fontsize=12)
ax.set_ylabel('Pathway', fontsize=12)

# Add diagnosis labels
diagnosis_labels = [f"{s}\n({samples[s]})" for s in top15.columns]
ax.set_xticklabels(diagnosis_labels, rotation=45, ha='right')

plt.tight_layout()
plt.savefig(figures_dir / 'pathway_heatmap.png', dpi=300, bbox_inches='tight')
print(f" Heatmap saved: {figures_dir / 'pathway_heatmap.png'}")

# Compare pathway diversity
pathway_richness = (all_pathways.drop('mean_abundance', axis=1) > 0).sum(axis=0)

fig2, ax2 = plt.subplots(figsize=(10, 6))

colors = {'nonIBD': 'green', 'CD': 'red', 'UC': 'orange'}
bar_colors = [colors[samples[s]] for s in pathway_richness.index]

ax2.bar(range(len(pathway_richness)), pathway_richness.values,
        color=bar_colors, alpha=0.7, edgecolor='black')
ax2.set_xticks(range(len(pathway_richness)))
ax2.set_xticklabels([f"{s}\n({samples[s]})" for s in pathway_richness.index],
                     rotation=45, ha='right')
ax2.set_ylabel('Number of Pathways Detected', fontsize=12)
ax2.set_title('Functional Pathway Richness', fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3, axis='y')

from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=colors[d], label=d, alpha=0.7) 
                   for d in ['nonIBD', 'CD', 'UC']]
ax2.legend(handles=legend_elements, loc='upper right')

plt.tight_layout()
plt.savefig(figures_dir / 'pathway_richness.png', dpi=300, bbox_inches='tight')
print(f" Pathway richness saved: {figures_dir / 'pathway_richness.png'}")

print("\n" + "="*60)
print("FUNCTIONAL ANALYSIS COMPLETE!")
print("="*60)
