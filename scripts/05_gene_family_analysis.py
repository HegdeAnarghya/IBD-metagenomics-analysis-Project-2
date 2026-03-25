import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import glob

print("="*60)
print("GENE FAMILY ANALYSIS")
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

print(f"\nAnalyzing gene families for {len(samples)} samples")

# Read gene family files
genefamily_data = {}

for sample_id, diagnosis in samples.items():
    # Find the gene family file
    filepath = data_dir / f"{sample_id}_genefamilies_cpm.tsv"
    
    if filepath.exists():
        print(f"\nReading {sample_id} ({diagnosis})...")
        
        # Read gene families
        df = pd.read_csv(filepath, sep='\t', comment='#')
        
        # First column is gene family, second is abundance
        df.columns = ['GeneFamily', 'Abundance_CPM']
        
        # Remove stratified (contains '|')
        df = df[~df['GeneFamily'].str.contains(r'\|', na=False)]
        
        # Remove unmapped/unintegrated
        df = df[~df['GeneFamily'].str.contains('UNMAPPED|UNINTEGRATED', na=False)]
        
        genefamily_data[sample_id] = {
            'data': df,
            'diagnosis': diagnosis
        }
        
        print(f"  {len(df)} gene families detected")
        print(f"  Total abundance: {df['Abundance_CPM'].sum():.1f} CPM")
    else:
        print(f"  ⚠️ File not found: {filepath}")

# Create combined gene family matrix
print("\n" + "="*60)
print("CREATING GENE FAMILY MATRIX")
print("="*60)

all_genefamilies = pd.DataFrame()

for sample_id, info in genefamily_data.items():
    temp = info['data'].copy()
    temp.columns = ['GeneFamily', sample_id]
    
    if all_genefamilies.empty:
        all_genefamilies = temp
    else:
        all_genefamilies = all_genefamilies.merge(temp, on='GeneFamily', how='outer')

all_genefamilies = all_genefamilies.fillna(0)
all_genefamilies = all_genefamilies.set_index('GeneFamily')

print(f"\nTotal unique gene families: {len(all_genefamilies)}")
print(f"Samples: {list(all_genefamilies.columns)}")

# Identify genes of interest
print("\n" + "="*60)
print("ANTIBIOTIC RESISTANCE GENES")
print("="*60)

# Search for AMR genes
amr_keywords = ['resist', 'beta.lactam', 'betalactam', 'mecA', 'vanA', 
                'tetracycline', 'aminoglycoside', 'macrolide', 'quinolone']

amr_genes = all_genefamilies[
    all_genefamilies.index.str.contains('|'.join(amr_keywords), case=False, na=False)
]

print(f"\nAntibiotic resistance genes found: {len(amr_genes)}")

if len(amr_genes) > 0:
    print("\nTop 10 AMR genes:")
    amr_genes['mean_abundance'] = amr_genes.mean(axis=1)
    top_amr = amr_genes.nlargest(10, 'mean_abundance')
    
    for gene in top_amr.index:
        mean_val = top_amr.loc[gene, 'mean_abundance']
        print(f"\n  {gene}")
        print(f"    Mean abundance: {mean_val:.2f} CPM")
        for sample in all_genefamilies.columns:
            val = all_genefamilies.loc[gene, sample]
            if val > 0:
                print(f"      {sample} ({samples[sample]}): {val:.2f} CPM")

# Butyrate synthesis genes
print("\n" + "="*60)
print("BUTYRATE SYNTHESIS GENES")
print("="*60)

butyrate_keywords = ['butyrate', 'butyryl', 'butyrat']

butyrate_genes = all_genefamilies[
    all_genefamilies.index.str.contains('|'.join(butyrate_keywords), case=False, na=False)
]

print(f"\nButyrate-related genes found: {len(butyrate_genes)}")

if len(butyrate_genes) > 0:
    print("\nButyrate synthesis genes:")
    butyrate_genes['mean_abundance'] = butyrate_genes.mean(axis=1)
    
    for gene in butyrate_genes.index:
        mean_val = butyrate_genes.loc[gene, 'mean_abundance']
        print(f"\n  {gene}")
        print(f"    Mean abundance: {mean_val:.2f} CPM")
        for sample in all_genefamilies.columns:
            val = all_genefamilies.loc[gene, sample]
            diagnosis = samples[sample]
            print(f"      {sample} ({diagnosis}): {val:.2f} CPM")

# Vitamin synthesis genes
print("\n" + "="*60)
print("VITAMIN SYNTHESIS GENES")
print("="*60)

vitamin_keywords = ['vitamin', 'cobalamin', 'folate', 'thiamin', 'riboflavin']

vitamin_genes = all_genefamilies[
    all_genefamilies.index.str.contains('|'.join(vitamin_keywords), case=False, na=False)
]

print(f"\nVitamin synthesis genes found: {len(vitamin_genes)}")

if len(vitamin_genes) > 0:
    vitamin_genes['mean_abundance'] = vitamin_genes.mean(axis=1)
    top_vit = vitamin_genes.nlargest(5, 'mean_abundance')
    
    print("\nTop 5 vitamin synthesis genes:")
    for gene in top_vit.index:
        mean_val = top_vit.loc[gene, 'mean_abundance']
        print(f"  {gene}: {mean_val:.2f} CPM (mean)")

# Top 20 most abundant gene families overall
print("\n" + "="*60)
print("TOP 20 MOST ABUNDANT GENE FAMILIES")
print("="*60)

all_genefamilies['mean_abundance'] = all_genefamilies.mean(axis=1)
top20_genes = all_genefamilies.nlargest(20, 'mean_abundance')

for idx, gene in enumerate(top20_genes.index, 1):
    mean_val = top20_genes.loc[gene, 'mean_abundance']
    print(f"\n{idx}. {gene}")
    print(f"   Mean: {mean_val:.2f} CPM")

# Save results
top20_genes.drop('mean_abundance', axis=1).to_csv(results_dir / 'top20_gene_families.csv')

if len(butyrate_genes) > 0:
    butyrate_genes.drop('mean_abundance', axis=1).to_csv(results_dir / 'butyrate_genes.csv')

print(f"\n Results saved to {results_dir}/")

# Visualization
print("\n" + "="*60)
print("CREATING VISUALIZATIONS")
print("="*60)

# Plot gene family richness
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Plot 1: Gene family richness
ax1 = axes[0]
gene_richness = (all_genefamilies.drop('mean_abundance', axis=1) > 0).sum(axis=0)

colors = {'nonIBD': 'green', 'CD': 'red', 'UC': 'orange'}
bar_colors = [colors[samples[s]] for s in gene_richness.index]

ax1.bar(range(len(gene_richness)), gene_richness.values,
        color=bar_colors, alpha=0.7, edgecolor='black')
ax1.set_xticks(range(len(gene_richness)))
ax1.set_xticklabels([f"{s}\n({samples[s]})" for s in gene_richness.index],
                     rotation=45, ha='right')
ax1.set_ylabel('Number of Gene Families Detected', fontsize=12)
ax1.set_title('Gene Family Richness', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3, axis='y')

from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=colors[d], label=d, alpha=0.7) 
                   for d in ['nonIBD', 'CD', 'UC']]
ax1.legend(handles=legend_elements, loc='upper right')

# Plot 2: Butyrate gene abundance comparison
ax2 = axes[1]

if len(butyrate_genes) > 0:
    butyrate_total = butyrate_genes.drop('mean_abundance', axis=1).sum(axis=0)
    
    ax2.bar(range(len(butyrate_total)), butyrate_total.values,
            color=bar_colors, alpha=0.7, edgecolor='black')
    ax2.set_xticks(range(len(butyrate_total)))
    ax2.set_xticklabels([f"{s}\n({samples[s]})" for s in butyrate_total.index],
                         rotation=45, ha='right')
    ax2.set_ylabel('Total Butyrate Gene Abundance (CPM)', fontsize=12)
    ax2.set_title('Butyrate Synthesis Gene Capacity', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Add value labels
    for i, val in enumerate(butyrate_total.values):
        ax2.text(i, val, f'{val:.0f}', ha='center', va='bottom', fontweight='bold')
else:
    ax2.text(0.5, 0.5, 'No butyrate genes detected', 
             ha='center', va='center', transform=ax2.transAxes, fontsize=14)
    ax2.axis('off')

plt.tight_layout()
plt.savefig(figures_dir / 'gene_family_analysis.png', dpi=300, bbox_inches='tight')
print(f" Figure saved: {figures_dir / 'gene_family_analysis.png'}")

print("\n" + "="*60)
print("GENE FAMILY ANALYSIS COMPLETE!")
print("="*60)
