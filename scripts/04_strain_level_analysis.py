import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

print("="*60)
print("STRAIN-LEVEL ANALYSIS")
print("="*60)

# File paths
data_dir = Path('../data')
results_dir = Path('../results')
figures_dir = Path('../figures')

# Sample information
samples = {
    'sample1_taxonomic.tsv': {'id': 'MSM5LLHV', 'diagnosis': 'UC'},
    'sample2_taxonomic.tsv': {'id': 'HSM7CZ2A', 'diagnosis': 'nonIBD'},
    'sample3_taxonomic.tsv': {'id': 'HSM6XRQE', 'diagnosis': 'UC'},
    'sample4_taxonomic.tsv': {'id': 'CSM5FZ4C', 'diagnosis': 'CD'},
    'sample5_taxonomic.tsv': {'id': 'CSM9X1ZO', 'diagnosis': 'UC'},
}

print(f"\nAnalyzing strain-level data for {len(samples)} samples")

# Read all taxonomic profiles
all_profiles = {}

for filename, info in samples.items():
    filepath = data_dir / filename
    
    # Read MetaPhlAn output
    with open(filepath) as f:
        lines = [l.strip() for l in f if not l.startswith('#')]
    
    # Parse data
    data = []
    for line in lines:
        parts = line.split('\t')
        if len(parts) >= 3:
            clade = parts[0]
            abundance = float(parts[2]) if len(parts) > 2 else 0.0
            data.append({'clade_name': clade, 'relative_abundance': abundance})
    
    df = pd.DataFrame(data)
    all_profiles[info['id']] = {
        'data': df,
        'diagnosis': info['diagnosis']
    }

# Extract strain-level data (contains 't__')
print("\n" + "="*60)
print("STRAIN-LEVEL TAXONOMIC RESOLUTION")
print("="*60)

strain_data = {}

for sample_id, profile in all_profiles.items():
    df = profile['data']
    
    # Filter for strain level (contains 't__')
    strains = df[df['clade_name'].str.contains(r'\|t__', na=False)].copy()
    
    # Extract strain name and species
    strains['strain'] = strains['clade_name'].str.extract(r't__([^|#]+)')
    strains['species'] = strains['clade_name'].str.extract(r's__([^|]+)')
    
    strains['strain'] = strains['strain'].str.strip()
    strains['species'] = strains['species'].str.strip()
    
    strain_data[sample_id] = strains[['species', 'strain', 'relative_abundance']]
    
    print(f"\n{sample_id} ({all_profiles[sample_id]['diagnosis']}): {len(strains)} strains detected")
    
    if len(strains) > 0:
        print(f"Top 5 abundant strains:")
        top5 = strains.nlargest(min(5, len(strains)), 'relative_abundance')
        for _, row in top5.iterrows():
            print(f"  {row['species']} | {row['strain']}: {row['relative_abundance']:.2f}%")
    else:
        print("  No strain-level resolution available")

# Identify strains of key species
print("\n" + "="*60)
print("KEY SPECIES: STRAIN VARIATION")
print("="*60)

# Combine all strain data
all_strains = pd.DataFrame()
for sample_id, data in strain_data.items():
    temp = data.copy()
    temp['sample'] = sample_id
    temp['diagnosis'] = all_profiles[sample_id]['diagnosis']
    all_strains = pd.concat([all_strains, temp], ignore_index=True)

# Focus on species present in multiple samples
species_counts = all_strains['species'].value_counts()
common_species = species_counts[species_counts >= 2].index.tolist()

print(f"\nSpecies with strain-level data in ≥2 samples: {len(common_species)}")

# Analyze strain diversity within species
print("\n" + "="*60)
print("INTRA-SPECIES STRAIN DIVERSITY")
print("="*60)

for species in common_species[:10]:  # Top 10 most common
    species_strains = all_strains[all_strains['species'] == species]
    
    print(f"\n{species}:")
    print(f"  Found in {len(species_strains)} sample-strain combinations")
    print(f"  Unique strains: {species_strains['strain'].nunique()}")
    
    # Show which samples have which strains
    for sample in species_strains['sample'].unique():
        sample_strains = species_strains[species_strains['sample'] == sample]
        diagnosis = sample_strains['diagnosis'].iloc[0]
        for _, row in sample_strains.iterrows():
            print(f"    {sample} ({diagnosis}): {row['strain']} - {row['relative_abundance']:.2f}%")

# Compare Bacteroides strains (most abundant genus)
print("\n" + "="*60)
print("BACTEROIDES STRAIN COMPARISON")
print("="*60)

bacteroides_strains = all_strains[all_strains['species'].str.contains('Bacteroides', na=False)]

if len(bacteroides_strains) > 0:
    print(f"\nTotal Bacteroides strain detections: {len(bacteroides_strains)}")
    print(f"Unique Bacteroides species: {bacteroides_strains['species'].nunique()}")
    print(f"Unique Bacteroides strains: {bacteroides_strains['strain'].nunique()}")
    
    # Create pivot table
    bact_pivot = bacteroides_strains.pivot_table(
        index='strain',
        columns='sample',
        values='relative_abundance',
        fill_value=0
    )
    
    print("\nBacteroides Strain Abundance Matrix:")
    print(bact_pivot.to_string())
    
    # Save
    bact_pivot.to_csv(results_dir / 'bacteroides_strains.csv')

# Save all strain data
all_strains.to_csv(results_dir / 'all_strains_detected.csv', index=False)
print(f"\n Strain data saved to {results_dir}/")

# Visualization: Strain diversity per sample
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Plot 1: Number of strains per sample
ax1 = axes[0]
strain_counts = all_strains.groupby('sample')['strain'].nunique()

colors = {'nonIBD': 'green', 'CD': 'red', 'UC': 'orange'}
sample_colors = [colors[all_profiles[s]['diagnosis']] for s in strain_counts.index]

ax1.bar(range(len(strain_counts)), strain_counts.values,
        color=sample_colors, alpha=0.7, edgecolor='black')
ax1.set_xticks(range(len(strain_counts)))
ax1.set_xticklabels([f"{s}\n({all_profiles[s]['diagnosis']})" 
                     for s in strain_counts.index], rotation=45, ha='right')
ax1.set_ylabel('Number of Strains Detected', fontsize=12)
ax1.set_title('Strain-Level Resolution by Sample', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3, axis='y')

from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=colors[d], label=d, alpha=0.7) 
                   for d in ['nonIBD', 'CD', 'UC']]
ax1.legend(handles=legend_elements, loc='upper right')

# Plot 2: Top species with strain diversity
ax2 = axes[1]

# Count strains per species
species_strain_counts = all_strains.groupby('species')['strain'].nunique().nlargest(10)

ax2.barh(range(len(species_strain_counts)), species_strain_counts.values,
         color='steelblue', alpha=0.7, edgecolor='black')
ax2.set_yticks(range(len(species_strain_counts)))
ax2.set_yticklabels([s.replace('_', ' ') for s in species_strain_counts.index], 
                     fontsize=9)
ax2.invert_yaxis()
ax2.set_xlabel('Number of Unique Strains', fontsize=12)
ax2.set_title('Top 10 Species by Strain Diversity', fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3, axis='x')

plt.tight_layout()
plt.savefig(figures_dir / 'strain_level_analysis.png', dpi=300, bbox_inches='tight')
print(f" Figure saved: {figures_dir / 'strain_level_analysis.png'}")

print("\n" + "="*60)
print("STRAIN ANALYSIS COMPLETE!")
print("="*60)
