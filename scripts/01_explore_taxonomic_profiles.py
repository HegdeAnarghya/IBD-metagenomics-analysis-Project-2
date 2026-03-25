import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

print("="*60)
print("METAGENOMIC TAXONOMIC PROFILE ANALYSIS")
print("="*60)

# File paths
data_dir = Path('../data')
results_dir = Path('../results')
figures_dir = Path('../figures')

# Create directories if needed
results_dir.mkdir(exist_ok=True)
figures_dir.mkdir(exist_ok=True)

# Sample information
samples = {
    'sample1_taxonomic.tsv': {'id': 'MSM5LLHV', 'diagnosis': 'UC'},
    'sample2_taxonomic.tsv': {'id': 'HSM7CZ2A', 'diagnosis': 'nonIBD'},
    'sample3_taxonomic.tsv': {'id': 'HSM6XRQE', 'diagnosis': 'UC'},
    'sample4_taxonomic.tsv': {'id': 'CSM5FZ4C', 'diagnosis': 'CD'},
    'sample5_taxonomic.tsv': {'id': 'CSM9X1ZO', 'diagnosis': 'UC'},
}

print(f"\nAnalyzing {len(samples)} samples:")
for fname, info in samples.items():
    print(f"  {info['id']}: {info['diagnosis']}")

# Read all taxonomic profiles
all_profiles = {}

for filename, info in samples.items():
    filepath = data_dir / filename
    
    # Read MetaPhlAn output - flexible column reading
    with open(filepath) as f:
        lines = [l.strip() for l in f if not l.startswith('#')]
    
    # Parse data
    data = []
    for line in lines:
        parts = line.split('\t')
        if len(parts) >= 3:
            clade = parts[0]
            ncbi_id = parts[1] if len(parts) > 1 else ''
            abundance = float(parts[2]) if len(parts) > 2 else 0.0
            data.append({'clade_name': clade, 'NCBI_tax_id': ncbi_id, 
                        'relative_abundance': abundance})
    
    df = pd.DataFrame(data)
    
    # Store
    all_profiles[info['id']] = {
        'data': df,
        'diagnosis': info['diagnosis']
    }
    
    print(f"\n{info['id']} ({info['diagnosis']}): {len(df)} clades detected")

# Extract species-level abundances
print("\n" + "="*60)
print("SPECIES-LEVEL ANALYSIS")
print("="*60)

species_data = {}

for sample_id, profile in all_profiles.items():
    df = profile['data']
    
    # Filter for species level (contains 's__' but not 't__')
    species = df[df['clade_name'].str.contains(r'\|s__', na=False) & 
                 ~df['clade_name'].str.contains(r'\|t__', na=False)].copy()
    
    # Extract species name
    species['species'] = species['clade_name'].str.extract(r's__([^|#]+)')
    species['species'] = species['species'].str.strip()
    
    species_data[sample_id] = species[['species', 'relative_abundance']]
    
    print(f"\n{sample_id}: {len(species)} species detected")
    print(f"Top 5 abundant species:")
    if len(species) > 0:
        top5 = species.nlargest(min(5, len(species)), 'relative_abundance')
        for _, row in top5.iterrows():
            print(f"  {row['species']}: {row['relative_abundance']:.2f}%")

# Create combined species abundance matrix
print("\n" + "="*60)
print("CREATING ABUNDANCE MATRIX")
print("="*60)

# Combine all species
all_species = pd.DataFrame()
for sample_id, data in species_data.items():
    temp = data.copy()
    temp.columns = ['species', sample_id]
    if all_species.empty:
        all_species = temp
    else:
        all_species = all_species.merge(temp, on='species', how='outer')

all_species = all_species.fillna(0)
all_species = all_species.set_index('species')

print(f"\nTotal unique species: {len(all_species)}")
print(f"Samples: {list(all_species.columns)}")

# Calculate alpha diversity (species richness)
richness = (all_species > 0).sum(axis=0)

# Shannon diversity
def shannon_diversity(abundances):
    # Convert to proportions
    proportions = abundances / abundances.sum()
    # Remove zeros
    proportions = proportions[proportions > 0]
    # Calculate Shannon
    return -np.sum(proportions * np.log(proportions))

shannon = all_species.apply(shannon_diversity, axis=0)

# Create diversity summary
diversity_df = pd.DataFrame({
    'Sample': all_species.columns,
    'Diagnosis': [all_profiles[sid]['diagnosis'] for sid in all_species.columns],
    'Species_Richness': richness.values,
    'Shannon_Diversity': shannon.values
})

print("\n" + "="*60)
print("ALPHA DIVERSITY (Species-level)")
print("="*60)
print(diversity_df.to_string(index=False))

# Save results
diversity_df.to_csv(results_dir / 'metagenomic_alpha_diversity.csv', index=False)
all_species.to_csv(results_dir / 'species_abundance_matrix.csv')

print(f"\n Results saved to {results_dir}/")

# Visualizations
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Plot 1: Species richness
ax1 = axes[0]
colors = {'nonIBD': 'green', 'CD': 'red', 'UC': 'orange'}
diagnosis_colors = [colors[d] for d in diversity_df['Diagnosis']]

ax1.bar(range(len(diversity_df)), diversity_df['Species_Richness'], 
        color=diagnosis_colors, alpha=0.7, edgecolor='black')
ax1.set_xticks(range(len(diversity_df)))
ax1.set_xticklabels(diversity_df['Sample'], rotation=45, ha='right')
ax1.set_ylabel('Number of Species', fontsize=12)
ax1.set_title('Species Richness (Metagenomic)', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3, axis='y')

# Add legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=colors[d], label=d, alpha=0.7) 
                   for d in ['nonIBD', 'CD', 'UC']]
ax1.legend(handles=legend_elements, loc='upper right')

# Plot 2: Shannon diversity
ax2 = axes[1]
ax2.bar(range(len(diversity_df)), diversity_df['Shannon_Diversity'], 
        color=diagnosis_colors, alpha=0.7, edgecolor='black')
ax2.set_xticks(range(len(diversity_df)))
ax2.set_xticklabels(diversity_df['Sample'], rotation=45, ha='right')
ax2.set_ylabel('Shannon Diversity Index', fontsize=12)
ax2.set_title('Shannon Diversity (Metagenomic)', fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(figures_dir / 'metagenomic_alpha_diversity.png', dpi=300, bbox_inches='tight')
print(f"Figure saved: {figures_dir / 'metagenomic_alpha_diversity.png'}")

print("\n" + "="*60)
print("ANALYSIS COMPLETE!")
print("="*60)
