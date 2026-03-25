import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

print("="*60)
print("GENE FUNCTION ANNOTATION")
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

print("\nAnalyzing EC (Enzyme Commission) annotated genes")
print("EC numbers provide functional classification of enzymes")

# Read EC-annotated gene files
ec_data = {}

for sample_id, diagnosis in samples.items():
    # Find the EC file
    filepath = data_dir / f"{sample_id}_level4ec.tsv"
    
    if filepath.exists():
        print(f"\nReading {sample_id} ({diagnosis})...")
        
        # Read EC data
        df = pd.read_csv(filepath, sep='\t', comment='#')
        df.columns = ['EC_number', 'Abundance_CPM']
        
        # Remove stratified and unmapped
        df = df[~df['EC_number'].str.contains(r'\|', na=False)]
        df = df[~df['EC_number'].str.contains('UNMAPPED|UNINTEGRATED', na=False)]
        
        ec_data[sample_id] = {
            'data': df,
            'diagnosis': diagnosis
        }
        
        print(f"  {len(df)} EC enzymes detected")
        print(f"  Total abundance: {df['Abundance_CPM'].sum():.1f} CPM")

# Create combined EC matrix
print("\n" + "="*60)
print("EC ENZYME CLASSIFICATION")
print("="*60)

all_ec = pd.DataFrame()

for sample_id, info in ec_data.items():
    temp = info['data'].copy()
    temp.columns = ['EC_number', sample_id]
    
    if all_ec.empty:
        all_ec = temp
    else:
        all_ec = all_ec.merge(temp, on='EC_number', how='outer')

all_ec = all_ec.fillna(0)
all_ec = all_ec.set_index('EC_number')

print(f"\nTotal unique EC enzymes: {len(all_ec)}")

# EC classification
print("\n" + "="*60)
print("EC CLASS DISTRIBUTION")
print("="*60)

print("\nEC Classes:")
print("  1.x.x.x - Oxidoreductases (electron transfer)")
print("  2.x.x.x - Transferases (group transfer)")
print("  3.x.x.x - Hydrolases (hydrolysis)")
print("  4.x.x.x - Lyases (bond cleavage)")
print("  5.x.x.x - Isomerases (isomerization)")
print("  6.x.x.x - Ligases (bond formation)")

# Extract EC class
all_ec_reset = all_ec.reset_index()
all_ec_reset['EC_class'] = all_ec_reset['EC_number'].str.extract(r'^(\d+)\.')

# Sum by class for each sample
for sample in samples.keys():
    print(f"\n{sample} ({samples[sample]}):")
    class_sums = all_ec_reset.groupby('EC_class')[sample].sum().sort_values(ascending=False)
    
    ec_names = {
        '1': 'Oxidoreductases',
        '2': 'Transferases',
        '3': 'Hydrolases',
        '4': 'Lyases',
        '5': 'Isomerases',
        '6': 'Ligases'
    }
    
    for ec_class, abundance in class_sums.items():
        if ec_class in ec_names:
            print(f"  EC {ec_class} ({ec_names[ec_class]}): {abundance:.1f} CPM")

# Search for specific functional enzymes
print("\n" + "="*60)
print("KEY ENZYMES OF INTEREST")
print("="*60)

# Butyrate kinase (EC 2.7.2.7)
print("\n1. BUTYRATE KINASE (EC 2.7.2.7)")
print("   Function: Final step of butyrate synthesis")

if '2.7.2.7' in all_ec.index:
    for sample in samples.keys():
        val = all_ec.loc['2.7.2.7', sample]
        print(f"   {sample} ({samples[sample]}): {val:.2f} CPM")
else:
    print("   Not detected")

# Butyryl-CoA dehydrogenase (EC 1.3.8.1)
print("\n2. BUTYRYL-CoA DEHYDROGENASE (EC 1.3.8.1)")
print("   Function: Butyrate synthesis pathway")

if '1.3.8.1' in all_ec.index:
    for sample in samples.keys():
        val = all_ec.loc['1.3.8.1', sample]
        print(f"   {sample} ({samples[sample]}): {val:.2f} CPM")
else:
    print("   Not detected")

# Beta-lactamase (EC 3.5.2.6)
print("\n3. BETA-LACTAMASE (EC 3.5.2.6)")
print("   Function: Antibiotic resistance (degrades penicillin)")

if '3.5.2.6' in all_ec.index:
    for sample in samples.keys():
        val = all_ec.loc['3.5.2.6', sample]
        print(f"   {sample} ({samples[sample]}): {val:.2f} CPM")
else:
    print("   Not detected (good - no resistance)")

# Lactate dehydrogenase (EC 1.1.1.27)
print("\n4. LACTATE DEHYDROGENASE (EC 1.1.1.27)")
print("   Function: Fermentation")

if '1.1.1.27' in all_ec.index:
    for sample in samples.keys():
        val = all_ec.loc['1.1.1.27', sample]
        print(f"   {sample} ({samples[sample]}): {val:.2f} CPM")
else:
    print("   Not detected")

# Top 20 most abundant enzymes
print("\n" + "="*60)
print("TOP 20 MOST ABUNDANT ENZYMES")
print("="*60)

all_ec['mean_abundance'] = all_ec.mean(axis=1)
top20_ec = all_ec.nlargest(20, 'mean_abundance')

for idx, ec in enumerate(top20_ec.index, 1):
    mean_val = top20_ec.loc[ec, 'mean_abundance']
    print(f"\n{idx}. EC {ec}")
    print(f"   Mean abundance: {mean_val:.2f} CPM")

# Save results
top20_ec.drop('mean_abundance', axis=1).to_csv(results_dir / 'top20_ec_enzymes.csv')
all_ec_reset.to_csv(results_dir / 'all_ec_enzymes.csv', index=False)

print(f"\n✅ Results saved to {results_dir}/")

# Visualization
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Plot 1: EC richness
ax1 = axes[0, 0]
ec_richness = (all_ec.drop('mean_abundance', axis=1) > 0).sum(axis=0)

colors = {'nonIBD': 'green', 'CD': 'red', 'UC': 'orange'}
bar_colors = [colors[samples[s]] for s in ec_richness.index]

ax1.bar(range(len(ec_richness)), ec_richness.values,
        color=bar_colors, alpha=0.7, edgecolor='black')
ax1.set_xticks(range(len(ec_richness)))
ax1.set_xticklabels([f"{s}\n({samples[s]})" for s in ec_richness.index],
                     rotation=45, ha='right')
ax1.set_ylabel('Number of EC Enzymes', fontsize=12)
ax1.set_title('Enzyme Richness by Sample', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3, axis='y')

# Plot 2: EC class distribution (stacked bar)
ax2 = axes[0, 1]

class_data = []
for sample in samples.keys():
    class_sums = all_ec_reset.groupby('EC_class')[sample].sum()
    class_data.append(class_sums)

class_df = pd.DataFrame(class_data, index=samples.keys()).fillna(0)
class_df.columns = [f"EC {c}" for c in class_df.columns]

class_df.plot(kind='bar', stacked=True, ax=ax2, 
              color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b'],
              alpha=0.8)
ax2.set_xticklabels([f"{s}\n({samples[s]})" for s in class_df.index],
                     rotation=45, ha='right')
ax2.set_ylabel('Total EC Abundance (CPM)', fontsize=12)
ax2.set_title('EC Class Distribution', fontsize=14, fontweight='bold')
ax2.legend(title='EC Class', bbox_to_anchor=(1.05, 1), loc='upper left')
ax2.grid(True, alpha=0.3, axis='y')

# Plot 3: Total functional capacity
ax3 = axes[1, 0]

total_capacity = all_ec.drop('mean_abundance', axis=1).sum(axis=0)

ax3.bar(range(len(total_capacity)), total_capacity.values,
        color=bar_colors, alpha=0.7, edgecolor='black')
ax3.set_xticks(range(len(total_capacity)))
ax3.set_xticklabels([f"{s}\n({samples[s]})" for s in total_capacity.index],
                     rotation=45, ha='right')
ax3.set_ylabel('Total EC Abundance (CPM)', fontsize=12)
ax3.set_title('Total Enzymatic Capacity', fontsize=14, fontweight='bold')
ax3.grid(True, alpha=0.3, axis='y')

# Add value labels
for i, val in enumerate(total_capacity.values):
    ax3.text(i, val, f'{val:.0f}', ha='center', va='bottom', fontweight='bold')

# Plot 4: Top 10 enzymes heatmap
ax4 = axes[1, 1]

top10_for_plot = all_ec.nlargest(10, 'mean_abundance').drop('mean_abundance', axis=1)

sns.heatmap(top10_for_plot, cmap='YlOrRd', annot=True, fmt='.0f',
            cbar_kws={'label': 'Abundance (CPM)'}, ax=ax4,
            yticklabels=[f"EC {ec}" for ec in top10_for_plot.index])
ax4.set_title('Top 10 Most Abundant Enzymes', fontsize=14, fontweight='bold')
ax4.set_xlabel('Sample', fontsize=12)
ax4.set_xticklabels([f"{s}\n({samples[s]})" for s in top10_for_plot.columns],
                     rotation=45, ha='right')

plt.tight_layout()
plt.savefig(figures_dir / 'ec_enzyme_analysis.png', dpi=300, bbox_inches='tight')
print(f" Figure saved: {figures_dir / 'ec_enzyme_analysis.png'}")

print("\n" + "="*60)
print("EC ANALYSIS COMPLETE!")
print("="*60)
