import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

print("="*60)
print("16S vs METAGENOMICS COMPARISON")
print("="*60)

# File paths
results_dir = Path('../results')
figures_dir = Path('../figures')

# Load metagenomic results
mg_diversity = pd.read_csv(results_dir / 'metagenomic_alpha_diversity.csv')

print("\n16S RESULTS (From Project 1):")
print("-" * 60)
print("Dataset: Gevers et al. 2014 - 1,233 samples")
print("Method: 16S rRNA V4 region sequencing")
print("\nAlpha Diversity (Shannon Index):")
print("  Healthy: 4.46 ± 0.82 (n=313)")
print("  CD: 4.23 ± 1.11 (n=661)")
print("  Difference: 5% reduction in CD")
print("  P-value: 0.009 (significant)")
print("\nKey Findings:")
print("  - 585 bacteria significantly different (FDR < 0.05)")
print("  - Enterococcus: 14× higher in CD")
print("  - Roseburia: 71% depleted in CD")
print("  - ML Prediction: 89.4% AUC")
print("\nLimitations:")
print("  - Genus-level resolution only")
print("  - No direct functional measurement")

print("\n" + "="*60)
print("METAGENOMIC RESULTS (Project 2):")
print("-" * 60)
print("Dataset: HMP2 (Lloyd-Price et al. 2019) - 5 samples")
print("Method: Whole genome shotgun sequencing")
print("\nAlpha Diversity (Species-level Shannon):")
for _, row in mg_diversity.iterrows():
    print(f"  {row['Sample']} ({row['Diagnosis']}): {row['Shannon_Diversity']:.2f} "
          f"({int(row['Species_Richness'])} species)")

print("\nKey Findings:")
print("  - 312 unique metabolic pathways detected")
print("  - HSM6XRQE: Catastrophic dysbiosis (9 species, 9 pathways)")
print("  - Butyrate pathways detected in most samples")
print("  - Functional capacity varies 1000-fold between samples!")

# Calculate stats
healthy_mg = mg_diversity[mg_diversity['Diagnosis'] == 'nonIBD']['Shannon_Diversity'].values
cd_mg = mg_diversity[mg_diversity['Diagnosis'] == 'CD']['Shannon_Diversity'].values
uc_mg = mg_diversity[mg_diversity['Diagnosis'] == 'UC']['Shannon_Diversity'].values

print("\n" + "="*60)
print("COMPARISON SUMMARY")
print("="*60)

print("\n1. RESOLUTION:")
print("  16S: Genus-level (8,513 OTUs)")
print("  Metagenomics: Species/strain-level (109 species)")
print("  → Metagenomics provides strain-level precision")

print("\n2. SHANNON DIVERSITY VALUES:")
print("  16S - Healthy: 4.46")
print(f"  Metagenomics - Healthy: {healthy_mg[0]:.2f}")
print("  → Lower in metagenomics (expected - stricter species definition)")

print("\n3. FUNCTIONAL INFORMATION:")
print("  16S: None (taxonomy only)")
print("  Metagenomics: MEASURED directly (actual genes & pathways)")
print("  → 312 pathways with CPM abundances")

print("\n4. DISEASE HETEROGENEITY:")
print("  16S: Average differences between groups")
print("  Metagenomics: Individual variation revealed")
print(f"  → UC Shannon range: {uc_mg.min():.2f} to {uc_mg.max():.2f}")
print("  → 2.3-fold variation within same diagnosis!")

print("\n5. SAMPLE SIZE vs DEPTH:")
print("  16S: 1,233 samples, genus-level")
print("  Metagenomics: 5 samples, strain + function")
print("  → Complementary approaches!")

# Create comparison visualization
fig = plt.figure(figsize=(18, 12))
gs = fig.add_gridspec(3, 2, hspace=0.4, wspace=0.35)

# Plot 1: Shannon diversity comparison
ax1 = fig.add_subplot(gs[0, :])

# 16S data (from Project 1 summary)
categories = ['Healthy\n(16S)', 'CD\n(16S)', 'Healthy\n(MGX)', 'CD\n(MGX)', 
              'UC\n(MGX)\nMild', 'UC\n(MGX)\nSevere']
values = [4.46, 4.23, healthy_mg[0], cd_mg[0], 
          uc_mg[uc_mg != uc_mg.min()].mean(), uc_mg.min()]
colors_list = ['green', 'red', 'green', 'red', 'orange', 'darkred']

bars = ax1.bar(range(len(categories)), values, color=colors_list, alpha=0.7, 
               edgecolor='black', linewidth=2)
ax1.set_xticks(range(len(categories)))
ax1.set_xticklabels(categories, fontsize=11)
ax1.set_ylabel('Shannon Diversity Index', fontsize=13, fontweight='bold')
ax1.set_title('Alpha Diversity: 16S vs Metagenomics', fontsize=15, fontweight='bold')
ax1.axhline(y=4.23, color='red', linestyle='--', alpha=0.3, linewidth=2, label='16S CD mean')
ax1.axhline(y=4.46, color='green', linestyle='--', alpha=0.3, linewidth=2, label='16S Healthy mean')
ax1.grid(True, alpha=0.3, axis='y')
ax1.legend(loc='upper right')

# Add value labels on bars
for i, (bar, val) in enumerate(zip(bars, values)):
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height,
             f'{val:.2f}', ha='center', va='bottom', fontweight='bold')

# Plot 2: Species vs OTU richness
ax2 = fig.add_subplot(gs[1, 0])

richness_data = pd.DataFrame({
    'Method': ['16S\n(avg)', 'Metagenomics\n(healthy)', 'Metagenomics\n(CD)', 
               'Metagenomics\n(UC avg)'],
    'Richness': [400, int(mg_diversity[mg_diversity['Diagnosis']=='nonIBD']['Species_Richness'].values[0]), 
                 int(mg_diversity[mg_diversity['Diagnosis']=='CD']['Species_Richness'].values[0]), 
                 int(mg_diversity[mg_diversity['Diagnosis']=='UC']['Species_Richness'].mean())],
    'Color': ['blue', 'green', 'red', 'orange']
})

ax2.bar(range(len(richness_data)), richness_data['Richness'], 
        color=richness_data['Color'], alpha=0.7, edgecolor='black')
ax2.set_xticks(range(len(richness_data)))
ax2.set_xticklabels(richness_data['Method'], fontsize=10)
ax2.set_ylabel('Feature Count', fontsize=11, fontweight='bold')
ax2.set_title('Taxonomic Resolution: OTUs vs Species', fontsize=13, fontweight='bold')
ax2.grid(True, alpha=0.3, axis='y')

# Add value labels
for i, val in enumerate(richness_data['Richness']):
    ax2.text(i, val, f'{val}', ha='center', va='bottom', fontweight='bold')

# Plot 3: Information content
ax3 = fig.add_subplot(gs[1, 1])

info_categories = ['Taxonomy\nOnly\n(16S)', 'Taxonomy\n+\nStrains\n(MGX)', 
                   'Functions\n+\nPathways\n(MGX)', 'Genes\n+\nAMR\n(MGX)']
info_values = [1, 2.5, 4, 5]
info_colors = ['lightblue', 'orange', 'darkorange', 'darkgreen']

bars3 = ax3.bar(range(len(info_categories)), info_values, 
                color=info_colors, alpha=0.7, edgecolor='black')
ax3.set_xticks(range(len(info_categories)))
ax3.set_xticklabels(info_categories, fontsize=9)
ax3.set_ylabel('Information Depth\n(Relative Scale)', fontsize=11, fontweight='bold')
ax3.set_title('Data Richness: 16S vs Metagenomics', fontsize=13, fontweight='bold')
ax3.set_ylim(0, 6)
ax3.grid(True, alpha=0.3, axis='y')

# Add value labels
for bar, val in zip(bars3, info_values):
    height = bar.get_height()
    ax3.text(bar.get_x() + bar.get_width()/2., height,
             f'{val:.1f}×', ha='center', va='bottom', fontweight='bold')

# Plot 4: Sample size vs depth trade-off
ax4 = fig.add_subplot(gs[2, 0])

methods = ['16S\nProject 1', 'Metagenomics\nProject 2']
sample_sizes = [1233, 5]
depth_scores = [3, 9]  # Arbitrary scale

x = np.arange(len(methods))
width = 0.35

bars1 = ax4.bar(x - width/2, np.log10(sample_sizes), width, 
                label='Sample Size (log10)', color='steelblue', alpha=0.7)
bars2 = ax4.bar(x + width/2, depth_scores, width, 
                label='Data Depth (1-10)', color='coral', alpha=0.7)

ax4.set_ylabel('Score', fontsize=11, fontweight='bold')
ax4.set_title('Sample Size vs Data Depth Trade-off', fontsize=13, fontweight='bold')
ax4.set_xticks(x)
ax4.set_xticklabels(methods, fontsize=11)
ax4.legend()
ax4.grid(True, alpha=0.3, axis='y')

# Add value labels
for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height,
                 f'{height:.1f}', ha='center', va='bottom', fontweight='bold', fontsize=9)

# Plot 5: Key strengths comparison
ax5 = fig.add_subplot(gs[2, 1])
ax5.axis('off')

strengths_text = """
16S STRENGTHS:
✓ Large sample size (n=1,233)
✓ Statistical power for detection
✓ Cost-effective ($50-100/sample)
✓ Standardized protocols
✓ Fast turnaround (2-3 days)
✓ Population-level insights
✓ Genus-level often sufficient

METAGENOMICS STRENGTHS:
✓ Strain-level resolution
✓ Direct functional measurement
✓ 312 pathways quantified
✓ Gene content analysis
✓ AMR gene detection
✓ Individual precision medicine
✓ Reveals individual variation

BEST PRACTICE:
- Use 16S for discovery (large n)
- Use metagenomics for validation
- Combine for comprehensive insights
"""

ax5.text(0.05, 0.95, strengths_text, transform=ax5.transAxes,
         fontsize=9.5, verticalalignment='top', fontfamily='monospace',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

plt.suptitle('Project 1 (16S) vs Project 2 (Metagenomics): Comprehensive Comparison',
             fontsize=16, fontweight='bold', y=0.99)

plt.savefig(figures_dir / '16S_vs_metagenomics_comparison.png', dpi=300, bbox_inches='tight')
print(f"\n Comparison figure saved: {figures_dir / '16S_vs_metagenomics_comparison.png'}")

print("\n" + "="*60)
print("COMPARISON COMPLETE!")
print("="*60)
