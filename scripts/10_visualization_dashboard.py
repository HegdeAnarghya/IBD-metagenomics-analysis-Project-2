import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

print("=== Session 7: Comprehensive Visualization Dashboard ===\n")

# Sample metadata
sample_labels = {
    'MSM5LLHV': 'UC',
    'HSM7CZ2A': 'nonIBD',
    'HSM6XRQE': 'UC',
    'CSM5FZ4C': 'CD',
    'CSM9X1ZO': 'UC'
}
sample_colors = {'nonIBD': '#2ecc71', 'UC': '#e74c3c', 'CD': '#e67e22'}

# Load all result files
print("Loading results...")
species = pd.read_csv('../results/species_abundance_matrix.csv', index_col=0)
species.columns = [c[:-2] if c.endswith('_P') else c for c in species.columns]

pathways = pd.read_csv('../results/top20_pathways.csv', index_col=0)
pathways.columns = [c[:-2] if c.endswith('_P') else c for c in pathways.columns]

diversity = pd.read_csv('../results/metagenomic_alpha_diversity.csv', index_col=0)
diversity.index = [i[:-2] if str(i).endswith('_P') else str(i) for i in diversity.index]

feat_imp = pd.read_csv('../results/pathway_feature_importances.csv')

sp_pw_corr = pd.read_csv('../results/species_pathway_correlation.csv', index_col=0)

print("All files loaded.")

# ============================================================
# BUILD DASHBOARD FIGURE
# ============================================================
fig = plt.figure(figsize=(20, 22))
fig.patch.set_facecolor('#f8f9fa')
gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.45, wspace=0.35)

fig.suptitle('IBD Metagenomics Analysis Dashboard\nHMP2 Data: 5 samples (1 nonIBD, 3 UC, 1 CD)',
             fontsize=16, fontweight='bold', y=0.98)

# ── PANEL 1: Alpha Diversity (top-left) ──────────────────────
ax1 = fig.add_subplot(gs[0, 0])
div_samples = [s for s in diversity.index if s in sample_labels]
div_col = [c for c in diversity.columns if 'shannon' in c.lower() or 'diversity' in c.lower() or 'Shannon' in c]
if not div_col:
    div_col = [diversity.columns[0]]
div_col = div_col[0]

div_vals = diversity.loc[div_samples, div_col]
bar_colors = [sample_colors[sample_labels[s]] for s in div_samples]
bars = ax1.bar(range(len(div_samples)), div_vals.values, color=bar_colors,
               edgecolor='white', linewidth=1.5, width=0.6)
ax1.set_xticks(range(len(div_samples)))
ax1.set_xticklabels([f"{s}\n({sample_labels[s]})" for s in div_samples], fontsize=8)
ax1.set_ylabel('Shannon Diversity Index', fontsize=10)
ax1.set_title('Panel 1: Alpha Diversity by Sample', fontweight='bold', fontsize=11)
ax1.set_facecolor('white')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
for bar, val in zip(bars, div_vals.values):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
             f'{val:.2f}', ha='center', va='bottom', fontsize=8, fontweight='bold')
# Legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=v, label=k) for k, v in sample_colors.items()]
ax1.legend(handles=legend_elements, loc='upper right', fontsize=8)
print("Panel 1 done.")

# ── PANEL 2: Top 10 Species Abundance (top-right) ────────────
ax2 = fig.add_subplot(gs[0, 1])
samples_ordered = ['HSM7CZ2A', 'MSM5LLHV', 'HSM6XRQE', 'CSM5FZ4C', 'CSM9X1ZO']
samples_ordered = [s for s in samples_ordered if s in species.columns]
top10_species = species[samples_ordered].mean(axis=1).nlargest(10).index
sp_sub = species.loc[top10_species, samples_ordered]

sp_sub_plot = sp_sub.copy()
sp_sub_plot.index = [i.split('|')[-1].replace('s__', '').replace('_', ' ')
                     if '|' in i else i.replace('s__', '').replace('_', ' ')
                     for i in sp_sub_plot.index]

colors_panel2 = plt.cm.Set3(np.linspace(0, 1, len(top10_species)))
bottom = np.zeros(len(samples_ordered))
for i, (idx, row) in enumerate(sp_sub_plot.iterrows()):
    ax2.bar(range(len(samples_ordered)), row.values, bottom=bottom,
            color=colors_panel2[i], label=idx[:20], edgecolor='white', linewidth=0.5)
    bottom += row.values

ax2.set_xticks(range(len(samples_ordered)))
ax2.set_xticklabels([f"{s}\n({sample_labels.get(s,'?')})" for s in samples_ordered], fontsize=8)
ax2.set_ylabel('Abundance (CPM)', fontsize=10)
ax2.set_title('Panel 2: Top 10 Species Abundance', fontweight='bold', fontsize=11)
ax2.set_facecolor('white')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.legend(loc='upper right', fontsize=6, ncol=1, bbox_to_anchor=(1.35, 1.0))
print("Panel 2 done.")

# ── PANEL 3: Top 10 Pathway Abundances heatmap (middle-left) ─
ax3 = fig.add_subplot(gs[1, 0])
pw_samples = [s for s in samples_ordered if s in pathways.columns]
pw_top10 = pathways[pw_samples].mean(axis=1).nlargest(10).index
pw_sub = pathways.loc[pw_top10, pw_samples]
pw_sub.index = [i.split(':')[0] for i in pw_sub.index]

pw_norm = pw_sub.div(pw_sub.max(axis=1) + 1e-10, axis=0)
sns.heatmap(pw_norm, ax=ax3, cmap='YlOrRd', linewidths=0.5,
            xticklabels=[f"{s}\n({sample_labels.get(s,'?')})" for s in pw_samples],
            yticklabels=pw_sub.index, cbar_kws={'label': 'Normalized\nAbundance'},
            annot=False)
ax3.set_title('Panel 3: Top 10 Pathway Abundances\n(Normalized)', fontweight='bold', fontsize=11)
ax3.tick_params(axis='x', labelsize=8)
ax3.tick_params(axis='y', labelsize=7)
print("Panel 3 done.")

# ── PANEL 4: ML Feature Importances (middle-right) ───────────
ax4 = fig.add_subplot(gs[1, 1])
top10_feat = feat_imp.head(10).copy()
top10_feat['label'] = [p.split(':')[0] for p in top10_feat['pathway']]
colors_feat = ['#e74c3c' if i < 5 else '#3498db' for i in range(len(top10_feat))]
bars4 = ax4.barh(range(len(top10_feat)), top10_feat['importance'].values,
                 color=colors_feat, edgecolor='white')
ax4.set_yticks(range(len(top10_feat)))
ax4.set_yticklabels(top10_feat['label'], fontsize=8)
ax4.invert_yaxis()
ax4.set_xlabel('Mean Feature Importance (RF)', fontsize=10)
ax4.set_title('Panel 4: Top ML-Predictive Pathways\n(Red=Top 5)', fontweight='bold', fontsize=11)
ax4.set_facecolor('white')
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.axvline(x=top10_feat['importance'].mean(), color='gray', linestyle='--', alpha=0.7)
print("Panel 4 done.")

# ── PANEL 5: Species-Pathway Correlation (bottom-left) ───────
ax5 = fig.add_subplot(gs[2, 0])
# Pick top 8 species and top 8 pathways by variance
sp_var = sp_pw_corr.var(axis=1).nlargest(8).index
pw_var = sp_pw_corr.var(axis=0).nlargest(8).index
corr_sub = sp_pw_corr.loc[sp_var, pw_var]
corr_sub.index = [i.split('|')[-1].replace('s__', '').replace('_', ' ')[:18]
                  if '|' in i else i.replace('s__', '').replace('_', ' ')[:18]
                  for i in corr_sub.index]
corr_sub.columns = [c.split(':')[0] for c in corr_sub.columns]

sns.heatmap(corr_sub, ax=ax5, cmap='RdBu_r', center=0,
            linewidths=0.5, annot=True, fmt='.2f', annot_kws={'size': 7},
            cbar_kws={'label': "Spearman r"})
ax5.set_title('Panel 5: Species-Pathway Correlations\n(Spearman r)', fontweight='bold', fontsize=11)
ax5.tick_params(axis='x', labelsize=7, rotation=45)
ax5.tick_params(axis='y', labelsize=7)
print("Panel 5 done.")

# ── PANEL 6: Project Summary Stats Table (bottom-right) ──────
ax6 = fig.add_subplot(gs[2, 1])
ax6.axis('off')

table_data = [
    ['Metric', 'Value'],
    ['Total samples', '5'],
    ['Diagnosis breakdown', '3 UC, 1 CD, 1 nonIBD'],
    ['Species detected', '109'],
    ['Pathways analyzed', '20'],
    ['EC enzymes decoded', '1,537'],
    ['Top species (IBD+)', 'Bacteroides vulgatus'],
    ['Top species (IBD−)', 'Faecalibacterium prausnitzii'],
    ['Top ML pathway', 'tRNA-CHARGING-PWY'],
    ['LOOCV AUC (n=5)', '0.000 (class imbalance)'],
    ['16S ML baseline AUC', '0.894 (larger dataset)'],
    ['Data source', 'NIH HMP2 / IBDMDB'],
]

table = ax6.table(cellText=table_data[1:], colLabels=table_data[0],
                  cellLoc='left', loc='center',
                  colWidths=[0.52, 0.48])
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 1.6)

for (row, col), cell in table.get_celld().items():
    if row == 0:
        cell.set_facecolor('#2c3e50')
        cell.set_text_props(color='white', fontweight='bold')
    elif row % 2 == 0:
        cell.set_facecolor('#ecf0f1')
    else:
        cell.set_facecolor('white')
    cell.set_edgecolor('#bdc3c7')

ax6.set_title('Panel 6: Project Summary', fontweight='bold', fontsize=11)
print("Panel 6 done.")

# Save
plt.savefig('../figures/session7_dashboard.png', dpi=150, bbox_inches='tight',
            facecolor='#f8f9fa')
plt.close()
print("\nSaved: session7_dashboard.png")
print("\n=== Session 7 Complete! ===")
print("Dashboard saved to figures/session7_dashboard.png")
