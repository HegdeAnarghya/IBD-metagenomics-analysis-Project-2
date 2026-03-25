#!/usr/bin/env python3
"""
=============================================================================
SESSION 5: Species-Function Correlation & Co-abundance Network Analysis
IBD Metagenomics Project
=============================================================================
WHAT THIS DOES:
  Part A: Which bacteria correlate with which metabolic pathways?
  Part B: Which bacteria co-occur or exclude each other?

SAMPLES:
  HSM7CZ2A -> nonIBD (Healthy)
  MSM5LLHV -> UC
  HSM6XRQE -> UC
  CSM9X1ZO -> UC
  CSM5FZ4C -> CD

NOTE: n=5 samples -> correlations are exploratory, not statistically powered.
=============================================================================
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

BASE    = '/mnt/e/IBD-metagenomics-analysis'
RESULTS = f'{BASE}/results'
FIGURES = f'{BASE}/figures'

SAMPLE_INFO = {
    'MSM5LLHV': 'UC',
    'HSM7CZ2A': 'nonIBD',
    'HSM6XRQE': 'UC',
    'CSM5FZ4C': 'CD',
    'CSM9X1ZO': 'UC'
}
COLORS = {'nonIBD': '#27ae60', 'UC': '#e74c3c', 'CD': '#8e44ad'}

def strip_P(df):
    df.columns = [c.replace('_P', '') for c in df.columns]
    return df

def corr_matrix(df_rows, df_cols):
    result = pd.DataFrame(index=df_rows.index, columns=df_cols.index, dtype=float)
    for a in df_rows.index:
        for b in df_cols.index:
            r, _ = spearmanr(df_rows.loc[a].values, df_cols.loc[b].values)
            result.loc[a, b] = r
    return result.astype(float)

# ── STEP 1: Load ──────────────────────────────────────────────────────────────
print("="*65)
print("  SESSION 5: SPECIES-FUNCTION CORRELATION ANALYSIS")
print("="*65)
print("\n[1/5] Loading data...")

species_df = strip_P(pd.read_csv(f'{RESULTS}/species_abundance_matrix.csv', index_col=0))
pathway_df = strip_P(pd.read_csv(f'{RESULTS}/top20_pathways.csv', index_col=0))
ec_df      = strip_P(pd.read_csv(f'{RESULTS}/top20_ec_enzymes.csv', index_col=0))
if 'EC_class' in ec_df.columns:
    ec_df = ec_df.drop(columns=['EC_class'])

print(f"  Species : {species_df.shape[0]} x {species_df.shape[1]}")
print(f"  Pathways: {pathway_df.shape[0]} x {pathway_df.shape[1]}")
print(f"  Enzymes : {ec_df.shape[0]} x {ec_df.shape[1]}")

# ── STEP 2: Align samples ─────────────────────────────────────────────────────
print("\n[2/5] Aligning samples...")
common = sorted(set(species_df.columns) & set(pathway_df.columns) & set(ec_df.columns))
print(f"  Common samples: {common}")
species_df = species_df[common]
pathway_df = pathway_df[common]
ec_df      = ec_df[common]

# ── STEP 3: Top 20 species ────────────────────────────────────────────────────
print("\n[3/5] Selecting top 20 species by mean abundance...")
top20     = species_df.mean(axis=1).nlargest(20).index
sp_top    = species_df.loc[top20]
print(f"  Top species: {list(top20[:3])} ...")

# ── STEP 4: Correlations ──────────────────────────────────────────────────────
print("\n[4/5] Computing Spearman correlations...")
print("  species x pathways ...")
sp_pw = corr_matrix(sp_top, pathway_df)
print("  species x EC enzymes ...")
sp_ec = corr_matrix(sp_top, ec_df)
print("  species x species (co-abundance) ...")
sp_sp = corr_matrix(sp_top, sp_top)

sp_pw.to_csv(f'{RESULTS}/species_pathway_correlation.csv')
sp_ec.to_csv(f'{RESULTS}/species_ec_correlation.csv')
sp_sp.to_csv(f'{RESULTS}/species_coabundance_correlation.csv')
print("  Saved 3 correlation CSVs.")

# ── STEP 5: Figures ───────────────────────────────────────────────────────────
print("\n[5/5] Generating figures...")
plt.rcParams.update({'font.family': 'DejaVu Sans'})

# Fig 1: Species-Pathway heatmap
print("  Fig 1: species-pathway heatmap...")
fig, ax = plt.subplots(figsize=(13, 8))
pw_labels = [p.split(':')[0] if ':' in p else p[:20] for p in sp_pw.columns]
sns.heatmap(sp_pw, ax=ax, cmap='RdBu_r', vmin=-1, vmax=1, center=0,
            xticklabels=pw_labels,
            yticklabels=[s.replace('_',' ') for s in sp_pw.index],
            linewidths=0.4, linecolor='#eeeeee',
            cbar_kws={'label': 'Spearman r', 'shrink': 0.75})
ax.set_title('Species-Pathway Correlations\n"Who is producing which metabolic pathway?"',
             fontsize=13, fontweight='bold', pad=15)
ax.set_xlabel('Metabolic Pathway', fontsize=10)
ax.set_ylabel('Bacterial Species', fontsize=10)
ax.tick_params(axis='x', rotation=50, labelsize=7)
ax.tick_params(axis='y', rotation=0,  labelsize=8)
plt.tight_layout()
plt.savefig(f'{FIGURES}/session5_fig1_species_pathway_heatmap.png', dpi=150, bbox_inches='tight')
plt.close()
print("    Saved: session5_fig1_species_pathway_heatmap.png")

# Fig 2: Species-EC heatmap
print("  Fig 2: species-EC heatmap...")
fig, ax = plt.subplots(figsize=(13, 8))
sns.heatmap(sp_ec, ax=ax, cmap='RdBu_r', vmin=-1, vmax=1, center=0,
            xticklabels=list(sp_ec.columns),
            yticklabels=[s.replace('_',' ') for s in sp_ec.index],
            linewidths=0.4, linecolor='#eeeeee',
            cbar_kws={'label': 'Spearman r', 'shrink': 0.75})
ax.set_title('Species-Enzyme (EC) Correlations\n"Which bacteria carry which enzymatic functions?"',
             fontsize=13, fontweight='bold', pad=15)
ax.set_xlabel('EC Enzyme', fontsize=10)
ax.set_ylabel('Bacterial Species', fontsize=10)
ax.tick_params(axis='x', rotation=45, labelsize=8)
ax.tick_params(axis='y', rotation=0,  labelsize=8)
plt.tight_layout()
plt.savefig(f'{FIGURES}/session5_fig2_species_ec_heatmap.png', dpi=150, bbox_inches='tight')
plt.close()
print("    Saved: session5_fig2_species_ec_heatmap.png")

# Fig 3: Co-abundance heatmap
print("  Fig 3: co-abundance heatmap...")
mask = np.eye(sp_sp.shape[0], dtype=bool)
fig, ax = plt.subplots(figsize=(11, 9))
sns.heatmap(sp_sp, ax=ax, cmap='RdBu_r', vmin=-1, vmax=1, center=0,
            mask=mask, square=True,
            xticklabels=[s.replace('_',' ') for s in sp_sp.columns],
            yticklabels=[s.replace('_',' ') for s in sp_sp.index],
            linewidths=0.3, linecolor='#eeeeee',
            cbar_kws={'label': 'Spearman r', 'shrink': 0.75})
ax.set_title('Species Co-abundance\n"Who lives with whom?" (Red=co-occur, Blue=exclude each other)',
             fontsize=13, fontweight='bold', pad=15)
ax.tick_params(axis='x', rotation=55, labelsize=7)
ax.tick_params(axis='y', rotation=0,  labelsize=7)
plt.tight_layout()
plt.savefig(f'{FIGURES}/session5_fig3_coabundance_heatmap.png', dpi=150, bbox_inches='tight')
plt.close()
print("    Saved: session5_fig3_coabundance_heatmap.png")

# Fig 4: Top associations bar chart
print("  Fig 4: top associations bar chart...")
pairs = []
for sp in sp_pw.index:
    for pw in sp_pw.columns:
        pairs.append({
            'Species':    sp.replace('_',' '),
            'Pathway':    pw.split(':')[0] if ':' in pw else pw[:22],
            'Spearman_r': sp_pw.loc[sp, pw]
        })
pairs_df = pd.DataFrame(pairs).dropna().sort_values('Spearman_r')
combined = pd.concat([pairs_df.head(10), pairs_df.tail(10)]).reset_index(drop=True)
combined['Label'] = combined['Species'] + '  ↔  ' + combined['Pathway']

fig, ax = plt.subplots(figsize=(11, 10))
bar_colors = ['#3498db' if r < 0 else '#e74c3c' for r in combined['Spearman_r']]
ax.barh(combined['Label'], combined['Spearman_r'],
        color=bar_colors, edgecolor='white', linewidth=0.5, height=0.7)
ax.axvline(x=0, color='black', linewidth=1.2)
ax.axvline(x=0.5,  color='gray', linewidth=0.7, linestyle='--', alpha=0.5)
ax.axvline(x=-0.5, color='gray', linewidth=0.7, linestyle='--', alpha=0.5)
ax.set_xlabel('Spearman Correlation (r)', fontsize=11)
ax.set_title('Top Species-Pathway Associations\nRed=positive (species likely produces pathway) | Blue=negative',
             fontsize=12, fontweight='bold', pad=12)
ax.set_xlim(-1.15, 1.15)
ax.tick_params(axis='y', labelsize=8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(f'{FIGURES}/session5_fig4_top_associations.png', dpi=150, bbox_inches='tight')
plt.close()
print("    Saved: session5_fig4_top_associations.png")

# ── SUMMARY ───────────────────────────────────────────────────────────────────
print("\n" + "="*65)
print("  SESSION 5 COMPLETE!")
print("="*65)
print(f"""
RESULTS: {RESULTS}/
  species_pathway_correlation.csv
  species_ec_correlation.csv
  species_coabundance_correlation.csv

FIGURES: {FIGURES}/
  session5_fig1_species_pathway_heatmap.png
  session5_fig2_species_ec_heatmap.png
  session5_fig3_coabundance_heatmap.png
  session5_fig4_top_associations.png

HOW TO READ YOUR RESULTS:
  Fig1: RED cell = that species correlates with that pathway being active
        BLUE cell = that species is absent when that pathway is active
  Fig3: RED pair = these species co-occur (potential mutualism/sharing)
        BLUE pair = these species exclude each other (competition)
  Fig4: The strongest individual species-pathway links — your key findings!

""")
