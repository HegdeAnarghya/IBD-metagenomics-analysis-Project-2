import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix
from sklearn.preprocessing import LabelEncoder
import warnings
warnings.filterwarnings('ignore')

print("=== Session 6: Machine Learning on Pathway Data ===\n")

# Load pathway abundance data
print("Loading data...")
pathways = pd.read_csv('../results/top20_pathways.csv', index_col=0)
print(f"Pathway matrix shape: {pathways.shape}")
print(f"Samples: {list(pathways.columns)}")
pathways.columns = [c[:-2] if c.endswith('_P') else c for c in pathways.columns]

# Define labels: IBD=1, nonIBD=0
sample_labels = {
    'MSM5LLHV': 1,  # UC
    'HSM7CZ2A': 0,  # nonIBD
    'HSM6XRQE': 1,  # UC
    'CSM5FZ4C': 1,  # CD
    'CSM9X1ZO': 1   # UC
}

# Align samples
samples = [s for s in pathways.columns if s in sample_labels]
X = pathways[samples].T  # shape: samples x features
y = np.array([sample_labels[s] for s in samples])
print(f"\nFeature matrix: {X.shape[0]} samples x {X.shape[1]} features")
print(f"Labels: {dict(zip(samples, y))}")
print(f"Class balance: {sum(y==0)} nonIBD, {sum(y==1)} IBD")

# Leave-One-Out Cross Validation
print("\nRunning Leave-One-Out Cross Validation...")
loo = LeaveOneOut()
rf = RandomForestClassifier(n_estimators=100, random_state=42, max_features='sqrt')

y_true = []
y_pred_proba = []
feature_importances = np.zeros(X.shape[1])

for train_idx, test_idx in loo.split(X):
    X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]
    
    rf.fit(X_train, y_train)
    proba = rf.predict_proba(X_test)
    
    # Handle case where only one class in training
    if len(rf.classes_) == 2:
        y_pred_proba.append(proba[0][1])
    else:
        # Only one class seen in training - predict that class's probability
        y_pred_proba.append(float(rf.classes_[0]))
    
    y_true.append(y_test[0])
    feature_importances += rf.feature_importances_

feature_importances /= len(samples)

y_true = np.array(y_true)
y_pred_proba = np.array(y_pred_proba)

print(f"True labels:       {y_true}")
print(f"Predicted probas:  {np.round(y_pred_proba, 3)}")

# Calculate AUC
try:
    auc = roc_auc_score(y_true, y_pred_proba)
    print(f"\nLOOCV AUC: {auc:.3f}")
    print(f"16S baseline AUC: 0.894")
    if auc > 0.894:
        print("Pathway ML OUTPERFORMS 16S baseline!")
    elif auc == 0.894:
        print("Pathway ML matches 16S baseline.")
    else:
        print("16S baseline outperforms pathway ML (expected with only 5 samples).")
except Exception as e:
    print(f"AUC calculation note: {e}")
    auc = 0.5

# Feature importances
feat_df = pd.DataFrame({
    'pathway': X.columns,
    'importance': feature_importances
}).sort_values('importance', ascending=False)
print(f"\nTop 5 most predictive pathways:")
print(feat_df.head())

# Save results
feat_df.to_csv('../results/pathway_feature_importances.csv', index=False)
print("\nSaved: pathway_feature_importances.csv")

# --- FIGURE 1: Feature Importances ---
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
fig.suptitle('Session 6: ML on Pathway Data (LOOCV, n=5)', fontsize=14, fontweight='bold')

# Plot 1: Feature importances
colors = ['#e74c3c' if i < 5 else '#3498db' for i in range(len(feat_df))]
axes[0].barh(range(len(feat_df)), feat_df['importance'].values, color=colors)
axes[0].set_yticks(range(len(feat_df)))
axes[0].set_yticklabels([p[:25] for p in feat_df['pathway']], fontsize=8)
axes[0].invert_yaxis()
axes[0].set_xlabel('Mean Feature Importance')
axes[0].set_title('Pathway Feature Importances\n(Red = Top 5 predictors)')
axes[0].axvline(x=feat_df['importance'].mean(), color='gray', linestyle='--', alpha=0.7, label='Mean')
axes[0].legend()

# Plot 2: ROC Curve
try:
    fpr, tpr, _ = roc_curve(y_true, y_pred_proba)
    axes[1].plot(fpr, tpr, color='#e74c3c', linewidth=2.5, label=f'Pathway ML (AUC={auc:.3f})')
    axes[1].plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Random (AUC=0.500)')
    axes[1].axhline(y=0.894, color='#3498db', linestyle=':', linewidth=2, label='16S Baseline (AUC=0.894)')
    axes[1].set_xlabel('False Positive Rate')
    axes[1].set_ylabel('True Positive Rate')
    axes[1].set_title('ROC Curve: Pathway ML vs 16S Baseline')
    axes[1].legend(loc='lower right')
    axes[1].set_xlim([0, 1])
    axes[1].set_ylim([0, 1.05])
except Exception as e:
    axes[1].text(0.5, 0.5, f'ROC not available\n(n=5 limitation)', 
                ha='center', va='center', transform=axes[1].transAxes, fontsize=12)
    axes[1].set_title('ROC Curve')

plt.tight_layout()
plt.savefig('../figures/session6_fig1_ml_results.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: session6_fig1_ml_results.png")

# --- FIGURE 2: Sample Prediction Summary ---
fig, ax = plt.subplots(figsize=(10, 5))
x_pos = range(len(samples))
bar_colors = ['#2ecc71' if l == 0 else '#e74c3c' for l in y_true]
bars = ax.bar(x_pos, y_pred_proba, color=bar_colors, alpha=0.8, edgecolor='black')
ax.axhline(y=0.5, color='black', linestyle='--', alpha=0.7, label='Decision boundary (0.5)')
ax.set_xticks(x_pos)
ax.set_xticklabels([f"{s}\n({'IBD' if l==1 else 'nonIBD'})" for s, l in zip(samples, y_true)])
ax.set_ylabel('Predicted IBD Probability')
ax.set_title('LOOCV Predictions per Sample\n(Green=nonIBD, Red=IBD actual label)')
ax.set_ylim([0, 1.1])
ax.legend()

for bar, prob in zip(bars, y_pred_proba):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02, 
            f'{prob:.2f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

plt.tight_layout()
plt.savefig('../figures/session6_fig2_sample_predictions.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: session6_fig2_sample_predictions.png")

print("\n=== Session 6 Complete! ===")
print(f"AUC: {auc:.3f} (pathway ML) vs 0.894 (16S baseline)")
print("Figures: session6_fig1_ml_results.png, session6_fig2_sample_predictions.png")
print("CSV: pathway_feature_importances.csv")
