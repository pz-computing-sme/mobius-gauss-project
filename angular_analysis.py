"""
Angular anisotropy analysis of the Möbius function on Gaussian integers.
Generates a publication‑ready PDF plot for the MPCPS article.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ============================================================
# 1. Load data
# ============================================================
df = pd.read_csv("angular_data.csv")

# Total sum of mu(alpha) over all bins (should equal S(10^8) = -1184)
S_total = df['sum_mu'].sum()
total_count = df['count'].sum()
sqrtX = 10000.0   # X^{1/2} for X = 10^8

# Expected isotropic contribution per bin
df['expected'] = (df['count'] / total_count) * S_total

# Normalised deviation
df['delta'] = (df['sum_mu'] - df['expected']) / sqrtX

# ============================================================
# 2. Set Seaborn style (professional, publication-ready)
# ============================================================
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.8, rc={"lines.linewidth": 2.5})

# ============================================================
# 3. Create figure
# ============================================================
fig, ax = plt.subplots(figsize=(10, 6))

# Main plot: delta vs angle
ax.plot(df['angle_deg'], df['delta'], 
        color='#2E86AB', linewidth=2.5, label=r'$\Delta_{\theta} / X^{1/2}$')

# Zero reference line
ax.axhline(y=0, color='black', linestyle='--', linewidth=1.5, alpha=0.6)

# Highlight axes (0°, 90°, 180°, 270°) and diagonals (45°, 135°, 225°, 315°)
for angle in [0, 90, 180, 270]:
    ax.axvline(x=angle, color='#A23B72', linestyle=':', linewidth=1.5, alpha=0.5, label='Axes' if angle == 0 else "")
for angle in [45, 135, 225, 315]:
    ax.axvline(x=angle, color='#F18F01', linestyle=':', linewidth=1.5, alpha=0.5, label='Diagonals' if angle == 45 else "")

# Labels and title
ax.set_xlabel(r'Angle $\theta$ (degrees)', fontsize=16)
ax.set_ylabel(r'$\Delta_{\theta} / X^{1/2}$', fontsize=16)
ax.set_title('Angular Equidistribution of the Möbius Function on $\mathbb{Z}[i]$', fontsize=18, pad=15)

# Legend (only show unique labels)
handles, labels = ax.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax.legend(by_label.values(), by_label.keys(), loc='upper right', fontsize=12)

# Grid styling
ax.grid(True, linestyle='--', alpha=0.4)

# ============================================================
# 4. Add a small statistics box (optional but elegant)
# ============================================================
stats_text = (
    f"$X = 10^8$\n"
    f"$S(X) = {S_total:.0f}$\n"
    f"$\\sigma(\\Delta) = {df['delta'].std():.4f}$\n"
    f"$\\Delta_{{\\max}} = {df['delta'].max():.4f}$\n"
    f"$\\Delta_{{\\min}} = {df['delta'].min():.4f}$"
)
ax.text(0.98, 0.02, stats_text, transform=ax.transAxes,
        fontsize=11, verticalalignment='bottom', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# ============================================================
# 5. Save figure
# ============================================================
plt.tight_layout()
plt.savefig('angular_equidistribution.pdf', dpi=300, bbox_inches='tight')
print("Plot saved as 'angular_equidistribution.pdf'")

# ============================================================
# 6. Optional: also show the plot on screen
# ============================================================
plt.show()