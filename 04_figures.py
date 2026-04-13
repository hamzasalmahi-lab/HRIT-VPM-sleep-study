"""
HRIT VPM Study — Step 4: Publication Figures
=============================================
Generates three figures for the empirical paper.

Figure 1: Group-level VPM demonstration
    Panel A: Mean RSI proxy (aperiodic slope) aligned to N3 transition
    Panel B: Mean RSI variance aligned to N3 transition
    Panel C: Overlay showing variance peak precedes mean threshold crossing

Figure 2: Individual subject lag distribution
    Panel A: Histogram of lags (positive = VPM confirmed)
    Panel B: Individual subject traces (spaghetti)

Figure 3: Power-law fit (H2, exploratory)
    β distribution vs HRIT prediction (0.546)
"""

import os
import glob
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

RESULTS_DIR    = 'results'
COMPLEXITY_DIR = 'complexity'
FIGURES_DIR    = 'figures'
os.makedirs(FIGURES_DIR, exist_ok=True)

BLUE   = '#2E75B6'
ORANGE = '#C55A11'
RED    = '#C00000'
GREEN  = '#538135'
GREY   = '#595959'
LBLUE  = '#BDD7EE'

plt.rcParams.update({
    'font.family':        ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size':          10,
    'axes.spines.top':    False,
    'axes.spines.right':  False,
    'axes.linewidth':     0.8,
    'figure.dpi':         300,
    'savefig.dpi':        300,
    'savefig.facecolor':  'white',
    'savefig.bbox':       'tight',
})

EPOCH_SEC   = 30
PRE_EPOCHS  = 20
POST_EPOCHS = 10

# ── Load results ──────────────────────────────────────────────────────────────

results_df = pd.read_csv(os.path.join(RESULTS_DIR, 'vpm_subject_results.csv'))
summary    = pd.read_csv(os.path.join(RESULTS_DIR, 'vpm_summary_stats.csv')).iloc[0]

print(f'Loaded {len(results_df)} subjects')

# ── Reconstruct time series ───────────────────────────────────────────────────
# Load individual subject time series
all_mean = []
all_var  = []
all_rel  = []

for _, row in results_df.iterrows():
    try:
        t     = np.array(json.loads(row['rsi_time']))
        m     = np.array(json.loads(row['rsi_mean']))
        v     = np.array(json.loads(row['rsi_var']))
        all_mean.append((t, m))
        all_var.append((t, v))
        all_rel.append(t)
    except Exception:
        pass

# Interpolate to common time grid for averaging
t_grid = np.arange(-PRE_EPOCHS, POST_EPOCHS + 1)
mean_mat = np.full((len(all_mean), len(t_grid)), np.nan)
var_mat  = np.full((len(all_var),  len(t_grid)), np.nan)

for i, ((t, m), (_, v)) in enumerate(zip(all_mean, all_var)):
    for j, tg in enumerate(t_grid):
        idx = np.where(t == tg)[0]
        if len(idx):
            mean_mat[i, j] = m[idx[0]]
            var_mat[i, j]  = v[idx[0]]

# Normalise each subject to 0-1 range for group average
for i in range(mean_mat.shape[0]):
    row_m = mean_mat[i]
    valid = ~np.isnan(row_m)
    if valid.sum() > 3:
        rng = np.nanmax(row_m) - np.nanmin(row_m)
        if rng > 0:
            mean_mat[i] = (row_m - np.nanmin(row_m)) / rng

    row_v = var_mat[i]
    valid = ~np.isnan(row_v)
    if valid.sum() > 3:
        rng = np.nanmax(row_v) - np.nanmin(row_v)
        if rng > 0:
            var_mat[i] = (row_v - np.nanmin(row_v)) / rng

t_sec   = t_grid * EPOCH_SEC
gm_mean = np.nanmean(mean_mat, axis=0)
gm_var  = np.nanmean(var_mat,  axis=0)
gm_mean_se = np.nanstd(mean_mat, axis=0) / np.sqrt(np.sum(~np.isnan(mean_mat), axis=0))
gm_var_se  = np.nanstd(var_mat,  axis=0) / np.sqrt(np.sum(~np.isnan(var_mat),  axis=0))

# ── Figure 1: Group-level VPM ─────────────────────────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.6))

ax1, ax2 = axes

# Panel A: RSI mean
ax1.fill_between(t_sec,
                 gm_mean - gm_mean_se,
                 gm_mean + gm_mean_se,
                 color=LBLUE, alpha=0.5)
ax1.plot(t_sec, gm_mean, color=BLUE, linewidth=2.0, label='RSI proxy (mean)')
ax1.axvline(0, color=RED, linewidth=1.0, linestyle='--', alpha=0.8)
ax1.axhline(0.33, color=GREY, linewidth=0.8, linestyle=':', alpha=0.7)
ax1.axvspan(0, t_sec[-1], alpha=0.05, color=GREY)
ax1.text(2 * EPOCH_SEC, 0.36, 'N3 threshold', fontsize=7.5,
         color=GREY, style='italic')
ax1.text(0 + 5, 0.05, 'N3 onset', fontsize=7.5,
         color=RED, rotation=90, va='bottom')
ax1.set_xlabel('Time relative to N3 onset (s)', fontsize=9.5)
ax1.set_ylabel('Normalised RSI proxy\n(aperiodic slope)', fontsize=9.5)
ax1.set_title('A   RSI mean: gradual decline', fontsize=9.5,
              fontweight='bold', loc='left', pad=4)
ax1.legend(fontsize=8, frameon=False)

# Panel B: RSI variance
ax2.fill_between(t_sec,
                 gm_var - gm_var_se,
                 gm_var + gm_var_se,
                 color='#F4B183', alpha=0.5)
ax2.plot(t_sec, gm_var, color=ORANGE, linewidth=2.0, label='RSI proxy (variance)')
ax2.axvline(0, color=RED, linewidth=1.0, linestyle='--', alpha=0.8)

# Mark predicted variance peak (pre-transition)
pre_mask = t_sec < 0
if pre_mask.sum() > 0:
    pre_var  = gm_var.copy()
    pre_var[~pre_mask] = np.nan
    peak_idx = np.nanargmax(pre_var)
    t_peak   = t_sec[peak_idx]
    ax2.axvline(t_peak, color=ORANGE, linewidth=1.2,
                linestyle='-.', alpha=0.9)
    ax2.annotate(f'Variance peak\n({t_peak:.0f}s before N3)',
                 xy=(t_peak, gm_var[peak_idx]),
                 xytext=(t_peak - 3 * EPOCH_SEC, 0.80),
                 fontsize=7.5, color=ORANGE,
                 arrowprops=dict(arrowstyle='->', color=ORANGE, lw=0.9))

ax2.set_xlabel('Time relative to N3 onset (s)', fontsize=9.5)
ax2.set_ylabel('Normalised RSI variance', fontsize=9.5)
ax2.set_title('B   Variance: rises then falls before threshold',
              fontsize=9.5, fontweight='bold', loc='left', pad=4)
ax2.legend(fontsize=8, frameon=False)
ax2.text(2 * EPOCH_SEC, 0.05, 'N3 onset', fontsize=7.5,
         color=RED, rotation=90, va='bottom')

fig.suptitle(
    'Figure 1. Variance-Precedes-Mean in EEG Complexity at the '
    'Wake→N3 Sleep Transition',
    fontsize=9.5, fontweight='bold', x=0.02, ha='left', y=1.02)

fig.text(0.01, -0.04,
         'Group-level normalised RSI proxy (aperiodic slope) mean ± SE across '
         f'N={len(results_df)} subjects. Time 0 = first sustained N3 epoch.\n'
         'Shaded region: pre-transition window (negative time). '
         'Dashed red line: N3 onset. Dotted grey: N3 mean threshold.',
         fontsize=7.5, color=GREY)

plt.tight_layout()
fig1_path = os.path.join(FIGURES_DIR, 'fig1_vpm_group.png')
plt.savefig(fig1_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f'✓ Figure 1 saved: {fig1_path}')
plt.close()

# ── Figure 2: Lag distribution ────────────────────────────────────────────────

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.5, 3.4))

lags = results_df['lag_seconds'].dropna()
n    = len(lags)

# Panel A: Histogram
ax1.axvline(0, color=RED, linewidth=1.2, linestyle='--',
            alpha=0.8, label='VPM boundary (lag=0)')
ax1.axvline(lags.mean(), color=BLUE, linewidth=1.5,
            linestyle='-', label=f'Mean lag: {lags.mean():+.0f}s')
ax1.hist(lags, bins=min(15, n // 2 + 1), color=LBLUE,
         edgecolor=BLUE, linewidth=0.8, alpha=0.85)

# Shade positive region
ylim = ax1.get_ylim()
ax1.axvspan(0, lags.max() + 20, alpha=0.06, color=GREEN)
ax1.text(lags.max() * 0.5 + 10, (ylim[1] - ylim[0]) * 0.8,
         f'{(lags>0).sum()}/{n}\nVPM confirmed',
         fontsize=8, color=GREEN, ha='center', style='italic')

p_one = summary['p_one_sided']
ax1.set_xlabel('Lag: variance peak → mean threshold crossing (s)',
               fontsize=9.5)
ax1.set_ylabel('Number of subjects', fontsize=9.5)
ax1.set_title(f'A   Lag distribution (p = {p_one:.3f}, one-sided)',
              fontsize=9.5, fontweight='bold', loc='left', pad=4)
ax1.legend(fontsize=8, frameon=False)

# Panel B: Individual spaghetti — variance traces
shown = 0
for i, ((t, m), (_, v)) in enumerate(zip(all_mean[:15], all_var[:15])):
    t_s = np.array(t) * EPOCH_SEC
    ax2.plot(t_s, v / (np.nanmax(v) + 1e-9),
             color=ORANGE, alpha=0.25, linewidth=0.8)
    shown += 1

ax2.plot(t_sec, gm_var, color=ORANGE, linewidth=2.2,
         label='Group mean (variance)')
ax2.axvline(0, color=RED, linewidth=1.0, linestyle='--', alpha=0.8)

# Mean lag annotation
mean_lag = lags.mean()
ax2.axvline(mean_lag, color=BLUE, linewidth=1.2, linestyle=':',
            label=f'Mean variance peak: {mean_lag:+.0f}s')
ax2.set_xlabel('Time relative to N3 onset (s)', fontsize=9.5)
ax2.set_ylabel('Normalised RSI variance', fontsize=9.5)
ax2.set_title(f'B   Individual traces (n={shown} shown)',
              fontsize=9.5, fontweight='bold', loc='left', pad=4)
ax2.legend(fontsize=8, frameon=False)

fig.suptitle(
    'Figure 2. Individual VPM Lags and Subject Variance Trajectories',
    fontsize=9.5, fontweight='bold', x=0.02, ha='left', y=1.02)
plt.tight_layout()

fig2_path = os.path.join(FIGURES_DIR, 'fig2_vpm_lags.png')
plt.savefig(fig2_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f'✓ Figure 2 saved: {fig2_path}')
plt.close()

# ── Figure 3: Power-law exponent (H2) ────────────────────────────────────────

betas = results_df['beta_exp'].dropna()
n_b   = len(betas)

fig, ax = plt.subplots(figsize=(5.5, 3.6))

ax.axvline(0.546, color=RED, linewidth=1.5, linestyle='--',
           label='HRIT prediction: β = 0.546 (Verification 02)')
ax.axvline(0.5, color=GREY, linewidth=0.8, linestyle=':',
           label='Exact exponent: β = 0.500 (Proof 08)')
ax.axvline(betas.mean(), color=ORANGE, linewidth=1.5,
           label=f'Observed mean: β = {betas.mean():.3f}')

ax.hist(betas, bins=min(12, n_b // 2 + 1),
        color=LBLUE, edgecolor=BLUE, linewidth=0.8, alpha=0.85)

ci_low  = betas.mean() - 1.96 * betas.std() / np.sqrt(n_b)
ci_high = betas.mean() + 1.96 * betas.std() / np.sqrt(n_b)
ax.axvspan(ci_low, ci_high, alpha=0.15, color=ORANGE,
           label=f'95% CI: [{ci_low:.3f}, {ci_high:.3f}]')

t_stat, p_stat = stats.ttest_1samp(betas, 0.546)

ax.set_xlabel('Power-law exponent β (variance ~ distance-to-threshold^β)',
              fontsize=9.5)
ax.set_ylabel('Number of subjects', fontsize=9.5)
ax.set_title(
    f'Figure 3. Power-Law Exponent Distribution\n'
    f't vs 0.546: t = {t_stat:.2f}, p = {p_stat:.3f}, N = {n_b}',
    fontsize=9.5, fontweight='bold', pad=6, loc='left')
ax.legend(fontsize=8, frameon=False, loc='upper right')

fig.text(0.01, -0.04,
         'HRIT Verification 02 predicts β ≈ 0.546 (consistent with 0.5 at p > 0.05).\n'
         'If 95% CI includes 0.546, H2 is not falsified. Deviations quantify\n'
         'precision of the VPM early-warning signal in sleep-stage EEG.',
         fontsize=7.5, color=GREY)

plt.tight_layout()
fig3_path = os.path.join(FIGURES_DIR, 'fig3_powerlaw_beta.png')
plt.savefig(fig3_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f'✓ Figure 3 saved: {fig3_path}')
plt.close()

print(f'\nAll figures saved to: {FIGURES_DIR}/')
print('Next step: run 05_write_results.py')
