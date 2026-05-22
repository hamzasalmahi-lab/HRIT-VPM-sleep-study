"""
HRIT VPM Study 1 — Step 4: Publication Figures
===============================================
Reads results/vpm_subject_results.csv and generates three figures.

Figure 1  Group-level VPM — mean ± SE of normalised aperiodic slope
          and its rolling variance, aligned to N3 onset.
Figure 2  Individual lags — histogram + spaghetti variance traces.
Figure 3  Power-law β — observed distribution vs HRIT prediction (0.546).

Output: figures/fig1_vpm_group.png
        figures/fig2_vpm_lags.png
        figures/fig3_powerlaw_beta.png
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

RESULTS_DIR = 'results'
FIGURES_DIR = 'figures'
os.makedirs(FIGURES_DIR, exist_ok=True)

EPOCH_SEC   = 30
PRE_EPOCHS  = 20
POST_EPOCHS = 10

NAVY   = '#1a3a5c'
BLUE   = '#2471a3'
ORANGE = '#c0672a'
RED    = '#c0392b'
GREEN  = '#1e8449'
GREY   = '#566573'
LGREY  = '#d5d8dc'
LBLUE  = '#aed6f1'

plt.rcParams.update({
    'font.family':       'DejaVu Serif',
    'font.size':         10,
    'axes.spines.top':   False,
    'axes.spines.right': False,
    'axes.linewidth':    0.8,
    'lines.linewidth':   1.6,
    'figure.dpi':        300,
    'savefig.dpi':       300,
    'savefig.bbox':      'tight',
    'savefig.facecolor': 'white',
})


def load_time_series(results_df):
    """Reconstruct rolling mean/variance time series from stored JSON."""
    t_grid   = np.arange(-PRE_EPOCHS, POST_EPOCHS + 1)
    n_subs   = len(results_df)
    mean_mat = np.full((n_subs, len(t_grid)), np.nan)
    var_mat  = np.full((n_subs, len(t_grid)), np.nan)
    raw_ts   = []

    for i, (_, row) in enumerate(results_df.iterrows()):
        t = np.array(json.loads(row['rsi_time']))
        m = np.array(json.loads(row['rsi_mean']))
        v = np.array(json.loads(row['rsi_var']))
        raw_ts.append((t, m, v))

        for j, tg in enumerate(t_grid):
            idx = np.where(t == tg)[0]
            if len(idx):
                mean_mat[i, j] = m[idx[0]]
                var_mat[i, j]  = v[idx[0]]

    # Normalise each subject 0–1 for group average
    for i in range(n_subs):
        for mat in [mean_mat, var_mat]:
            row = mat[i]
            rng = np.nanmax(row) - np.nanmin(row)
            if rng > 0:
                mat[i] = (row - np.nanmin(row)) / rng

    t_sec = t_grid * EPOCH_SEC
    gm    = np.nanmean(mean_mat, axis=0)
    gv    = np.nanmean(var_mat, axis=0)
    se_m  = np.nanstd(mean_mat, axis=0) / np.sqrt(np.sum(~np.isnan(mean_mat), axis=0))
    se_v  = np.nanstd(var_mat, axis=0)  / np.sqrt(np.sum(~np.isnan(var_mat),  axis=0))

    return t_sec, gm, gv, se_m, se_v, raw_ts


# ── Load ──────────────────────────────────────────────────────────────────────

results_df = pd.read_csv(os.path.join(RESULTS_DIR, 'vpm_subject_results.csv'))
summary    = pd.read_csv(os.path.join(RESULTS_DIR, 'vpm_summary_stats.csv')).iloc[0]
n          = len(results_df)
lags       = results_df['lag_seconds'].dropna()
print(f"Loaded {n} subjects")

t_sec, gm, gv, se_m, se_v, raw_ts = load_time_series(results_df)

# ── Figure 1: Group-level VPM ─────────────────────────────────────────────────

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.5, 3.8))

# Panel A — rolling mean
ax1.fill_between(t_sec, gm - se_m, gm + se_m, color=LBLUE, alpha=0.55)
ax1.plot(t_sec, gm, color=BLUE, lw=2.0, label='Aperiodic slope (mean)')
ax1.axvline(0, color=RED, lw=1.0, ls='--', alpha=0.8)
ax1.axvspan(0, t_sec[-1], alpha=0.04, color=GREY)
ax1.text(5, 0.04, 'N3 onset', color=RED, fontsize=8, rotation=90, va='bottom')
ax1.set_xlabel('Time relative to N3 onset (s)')
ax1.set_ylabel('Normalised aperiodic slope')
ax1.set_title('A.  Mean: rises into N3', loc='left', fontweight='bold', pad=5)
ax1.legend(fontsize=8, frameon=False)

# Panel B — rolling variance with annotated peak
ax2.fill_between(t_sec, gv - se_v, gv + se_v, color='#f9c48c', alpha=0.55)
ax2.plot(t_sec, gv, color=ORANGE, lw=2.0, label='Aperiodic slope (variance)')
ax2.axvline(0, color=RED, lw=1.0, ls='--', alpha=0.8)

pre_gv        = gv.copy()
pre_gv[t_sec >= 0] = np.nan
peak_idx      = np.nanargmax(pre_gv)
t_peak        = t_sec[peak_idx]
ax2.axvline(t_peak, color=ORANGE, lw=1.2, ls='-.', alpha=0.85)
ax2.annotate(f'Variance peak\n({abs(t_peak):.0f}s before N3)',
             xy=(t_peak, gv[peak_idx]),
             xytext=(t_peak - 4 * EPOCH_SEC, 0.78),
             fontsize=8, color=ORANGE,
             arrowprops=dict(arrowstyle='->', color=ORANGE, lw=0.9))

ax2.text(5, 0.04, 'N3 onset', color=RED, fontsize=8, rotation=90, va='bottom')
ax2.set_xlabel('Time relative to N3 onset (s)')
ax2.set_ylabel('Normalised variance')
ax2.set_title('B.  Variance: peaks before mean crosses threshold',
              loc='left', fontweight='bold', pad=5)
ax2.legend(fontsize=8, frameon=False)

fig.suptitle(
    'Figure 1.  Variance-Precedes-Mean at the Wake→N3 Sleep Transition\n'
    f'Group mean ± SE, N = {n} subjects. Aperiodic slope proxy (Study 1).',
    fontsize=9.5, fontweight='bold', x=0.02, ha='left')
plt.tight_layout(rect=[0, 0, 1, 0.93])
path1 = os.path.join(FIGURES_DIR, 'fig1_vpm_group.png')
plt.savefig(path1)
print(f'✓  {path1}')
plt.close()

# ── Figure 2: Individual lags ─────────────────────────────────────────────────

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.5, 3.8))

# Panel A — lag histogram
ax1.axvline(0, color=RED, lw=1.2, ls='--', alpha=0.8, label='VPM boundary (lag = 0)')
ax1.axvline(lags.mean(), color=BLUE, lw=1.5,
            label=f'Mean = {lags.mean():+.0f}s ({lags.mean()/60:.1f} min)')
ax1.hist(lags, bins=max(6, len(lags) // 3),
         color=LBLUE, edgecolor=BLUE, lw=0.8, alpha=0.85)
ax1.axvspan(0, lags.max() + 50, alpha=0.06, color=GREEN)
ax1.text(lags.max() * 0.55, ax1.get_ylim()[1] * 0.7,
         f'{(lags>0).sum()}/{len(lags)}\nVPM ✓',
         fontsize=9, color=GREEN, ha='center', style='italic')
ax1.set_xlabel('Lag: variance peak → N3 threshold crossing (s)')
ax1.set_ylabel('N subjects')
ax1.set_title(f'A.  Lag distribution\n'
              f'W = {summary["wilcoxon_W"]:.0f},  '
              f'p < .001,  r = {summary["effect_r"]:.3f}',
              loc='left', fontweight='bold', pad=5)
ax1.legend(fontsize=8, frameon=False)

# Panel B — spaghetti variance traces
for t, m, v in raw_ts:
    t_s = t * EPOCH_SEC
    vn  = v / (np.nanmax(v) + 1e-9)
    ax2.plot(t_s, vn, color=ORANGE, alpha=0.20, lw=0.7)

ax2.plot(t_sec, gv, color=ORANGE, lw=2.2, label='Group mean (variance)')
ax2.axvline(0, color=RED, lw=1.0, ls='--', alpha=0.8, label='N3 onset')
ax2.axvline(lags.mean(), color=BLUE, lw=1.2, ls=':',
            label=f'Mean lag ({lags.mean():+.0f}s)')
ax2.set_xlabel('Time relative to N3 onset (s)')
ax2.set_ylabel('Normalised variance')
ax2.set_title(f'B.  Individual subject traces (N = {len(raw_ts)})',
              loc='left', fontweight='bold', pad=5)
ax2.legend(fontsize=8, frameon=False)

fig.suptitle('Figure 2.  Individual VPM Lags and Variance Trajectories\n'
             'Aperiodic slope proxy (Study 1).  Each grey trace = one subject.',
             fontsize=9.5, fontweight='bold', x=0.02, ha='left')
plt.tight_layout(rect=[0, 0, 1, 0.93])
path2 = os.path.join(FIGURES_DIR, 'fig2_vpm_lags.png')
plt.savefig(path2)
print(f'✓  {path2}')
plt.close()

# ── Figure 3: Power-law β (H2, exploratory) ──────────────────────────────────

betas  = results_df['beta_exp'].dropna()
n_b    = len(betas)

if n_b < 3:
    print("⚠  Fewer than 3 β estimates — skipping Figure 3")
else:
    fig, ax = plt.subplots(figsize=(5.5, 3.8))

    ci_lo = betas.mean() - 1.96 * betas.std() / np.sqrt(n_b)
    ci_hi = betas.mean() + 1.96 * betas.std() / np.sqrt(n_b)
    t_b, p_b = stats.ttest_1samp(betas, 0.546)

    ax.axvline(0.546, color=RED,  lw=1.5, ls='--',
               label='HRIT predicted: β = 0.546')
    ax.axvline(0.500, color=GREY, lw=0.8, ls=':',
               label='Exact exponent: β = 0.500')
    ax.axvline(betas.mean(), color=ORANGE, lw=1.5,
               label=f'Observed mean: β = {betas.mean():.3f}')
    ax.axvspan(ci_lo, ci_hi, alpha=0.15, color=ORANGE,
               label=f'95% CI: [{ci_lo:.3f}, {ci_hi:.3f}]')
    ax.hist(betas, bins=max(5, n_b // 2),
            color=LBLUE, edgecolor=BLUE, lw=0.8, alpha=0.85)

    ax.set_xlabel('Power-law exponent β  (variance ~ distance-to-threshold^β)')
    ax.set_ylabel('N subjects')
    ax.set_title(f'Figure 3.  H2 — Power-Law Exponent (Exploratory)\n'
                 f'Observed β vs HRIT prediction 0.546 '
                 f'(t = {t_b:.2f}, p = {p_b:.3f}, N = {n_b})',
                 loc='left', fontweight='bold', pad=5)
    ax.legend(fontsize=8, frameon=False, loc='upper left')

    path3 = os.path.join(FIGURES_DIR, 'fig3_powerlaw_beta.png')
    plt.tight_layout()
    plt.savefig(path3)
    print(f'✓  {path3}')
    plt.close()

print(f'\nAll figures saved to:  {FIGURES_DIR}/')
print('Run next:  python 05_write_results.py')
