"""
HRIT VPM Study 1 — Step 3 (FIXED): Core VPM Analysis
======================================================
Pre-registration: https://doi.org/10.17605/OSF.IO/X5ZPU

HRIT Prediction (H1):
    In the window preceding the Wake→N3 EEG complexity threshold
    crossing, rolling variance of the complexity proxy reaches its
    maximum BEFORE the rolling mean crosses the N3 threshold.
    Lag = t_cross_mean − t_peak_var > 0 (Wilcoxon, one-sided).

HRIT Prediction (H2, exploratory):
    Variance ~ (z_c − z)^β with predicted β ≈ 0.546.

COMPLEXITY MEASURE: Permutation entropy (order=6, delay=1)
    — decreases monotonically from Wake to N3 (Wake > N1 > N2 > N3)
    — correct direction for the threshold-crossing logic below

THRESHOLD CONVENTION:
    threshold = 75th percentile of N3 PE epochs (upper boundary of N3)
    t_cross_mean = first epoch where rolling mean drops BELOW threshold
    (mean falls from Wake level ~0.85 into N3 territory ~0.58)
    lag > 0 ↔ variance peaked BEFORE mean crossed ↔ VPM confirmed

BUG FIXED (vs prior version):
    Prior code used aperiodic slope, which INCREASES Wake→N3 in
    Sleep-EDF (N3 > Wake). This inverted the threshold logic and
    caused t_cross_mean to fire at the window start for every subject.
    Permutation entropy correctly decreases Wake→N3.
"""

import os
import glob
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import norm as norm_dist
import warnings
warnings.filterwarnings('ignore')

COMPLEXITY_DIR = 'complexity'
RESULTS_DIR    = 'results'
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── Analysis parameters ───────────────────────────────────────────────────────
EPOCH_SEC      = 30
WINDOW_EPOCHS  = 6     # rolling window: 6 × 30s = 3 min
PRE_EPOCHS     = 20    # look-back window: 20 × 30s = 10 min before N3
POST_EPOCHS    = 10    # look-ahead: 5 min after N3 onset
STAGE_N3       = 3
COMPLEXITY_COL = 'perm_entropy'   # ← FIXED: use PE, not aperiodic


def find_first_n3_transition(stages, min_consecutive=3):
    """First epoch of the first sustained Wake→N3 transition."""
    n = len(stages)
    for i in range(n - min_consecutive):
        if stages[i] < 3 and all(stages[i+j] == 3
                                   for j in range(1, min_consecutive + 1)):
            return i + 1
    return None


def compute_rolling_stats(values, window):
    s = pd.Series(values)
    return (s.rolling(window, center=True, min_periods=1).mean().values,
            s.rolling(window, center=True, min_periods=1).var().values)


def analyse_subject(df_sub):
    df_sub = df_sub.dropna(subset=[COMPLEXITY_COL]).reset_index(drop=True)
    if len(df_sub) < PRE_EPOCHS + POST_EPOCHS:
        return None

    transition_idx = find_first_n3_transition(df_sub['stage'].values)
    if transition_idx is None:
        return None

    # Analysis window
    start    = max(0, transition_idx - PRE_EPOCHS)
    end      = min(len(df_sub), transition_idx + POST_EPOCHS)
    win      = df_sub.iloc[start:end].copy()
    win['epoch_rel'] = np.arange(start, end) - transition_idx

    rsi_proxy  = win[COMPLEXITY_COL].values
    rel_pos    = win['epoch_rel'].values

    # Per-subject N3 threshold: PE 75th percentile of N3 epochs
    # (upper boundary of N3 = first level PE drops into from above)
    n3_ep   = df_sub[df_sub.stage == 3][COMPLEXITY_COL]
    wake_ep = df_sub[df_sub.stage == 0][COMPLEXITY_COL]
    if len(n3_ep) < 5 or len(wake_ep) < 5:
        return None

    threshold = float(np.percentile(n3_ep, 75))

    # Sanity check: Wake PE should be above threshold
    if wake_ep.mean() < threshold:
        # Wake and N3 PE are inverted for this subject — skip
        print(f"    ⚠ Inverted Wake/N3 PE direction — subject skipped")
        return None

    roll_mean, roll_var = compute_rolling_stats(rsi_proxy, WINDOW_EPOCHS)

    pre_mask = rel_pos < 0
    if pre_mask.sum() < 5:
        return None

    pre_rel  = rel_pos[pre_mask]
    pre_var  = roll_var[pre_mask]
    pre_mean = roll_mean[pre_mask]

    # ── t_peak_var ────────────────────────────────────────────────────────────
    t_peak_var = pre_rel[np.argmax(pre_var)]

    # ── t_cross_mean: first epoch mean DROPS BELOW threshold ─────────────────
    # PE falls from Wake (~0.85) into N3 territory (below threshold ~0.62)
    cross_idx = np.where(pre_mean <= threshold)[0]
    if len(cross_idx) == 0:
        # Mean never reached N3 level in pre-window — assign to post-window
        all_cross = np.where(roll_mean <= threshold)[0]
        t_cross_mean = rel_pos[all_cross[0]] if len(all_cross) else POST_EPOCHS
    else:
        t_cross_mean = pre_rel[cross_idx[0]]

    lag_epochs  = t_cross_mean - t_peak_var
    lag_seconds = lag_epochs * EPOCH_SEC

    # ── Power-law exponent β (H2) ─────────────────────────────────────────────
    t_to_thresh = t_cross_mean - pre_rel
    t_to_thresh = np.maximum(t_to_thresh, 0.5)
    beta_exp = np.nan
    if len(t_to_thresh) >= 5 and pre_var.max() > 0:
        try:
            from scipy.stats import linregress
            slope, *_ = linregress(np.log(t_to_thresh),
                                   np.log(pre_var + 1e-30))
            beta_exp = slope
        except Exception:
            pass

    # ── Variance-mean correlation (should be negative: var ↑ as mean ↓) ──────
    r_vm, p_vm = (np.nan, np.nan)
    if len(pre_mean) > 3:
        r_vm, p_vm = stats.pearsonr(pre_var, pre_mean)

    return {
        'transition_idx': transition_idx,
        'window_start':   start,
        'window_end':     end,
        't_peak_var':     int(t_peak_var),
        't_cross_mean':   int(t_cross_mean),
        'lag_epochs':     int(lag_epochs),
        'lag_seconds':    float(lag_seconds),
        'threshold_used': float(threshold),
        'n3_pe_mean':     float(n3_ep.mean()),
        'wake_pe_mean':   float(wake_ep.mean()),
        'r_var_mean':     float(r_vm),
        'p_var_mean':     float(p_vm),
        'beta_exp':       float(beta_exp) if not np.isnan(beta_exp) else None,
        'n_pre_epochs':   int(pre_mask.sum()),
        'n_n3_epochs':    int(len(n3_ep)),
        'rsi_time':       rel_pos.tolist(),
        'rsi_mean':       roll_mean.tolist(),
        'rsi_var':        roll_var.tolist(),
    }


# ── Main ─────────────────────────────────────────────────────────────────────

complexity_files = sorted(glob.glob(os.path.join(COMPLEXITY_DIR, 'SC*_complexity.csv')))
print(f'Found {len(complexity_files)} complexity files')
print(f'Complexity measure: {COMPLEXITY_COL} (decreases Wake→N3 ✓)\n')

subject_results = []

for fpath in complexity_files:
    sid    = os.path.basename(fpath).replace('_complexity.csv', '')
    df_sub = pd.read_csv(fpath)

    # Sanity check PE direction
    wake_pe = df_sub[df_sub.stage==0][COMPLEXITY_COL].mean()
    n3_pe   = df_sub[df_sub.stage==3][COMPLEXITY_COL].mean()
    direction = '✓' if wake_pe > n3_pe else '⚠ INVERTED'

    print(f'{sid}: Wake PE={wake_pe:.3f}  N3 PE={n3_pe:.3f}  {direction}', end='  ')

    result = analyse_subject(df_sub)
    if result is None:
        print('✗ no usable transition')
        continue

    result['subject'] = sid
    subject_results.append(result)

    lag_s = result['lag_seconds']
    beta  = result['beta_exp']
    vpm   = 'VPM ✓' if lag_s > 0 else ('tie' if lag_s == 0 else 'NOT VPM ✗')
    beta_str = f'{beta:.3f}' if beta is not None else 'N/A'
    print(f'lag={lag_s:+.0f}s  β={beta_str}  {vpm}')

if not subject_results:
    print('\nNo subjects had usable transitions.')
    raise SystemExit(1)

# ── Statistics ───────────────────────────────────────────────────────────────

results_df = pd.DataFrame(subject_results)
results_df.to_csv(os.path.join(RESULTS_DIR, 'vpm_subject_results.csv'), index=False)

lags = results_df['lag_seconds'].dropna()
n    = len(lags)

print(f'\n{"="*60}')
print('VPM ANALYSIS — STUDY 1 RESULTS')
print(f'{"="*60}')
print(f'N subjects with valid transitions: {n}')
print(f'Complexity measure: {COMPLEXITY_COL}')
print(f'\nLag statistics:')
print(f'  Mean:     {lags.mean():+.1f}s ({lags.mean()/60:.1f} min)')
print(f'  Median:   {lags.median():+.1f}s')
print(f'  SD:       {lags.std():.1f}s')
print(f'  Range:    [{lags.min():+.1f}, {lags.max():+.1f}]s')
print(f'  VPM confirmed: {(lags>0).sum()}/{n} ({100*(lags>0).mean():.0f}%)')

# H1
stat, p_one = stats.wilcoxon(lags, alternative='greater')
_, p_two    = stats.wilcoxon(lags, alternative='two-sided')
z           = norm_dist.ppf(1 - p_one)
r_eff       = z / np.sqrt(n)

print(f'\nH1 — Wilcoxon signed-rank (one-sided, lag > 0):')
print(f'  W = {stat:.0f}')
print(f'  p = {p_one:.4f}')
print(f'  r = {r_eff:.3f}')
print(f'  → {"CONFIRMED ✓" if p_one < 0.05 else "NOT CONFIRMED ✗"}')

# H2
betas = results_df['beta_exp'].dropna()
if len(betas) >= 3:
    ci_lo = betas.mean() - 1.96*betas.std()/np.sqrt(len(betas))
    ci_hi = betas.mean() + 1.96*betas.std()/np.sqrt(len(betas))
    t_b, p_b = stats.ttest_1samp(betas, 0.546)
    print(f'\nH2 — Power-law exponent β:')
    print(f'  Mean β = {betas.mean():.3f} ± {betas.std():.3f}')
    print(f'  95% CI = [{ci_lo:.3f}, {ci_hi:.3f}]')
    print(f'  t vs 0.546: t={t_b:.2f}, p={p_b:.3f}')
    pred_in_ci = ci_lo <= 0.546 <= ci_hi
    print(f'  Predicted β=0.546 {"inside" if pred_in_ci else "outside"} 95% CI')
    print(f'  → {"CONSISTENT ✓" if p_b > 0.05 else "DIFFERS ✗"} with HRIT prediction')
else:
    ci_lo = ci_hi = p_b = np.nan
    betas_mean = betas_sd = np.nan

summary = {
    'n_subjects':          n,
    'complexity_measure':  COMPLEXITY_COL,
    'lag_mean_sec':        lags.mean(),
    'lag_median_sec':      lags.median(),
    'lag_sd_sec':          lags.std(),
    'lag_mean_min':        lags.mean() / 60,
    'pct_positive_lag':    100*(lags>0).mean(),
    'n_positive_lag':      int((lags>0).sum()),
    'wilcoxon_W':          stat,
    'p_one_sided':         p_one,
    'p_two_sided':         p_two,
    'effect_r':            r_eff,
    'beta_mean':           betas.mean() if len(betas) >= 3 else np.nan,
    'beta_sd':             betas.std()  if len(betas) >= 3 else np.nan,
    'beta_ci_lo':          ci_lo if len(betas) >= 3 else np.nan,
    'beta_ci_hi':          ci_hi if len(betas) >= 3 else np.nan,
    'beta_vs_predicted_p': p_b   if len(betas) >= 3 else np.nan,
}
pd.DataFrame([summary]).to_csv(
    os.path.join(RESULTS_DIR, 'vpm_summary_stats.csv'), index=False)

print(f'\nResults saved to: {RESULTS_DIR}/')
print('Next step: run 04_figures.py')
