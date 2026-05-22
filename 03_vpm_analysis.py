"""
HRIT VPM Study 1 — Step 3 (PRE-REGISTERED): Core VPM Analysis
==============================================================
Pre-registration: https://doi.org/10.17605/OSF.IO/X5ZPU
Date registered:  13 Apr 2026

COMPLEXITY MEASURE: Aperiodic 1/f slope (−β from log-log PSD, Welch method,
    4-s windows, 50% overlap, 2–40 Hz, alpha-excluded 9–12 Hz)
    Computed in 02_compute_complexity.py → 'aperiodic' column
    Direction: INCREASES Wake→N3 (N3 > Wake — steeper 1/f during NREM)

N3 THRESHOLD (pre-registered Step 4):
    25th percentile of that subject's own N3 aperiodic slopes
    = lower bound of N3 = first level the rising mean enters N3 territory

CROSSING DIRECTION (pre-registered Step 5):
    t_cross_mean = first epoch where rolling mean EXCEEDS (>=) N3 threshold
    Mean rises from Wake level and crosses upward into N3 territory

LAG SIGN CONVENTION (pre-registered Step 6):
    lag = t_cross_mean − t_peak_var
    Positive lag → variance peaked BEFORE mean crossed → VPM confirmed ✓

VPM PREDICTIONS:
    H1 (primary):    lag > 0 (Wilcoxon one-sided, α = .05)
    H2 (exploratory): variance power-law β ≈ 0.546

HISTORICAL BUG (fixed in this version):
    Original code used 75th percentile + <= threshold — both calibrated
    for a DECREASING measure. For aperiodic slope (rising), this set the
    threshold above nearly all window values, forcing t_cross_mean = −20
    for 16/19 subjects and making the test meaningless.
    Fixed: 25th pctile + >= direction, exactly matching the pre-registration.
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

# ── Pre-registered parameters (do not change) ─────────────────────────────────
EPOCH_SEC      = 30
WINDOW_EPOCHS  = 6      # 6 × 30s = 3-minute rolling window (centred)
PRE_EPOCHS     = 20     # 20 × 30s = 10-minute pre-transition look-back
POST_EPOCHS    = 10     # 10 × 30s = 5-minute post-onset
MIN_PRE_EPOCHS = 5      # exclusion criterion (pre-registered)
MIN_N3_EPOCHS  = 3      # exclusion criterion (pre-registered)
COMPLEXITY_COL = 'aperiodic'   # pre-registered proxy


def find_first_n3_transition(stages, min_consecutive=3):
    """First epoch of first sustained Wake→N3 transition (≥3 consecutive N3)."""
    n = len(stages)
    for i in range(n - min_consecutive):
        if stages[i] < 3 and all(stages[i+j] == 3
                                   for j in range(1, min_consecutive + 1)):
            return i + 1
    return None


def rolling_stats(values, window):
    s = pd.Series(values)
    return (s.rolling(window, center=True, min_periods=1).mean().values,
            s.rolling(window, center=True, min_periods=1).var().values)


def analyse_subject(df_sub):
    df_sub = df_sub.dropna(subset=[COMPLEXITY_COL]).reset_index(drop=True)

    n3_ep   = df_sub[df_sub.stage == 3][COMPLEXITY_COL]
    wake_ep = df_sub[df_sub.stage == 0][COMPLEXITY_COL]

    # Pre-registered exclusion criteria
    if len(n3_ep) < MIN_N3_EPOCHS:
        return None, 'fewer than 3 N3 epochs'

    transition_idx = find_first_n3_transition(df_sub['stage'].values)
    if transition_idx is None:
        return None, 'no sustained Wake→N3 transition'

    start = max(0, transition_idx - PRE_EPOCHS)
    end   = min(len(df_sub), transition_idx + POST_EPOCHS)
    win   = df_sub.iloc[start:end].copy()
    win['epoch_rel'] = np.arange(start, end) - transition_idx

    rsi_proxy = win[COMPLEXITY_COL].values
    rel_pos   = win['epoch_rel'].values

    pre_mask = rel_pos < 0
    if pre_mask.sum() < MIN_PRE_EPOCHS:
        return None, 'fewer than 5 pre-transition epochs'

    pre_rel  = rel_pos[pre_mask]
    roll_mean, roll_var = rolling_stats(rsi_proxy, WINDOW_EPOCHS)
    pre_mean = roll_mean[pre_mask]
    pre_var  = roll_var[pre_mask]

    # ── Pre-registered Step 4: threshold = 25th pctile of N3 aperiodic ───────
    threshold = float(np.percentile(n3_ep, 25))

    # ── Pre-registered Step 5: t_cross_mean ──────────────────────────────────
    # Aperiodic slope RISES Wake→N3; mean first EXCEEDS (>=) 25th pctile of N3
    cross_pre = np.where(pre_mean >= threshold)[0]
    if len(cross_pre) == 0:
        cross_all = np.where(roll_mean >= threshold)[0]
        if len(cross_all) == 0:
            return None, 'mean never reached N3 threshold'
        t_cross_mean = int(rel_pos[cross_all[0]])
    else:
        t_cross_mean = int(pre_rel[cross_pre[0]])

    # ── Pre-registered Step 5: t_peak_var ────────────────────────────────────
    t_peak_var = int(pre_rel[np.argmax(pre_var)])

    # ── Pre-registered Step 6: lag ────────────────────────────────────────────
    lag_epochs  = t_cross_mean - t_peak_var
    lag_seconds = float(lag_epochs * EPOCH_SEC)

    # ── H2 (exploratory): power-law β ────────────────────────────────────────
    t_to_thresh = t_cross_mean - pre_rel
    t_to_thresh = np.maximum(t_to_thresh, 0.5)
    beta_exp = np.nan
    if len(t_to_thresh) >= 5 and pre_var.max() > 0:
        try:
            slope, *_ = stats.linregress(np.log(t_to_thresh),
                                         np.log(pre_var + 1e-30))
            beta_exp = float(slope)
        except Exception:
            pass

    r_vm, p_vm = (np.nan, np.nan)
    if len(pre_mean) > 3:
        r_vm, p_vm = stats.pearsonr(pre_var, pre_mean)

    return {
        'transition_idx': transition_idx,
        't_peak_var':     t_peak_var,
        't_cross_mean':   t_cross_mean,
        'lag_epochs':     lag_epochs,
        'lag_seconds':    lag_seconds,
        'threshold_used': threshold,
        'n3_ap_mean':     float(n3_ep.mean()),
        'wake_ap_mean':   float(wake_ep.mean()) if len(wake_ep) else np.nan,
        'r_var_mean':     float(r_vm),
        'p_var_mean':     float(p_vm),
        'beta_exp':       beta_exp if not np.isnan(beta_exp) else None,
        'n_pre_epochs':   int(pre_mask.sum()),
        'n_n3_epochs':    int(len(n3_ep)),
        'rsi_time':       rel_pos.tolist(),
        'rsi_mean':       roll_mean.tolist(),
        'rsi_var':        roll_var.tolist(),
    }, None


# ── Main ─────────────────────────────────────────────────────────────────────

files = sorted(glob.glob(os.path.join(COMPLEXITY_DIR, 'SC*_complexity.csv')))
print(f'VPM Analysis — Study 1')
print(f'Pre-registration: doi.org/10.17605/OSF.IO/X5ZPU')
print(f'Complexity: {COMPLEXITY_COL} | Threshold: 25th pctile N3 | Direction: >=')
print(f'Subjects found: {len(files)}\n')

results, excluded = [], []

for fpath in files:
    sid    = os.path.basename(fpath).replace('_complexity.csv', '')
    df_sub = pd.read_csv(fpath)

    wake_ap = df_sub[df_sub.stage==0][COMPLEXITY_COL].mean()
    n3_ap   = df_sub[df_sub.stage==3][COMPLEXITY_COL].mean()

    result, reason = analyse_subject(df_sub)
    if result is None:
        print(f'{sid}: EXCLUDED — {reason}')
        excluded.append({'subject': sid, 'reason': reason})
        continue

    result['subject'] = sid
    results.append(result)

    lag_s  = result['lag_seconds']
    vpm    = 'VPM ✓' if lag_s > 0 else ('tie' if lag_s == 0 else '✗')
    beta   = result['beta_exp']
    beta_s = f'{beta:.3f}' if beta is not None else 'N/A'
    print(f'{sid}: Wake={wake_ap:.3f} N3={n3_ap:.3f} '
          f'thresh={result["threshold_used"]:.3f} '
          f't_peak={result["t_peak_var"]:+d} t_cross={result["t_cross_mean"]:+d} '
          f'lag={lag_s:+.0f}s β={beta_s}  {vpm}')

if not results:
    print('\nNo subjects passed exclusion criteria.')
    raise SystemExit(1)

df_out = pd.DataFrame(results)
df_out.to_csv(os.path.join(RESULTS_DIR, 'vpm_subject_results.csv'), index=False)
if excluded:
    pd.DataFrame(excluded).to_csv(
        os.path.join(RESULTS_DIR, 'excluded_subjects.csv'), index=False)

lags = df_out['lag_seconds'].dropna()
n    = len(lags)

print(f'\n{"="*60}')
print(f'STUDY 1 — PRE-REGISTERED RESULTS')
print(f'Pre-registration: doi.org/10.17605/OSF.IO/X5ZPU')
print(f'{"="*60}')
print(f'N included: {n}  |  N excluded: {len(excluded)}')
print(f'\nLag (variance peak → mean threshold crossing):')
print(f'  Mean:   {lags.mean():+.1f}s ({lags.mean()/60:.1f} min)')
print(f'  Median: {lags.median():+.1f}s ({lags.median()/60:.1f} min)')
print(f'  SD:     {lags.std():.1f}s')
print(f'  Range:  [{lags.min():+.1f}, {lags.max():+.1f}]s')
print(f'  VPM confirmed: {(lags>0).sum()}/{n} ({100*(lags>0).mean():.0f}%)')

stat, p_one = stats.wilcoxon(lags, alternative='greater')
_, p_two    = stats.wilcoxon(lags, alternative='two-sided')
z = norm_dist.ppf(1-p_one)
r = z / np.sqrt(n)

print(f'\nH1 (primary) — Wilcoxon signed-rank (one-sided, lag > 0):')
print(f'  W = {stat:.0f}')
print(f'  p = {p_one:.4f}')
print(f'  r = {r:.3f}')
print(f'  → {"CONFIRMED ✓" if p_one < 0.05 else "NOT CONFIRMED ✗"}')

betas = pd.Series([r['beta_exp'] for r in results
                   if r['beta_exp'] is not None]).dropna()
beta_mean = beta_sd = ci_lo = ci_hi = beta_p = np.nan
if len(betas) >= 3:
    ci_lo = betas.mean() - 1.96*betas.std()/np.sqrt(len(betas))
    ci_hi = betas.mean() + 1.96*betas.std()/np.sqrt(len(betas))
    t_b, beta_p = stats.ttest_1samp(betas, 0.546)
    beta_mean, beta_sd = betas.mean(), betas.std()
    print(f'\nH2 (exploratory) — Power-law β:')
    print(f'  Mean β = {beta_mean:.3f} ± {beta_sd:.3f}')
    print(f'  95% CI = [{ci_lo:.3f}, {ci_hi:.3f}]')
    print(f'  t vs predicted 0.546: t={t_b:.2f}, p={beta_p:.3f}')
    print(f'  β=0.546 {"inside" if ci_lo<=0.546<=ci_hi else "outside"} 95% CI '
          f'→ {"CONSISTENT ✓" if beta_p > 0.05 else "DIFFERS ✗"}')

summary = {
    'n_subjects': n, 'n_excluded': len(excluded),
    'complexity_measure': COMPLEXITY_COL,
    'threshold_definition': '25th_pctile_N3_aperiodic',
    'crossing_direction': '>=',
    'lag_mean_sec': lags.mean(), 'lag_mean_min': lags.mean()/60,
    'lag_median_sec': lags.median(), 'lag_sd_sec': lags.std(),
    'pct_positive_lag': 100*(lags>0).mean(),
    'n_positive_lag': int((lags>0).sum()),
    'wilcoxon_W': stat, 'p_one_sided': p_one,
    'p_two_sided': p_two, 'effect_r': r,
    'beta_mean': beta_mean, 'beta_sd': beta_sd,
    'beta_ci_lo': ci_lo, 'beta_ci_hi': ci_hi,
    'beta_vs_predicted_p': beta_p,
    'preregistration_doi': '10.17605/OSF.IO/X5ZPU',
}
pd.DataFrame([summary]).to_csv(
    os.path.join(RESULTS_DIR, 'vpm_summary_stats.csv'), index=False)
print(f'\nSaved → {RESULTS_DIR}/')
print('Next: python 04_figures.py')
