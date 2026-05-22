"""
HRIT VPM Study 1b — Step 1: PEWS Analysis (Permutation Entropy Proxy)
=======================================================================
Pre-registration: https://doi.org/10.17605/OSF.IO/EYMSF

Convergent replication of Study 1 using permutation entropy (PE) as an
independent I_Sim proxy applied to the same Sleep-EDF Expanded dataset
with a completely independent analysis pipeline.

COMPLEXITY MEASURE: Permutation entropy (order=6, delay=1)
  — DECREASES Wake→N3 (Wake > N1 > N2 > N3)
  — Correct direction for standard falling-threshold VPM logic

N3 THRESHOLD (pre-registered):
  75th percentile of subject's own N3 PE epochs
  = upper boundary of N3 = first level PE falls into from above

CROSSING DIRECTION:
  t_cross_mean = first epoch where rolling mean DROPS BELOW (<=) threshold
  Mean falls from Wake level (~0.85) into N3 territory (~0.65)

LAG:
  lag = t_cross_mean − t_peak_var
  Positive → variance peaked BEFORE mean crossed → PEWS confirmed ✓

PRE-REGISTERED HYPOTHESES:
  H1: lag > 0 (Wilcoxon one-sided, α = .05)
  H2: r ≥ .633  (effect size replicates Study 1 threshold)
  H3: β ≈ 0.546 (power-law variance exponent, exploratory)

INPUTS:  ../complexity/<SID>_complexity.csv  (perm_entropy column)
         Run ../02_compute_complexity.py first if not already done.

OUTPUTS: results/pews_subject_results.csv
         results/pews_summary_stats.csv
         results/excluded_subjects.csv
"""

import os
import sys
import glob
import json
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import norm as norm_dist
import warnings
warnings.filterwarnings('ignore')

# ── Paths ─────────────────────────────────────────────────────────────────────
COMPLEXITY_DIR = os.path.join('..', 'complexity')   # shared with Study 1
RESULTS_DIR    = 'results'
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── Pre-registered parameters ─────────────────────────────────────────────────
EPOCH_SEC       = 30
WINDOW_EPOCHS   = 6       # 6 × 30s = 3-minute centred rolling window
PRE_EPOCHS      = 20      # 10-minute look-back before N3 onset
POST_EPOCHS     = 10      # 5-minute look-ahead after N3 onset
MIN_PRE_EPOCHS  = 5       # exclusion criterion
MIN_N3_EPOCHS   = 3       # exclusion criterion
COMPLEXITY_COL  = 'perm_entropy'   # PE decreases Wake→N3


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

    # Exclusion criteria
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

    # ── Threshold: 75th pctile of N3 PE (upper boundary of N3) ───────────────
    # PE falls from Wake (~0.85) into N3 territory; 75th pctile is the highest
    # N3 PE value — where the falling mean first enters N3 from above.
    threshold = float(np.percentile(n3_ep, 75))

    # Sanity check: Wake PE should be above threshold
    if len(wake_ep) and wake_ep.mean() <= threshold:
        return None, f'Wake PE ({wake_ep.mean():.3f}) not above N3 75th pctile ({threshold:.3f})'

    # ── t_cross_mean: first epoch mean DROPS BELOW threshold ─────────────────
    cross_pre = np.where(pre_mean <= threshold)[0]
    if len(cross_pre) == 0:
        cross_all = np.where(roll_mean <= threshold)[0]
        if len(cross_all) == 0:
            return None, 'mean never dropped to N3 PE level'
        t_cross_mean = int(rel_pos[cross_all[0]])
    else:
        t_cross_mean = int(pre_rel[cross_pre[0]])

    # ── t_peak_var: epoch of max variance in pre-window ───────────────────────
    t_peak_var = int(pre_rel[np.argmax(pre_var)])

    # ── Lag ───────────────────────────────────────────────────────────────────
    lag_epochs  = t_cross_mean - t_peak_var
    lag_seconds = float(lag_epochs * EPOCH_SEC)

    # ── H3 (exploratory): power-law β ────────────────────────────────────────
    t_to_thresh = t_cross_mean - pre_rel
    t_to_thresh = np.maximum(t_to_thresh, 0.5)
    beta_exp = None
    if len(t_to_thresh) >= 5 and pre_var.max() > 0:
        try:
            slope, *_ = stats.linregress(np.log(t_to_thresh),
                                         np.log(pre_var + 1e-30))
            beta_exp = float(slope)
        except Exception:
            pass

    r_vm, p_vm = np.nan, np.nan
    if len(pre_mean) > 3:
        r_vm, p_vm = stats.pearsonr(pre_var, pre_mean)

    return {
        'transition_idx': transition_idx,
        't_peak_var':     t_peak_var,
        't_cross_mean':   t_cross_mean,
        'lag_epochs':     lag_epochs,
        'lag_seconds':    lag_seconds,
        'threshold_used': threshold,
        'n3_pe_mean':     float(n3_ep.mean()),
        'wake_pe_mean':   float(wake_ep.mean()) if len(wake_ep) else np.nan,
        'r_var_mean':     float(r_vm),
        'p_var_mean':     float(p_vm),
        'beta_exp':       beta_exp,
        'n_pre_epochs':   int(pre_mask.sum()),
        'n_n3_epochs':    int(len(n3_ep)),
        'rsi_time':       rel_pos.tolist(),
        'rsi_mean':       roll_mean.tolist(),
        'rsi_var':        roll_var.tolist(),
    }, None


# ── Main ──────────────────────────────────────────────────────────────────────

files = sorted(glob.glob(os.path.join(COMPLEXITY_DIR, 'SC*_complexity.csv')))

if not files:
    sys.exit(
        f"No complexity files found in {COMPLEXITY_DIR}/\n"
        "Run ../02_compute_complexity.py first to generate them."
    )

print('PEWS Analysis — Study 1b')
print('Pre-registration: doi.org/10.17605/OSF.IO/EYMSF')
print(f'Proxy: {COMPLEXITY_COL}  |  Threshold: 75th pctile N3  |  Direction: <=')
print(f'Subjects found: {len(files)}\n')

results, excluded = [], []

for fpath in files:
    sid    = os.path.basename(fpath).replace('_complexity.csv', '')
    df_sub = pd.read_csv(fpath)

    wake_pe = df_sub[df_sub.stage == 0][COMPLEXITY_COL].mean()
    n3_pe   = df_sub[df_sub.stage == 3][COMPLEXITY_COL].mean()

    result, reason = analyse_subject(df_sub)
    if result is None:
        print(f'{sid}: EXCLUDED — {reason}')
        excluded.append({'subject': sid, 'reason': reason})
        continue

    result['subject'] = sid
    results.append(result)

    lag_s  = result['lag_seconds']
    vpm    = 'PEWS ✓' if lag_s > 0 else ('tie' if lag_s == 0 else '✗')
    beta   = result['beta_exp']
    beta_s = f'{beta:.3f}' if beta is not None else 'N/A'
    print(f'{sid}: Wake={wake_pe:.4f} N3={n3_pe:.4f} '
          f'thresh={result["threshold_used"]:.4f} '
          f't_peak={result["t_peak_var"]:+d} t_cross={result["t_cross_mean"]:+d} '
          f'lag={lag_s:+.0f}s ({lag_s/60:+.1f}min) β={beta_s}  {vpm}')

if not results:
    sys.exit('No subjects passed exclusion criteria.')

df_out = pd.DataFrame(results)
df_out.to_csv(os.path.join(RESULTS_DIR, 'pews_subject_results.csv'), index=False)
if excluded:
    pd.DataFrame(excluded).to_csv(
        os.path.join(RESULTS_DIR, 'excluded_subjects.csv'), index=False)

lags  = df_out['lag_seconds'].dropna()
betas = pd.Series([r['beta_exp'] for r in results
                   if r['beta_exp'] is not None]).dropna()
n     = len(lags)

print(f'\n{"="*60}')
print('STUDY 1b — PRE-REGISTERED RESULTS')
print('Pre-registration: doi.org/10.17605/OSF.IO/EYMSF')
print(f'{"="*60}')
print(f'N included: {n}  |  N excluded: {len(excluded)}')
print(f'\nLag (variance peak → N3 threshold crossing):')
print(f'  Mean:   {lags.mean():+.1f}s ({lags.mean()/60:.1f} min)')
print(f'  Median: {lags.median():+.1f}s ({lags.median()/60:.1f} min)')
print(f'  SD:     {lags.std():.1f}s')
print(f'  Range:  [{lags.min():+.1f}, {lags.max():+.1f}]s')
print(f'  PEWS confirmed: {(lags>0).sum()}/{n} ({100*(lags>0).mean():.0f}%)')

# H1: Wilcoxon
stat, p_one = stats.wilcoxon(lags, alternative='greater')
_, p_two    = stats.wilcoxon(lags, alternative='two-sided')
z = norm_dist.ppf(1 - p_one)
r = z / np.sqrt(n)
print(f'\nH1 — Wilcoxon (one-sided, lag > 0):')
print(f'  W = {stat:.0f}')
print(f'  p = {p_one:.4f}')
print(f'  r = {r:.3f}')
print(f'  → {"CONFIRMED ✓" if p_one < 0.05 else "NOT CONFIRMED ✗"}')

# H2: r >= .633
h2 = r >= 0.633
print(f'\nH2 — Effect size r ≥ .633 (Study 1 replication threshold):')
print(f'  r = {r:.3f}')
print(f'  → {"CONFIRMED ✓" if h2 else "NOT CONFIRMED ✗"} ({r:.3f} {"≥" if h2 else "<"} .633)')

# H3: beta ~ 0.546
beta_mean = beta_sd = ci_lo = ci_hi = beta_p = np.nan
if len(betas) >= 3:
    ci_lo = betas.mean() - 1.96 * betas.std() / np.sqrt(len(betas))
    ci_hi = betas.mean() + 1.96 * betas.std() / np.sqrt(len(betas))
    t_b, beta_p = stats.ttest_1samp(betas, 0.546)
    beta_mean, beta_sd = betas.mean(), betas.std()
    print(f'\nH3 — Power-law β ≈ 0.546 (exploratory):')
    print(f'  Mean β = {beta_mean:.3f} ± {beta_sd:.3f}')
    print(f'  95% CI = [{ci_lo:.3f}, {ci_hi:.3f}]')
    print(f'  t vs 0.546: t={t_b:.2f}, p={beta_p:.3f}')
    pred_in = ci_lo <= 0.546 <= ci_hi
    print(f'  β=0.546 {"inside" if pred_in else "outside"} 95% CI '
          f'→ {"CONSISTENT ✓" if beta_p > 0.05 else "NOT CONFIRMED ✗"}')

summary = {
    'n_subjects': n, 'n_excluded': len(excluded),
    'complexity_measure': COMPLEXITY_COL,
    'threshold_definition': '75th_pctile_N3_PE',
    'crossing_direction': '<=',
    'lag_mean_sec': lags.mean(), 'lag_mean_min': lags.mean() / 60,
    'lag_median_sec': lags.median(), 'lag_sd_sec': lags.std(),
    'pct_positive_lag': 100 * (lags > 0).mean(),
    'n_positive_lag': int((lags > 0).sum()),
    'wilcoxon_W': stat, 'p_one_sided': p_one,
    'p_two_sided': p_two, 'effect_r': r,
    'h1_confirmed': p_one < 0.05,
    'h2_confirmed': h2,
    'beta_mean': beta_mean, 'beta_sd': beta_sd,
    'beta_ci_lo': ci_lo, 'beta_ci_hi': ci_hi,
    'beta_vs_predicted_p': beta_p,
    'h3_confirmed': beta_p > 0.05 if not np.isnan(beta_p) else False,
    'preregistration_doi': '10.17605/OSF.IO/EYMSF',
}
pd.DataFrame([summary]).to_csv(
    os.path.join(RESULTS_DIR, 'pews_summary_stats.csv'), index=False)

print(f'\nSaved → {RESULTS_DIR}/')
print('Run next:  python 02_figures.py')
