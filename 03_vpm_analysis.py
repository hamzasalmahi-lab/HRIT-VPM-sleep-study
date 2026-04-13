"""
HRIT VPM Study — Step 3: Core VPM Analysis
============================================
The pre-registered analysis testing the VPM prediction.

HRIT Prediction (H1):
    In the 5-minute window preceding EEG complexity threshold
    crossing at the Wake→NREM2→NREM3 transition, rolling variance
    of the complexity proxy reaches its maximum BEFORE the rolling
    mean crosses the clinical threshold.

    Operationalised: temporal lag between variance peak and mean
    threshold crossing is significantly > 0.

HRIT Prediction (H2, exploratory):
    Variance ~ (t_threshold - t)^β with β ≈ 0.546
    Mean ~ (t_threshold - t)^γ with γ tested against 0.5

Analysis steps:
    1. For each subject, identify the first sustained W→N3 transition
    2. Define a 10-minute analysis window centred on transition
    3. Compute rolling mean and rolling variance of RSI proxy
    4. Identify: t_peak_var (time of variance maximum)
    5. Identify: t_cross_mean (time mean crosses N3 threshold)
    6. Compute lag = t_cross_mean - t_peak_var
    7. Test H1: lag > 0 (Wilcoxon signed-rank test, one-sided)
    8. Fit power law for H2
"""

import os
import glob
import numpy as np
import pandas as pd
from scipy import stats, optimize
import warnings
warnings.filterwarnings('ignore')

COMPLEXITY_DIR = 'complexity'
RESULTS_DIR    = 'results'
os.makedirs(RESULTS_DIR, exist_ok=True)

# Analysis parameters
EPOCH_SEC       = 30      # seconds per epoch
WINDOW_EPOCHS   = 6       # rolling window: 6 × 30s = 3 min
PRE_EPOCHS      = 20      # epochs before transition to include (10 min)
POST_EPOCHS     = 10      # epochs after transition to include (5 min)

# Threshold: the N3 threshold for mean RSI proxy
# We operationalise this as the mean aperiodic slope at N2/N3 boundary
# Defined per-subject as the 33rd percentile of their N3 epochs
STAGE_N3 = 3


def find_first_n3_transition(df):
    """
    Find the first sustained transition from Wake/N1/N2 to N3.
    'Sustained' = ≥ 3 consecutive N3 epochs.
    Returns the index of the first N3 epoch in the transition.
    """
    stages = df.stage.values
    n      = len(stages)

    for i in range(n - 3):
        # Check if current epoch is not N3 and next 3 are N3
        if stages[i] < 3 and all(stages[i+1:i+4] == 3):
            return i + 1  # First N3 epoch
    return None


def compute_rolling_stats(values, window):
    """Compute rolling mean and variance."""
    series = pd.Series(values)
    roll_mean = series.rolling(window, center=True, min_periods=1).mean()
    roll_var  = series.rolling(window, center=True, min_periods=1).var()
    return roll_mean.values, roll_var.values


def analyse_subject(df_sub):
    """
    Analyse one subject's complexity time series for VPM prediction.
    Returns dict with results, or None if transition not found.
    """
    df_sub = df_sub.dropna(subset=['aperiodic']).reset_index(drop=True)
    if len(df_sub) < PRE_EPOCHS + POST_EPOCHS:
        return None

    # Find the transition
    transition_idx = find_first_n3_transition(df_sub)
    if transition_idx is None:
        return None

    # Define analysis window
    start = max(0, transition_idx - PRE_EPOCHS)
    end   = min(len(df_sub), transition_idx + POST_EPOCHS)
    window_df = df_sub.iloc[start:end].copy()
    window_df['epoch_rel'] = window_df.index - transition_idx  # relative to transition

    # RSI proxy values in window
    rsi_proxy = window_df['aperiodic'].values

    # Per-subject N3 threshold = 25th percentile of all N3 epochs
    n3_epochs   = df_sub[df_sub.stage == 3]['aperiodic']
    wake_epochs  = df_sub[df_sub.stage == 0]['aperiodic']

    if len(n3_epochs) < 5 or len(wake_epochs) < 5:
        return None

    threshold = np.percentile(n3_epochs, 75)  # upper quartile of N3 = lower bound of N3 region

    # Rolling stats
    roll_mean, roll_var = compute_rolling_stats(rsi_proxy, WINDOW_EPOCHS)

    # Epoch positions relative to window start
    rel_positions = window_df['epoch_rel'].values  # epochs relative to transition

    # Find pre-transition epochs only (for VPM test)
    pre_mask = rel_positions < 0  # epochs before transition

    if pre_mask.sum() < 5:
        return None

    # ── VPM Core Prediction ──────────────────────────────────────────────────

    # t_peak_var: epoch (relative to transition) when variance is maximum
    #             in the pre-transition window
    pre_var  = roll_var[pre_mask]
    pre_rel  = rel_positions[pre_mask]
    peak_var_idx = np.argmax(pre_var)
    t_peak_var   = pre_rel[peak_var_idx]  # in epochs (negative = before transition)

    # t_cross_mean: epoch when rolling mean first crosses below threshold
    #              (in the full window, searching from start)
    cross_idx    = np.where(roll_mean <= threshold)[0]
    if len(cross_idx) == 0:
        t_cross_mean = POST_EPOCHS  # mean never crossed — assign max
    else:
        t_cross_mean = rel_positions[cross_idx[0]]

    # Lag: positive means variance peaked BEFORE mean crossed (VPM prediction)
    lag_epochs = t_cross_mean - t_peak_var
    lag_seconds = lag_epochs * EPOCH_SEC

    # ── Effect size (Cohen's d equivalent) ───────────────────────────────────
    pre_mean_values = roll_mean[pre_mask]
    pre_var_values  = roll_var[pre_mask]

    # Normalise to 0-1 range for comparison
    mean_range = pre_mean_values.max() - pre_mean_values.min() + 1e-10
    var_range  = pre_var_values.max()  - pre_var_values.min()  + 1e-10
    norm_mean  = (pre_mean_values - pre_mean_values.min()) / mean_range
    norm_var   = (pre_var_values  - pre_var_values.min())  / var_range

    # Correlation of normed var and normed mean (should be negative: var up, mean down)
    if len(norm_mean) > 3:
        r_vm, p_vm = stats.pearsonr(norm_var, norm_mean)
    else:
        r_vm, p_vm = np.nan, np.nan

    # ── Power-law fit (H2, exploratory) ──────────────────────────────────────
    # t_to_threshold = t_cross_mean - t (epochs remaining to threshold)
    # Var ~ t_to_threshold^β
    t_to_thresh = t_cross_mean - rel_positions[pre_mask]
    t_to_thresh = np.maximum(t_to_thresh, 0.5)  # avoid log(0)
    var_vals    = pre_var_values

    # Fit in log-log space
    beta_exp = np.nan
    if len(t_to_thresh) >= 5 and var_vals.max() > 0:
        try:
            log_t   = np.log(t_to_thresh)
            log_var = np.log(var_vals + 1e-30)
            slope, intercept, r2, p, se = stats.linregress(log_t, log_var)
            beta_exp = slope  # HRIT predicts β ≈ 0.546
        except Exception:
            pass

    return {
        'transition_idx':  transition_idx,
        'window_start':    start,
        'window_end':      end,
        't_peak_var':      t_peak_var,      # epochs before transition (negative)
        't_cross_mean':    t_cross_mean,    # epochs from window start
        'lag_epochs':      lag_epochs,      # VPM lag (positive = VPM confirmed)
        'lag_seconds':     lag_seconds,
        'threshold_used':  threshold,
        'n3_mean':         n3_epochs.mean(),
        'wake_mean':       wake_epochs.mean(),
        'r_var_mean':      r_vm,            # variance-mean correlation (should be neg)
        'p_var_mean':      p_vm,
        'beta_exp':        beta_exp,        # H2: power-law exponent (expect 0.546)
        'n_pre_epochs':    pre_mask.sum(),
        'n_n3_epochs':     len(n3_epochs),
        'rsi_time':        rel_positions.tolist(),
        'rsi_mean':        roll_mean.tolist(),
        'rsi_var':         roll_var.tolist(),
    }


# ── Main ─────────────────────────────────────────────────────────────────────

complexity_files = sorted(glob.glob(
    os.path.join(COMPLEXITY_DIR, 'SC*_complexity.csv')))

print(f'Found {len(complexity_files)} complexity files\n')

subject_results = []

for fpath in complexity_files:
    sid = os.path.basename(fpath).replace('_complexity.csv', '')
    df_sub = pd.read_csv(fpath)
    print(f'Analysing {sid}...')

    result = analyse_subject(df_sub)
    if result is None:
        print(f'  ✗ No usable Wake→N3 transition found')
        continue

    result['subject'] = sid
    subject_results.append(result)

    lag_s = result['lag_seconds']
    beta  = result['beta_exp']
    print(f'  ✓ Lag: {lag_s:+.0f}s | beta_exp: '
          f'{beta:.3f}' if not np.isnan(beta) else f'  ✓ Lag: {lag_s:+.0f}s | beta_exp: N/A')

if not subject_results:
    print('\nNo subjects had usable transitions.')
    print('Check that data was downloaded and complexity computed.')
    exit(1)

# ── Statistics ───────────────────────────────────────────────────────────────

results_df = pd.DataFrame(subject_results)
results_df.to_csv(os.path.join(RESULTS_DIR, 'vpm_subject_results.csv'),
                  index=False)

lags = results_df['lag_seconds'].dropna()
n    = len(lags)

print(f'\n{"="*55}')
print('VPM ANALYSIS RESULTS')
print(f'{"="*55}')
print(f'N subjects with valid transitions: {n}')
print(f'\nLag (variance peak before mean threshold crossing):')
print(f'  Mean:   {lags.mean():+.1f}s')
print(f'  Median: {lags.median():+.1f}s')
print(f'  SD:     {lags.std():.1f}s')
print(f'  Range:  [{lags.min():+.1f}, {lags.max():+.1f}]s')
print(f'  Positive lags (VPM confirmed): '
      f'{(lags > 0).sum()}/{n} ({100*(lags>0).mean():.0f}%)')

# H1: Wilcoxon signed-rank test (lag > 0, one-sided)
stat, p_two = stats.wilcoxon(lags, alternative='two-sided')
_, p_one    = stats.wilcoxon(lags, alternative='greater')

print(f'\nH1 (lag > 0): Wilcoxon signed-rank test')
print(f'  W = {stat:.1f}')
print(f'  p (one-sided, greater) = {p_one:.4f}')
print(f'  p (two-sided)          = {p_two:.4f}')
if p_one < 0.05:
    print(f'  → SUPPORTED: variance peaks before mean threshold crossing')
else:
    print(f'  → NOT SUPPORTED at α=0.05')

# Effect size: r = Z / sqrt(N)
from scipy.stats import norm as norm_dist
z_stat = norm_dist.ppf(1 - p_one)
r_eff  = z_stat / np.sqrt(n)
print(f'  Effect size r = {r_eff:.3f}')

# H2: Power-law exponent
betas = results_df['beta_exp'].dropna()
print(f'\nH2 (exploratory): Power-law exponent β')
print(f'  Mean β: {betas.mean():.3f} ± {betas.std():.3f}')
print(f'  HRIT prediction: β ≈ 0.546 (Verification 02)')
print(f'  95% CI: [{betas.mean() - 1.96*betas.std()/np.sqrt(len(betas)):.3f}, '
      f'{betas.mean() + 1.96*betas.std()/np.sqrt(len(betas)):.3f}]')
t_beta, p_beta = stats.ttest_1samp(betas.dropna(), 0.546)
print(f'  t-test vs 0.546: t={t_beta:.2f}, p={p_beta:.4f}')
if p_beta > 0.05:
    print(f'  → Consistent with predicted exponent (not significantly different from 0.546)')
else:
    print(f'  → Differs from predicted exponent 0.546')

# ── Save summary ─────────────────────────────────────────────────────────────
summary = {
    'n_subjects':       n,
    'lag_mean_sec':     lags.mean(),
    'lag_median_sec':   lags.median(),
    'lag_sd_sec':       lags.std(),
    'pct_positive_lag': 100 * (lags > 0).mean(),
    'wilcoxon_W':       stat,
    'p_one_sided':      p_one,
    'p_two_sided':      p_two,
    'effect_r':         r_eff,
    'beta_mean':        betas.mean() if len(betas) else np.nan,
    'beta_sd':          betas.std()  if len(betas) else np.nan,
    'beta_vs_predicted_p': p_beta    if len(betas) else np.nan,
}

pd.DataFrame([summary]).to_csv(
    os.path.join(RESULTS_DIR, 'vpm_summary_stats.csv'), index=False)

print(f'\nResults saved to: {RESULTS_DIR}/')
print('Next step: run 04_figures.py')
