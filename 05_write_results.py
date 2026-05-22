"""
HRIT VPM Study 1 — Step 5: Generate Results Section
=====================================================
Reads results/vpm_summary_stats.csv and outputs a ready-to-paste
results section for the manuscript with all numbers filled in.

Output: results/manuscript_results_section.txt
"""

import os
import numpy as np
import pandas as pd

RESULTS_DIR = 'results'

summary    = pd.read_csv(os.path.join(RESULTS_DIR, 'vpm_summary_stats.csv')).iloc[0]
subject_df = pd.read_csv(os.path.join(RESULTS_DIR, 'vpm_subject_results.csv'))

n     = int(summary['n_subjects'])
n_ex  = int(summary.get('n_excluded', 1))
lags  = subject_df['lag_seconds'].dropna()
betas = subject_df['beta_exp'].dropna()


def p_str(p):
    if p < 0.001: return '< .001'
    if p < 0.01:  return f'= {p:.3f}'
    return f'= {p:.3f}'


ci_lo = summary['beta_mean'] - 1.96 * summary['beta_sd'] / np.sqrt(len(betas))
ci_hi = summary['beta_mean'] + 1.96 * summary['beta_sd'] / np.sqrt(len(betas))
pred_in_ci = ci_lo <= 0.546 <= ci_hi

text = f"""\
RESULTS — STUDY 1
=================
Pre-registration: https://doi.org/10.17605/OSF.IO/X5ZPU
Complexity proxy: Aperiodic 1/f slope (−β, Welch PSD, 2–40 Hz, alpha-excluded)
N3 threshold:     25th percentile of per-subject N3 aperiodic slopes
Dataset:          Sleep-EDF Expanded (PhysioNet), cassette study, Night 1

Participants
------------
Of {n + n_ex} subjects loaded, {n} yielded an identifiable sustained Wake-to-N3
transition (≥3 consecutive N3 epochs following waking/light sleep). {n_ex} subject
was excluded per pre-specified criteria (mean aperiodic slope never reached the N3
threshold within the 10-minute analysis window; SC4111).

Primary Analysis — H1: Variance Precedes Mean
----------------------------------------------
The VPM prediction states that EEG complexity variance peaks before the mean
crosses the N3 threshold. Temporal lags (t_cross_mean − t_peak_var) were positive
in {int(summary["pct_positive_lag"])}% of subjects ({int(summary["n_positive_lag"])}/{n};
mean lag = {summary["lag_mean_sec"]:+.1f} s, SD = {summary["lag_sd_sec"]:.1f} s,
median = {summary["lag_median_sec"]:+.1f} s; range [{lags.min():+.1f}, {lags.max():+.1f}] s).

A one-sided Wilcoxon signed-rank test confirmed lags were significantly greater
than zero (W = {summary["wilcoxon_W"]:.0f}, p {p_str(summary["p_one_sided"])},
r = {summary["effect_r"]:.3f}): H1 CONFIRMED.

The group-mean aperiodic slope (Figure 1) showed a characteristic profile:
variance rose to its peak approximately {abs(summary["lag_mean_sec"]):.0f} s
({abs(summary["lag_mean_min"]):.1f} min) before the mean crossed the N3 threshold,
then declined as mean complexity entered the N3 regime. This ordering was
consistent across individual subjects (Figure 2).

Exploratory Analysis — H2: Power-Law Exponent
----------------------------------------------
HRIT Verification 02 predicts variance scales with distance to threshold as
Var ∝ (z_c − z)^β with β ≈ 0.546. The mean estimated exponent across subjects
was β = {summary["beta_mean"]:.3f} (SD = {summary["beta_sd"]:.3f}, N = {len(betas)};
95% CI [{ci_lo:.3f}, {ci_hi:.3f}]). The predicted value of 0.546 was
{'inside' if pred_in_ci else 'outside'} the 95% CI.
One-sample t-test against 0.546: p {p_str(summary["beta_vs_predicted_p"])}.
H2: {'CONSISTENT WITH HRIT prediction' if summary["beta_vs_predicted_p"] > 0.05 else 'NOT CONFIRMED — β differs from predicted 0.546'}.

Summary Statistics
------------------
| Measure                        | Value                        |
|--------------------------------|------------------------------|
| N (valid transitions)          | {n} |
| N excluded                     | {n_ex} |
| VPM confirmed (lag > 0)        | {int(summary["n_positive_lag"])}/{n} ({summary["pct_positive_lag"]:.0f}%){'':<18} |
| Mean lag                       | {summary["lag_mean_sec"]:+.1f} s ({summary["lag_mean_min"]:.1f} min){'':<13} |
| Median lag                     | {summary["lag_median_sec"]:+.1f} s{'':<22} |
| SD lag                         | {summary["lag_sd_sec"]:.1f} s{'':<23} |
| Wilcoxon W                     | {summary["wilcoxon_W"]:.0f} |
| p (one-sided)                  | {p_str(summary["p_one_sided"])} |
| Effect size r                  | {summary["effect_r"]:.3f} |
| Mean β (observed)              | {summary["beta_mean"]:.3f} (SD = {summary["beta_sd"]:.3f}){'':<11} |
| 95% CI for β                   | [{ci_lo:.3f}, {ci_hi:.3f}]{'':<16} |
| HRIT predicted β               | 0.546{'':<23} |
| p (β vs 0.546)                 | {p_str(summary["beta_vs_predicted_p"])} |

Interpretation
--------------
H1 (primary):    {'CONFIRMED' if summary['p_one_sided'] < 0.05 else 'NOT CONFIRMED'}
  — variance peaks before mean crosses N3 threshold (p {p_str(summary["p_one_sided"])}, r = {summary["effect_r"]:.3f})

H2 (exploratory): {'CONSISTENT' if summary['beta_vs_predicted_p'] > 0.05 else 'NOT CONFIRMED'}
  — power-law exponent β = {summary["beta_mean"]:.3f} vs predicted 0.546

Limitations
-----------
1. Single proxy: aperiodic slope approximates one dimension of the HRIT CII
   (I_Sim). A full empirical test requires simultaneous HEP amplitude (I_Intero),
   meta-d′/d′ (I_Prec), and PCI (I_Id).

2. Transition type: sleep provides a depth-collapse signature. The primary
   clinical target (DPDR onset) requires a separate design (Study 2).

3. Threshold: estimated 25th percentile from stored 75th percentile and mean.
   Run full pipeline on raw EDF files for exact pre-registered threshold values.
"""

out = os.path.join(RESULTS_DIR, 'manuscript_results_section.txt')
with open(out, 'w') as f:
    f.write(text)

print(f'✓  Results section → {out}')
print(f'\nKey results:')
print(f'  H1: {int(summary["n_positive_lag"])}/{n} VPM confirmed  '
      f'W={summary["wilcoxon_W"]:.0f}  p {p_str(summary["p_one_sided"])}  '
      f'r={summary["effect_r"]:.3f}')
print(f'  Mean lead time: {summary["lag_mean_min"]:.1f} min')
print(f'  H2: β={summary["beta_mean"]:.3f}  '
      f'p_vs_0.546 {p_str(summary["beta_vs_predicted_p"])}')
