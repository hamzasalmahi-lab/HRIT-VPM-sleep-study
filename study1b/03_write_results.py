"""
HRIT VPM Study 1b — Step 3: Generate Results Section
======================================================
Reads results/pews_summary_stats.csv and outputs a ready-to-paste
results section for the manuscript.

Output: results/manuscript_results_section.txt
"""

import os
import numpy as np
import pandas as pd

RESULTS_DIR = 'results'

summary    = pd.read_csv(os.path.join(RESULTS_DIR, 'pews_summary_stats.csv')).iloc[0]
subject_df = pd.read_csv(os.path.join(RESULTS_DIR, 'pews_subject_results.csv'))

n     = int(summary['n_subjects'])
n_ex  = int(summary.get('n_excluded', 0))
lags  = subject_df['lag_seconds'].dropna()
betas = subject_df['beta_exp'].dropna()


def p_str(p):
    if np.isnan(p):  return 'N/A'
    if p < 0.001:    return '< .001'
    if p < 0.01:     return f'= {p:.3f}'
    return f'= {p:.3f}'


ci_lo = summary['beta_mean'] - 1.96 * summary['beta_sd'] / np.sqrt(len(betas))
ci_hi = summary['beta_mean'] + 1.96 * summary['beta_sd'] / np.sqrt(len(betas))
pred_in_ci = ci_lo <= 0.546 <= ci_hi

text = f"""\
RESULTS — STUDY 1b
==================
Pre-registration: https://doi.org/10.17605/OSF.IO/EYMSF
Complexity proxy: Permutation entropy (order=6, delay=1) — decreases Wake→N3
N3 threshold:     75th percentile of per-subject N3 PE epochs
Dataset:          Sleep-EDF Expanded (PhysioNet), cassette study, Night 1, N=20
Relationship:     Independent replication of Study 1 using PE instead of aperiodic slope

Participants
------------
All {n + n_ex} subjects yielded an identifiable sustained Wake-to-N3 transition;
{n} subjects included, {n_ex} excluded per pre-specified criteria.

Primary Analysis — H1: Variance Precedes Mean
----------------------------------------------
Permutation entropy (PE) variance rose to its peak before the rolling mean
PE crossed the N3 threshold in {int(summary["n_positive_lag"])}/{n} subjects
({int(summary["pct_positive_lag"])}%).

Mean lag = {summary["lag_mean_sec"]:+.1f} s ({summary["lag_mean_min"]:.1f} min),
SD = {summary["lag_sd_sec"]:.1f} s, median = {summary["lag_median_sec"]:+.1f} s,
range [{lags.min():+.1f}, {lags.max():+.1f}] s.

One-sided Wilcoxon signed-rank: W = {summary["wilcoxon_W"]:.0f},
p {p_str(summary["p_one_sided"])}, r = {summary["effect_r"]:.3f}.
H1: {'CONFIRMED' if summary['h1_confirmed'] else 'NOT CONFIRMED'}

H2: Effect Size Replication
----------------------------
Pre-registered threshold: r ≥ .633 (Study 1 effect size).
Observed: r = {summary["effect_r"]:.3f}.
H2: {'CONFIRMED' if summary["h2_confirmed"] else 'NOT CONFIRMED'} \
({summary["effect_r"]:.3f} {"≥" if summary["h2_confirmed"] else "<"} .633)

H3 (Exploratory): Power-Law Exponent
-------------------------------------
HRIT Verification 02 predicts β ≈ 0.546.
Observed: β = {summary["beta_mean"]:.3f} (SD = {summary["beta_sd"]:.3f}, N = {len(betas)};
95% CI [{ci_lo:.3f}, {ci_hi:.3f}]).
β = 0.546 is {'inside' if pred_in_ci else 'outside'} the 95% CI.
t-test vs 0.546: p {p_str(summary["beta_vs_predicted_p"])}.
H3: {'CONSISTENT WITH HRIT prediction ✓' if summary["h3_confirmed"] else 'NOT CONFIRMED'}

Summary Statistics
------------------
| Measure                        | Value                    |
|--------------------------------|--------------------------|
| N included                     | {n}                        |
| N excluded                     | {n_ex}                        |
| PEWS confirmed (lag > 0)       | {int(summary["n_positive_lag"])}/{n} ({summary["pct_positive_lag"]:.0f}%)              |
| Mean lag                       | {summary["lag_mean_sec"]:+.1f} s ({summary["lag_mean_min"]:.1f} min)        |
| Median lag                     | {summary["lag_median_sec"]:+.1f} s                   |
| SD lag                         | {summary["lag_sd_sec"]:.1f} s                    |
| Wilcoxon W                     | {summary["wilcoxon_W"]:.0f}                      |
| p (one-sided)                  | {p_str(summary["p_one_sided"])}                  |
| r                              | {summary["effect_r"]:.3f}                    |
| H1 (lag > 0)                   | {'CONFIRMED' if summary['h1_confirmed'] else 'NOT CONFIRMED'}                |
| H2 (r ≥ .633)                  | {'CONFIRMED' if summary['h2_confirmed'] else 'NOT CONFIRMED'}                |
| Mean β (observed)              | {summary["beta_mean"]:.3f} (SD = {summary["beta_sd"]:.3f})          |
| 95% CI for β                   | [{ci_lo:.3f}, {ci_hi:.3f}]         |
| HRIT predicted β               | 0.546                    |
| p (β vs 0.546)                 | {p_str(summary["beta_vs_predicted_p"])}                  |
| H3 (β ≈ 0.546)                 | {'CONSISTENT' if summary['h3_confirmed'] else 'NOT CONFIRMED'}                |

Comparison with Study 1 (Aperiodic Slope)
------------------------------------------
| Measure        | Study 1 (Aperiodic)  | Study 1b (PE)        |
|----------------|----------------------|----------------------|
| N confirmed    | 18/18 (100%)         | {int(summary["n_positive_lag"])}/{n} ({int(summary["pct_positive_lag"])}%)        |
| Mean lead time | 3.6 min              | {summary["lag_mean_min"]:.1f} min               |
| W              | 171                  | {summary["wilcoxon_W"]:.0f}                   |
| p              | < .001               | {p_str(summary["p_one_sided"])}              |
| r              | .882                 | {summary["effect_r"]:.3f}                |
| H1             | CONFIRMED            | {'CONFIRMED' if summary['h1_confirmed'] else 'NOT CONFIRMED'}            |
| H2 (r ≥ .633)  | N/A                  | {'CONFIRMED' if summary['h2_confirmed'] else 'NOT CONFIRMED'}            |
| H3 (β = 0.546) | NOT CONFIRMED        | {'CONFIRMED' if summary['h3_confirmed'] else 'NOT CONFIRMED'}            |

Limitations
-----------
1. Threshold: 75th percentile of N3 PE is an approximation of the theoretical
   consciousness threshold; the exact HRIT threshold requires SPM25 simulation.
2. Single channel: EEG Fpz-Cz only. Multi-channel HEP and PCI require further studies.
3. Lead time shorter than the pre-registered 6.6 min — re-run on full raw EDF
   files with exact pre-registered parameters for definitive values.
"""

out = os.path.join(RESULTS_DIR, 'manuscript_results_section.txt')
with open(out, 'w') as f:
    f.write(text)

print(f'✓  Results section → {out}')
print(f'\nKey results:')
print(f'  H1: {int(summary["n_positive_lag"])}/{n} PEWS confirmed  '
      f'W={summary["wilcoxon_W"]:.0f}  p {p_str(summary["p_one_sided"])}  '
      f'r={summary["effect_r"]:.3f}')
print(f'  H2: r={summary["effect_r"]:.3f} {"≥" if summary["h2_confirmed"] else "<"} .633  '
      f'→ {"CONFIRMED" if summary["h2_confirmed"] else "NOT CONFIRMED"}')
print(f'  H3: β={summary["beta_mean"]:.3f}  '
      f'p_vs_0.546 {p_str(summary["beta_vs_predicted_p"])}  '
      f'→ {"CONFIRMED" if summary["h3_confirmed"] else "NOT CONFIRMED"}')
