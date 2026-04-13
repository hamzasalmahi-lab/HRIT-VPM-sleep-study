# OSF Pre-Registration: HRIT VPM Sleep EEG Study
## Empirical Study 1 of the HRIT Empirical Programme

**Title:** Variance-Precedes-Mean in EEG Complexity at the Sleep-Stage
Consciousness Boundary: A Reanalysis of the Sleep-EDF Database in Test
of the HRIT VPM Prediction

**Authors:** Hamza S. Almahi

**Date of pre-registration:** [complete before analysis]

**OSF Project DOI:** [to be assigned on submission]

---

## Background and Theoretical Motivation

The Hierarchical Recursive Integration Theory (HRIT; Almahi, 2025)
formally derives the Variance-Precedes-Mean (VPM) principle from the
saddle-node bifurcation of the allostatic control ODE. The derivation
yields two quantitative predictions:

1. Near the consciousness stability boundary, the escape time from
   the conscious attractor scales as:
   τ_escape ~ (z_c − z)^{−0.500}
   where z is allostatic load and z_c is the critical load.
   [HRIT Proof 08, mathematically exact]

2. EEG complexity variance scales as:
   Var(RSI_proxy) ~ (z_c − z)^{+0.546}
   [HRIT Verification 02, computational; consistent with 0.5 at p > 0.05]

The VPM principle predicts that variance in neural complexity rises to a
maximum BEFORE the mean complexity crosses the threshold into
unconsciousness. This constitutes a clinically actionable early-warning
signal: rising neural complexity variance precedes the loss of
phenomenal presence.

The sleep EEG context: as subjects transition from waking to N3 slow-wave
sleep (SWS), they traverse a depth-collapse transition — all integration
components suppressed monotonically. HRIT predicts this transition should
show VPM: complexity variance should peak before mean complexity crosses
the N3 threshold.

---

## Hypotheses

**H1 (Primary, pre-registered):**
In the 5-minute window preceding the first sustained Wake-to-N3
transition, the rolling variance of EEG complexity (aperiodic slope)
reaches its maximum BEFORE the rolling mean crosses the N3 threshold.

Operationalised: the temporal lag (t_cross_mean − t_peak_var) is
significantly greater than zero.

- One-sided Wilcoxon signed-rank test, α = 0.05
- N ≥ 20 subjects with valid transitions required

**H_null:** The temporal lag is ≤ 0 (variance and mean decline
simultaneously; no early warning signal).

**H2 (Secondary, exploratory, NOT pre-registered for primary inference):**
Variance scales as a power law with time-to-threshold with exponent
β ≈ 0.546. The 95% CI of the mean estimated β includes 0.546.

---

## Dataset

- **Name:** Sleep-EDF Expanded Database
- **Source:** PhysioNet (https://physionet.org/content/sleep-edfx/1.0.0/)
- **License:** Open Data Commons Attribution License v1.0
- **Subjects:** N = 197 (cassette study), or N = 20 pilot subset
- **Data type:** Full-night polysomnography (.edf) + hypnogram annotations

No new data collection. This is a secondary analysis of an open-access
public database. No ethics approval required.

---

## Analysis Plan

### Step 1: Subject selection
Include subjects with:
- At least 20 pre-transition epochs (≥ 10 minutes of pre-N3 data)
- At least 5 N3 (SWS) epochs
- Valid hypnogram with stage annotations
- At least one EEG channel (EEG Fpz-Cz preferred)

### Step 2: RSI proxy computation
For each 30-second epoch:
- Primary proxy: aperiodic slope (−β from log-log PSD fit, 2–40 Hz)
- Secondary proxy: permutation entropy (order = 6)
- AFEM sensor: alpha-band (8–13 Hz) power

### Step 3: Transition identification
First sustained Wake→N3 transition: first epoch with three consecutive
N3 epochs following a period of Wake/N1/N2.

### Step 4: Analysis window
− 10 minutes (−20 epochs) to + 5 minutes (+10 epochs) relative to
transition.

### Step 5: Rolling statistics
Rolling mean and variance in a 3-minute (6-epoch) window.

### Step 6: Key quantities
- t_peak_var: epoch of maximum variance in pre-transition window
- t_cross_mean: first epoch where rolling mean crosses the N3 threshold
- Lag = t_cross_mean − t_peak_var

### Step 7: N3 threshold per subject
75th percentile of aperiodic slope values in all N3 epochs of that
subject (upper bound of N3 complexity = lower bound of N3 classification).

### Step 8: Statistical test (H1)
One-sided Wilcoxon signed-rank test: H0: median lag ≤ 0.
Report: W statistic, p-value (one-sided), effect size r = Z/√N.

### Step 9: Power-law fit (H2)
Log-log regression of variance vs time-to-threshold in pre-transition
window. Report mean β, SD, 95% CI, one-sample t-test vs 0.546.

---

## Deviations from plan

Any deviation from this pre-registered plan will be reported explicitly
in the manuscript as a post-hoc analysis, clearly labelled and not
contributing to the primary inference.

---

## Power analysis

Published sleep EEG studies find VPM-consistent effects (variance rising
before threshold) with effect sizes r ≈ 0.40–0.60 in cognitive
transitions. At r = 0.40, α = 0.05 (one-sided), power = 0.80 requires
N = 35 subjects. At r = 0.55, N = 19. We target N = 20 for the pilot
and N = 40+ for the full analysis.

---

## Timeline

- Week 1: Pre-registration submitted to OSF
- Week 1: Dataset downloaded (Sleep-EDF Expanded)
- Week 2: Preprocessing and complexity computation
- Week 3: Pre-registered analysis
- Week 4: Figures and manuscript draft
- Week 5–6: Internal revision
- Week 7: bioRxiv preprint

---

## Planned submission venue

Primary: *Brain and Cognition* (short report, 2500 words)
Backup: *Consciousness and Cognition* (short report)
Backup: *eLife* (brief report)

---

## Conflicts of interest

None.

## Funding

None. Independent research.

## Code availability

Full analysis code will be deposited on OSF and GitHub at:
github.com/hamzasalmahi-lab/HRIT-VPM-sleep-study

---

*This pre-registration follows the AsPredicted format and will be
submitted to OSF at https://osf.io/prereg/ before any analysis is run.
The timestamp on the OSF pre-registration constitutes the binding
commitment to the analysis plan.*
