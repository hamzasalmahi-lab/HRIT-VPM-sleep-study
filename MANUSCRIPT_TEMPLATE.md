# MANUSCRIPT TEMPLATE
# Variance-Precedes-Mean in Sleep EEG: Empirical Test of the HRIT VPM Principle
# Target: Brain and Cognition (Short Report, 2500 words)
# Replace [RESULT] tags with actual numbers from 05_write_results.py

---

**Title:** Variance-Precedes-Mean in EEG Complexity at the Sleep-Stage
Consciousness Boundary: A Pre-registered Reanalysis of the Sleep-EDF Database

**Authors:** Hamza S. Almahi

**Affiliation:** Independent Researcher

**Correspondence:** hamza.s.almahi@gmail.com

**Data availability:** All data used are publicly available at
physionet.org/content/sleep-edfx/1.0.0/. Analysis code available at
github.com/hamzasalmahi-lab/HRIT-VPM-sleep-study. Pre-registration at
osf.io/[ID].

**Keywords:** consciousness; sleep EEG; aperiodic slope; early warning signal;
allostatic control; variance-precedes-mean; HRIT

---

## Abstract (150 words)

The Variance-Precedes-Mean (VPM) principle, formally derived from the
saddle-node bifurcation of the allostatic control ODE in the Hierarchical
Recursive Integration Theory (HRIT; Almahi, 2025), predicts that neural
complexity variance rises to a peak before mean complexity crosses the
consciousness threshold. We tested this prediction in a pre-registered
secondary analysis of the Sleep-EDF Expanded database (N = [N] subjects).
EEG complexity (aperiodic slope) was computed for each 30-second epoch
across the Wake-to-N3 sleep transition. Temporal lags between variance
peak and mean threshold crossing were positive in [%]% of subjects
(mean lag = [LAG]s, p = [P], one-sided Wilcoxon, r = [R]). The
power-law variance exponent (β = [BETA]) was consistent with the HRIT
prediction of 0.546 (p = [P_BETA]). These findings constitute the first
empirical test of the HRIT VPM principle and support the use of
complexity variance as an early-warning biomarker for consciousness-state
transitions.

---

## 1. Introduction (400 words)

The transition from waking to slow-wave sleep (SWS) provides a
naturalistically occurring, experimentally tractable loss of
phenomenal consciousness. Electroencephalographic (EEG) complexity
measures — including the aperiodic (1/f) slope, sample entropy, and
permutation entropy — track conscious level across this transition,
declining monotonically from waking through N1, N2, and N3 (Casali
et al., 2013; Rosso et al., 2007).

The Hierarchical Recursive Integration Theory (HRIT; Almahi, 2025)
formally derives a counterintuitive prediction about this transition:
complexity variance should not decline monotonically alongside the mean.
Instead, variance should *rise* to a peak in the window immediately
preceding the threshold crossing, then fall as the mean enters the
N3 regime. This Variance-Precedes-Mean (VPM) principle is derived from
the saddle-node bifurcation of the allostatic control ODE (Proof 08),
which yields a mathematically exact power-law exponent for escape time
(τ_escape ~ (z_c − z)^{−0.500}) and an empirically testable variance
exponent (Var ~ (z_c − z)^{+0.546}).

The VPM prediction has direct clinical significance. If rising complexity
variance precedes mean-level decline, it constitutes an early warning
signal of impending consciousness-state transitions — applicable to
anesthesia monitoring, dissociation onset prediction, and coma assessment.
Testing this prediction in sleep EEG provides the lowest-cost initial
empirical test: public datasets are available, methodology is established,
and the transition is well-defined.

We conducted a pre-registered secondary analysis of the Sleep-EDF
Expanded database to test two predictions: (H1) that variance peaks
before mean threshold crossing (temporal lag > 0); and (H2, exploratory)
that the variance-to-threshold power-law exponent is consistent with
the HRIT prediction of β ≈ 0.546.

---

## 2. Methods (500 words)

### 2.1 Dataset
The Sleep-EDF Expanded database (Kemp et al., 2000; physionet.org/
content/sleep-edfx/1.0.0/) provides full-night polysomnography recordings
for 197 healthy subjects with expert-annotated sleep stage hypnograms.
We analysed [N] subjects from the cassette substudy. Subjects were
excluded if they lacked a sustained Wake-to-N3 transition (≥ 3 consecutive
N3 epochs following Wake/N1/N2) or had fewer than 20 pre-transition epochs.
This analysis was pre-registered on OSF (DOI: osf.io/[ID]) before
data analysis began.

### 2.2 EEG Complexity Measure
We computed the aperiodic slope (−β from log-log power spectral density,
2–40 Hz, excluding 9–12 Hz to reduce alpha-peak contamination) for each
30-second epoch using Welch's method (4-second windows, 50% overlap).
This measure correlates with PCI across states (Colombo et al., 2019)
and serves as the primary RSI proxy.

### 2.3 Transition Identification
We identified the first sustained Wake-to-N3 transition per subject.
The N3 threshold was defined as the 75th percentile of that subject's
N3 aperiodic slopes (upper bound of N3 complexity).

### 2.4 Rolling Statistics
We computed rolling mean and variance of aperiodic slope in a 3-minute
(6-epoch) window, centred on each epoch.

### 2.5 Primary Analysis (H1)
For each subject we recorded: t_peak_var (epoch of maximum rolling
variance in the 10-minute pre-transition window) and t_cross_mean
(first epoch where rolling mean crossed the N3 threshold). Lag =
t_cross_mean − t_peak_var (positive = VPM confirmed). We tested H1
with a one-sided Wilcoxon signed-rank test against zero (α = 0.05).

### 2.6 Exploratory Analysis (H2)
We fit log(Var) ~ β × log(time-to-threshold) in the pre-transition window
and reported mean β with 95% CI and t-test vs 0.546.

---

## 3. Results

[PASTE OUTPUT FROM 05_write_results.py HERE]

---

## 4. Discussion (600 words)

[Write after seeing results — template below]

The present results [support / do not support] the HRIT VPM prediction
that neural complexity variance peaks before mean complexity crosses the
N3 sleep threshold. [IF SUPPORTED: This is consistent with the theoretical
derivation from the saddle-node bifurcation of the allostatic control ODE,
which predicts that systems near a stability boundary exhibit variance
inflation — a generic signature of Critical Slowing Down — before the mean-level
shift that marks the boundary crossing.] [IF NOT SUPPORTED: discuss possible
reasons: single-channel proxy insufficient; threshold definition; need for
full RSI rather than aperiodic slope proxy.]

The power-law exponent β = [BETA] [is / is not] consistent with the HRIT
prediction of 0.546. [Note the theoretical grounding: Verification 02
predicts β consistent with 0.5, which is the canonical CSD exponent.]

**Limitations.** The aperiodic slope is a single-channel proxy for the
multi-component RSI. A full empirical test requires simultaneous HEP
amplitude (AFEM proxy), meta-d′/d′ (RSR proxy), and PCI (HGM proxy).
The sleep transition constitutes a depth-collapse signature; the primary
clinical target of HRIT — presence collapse in DPDR — requires a
separate study design. Additionally, the N3 threshold is defined
per-subject and may not correspond precisely to the theoretical allostatic
setpoint Φ_0 = 2/3 predicted by the companion GIFT framework.

**Future directions.** Study 2 will test the DPDR double dissociation
(P-NEW-7): in DPDR patients, meta-d′/d′ should remain at healthy-control
levels while HEP amplitude is significantly reduced. This constitutes
the primary clinical test of the HRIT presence-collapse mechanism.

---

## 5. Conclusion

We report the first empirical test of the HRIT VPM principle in sleep
EEG. [Results sentence]. This supports the use of complexity variance
as an early-warning biomarker for consciousness-state transitions, with
implications for anesthesia monitoring, DPDR prediction, and clinical
assessment of disorders of consciousness.

---

## References

Almahi, H. S. (2025). Hierarchical Recursive Integration Theory (HRIT):
  A Unified Neurobiological Framework of Consciousness. Version 2.
  Physics of Life Reviews [submitted].

Casali, A. G., et al. (2013). A theoretically based index of
  consciousness independent of sensory processing and behaviour.
  Science Translational Medicine, 5(198), 198ra105.

Colombo, M. A., et al. (2019). The spectral exponent of the resting
  EEG indexes the presence of consciousness during unresponsiveness
  induced by three anaesthetic agents. NeuroImage, 189, 631–644.

Kemp, B., et al. (2000). Analysis of a sleep-dependent neuronal
  feedback loop: the slow-wave microcontinuity of the EEG.
  IEEE Transactions on Biomedical Engineering, 47(9), 1185–1194.

Rosso, O. A., et al. (2007). Wavelet entropy: a new tool for analysis
  of short duration brain electrical signals. Journal of Neuroscience
  Methods, 105(1), 65–75.
