# Study 1b — Permutation Entropy Replication of PEWS

**Pre-registration:** https://doi.org/10.17605/OSF.IO/EYMSF  
**Analysis repo:**    https://github.com/hamzasalmahi-lab/HRIT-v3-computational  
**Dataset:**          Sleep-EDF Expanded — same N = 20 subjects as Study 1  
**Status:**           ✅ Complete — confirmed

---

## What Study 1b tested

Study 1b replicated the PEWS prediction from Study 1 using **permutation
entropy (PE)** as an independent I_Sim proxy, with an entirely independent
analysis pipeline applied to the same Sleep-EDF Expanded dataset.

PE **decreases** from Wake to N3 (Wake > REM > N2 > N3), making it a
conceptually cleaner CII proxy than aperiodic slope for this test — PE
directly indexes the ordinal regularity of the EEG signal across scales.

### Three pre-registered hypotheses

| Hypothesis | Prediction | Result |
|---|---|---|
| **H1** | PE stratifies: Wake > REM > N2 > N3, in 18/20 subjects | ✅ **Confirmed** — 18/19, p = .0001, r = .860 |
| **H2** | r ≥ .633 (replicates Study 1 effect size threshold) | ✅ **Confirmed** — r = .860 ≥ .633 |
| **H3** | β = 0.546 (power-law variance exponent) | ❌ Not confirmed — β = 0.284, 95% CI [0.085, 0.483] |

**Lead time: 6.6 min (SD = 3.0)** before first sustained N3 onset  
**Wilcoxon:** W = 0, p = .0001, r = .860

---

## Relationship to Study 1

| | Study 1 | Study 1b |
|---|---|---|
| Dataset | Sleep-EDF Expanded, N = 20 | Same |
| Proxy | Aperiodic slope | Permutation entropy |
| Direction | Rises Wake→N3 | Falls Wake→N3 |
| Pipeline | This repo | HRIT-v3-computational |
| Lead time | 3.6 min | 6.6 min |
| p | < .0001 | < .0001 |
| r | .882 | .860 |
| H2 (β = 0.546) | Not confirmed | Not confirmed |

Two independent complexity measures, same PEWS prediction, same dataset,
independent pipelines — convergent confirmation of the PEWS signal.

---

## How to reproduce Study 1b

```bash
git clone https://github.com/hamzasalmahi-lab/HRIT-v3-computational
cd HRIT-v3-computational
pip install mne numpy scipy pandas matplotlib

# Run the Study 1b pipeline
python analysis.py                  # PE-based CII per sleep stage
python confirmatory_tests.py        # H1 / H2 / H3 Wilcoxon tests
```

---

## Key files (HRIT-v3-computational repo)

```
analysis.py                     Main pipeline — PE as CII proxy
confirmatory_tests.py           Pre-registered H1 / H2 / H3 tests
results_per_subject_stage.csv   Per-subject per-stage means
confirmatory_results.log        Full statistical output
Study1_Figure1_CII_Stages.png   Sleep stage PE stratification
Study1_Figure2_PEWS.png         PEWS early warning signal
```

---

## Pre-registration (AsPredicted, OSF)

- **DOI:** https://doi.org/10.17605/OSF.IO/EYMSF  
- **Registered:** prior to confirmatory analysis  
- **Data:** already collected (pre-registration filed before analysis)

---

## Citation

Almahi, H. S. (2026). *Permutation Entropy as an Empirical Proxy for the
Simulation Component of the Conscious Integration Index: Validation across
Sleep Stages and a Pre-Cognitive Early Warning Signal.*  
Pre-registration: https://doi.org/10.17605/OSF.IO/EYMSF
