# Study 1b — Convergent Replication of PEWS Using Permutation Entropy

**Pre-registration:** https://doi.org/10.17605/OSF.IO/EYMSF  
**Analysis repo:** https://github.com/hamzasalmahi-lab/HRIT-v3-computational  
**Dataset:** Sleep-EDF Expanded, PhysioNet — same N=20 subjects as Study 1  
**Pre-registered:** Prior to analysis, timestamped on OSF

---

## What Study 1b tested

Study 1b replicated the PEWS prediction from Study 1 using permutation
entropy (PE) as an independent I_Sim proxy, applied to the same Sleep-EDF
Expanded dataset with an independent analysis pipeline.

Three pre-registered hypotheses:

| Hypothesis | Prediction | Result |
|------------|-----------|--------|
| H1 | PE stratifies: Wake > REM > N2 > N3 | **Confirmed** — 18/19 subjects, p = .0001, r = .860 |
| H2 | r ≥ .633 (replicates Study 1 effect) | **Confirmed** — r = .860 ≥ .633 |
| H3 | β = 0.546 (power-law exponent) | Not confirmed — β = 0.284, 95% CI [0.085, 0.483] |

Mean lead time: **6.6 minutes** (SD = 3.0) before N3 onset  
Wilcoxon: W = 0, p = .0001, r = .860

---

## Relationship to Study 1

Study 1 tested the PEWS prediction using aperiodic slope as the I_Sim
proxy. Study 1b replicated it with permutation entropy using an independent
pipeline. The two studies use:

- Same dataset (Sleep-EDF Expanded, N=20 cassette subjects)
- Same pre-registered prediction (variance rises before threshold crossing)  
- Different complexity proxies (aperiodic slope vs permutation entropy)
- Different analysis pipelines (independent code)

Convergent confirmation across two independent complexity measures
strengthens the PEWS empirical evidence base.

---

## Key files (in HRIT-v3-computational repo)

```
analysis.py                    Main pipeline (PE-based CII estimation)
confirmatory_tests.py          Pre-registered H1/H2/H3 tests
results_per_subject_stage.csv  Per-subject per-stage PE and HEP means
confirmatory_results.log       Full statistical output
Study1_Figure1_CII_Stages.png  Sleep stage PE stratification
Study1_Figure2_PEWS.png        PEWS early warning signal
```

---

## Citation

Almahi, H. S. (2026). Permutation Entropy as an Empirical Proxy for the
Simulation Component of the Conscious Integration Index: Validation across
Sleep Stages and a Pre-Cognitive Early Warning Signal.
Pre-registration: https://doi.org/10.17605/OSF.IO/EYMSF
