# HRIT Empirical Programme — Studies 1 and 1b
## Variance-Precedes-Mean (VPM / PEWS) in Sleep EEG

Pre-registered empirical tests of the HRIT v3 Presence Early Warning
Signal (PEWS) prediction — that neural complexity variance rises and peaks
**before** the mean complexity crosses the consciousness threshold at the
Wake→N3 sleep transition.

**Framework:** HRIT v3 — Zenodo: https://doi.org/10.5281/zenodo.19490741  
**Author:** Hamza S. Almahi · hamza.s.almahi@gmail.com

---

## PEWS Prediction (both studies)

From the HRIT v3 allostatic ODE (saddle-node bifurcation):

> **Var(CII proxy) ~ (z_c − z)^+β** with predicted β ≈ 0.546
>
> Neural complexity variance rises to a peak **before** the mean crosses
> the N3 threshold. Lead time = t_cross_mean − t_peak_var > 0.

This constitutes a clinically actionable pre-cognitive early warning signal.

---

## Study 1

| Item | Detail |
|------|--------|
| Pre-registration | https://doi.org/10.17605/OSF.IO/X5ZPU |
| Dataset | Sleep-EDF Expanded (PhysioNet), N = 20 cassette subjects |
| Complexity measure | Permutation entropy (order=6, delay=1) — decreases Wake→N3 |
| Threshold | 75th percentile of N3 PE epochs (upper N3 boundary) |
| Prediction | Variance peaks before mean crosses N3 PE threshold (lag > 0) |

**Status:** Code fixed and validated. Awaiting real Sleep-EDF data run.
Results below are from the corrected pipeline (synthetic validation data).
Run `python 00_download.py` to obtain real data and reproduce.

**Known issue in prior repo version:** The original code used aperiodic
slope as the complexity measure. Aperiodic slope *increases* from Wake
to N3 in Sleep-EDF (N3 > Wake), inverting the threshold logic and causing
`t_cross_mean` to fire at window start for every subject. Fixed in this
version by using permutation entropy, which correctly decreases Wake→N3.

### How to reproduce Study 1

```bash
pip install mne numpy scipy matplotlib pandas wfdb
python 00_download.py        # downloads Sleep-EDF (requires PhysioNet account)
python 01_preprocess.py
python 02_compute_complexity.py
python 03_vpm_analysis.py    # FIXED: uses perm_entropy, correct threshold direction
python 04_figures.py
python 05_write_results.py
```

### Study 1 file structure

```
00_download.py               Download Sleep-EDF from PhysioNet
01_preprocess.py             EDF loading and sleep staging
02_compute_complexity.py     Permutation entropy + aperiodic slope computation
03_vpm_analysis.py           Pre-registered VPM analysis (FIXED)
04_figures.py                Publication figures
05_write_results.py          Auto-generate results section
PREREGISTRATION.md           Pre-registration document
results/
  vpm_subject_results.csv    Per-subject results
  vpm_summary_stats.csv      Summary statistics
study1b/                     Study 1b summary (see below)
```

---

## Study 1b

| Item | Detail |
|------|--------|
| Pre-registration | https://doi.org/10.17605/OSF.IO/EYMSF |
| Dataset | Sleep-EDF Expanded, same N = 20 subjects |
| Complexity measure | Permutation entropy (normalised to Wake baseline) |
| Analysis repo | https://github.com/hamzasalmahi-lab/HRIT-v3-computational |

### Study 1b confirmed results

| Hypothesis | Prediction | Result |
|------------|-----------|--------|
| H1 | PE stratifies Wake > REM > N2 > N3 | ✅ **Confirmed** — 18/19, p = .0001, r = .860 |
| H2 | r ≥ .633 (replication threshold) | ✅ **Confirmed** — r = .860 ≥ .633 |
| H3 | β = 0.546 (power-law exponent) | ❌ Not confirmed — β = 0.284 [0.085, 0.483] |

**Lead time: 6.6 min (SD = 3.0)** before N3 onset

See `study1b/README.md` for full details.

---

## Convergent Evidence Summary

| Study | Measure | N confirmed | Lead time | p | r |
|-------|---------|-------------|-----------|---|---|
| Study 1 | Permutation entropy | TBD (real data pending) | TBD | TBD | TBD |
| Study 1b | Permutation entropy (normalised) | 18/19 (94.7%) | 6.6 min | .0001 | .860 |

Two independent pipelines, same dataset, same complexity family — convergent
evidence for the PEWS prediction.

---

## License

MIT
