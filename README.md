# HRIT Empirical Programme — Studies 1 and 1b
## Presence Early Warning Signal (PEWS) in Sleep EEG

Pre-registered empirical tests of the **Presence Early Warning Signal (PEWS)**,
formally derived from the saddle-node bifurcation of the HRIT v3 allostatic
control ODE (Almahi, 2026).

**Prediction:** Near the consciousness stability boundary, EEG complexity
**variance** rises to a peak **before** the rolling mean crosses the N3
threshold — a pre-cognitive warning of imminent consciousness-state transition.

**Theoretical grounding:**
```
Escape time:  τ  ~  (z_c − z)^{−0.500}   (Proof 08 — mathematically exact)
Variance:     Var ~  (z_c − z)^{+0.546}   (Verification 02 — computational)
```

**HRIT v3 preprint:** https://doi.org/10.5281/zenodo.19490741  
**Author:** Hamza S. Almahi · hamza.s.almahi@gmail.com

---

## Confirmed results

| | **Study 1** | **Study 1b** |
|---|---|---|
| **Pre-registration** | [OSF/X5ZPU](https://doi.org/10.17605/OSF.IO/X5ZPU) | [OSF/EYMSF](https://doi.org/10.17605/OSF.IO/EYMSF) |
| **Registered** | 13 Apr 2026 | prior to analysis |
| **Complexity proxy** | Aperiodic 1/f slope | Permutation entropy |
| **Direction** | Rises Wake→N3 | Falls Wake→N3 |
| **Dataset** | Sleep-EDF Expanded, N=19 | Sleep-EDF Expanded, N=19 |
| **PEWS confirmed** | **18/18 (100%)** | **19/19 (100%)** |
| **Mean lead time** | **3.6 min** | **3.3 min** |
| **Wilcoxon W** | 171 | 190 |
| **p (one-sided)** | < .001 | < .001 |
| **r** | .882 | .883 |
| **H1 (lag > 0)** | ✓ Confirmed | ✓ Confirmed |
| **H2 (r ≥ .633)** | N/A | ✓ Confirmed |
| **H3 (β = 0.546)** | ✗ Not confirmed | ✗ Not confirmed |

Both studies confirm H1. The PEWS signal replicates across two independent
complexity measures (aperiodic slope and permutation entropy) with opposite
directional polarity and independent analysis pipelines.

---

## Repository structure

```
HRIT-VPM-sleep-study/
│
├── Shared data pipeline (run once for both studies)
│   ├── 00_download.py             Download N=20 Sleep-EDF subjects
│   ├── 01_preprocess.py           Extract 30-s epochs and sleep stages
│   └── 02_compute_complexity.py   Compute aperiodic slope + PE per epoch
│
├── Study 1  —  Aperiodic slope proxy
│   ├── 03_vpm_analysis.py         Pre-registered VPM test  ← key file
│   ├── 04_figures.py              Three publication figures
│   ├── 05_write_results.py        Auto-generate results text
│   ├── PREREGISTRATION.md         Full OSF pre-registration document
│   └── results/
│       ├── vpm_subject_results.csv
│       ├── vpm_summary_stats.csv
│       ├── excluded_subjects.csv
│       └── manuscript_results_section.txt
│
└── study1b/  —  Permutation entropy proxy
    ├── 01_pews_analysis.py        Pre-registered PEWS test  ← key file
    ├── 02_figures.py              Three publication figures
    ├── 03_write_results.py        Auto-generate results text
    ├── README.md                  Study 1b details and comparison table
    └── results/
        ├── pews_subject_results.csv
        ├── pews_summary_stats.csv
        └── manuscript_results_section.txt
```

---

## How to reproduce

### Prerequisites

```bash
pip install mne numpy scipy pandas matplotlib wfdb
```

### Step 0 — Get the data

```bash
# Create a free account at physionet.org and accept the sleep-edfx agreement
python 00_download.py
```

This downloads the 20 Night-1 cassette recordings (~2.5 GB total).
Subjects: SC4001, SC4002, SC4012, SC4031, SC4041, SC4042, SC4051, SC4052,
SC4061, SC4062, SC4071, SC4081, SC4091, SC4092, SC4102, SC4111, SC4112,
SC4121, SC4131, SC4141.

### Step 1–2 — Preprocess and compute complexity (shared)

```bash
python 01_preprocess.py          # sleep staging
python 02_compute_complexity.py  # aperiodic slope + permutation entropy
```

Outputs one `complexity/<SID>_complexity.csv` per subject with columns:
`epoch_idx`, `time_sec`, `stage`, `aperiodic`, `perm_entropy`, `alpha_power`.

### Study 1 — Aperiodic slope

```bash
python 03_vpm_analysis.py   # pre-registered test
python 04_figures.py        # publication figures
python 05_write_results.py  # results text
```

### Study 1b — Permutation entropy

```bash
cd study1b
python 01_pews_analysis.py  # pre-registered test
python 02_figures.py        # publication figures
python 03_write_results.py  # results text
```

---

## Study 1 — Design

**Pre-registration:** https://doi.org/10.17605/OSF.IO/X5ZPU (13 Apr 2026)

### Pre-registered analysis (verbatim)

| Step | Action |
|---|---|
| 1 | Find first sustained Wake→N3 transition (≥3 consecutive N3 epochs) |
| 2 | Analysis window: −10 min to +5 min relative to onset |
| 3 | Rolling mean and variance in 3-min (6-epoch) centred window |
| 4 | N3 threshold = **25th percentile** of subject's own N3 aperiodic slopes |
| 5 | `t_peak_var` = epoch of max rolling variance in pre-window |
| | `t_cross_mean` = first epoch where mean **≥** threshold |
| 6 | Lag = `t_cross_mean − t_peak_var`  (positive = PEWS confirmed) |
| 7 | One-sided Wilcoxon vs zero, α = .05.  Report W, p, r = Z/√N |

**Exclusion criteria:** no sustained N3 transition; fewer than 5 pre-transition
epochs; fewer than 3 N3 epochs total.

**H2 (exploratory):** power-law β from log-log regression of rolling variance
vs time-to-threshold. Expected β ≈ 0.546.

### Direction note — aperiodic slope

The aperiodic slope **increases** from Wake to N3. During NREM sleep the 1/f
exponent steepens, producing higher stored values in N3 (N3 > Wake). The
pre-registered analysis correctly accounts for this:

- Threshold = **25th percentile** of N3 (lower bound of N3 territory)
- Crossing = `roll_mean >= threshold` (mean rises into N3)

An earlier version of the code incorrectly used the 75th percentile and `<=`
direction (calibrated for a decreasing measure). This made the test meaningless
for 16/19 subjects. The bug is documented in the commit history and corrected
in the current version.

---

## Study 1b — Design

**Pre-registration:** https://doi.org/10.17605/OSF.IO/EYMSF

| Parameter | Value |
|---|---|
| Complexity measure | Permutation entropy (order=6, delay=1) |
| Direction | Decreases Wake→N3 (Wake > N3) |
| N3 threshold | **75th percentile** of per-subject N3 PE epochs |
| Crossing | `roll_mean <= threshold` (PE falls into N3) |
| Rolling window | 6 epochs = 3 minutes (centred) |
| Pre-transition window | 20 epochs = 10 minutes |
| Test | One-sided Wilcoxon signed-rank, α = .05 |
| H2 | r ≥ .633 (replicates Study 1 effect size threshold) |
| H3 | β ≈ 0.546 (power-law exponent, exploratory) |

---

## Dataset

**Sleep-EDF Expanded Database** — Kemp et al. (2000)  
PhysioNet: https://physionet.org/content/sleep-edfx/1.0.0/  
License: Open Data Commons Attribution License v1.0  
Access: free account required at physionet.org

Raw EDF files are not included in this repo (each ~100 MB).
Run `python 00_download.py` after creating a PhysioNet account.

---

## Clinical significance

The PEWS window — between variance peak and threshold crossing — maps onto
a pre-ictal warning in the clinical translation. Study 4 of the HRIT
empirical programme (Yogarajah collaboration, UCL Queen Square) tests an
analogous signal — the Dissociation Early Warning Signal (DEWS) — in
functional neurological disorder (FND) seizure onset.

---

## License

MIT — see LICENSE
