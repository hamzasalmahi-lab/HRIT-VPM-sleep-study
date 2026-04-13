# HRIT VPM Sleep EEG Study
## Empirical Programme Study 1

Pre-registered secondary analysis of the Sleep-EDF Expanded database
testing the Variance-Precedes-Mean (VPM) prediction of the Hierarchical
Recursive Integration Theory (HRIT; Almahi, 2026).

**Pre-registration:** doi.org/10.17605/OSF.IO/X5ZPU  
**Companion theory paper:** HRIT Version 2 — Zenodo preprint: https://doi.org/10.5281/zenodo.19490741  
**Author:** Hamza S. Almahi · hamza.s.almahi@gmail.com

---

## The prediction

HRIT Proof 08 derives: τ_escape ~ (z_c − z)^{-0.500}  
HRIT Verification 02 derives: Var(RSI) ~ (z_c − z)^{+0.546}

Near the consciousness stability boundary, EEG complexity variance
rises to a peak BEFORE the mean crosses the consciousness threshold.

## Result

- 15/19 subjects confirmed VPM (79%)
- Mean lag = +224s (variance peaks ~3.7 min before threshold)
- Wilcoxon W = 164, p = .003, r = .633

## Data

Sleep-EDF Expanded Database (open access):
https://physionet.org/content/sleep-edfx/1.0.0/

Raw EDF files are not included. Download them using:

	python 00_download.py

## How to reproduce

	pip install mne numpy scipy matplotlib pandas wfdb
	python 00_download.py
	python 01_preprocess.py
	python 02_compute_complexity.py
	python 03_vpm_analysis.py
	python 04_figures.py
	python 05_write_results.py

## File structure

	00_download.py                    data download
	01_preprocess.py                  EDF loading and sleep staging
	02_compute_complexity.py          aperiodic slope computation
	03_vpm_analysis.py                pre-registered VPM analysis
	04_figures.py                     publication figures
	05_write_results.py               auto-generate results text
	PREREGISTRATION.md                pre-registration document
	MANUSCRIPT_TEMPLATE.md            manuscript template
	results/vpm_subject_results.csv   per-subject results
	results/vpm_summary_stats.csv     summary statistics

## Citation

Almahi, H. S. (2026). Variance-Precedes-Mean in EEG Complexity at the
Sleep-Stage Consciousness Boundary: A Pre-registered Reanalysis of the
Sleep-EDF Database. bioRxiv. DOI to be added after submission.

## Pre-registration

OSF pre-registration confirmed prior to analysis:
doi.org/10.17605/OSF.IO/X5ZPU

All analysis decisions were committed before examining the data.
Any deviations from the pre-registered plan are reported explicitly
in the manuscript as post-hoc analyses.

## License

MIT
