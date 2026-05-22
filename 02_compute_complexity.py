"""
HRIT VPM Study 1 — Step 2: Compute EEG Complexity
===================================================
For each 30-second epoch computes:

  aperiodic  — 1/f slope (−β from log-log PSD, Welch, 4-s windows,
                50% overlap, 2–40 Hz, alpha-excluded 9–12 Hz)
                PRE-REGISTERED primary proxy for Study 1.
                INCREASES Wake→N3 (steeper 1/f during NREM; N3 > Wake).

  perm_entropy — Permutation entropy (order=6, delay=1), normalised to
                 [0, log₂(6!)]. DECREASES Wake→N3 (Wake > N3).
                 Primary proxy for Study 1b.

  alpha_power  — Mean power in 8–13 Hz band.

Direction note: the aperiodic slope stored here is -slope from log-log
regression, so higher values = steeper 1/f = more N3-like. This is correct
and expected. The VPM analysis (03_vpm_analysis.py) is written for this
direction using a 25th-percentile N3 threshold and >= crossing.

Output: complexity/<SID>_complexity.csv  (one file per subject)
        complexity/all_complexity.csv    (combined)
"""

import os
import glob
import math
import warnings
import numpy as np
import pandas as pd
from scipy import signal
from scipy.stats import linregress
from collections import Counter

warnings.filterwarnings('ignore')

try:
    import mne
except ImportError:
    raise SystemExit("MNE not installed. Run:  pip install mne")

DATA_DIR    = 'data'
OUTPUT_DIR  = 'complexity'
EPOCH_SEC   = 30
SFREQ       = 100
CHANNELS    = ['EEG Fpz-Cz', 'EEG Pz-Oz']

STAGE_MAP = {
    'Sleep stage W': 0, 'Sleep stage 1': 1, 'Sleep stage 2': 2,
    'Sleep stage 3': 3, 'Sleep stage 4': 3,
    'Sleep stage R': 4, 'Sleep stage ?': -1, 'Movement time': -1,
}

os.makedirs(OUTPUT_DIR, exist_ok=True)


def aperiodic_slope(epoch, sfreq, freq_range=(2, 40), alpha_excl=(9, 12)):
    """
    Aperiodic 1/f exponent via log-log PSD regression.
    Returns −slope (positive; higher = steeper 1/f = more NREM-like).
    """
    freqs, psd = signal.welch(epoch, fs=sfreq,
                               nperseg=int(sfreq * 4),
                               noverlap=int(sfreq * 2))
    mask  = (freqs >= freq_range[0]) & (freqs <= freq_range[1])
    mask &= ~((freqs >= alpha_excl[0]) & (freqs <= alpha_excl[1]))
    if mask.sum() < 5:
        return np.nan
    slope, *_ = linregress(np.log10(freqs[mask]), np.log10(psd[mask] + 1e-30))
    return float(-slope)


def permutation_entropy(epoch, order=6, delay=1):
    """
    Normalised permutation entropy (Bandt & Pompe, 2002).
    Returns value in [0, 1]; higher = more random = more Wake-like.
    """
    n = len(epoch)
    if n < order * delay:
        return np.nan
    motifs = [tuple(np.argsort(epoch[i:i + order * delay:delay]))
              for i in range(n - (order - 1) * delay)]
    counts = Counter(motifs)
    total  = len(motifs)
    probs  = np.array([v / total for v in counts.values()])
    pe     = -np.sum(probs * np.log2(probs + 1e-30))
    return float(pe / np.log2(math.factorial(order)))


def alpha_power(epoch, sfreq, band=(8, 13)):
    freqs, psd = signal.welch(epoch, fs=sfreq,
                               nperseg=int(sfreq * 4),
                               noverlap=int(sfreq * 2))
    mask = (freqs >= band[0]) & (freqs <= band[1])
    return float(np.mean(psd[mask])) if mask.sum() else np.nan


def process_subject(psg_path, hyp_path):
    raw = mne.io.read_raw_edf(psg_path, preload=True, verbose=False)
    available = [ch for ch in CHANNELS if ch in raw.ch_names]
    if not available:
        return []
    raw.pick(available, verbose=False)
    if raw.info['sfreq'] != SFREQ:
        raw.resample(SFREQ, verbose=False)

    annot  = mne.read_annotations(hyp_path)
    stages = []
    for dur, label in zip(annot.duration, annot.description):
        code = STAGE_MAP.get(label, -1)
        stages.extend([code] * max(1, int(round(dur / EPOCH_SEC))))

    sfreq_used = raw.info['sfreq']
    epoch_samp = int(EPOCH_SEC * sfreq_used)
    data       = raw.get_data()[0]   # first channel (Fpz-Cz preferred)
    n_epochs   = min(len(stages), len(data) // epoch_samp)
    records    = []

    for i in range(n_epochs):
        stage = stages[i] if i < len(stages) else -1
        if stage == -1:
            continue
        ep = data[i * epoch_samp:(i + 1) * epoch_samp]
        records.append({
            'epoch_idx':    i,
            'time_sec':     i * EPOCH_SEC,
            'stage':        stage,
            'aperiodic':    aperiodic_slope(ep, sfreq_used),
            'perm_entropy': permutation_entropy(ep),
            'alpha_power':  alpha_power(ep, sfreq_used),
        })
    return records


# ── Main ──────────────────────────────────────────────────────────────────────

psg_files = sorted(glob.glob(os.path.join(DATA_DIR, '*E0-PSG.edf')))
print(f"Computing EEG complexity for {len(psg_files)} subjects...\n")
print("Measures: aperiodic slope (Study 1 proxy) + permutation entropy (Study 1b)")
print()

if not psg_files:
    raise SystemExit(f"No EDF files in {DATA_DIR}/. Run 00_download.py first.")

all_dfs = []

for psg_path in psg_files:
    hyp_path = psg_path.replace('E0-PSG.edf', 'EC-Hypnogram.edf')
    if not os.path.exists(hyp_path):
        continue

    sid = os.path.basename(psg_path).replace('E0-PSG.edf', '')
    print(f"  {sid}...", end=' ', flush=True)

    try:
        records = process_subject(psg_path, hyp_path)
        if not records:
            print("✗  no valid epochs")
            continue

        df      = pd.DataFrame(records)
        df['subject'] = sid

        wake_ap = df[df.stage == 0]['aperiodic'].mean()
        n3_ap   = df[df.stage == 3]['aperiodic'].mean()
        wake_pe = df[df.stage == 0]['perm_entropy'].mean()
        n3_pe   = df[df.stage == 3]['perm_entropy'].mean()

        # N3 > Wake for aperiodic (expected); Wake > N3 for PE (expected)
        ap_ok = '✓' if n3_ap > wake_ap else '⚠'
        pe_ok = '✓' if wake_pe > n3_pe else '⚠'

        print(f"✓  {len(records)} epochs | "
              f"Aperiodic Wake={wake_ap:.2f} N3={n3_ap:.2f} {ap_ok} | "
              f"PE Wake={wake_pe:.3f} N3={n3_pe:.3f} {pe_ok}")

        df.to_csv(os.path.join(OUTPUT_DIR, f'{sid}_complexity.csv'), index=False)
        all_dfs.append(df)

    except Exception as e:
        print(f"✗  {e}")

if all_dfs:
    combined = pd.concat(all_dfs, ignore_index=True)
    combined.to_csv(os.path.join(OUTPUT_DIR, 'all_complexity.csv'), index=False)
    print(f"\n✓  {len(all_dfs)} subjects → complexity/<SID>_complexity.csv")
    print("   Run next:  python 03_vpm_analysis.py")
else:
    raise SystemExit("No subjects processed. Check EDF files.")
