"""
HRIT VPM Study — Step 2 (FIXED): Compute EEG Complexity Metrics
================================================================
Fixes: numpy 2.0 removed np.math — use stdlib math module instead.
Also fixes: legacy pick_channels() warning → use inst.pick().

Run:
    python 02_compute_complexity_fixed.py
"""

import os
import glob
import math          # ← stdlib math, not numpy.math
import numpy as np
import pandas as pd
from scipy import signal
from scipy.stats import linregress
import warnings
warnings.filterwarnings('ignore')

try:
    import mne
    MNE_AVAILABLE = True
except ImportError:
    MNE_AVAILABLE = False
    print("ERROR: MNE not installed. Run: pip install mne")

DATA_DIR    = 'data'
OUTPUT_DIR  = 'complexity'
os.makedirs(OUTPUT_DIR, exist_ok=True)

SFREQ_TARGET = 100
EPOCH_SEC    = 30
ALPHA_BAND   = (8, 13)
APERIODIC_RANGE = (2, 40)

STAGE_MAP = {
    'Sleep stage W': 0, 'Sleep stage 1': 1,
    'Sleep stage 2': 2, 'Sleep stage 3': 3,
    'Sleep stage 4': 3, 'Sleep stage R': 4,
}


def compute_aperiodic_slope(epoch_data, sfreq):
    """1/f exponent via log-log PSD regression (2–40 Hz, alpha-excluded)."""
    freqs, psd = signal.welch(epoch_data, fs=sfreq,
                               nperseg=int(sfreq * 4),
                               noverlap=int(sfreq * 2))
    mask       = (freqs >= APERIODIC_RANGE[0]) & (freqs <= APERIODIC_RANGE[1])
    alpha_mask = (freqs >= 9) & (freqs <= 12)
    mask       = mask & ~alpha_mask

    f_sel, psd_sel = freqs[mask], psd[mask]
    if len(f_sel) < 5:
        return np.nan

    slope, *_ = linregress(np.log10(f_sel), np.log10(psd_sel + 1e-30))
    return -slope   # higher = steeper = more complex


def compute_permutation_entropy(epoch_data, order=6, delay=1):
    """Normalised permutation entropy (Bandt & Pompe 2002)."""
    n = len(epoch_data)
    if n < order * delay:
        return np.nan

    from collections import Counter
    motifs = []
    for i in range(n - (order - 1) * delay):
        pattern = epoch_data[i:i + order * delay:delay]
        motifs.append(tuple(np.argsort(pattern)))

    counts = Counter(motifs)
    total  = len(motifs)
    probs  = np.array([v / total for v in counts.values()])
    pe     = -np.sum(probs * np.log2(probs + 1e-30))
    pe    /= np.log2(math.factorial(order))   # ← stdlib math, not np.math
    return pe


def compute_alpha_power(epoch_data, sfreq):
    """Alpha-band (8–13 Hz) mean power."""
    freqs, psd = signal.welch(epoch_data, fs=sfreq,
                               nperseg=int(sfreq * 4),
                               noverlap=int(sfreq * 2))
    mask = (freqs >= ALPHA_BAND[0]) & (freqs <= ALPHA_BAND[1])
    return np.mean(psd[mask]) if mask.sum() else np.nan


def process_subject(psg_path, hyp_path):
    """Return list of per-epoch complexity records for one subject."""
    if not MNE_AVAILABLE:
        return []

    try:
        raw = mne.io.read_raw_edf(psg_path, preload=True, verbose=False)

        # Use inst.pick() — avoids legacy warning
        available = [ch for ch in ['EEG Fpz-Cz', 'EEG Pz-Oz']
                     if ch in raw.ch_names]
        if not available:
            print(f"  No standard EEG channels found in {psg_path}")
            return []
        raw.pick(available)   # ← modern API, no warning

        if raw.info['sfreq'] != SFREQ_TARGET:
            raw.resample(SFREQ_TARGET, verbose=False)

        annot = mne.read_annotations(hyp_path)

        data          = raw.get_data()[0]   # first channel only
        sfreq         = raw.info['sfreq']
        epoch_samples = int(EPOCH_SEC * sfreq)

        # Build stage list from hypnogram
        stages = []
        for onset, dur, label in zip(annot.onset, annot.duration,
                                     annot.description):
            stage = STAGE_MAP.get(label, -1)
            n_ep  = max(1, int(round(dur / EPOCH_SEC)))
            stages.extend([stage] * n_ep)

        n_epochs = min(len(stages), len(data) // epoch_samples)
        records  = []

        for i in range(n_epochs):
            stage = stages[i] if i < len(stages) else -1
            if stage == -1:
                continue
            start      = i * epoch_samples
            epoch_data = data[start:start + epoch_samples]

            records.append({
                'epoch_idx':    i,
                'time_sec':     i * EPOCH_SEC,
                'stage':        stage,
                'aperiodic':    compute_aperiodic_slope(epoch_data, sfreq),
                'perm_entropy': compute_permutation_entropy(epoch_data),
                'alpha_power':  compute_alpha_power(epoch_data, sfreq),
            })

        return records

    except Exception as e:
        print(f"  ERROR: {e}")
        return []


# ── Main ─────────────────────────────────────────────────────────────────────

psg_files = sorted(glob.glob(os.path.join(DATA_DIR, '*E0-PSG.edf')))
print(f'Computing complexity for {len(psg_files)} subjects...\n')

if not psg_files:
    print("No EDF files found. Run download script first.")
    raise SystemExit(1)

all_dfs = []

for psg_path in psg_files:
    hyp_path = psg_path.replace('E0-PSG.edf', 'EC-Hypnogram.edf')
    if not os.path.exists(hyp_path):
        print(f'Skipping {psg_path} — no hypnogram')
        continue

    sid = os.path.basename(psg_path).replace('E0-PSG.edf', '')
    print(f'Processing {sid}...', end=' ', flush=True)

    records = process_subject(psg_path, hyp_path)

    if not records:
        print('✗ no records')
        continue

    df = pd.DataFrame(records)
    df['subject'] = sid

    # Quick sanity check: aperiodic should be highest in Wake
    by_stage = df.groupby('stage')['aperiodic'].mean()
    wake_ap  = by_stage.get(0, np.nan)
    n3_ap    = by_stage.get(3, np.nan)
    direction = '✓' if wake_ap > n3_ap else '⚠ INVERTED'

    print(f'✓ {len(records)} epochs | '
          f'Wake={wake_ap:.2f} N3={n3_ap:.2f} {direction}')

    out = os.path.join(OUTPUT_DIR, f'{sid}_complexity.csv')
    df.to_csv(out, index=False)
    all_dfs.append(df)

if all_dfs:
    combined = pd.concat(all_dfs, ignore_index=True)
    combined.to_csv(os.path.join(OUTPUT_DIR, 'all_complexity.csv'), index=False)

    print(f'\n{"="*55}')
    print(f'Done: {len(all_dfs)} subjects, {len(combined)} epochs')
    print(f'\nAperiodic slope by stage (should decrease Wake→N3):')
    for stage, name in [(0,'Wake'),(1,'N1'),(2,'N2'),(3,'N3'),(4,'REM')]:
        sub = combined[combined.stage == stage]['aperiodic'].dropna()
        if len(sub):
            print(f'  {name}: {sub.mean():.3f} ± {sub.std():.3f}')

    print(f'\nOutput: {OUTPUT_DIR}/all_complexity.csv')
    print('Next step: upload all_complexity.csv here '
          'and I will run the VPM analysis.')
else:
    print('\nNo subjects processed successfully.')
