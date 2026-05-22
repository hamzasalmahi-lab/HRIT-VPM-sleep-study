"""
HRIT VPM Study 1 — Step 1: Preprocess Sleep-EDF Recordings
============================================================
Loads raw EDF files and hypnogram annotations.
Extracts 30-second sleep epochs with stage labels.
Saves per-subject epoch metadata (no raw signals — too large for git).

Sleep stages used:
    0 = Wake      3 = N3 / SWS (Stages 3+4 merged)
    1 = N1        4 = REM
    2 = N2       -1 = Movement / unscored (excluded)

Channels:
    Primary:   EEG Fpz-Cz   (frontal — most complete across subjects)
    Secondary: EEG Pz-Oz    (parietal)

Output: preprocessed/epoch_metadata.csv
"""

import os
import glob
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')

try:
    import mne
except ImportError:
    raise SystemExit("MNE not installed. Run:  pip install mne")

DATA_DIR    = 'data'
OUTPUT_DIR  = 'preprocessed'
EPOCH_SEC   = 30
SFREQ       = 100   # downsample target (Hz)
CHANNELS    = ['EEG Fpz-Cz', 'EEG Pz-Oz']

STAGE_MAP = {
    'Sleep stage W': 0, 'Sleep stage 1': 1, 'Sleep stage 2': 2,
    'Sleep stage 3': 3, 'Sleep stage 4': 3,   # merge SWS stages
    'Sleep stage R': 4, 'Sleep stage ?': -1, 'Movement time': -1,
}

os.makedirs(OUTPUT_DIR, exist_ok=True)


def process_subject(psg_path, hyp_path):
    sid  = os.path.basename(psg_path).replace('E0-PSG.edf', '')
    raw  = mne.io.read_raw_edf(psg_path, preload=True, verbose=False)

    available = [ch for ch in CHANNELS if ch in raw.ch_names]
    if not available:
        print(f"  {sid}: no standard EEG channels — skipping")
        return []
    raw.pick(available, verbose=False)
    if raw.info['sfreq'] != SFREQ:
        raw.resample(SFREQ, verbose=False)

    annot = mne.read_annotations(hyp_path)
    stages = []
    for dur, label in zip(annot.duration, annot.description):
        code = STAGE_MAP.get(label, -1)
        stages.extend([code] * max(1, int(round(dur / EPOCH_SEC))))

    sfreq       = raw.info['sfreq']
    epoch_samp  = int(EPOCH_SEC * sfreq)
    n_epochs    = min(len(stages), raw.get_data().shape[1] // epoch_samp)
    records     = []

    for i in range(n_epochs):
        stage = stages[i] if i < len(stages) else -1
        if stage == -1:
            continue
        records.append({
            'subject': sid, 'epoch_idx': i,
            'time_sec': i * EPOCH_SEC, 'stage': stage,
        })
    return records


# ── Main ──────────────────────────────────────────────────────────────────────

psg_files = sorted(glob.glob(os.path.join(DATA_DIR, '*E0-PSG.edf')))
print(f"Preprocessing {len(psg_files)} subjects...\n")

if not psg_files:
    raise SystemExit(f"No EDF files in {DATA_DIR}/. Run 00_download.py first.")

all_records = []

for psg_path in psg_files:
    hyp_path = psg_path.replace('E0-PSG.edf', 'EC-Hypnogram.edf')
    if not os.path.exists(hyp_path):
        print(f"  SKIP {os.path.basename(psg_path)} — no hypnogram")
        continue

    sid = os.path.basename(psg_path).replace('E0-PSG.edf', '')
    print(f"  {sid}...", end=' ', flush=True)
    try:
        records = process_subject(psg_path, hyp_path)
        if records:
            all_records.extend(records)
            stage_counts = pd.Series([r['stage'] for r in records]).value_counts()
            print(f"✓  {len(records)} epochs  "
                  f"(W:{stage_counts.get(0,0)} N1:{stage_counts.get(1,0)} "
                  f"N2:{stage_counts.get(2,0)} N3:{stage_counts.get(3,0)} "
                  f"REM:{stage_counts.get(4,0)})")
    except Exception as e:
        print(f"✗  {e}")

meta = pd.DataFrame(all_records)
meta.to_csv(os.path.join(OUTPUT_DIR, 'epoch_metadata.csv'), index=False)

n_subs = meta['subject'].nunique()
print(f"\n✓  {n_subs} subjects, {len(meta)} epochs → preprocessed/epoch_metadata.csv")
print("   Run next:  python 02_compute_complexity.py")
