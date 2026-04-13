"""
HRIT VPM Study — Step 1: Preprocessing
========================================
Loads Sleep-EDF .edf files and hypnogram annotations.
Extracts 30-second epochs with sleep stage labels.
Outputs a clean per-subject dataframe.

Sleep stages (Sleep-EDF annotation codes):
    W  = Waking (0)
    N1 = Stage 1 NREM (1)
    N2 = Stage 2 NREM (2)
    N3 = Stage 3 NREM / SWS (3)
    R  = REM (4)
    M  = Movement (5, excluded)

Channels used:
    EEG Fpz-Cz  — frontal, most available across subjects
    EEG Pz-Oz   — parietal

Requirements:
    pip install mne numpy scipy pandas
"""

import os
import glob
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

try:
    import mne
    MNE_AVAILABLE = True
except ImportError:
    MNE_AVAILABLE = False
    print("WARNING: MNE not found. Using fallback EDF reader.")
    print("Install with: pip install mne")

DATA_DIR   = 'data'
OUTPUT_DIR = 'preprocessed'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Epoch parameters
EPOCH_SEC    = 30      # standard PSG epoch
SFREQ_TARGET = 100     # downsample to 100 Hz
STAGE_MAP = {'Sleep stage W': 0, 'Sleep stage 1': 1,
             'Sleep stage 2': 2, 'Sleep stage 3': 3,
             'Sleep stage 4': 3,  # Stage 4 = SWS, merge with N3
             'Sleep stage R': 4,
             'Sleep stage ?': -1, 'Movement time': -1}

CHANNELS = ['EEG Fpz-Cz', 'EEG Pz-Oz']


def load_edf_mne(psg_path, hyp_path):
    """Load PSG and hypnogram using MNE."""
    # Load PSG
    raw = mne.io.read_raw_edf(psg_path, preload=True, verbose=False)
    raw.pick_channels([ch for ch in CHANNELS if ch in raw.ch_names])
    if raw.info['sfreq'] != SFREQ_TARGET:
        raw.resample(SFREQ_TARGET, verbose=False)

    # Load hypnogram
    annot = mne.read_annotations(hyp_path)
    raw.set_annotations(annot, verbose=False)

    # Extract events from annotations
    events, event_id = mne.events_from_annotations(
        raw, chunk_duration=EPOCH_SEC, verbose=False)

    # Map annotations to stage codes
    stage_labels = []
    for onset, dur, label in zip(annot.onset, annot.duration, annot.description):
        stage = STAGE_MAP.get(label, -1)
        n_epochs_in_annotation = max(1, int(round(dur / EPOCH_SEC)))
        stage_labels.extend([stage] * n_epochs_in_annotation)

    return raw, stage_labels


def load_edf_fallback(psg_path, hyp_path):
    """Minimal EDF reader without MNE (single-channel)."""
    # Read hypnogram as text annotations (Sleep-EDF format)
    stages = []
    try:
        with open(hyp_path, 'rb') as f:
            content = f.read().decode('utf-8', errors='ignore')
        for line in content.split('\n'):
            for stage_name, stage_code in STAGE_MAP.items():
                if stage_name in line:
                    stages.append(stage_code)
    except Exception:
        pass

    # Read EDF header for basic info
    # Return None for signal — user should install MNE for full analysis
    print(f"  Fallback: extracted {len(stages)} stage annotations")
    print(f"  Install MNE for full signal analysis: pip install mne")
    return None, stages


def process_subject(psg_path, hyp_path):
    """Process a single subject. Returns DataFrame of epochs."""
    subject_id = os.path.basename(psg_path).replace('E0-PSG.edf', '')

    if MNE_AVAILABLE:
        raw, stage_labels = load_edf_mne(psg_path, hyp_path)
        if raw is None:
            return None

        sfreq = raw.info['sfreq']
        epoch_samples = int(EPOCH_SEC * sfreq)
        data = raw.get_data()  # shape: (n_channels, n_samples)

        n_epochs = min(len(stage_labels),
                       data.shape[1] // epoch_samples)

        records = []
        for i in range(n_epochs):
            stage = stage_labels[i] if i < len(stage_labels) else -1
            if stage == -1:
                continue
            start = i * epoch_samples
            end   = start + epoch_samples
            epoch_data = data[:, start:end]

            records.append({
                'subject':    subject_id,
                'epoch_idx':  i,
                'stage':      stage,
                'time_sec':   i * EPOCH_SEC,
                'eeg_Fpz':    epoch_data[0] if epoch_data.shape[0] > 0 else None,
                'eeg_Pz':     epoch_data[1] if epoch_data.shape[0] > 1 else None,
                'sfreq':      sfreq,
            })

        return records
    else:
        _, stage_labels = load_edf_fallback(psg_path, hyp_path)
        records = []
        for i, stage in enumerate(stage_labels):
            if stage == -1:
                continue
            records.append({
                'subject':   subject_id,
                'epoch_idx': i,
                'stage':     stage,
                'time_sec':  i * EPOCH_SEC,
                'eeg_Fpz':   None,
                'eeg_Pz':    None,
                'sfreq':     None,
            })
        return records


# ── Main ─────────────────────────────────────────────────────────────────────

psg_files = sorted(glob.glob(os.path.join(DATA_DIR, '*E0-PSG.edf')))
print(f'Found {len(psg_files)} PSG files in {DATA_DIR}/')
print()

if len(psg_files) == 0:
    print("No EDF files found. Run 00_download_data.py first.")
    print(f"Expected files in: {os.path.abspath(DATA_DIR)}/")
    exit(1)

all_records   = []
subjects_done = []

for psg_path in psg_files:
    hyp_path = psg_path.replace('E0-PSG.edf', 'EC-Hypnogram.edf')
    if not os.path.exists(hyp_path):
        print(f'  Skipping {psg_path} — no hypnogram found')
        continue

    sid = os.path.basename(psg_path).replace('E0-PSG.edf', '')
    print(f'Processing {sid}...')

    try:
        records = process_subject(psg_path, hyp_path)
        if records:
            all_records.extend(records)
            subjects_done.append(sid)
            print(f'  ✓ {len(records)} epochs')
    except Exception as e:
        print(f'  ✗ ERROR: {e}')

# Save stage sequence (without raw signal data — too large)
meta_records = [{k: v for k, v in r.items()
                 if k not in ('eeg_Fpz', 'eeg_Pz')}
                for r in all_records]

meta_df = pd.DataFrame(meta_records)
meta_df.to_csv(os.path.join(OUTPUT_DIR, 'epoch_metadata.csv'), index=False)

print(f'\nDone: {len(subjects_done)} subjects, {len(all_records)} epochs')
print(f'Stage distribution:')
for stage, name in [(0,'Wake'), (1,'N1'), (2,'N2'), (3,'N3/SWS'), (4,'REM')]:
    n = sum(1 for r in all_records if r['stage'] == stage)
    print(f'  {name}: {n} epochs ({n * EPOCH_SEC / 3600:.1f} hours)')

print(f'\nOutput saved to: {OUTPUT_DIR}/epoch_metadata.csv')
print('Next step: run 02_compute_complexity.py')
