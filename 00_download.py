"""
HRIT VPM Study 1 — Step 0: Download Sleep-EDF Data
====================================================
Downloads N=20 Night-1 cassette recordings from the Sleep-EDF Expanded
dataset (PhysioNet). Requires a free PhysioNet account.

  1. Create account at: https://physionet.org/register/
  2. Accept the data use agreement for sleep-edfx
  3. Run: python 00_download.py

Each subject = one PSG (.edf) + one hypnogram (.edf) ≈ 100–150 MB.
Total download: ~2.5 GB for N=20.

Pre-registration: https://doi.org/10.17605/OSF.IO/X5ZPU
Dataset:          https://physionet.org/content/sleep-edfx/1.0.0/
"""

import os
import sys
import shutil
import urllib.request

try:
    import wfdb
except ImportError:
    sys.exit("wfdb not installed. Run:  pip install wfdb mne numpy scipy pandas matplotlib")

DATA_DIR = 'data'
DATABASE = 'sleep-edfx'
N_TARGET = 20

# Study 1 subjects (pre-registered N=20 cassette study, Night 1)
STUDY1_SUBJECTS = [
    'SC4001', 'SC4002', 'SC4012', 'SC4031', 'SC4041', 'SC4042',
    'SC4051', 'SC4052', 'SC4061', 'SC4062', 'SC4071', 'SC4081',
    'SC4091', 'SC4092', 'SC4102', 'SC4111', 'SC4112', 'SC4121', 'SC4131',
    'SC4141',
]

os.makedirs(DATA_DIR, exist_ok=True)


def is_complete(sid):
    psg = os.path.join(DATA_DIR, f'{sid}E0-PSG.edf')
    hyp = os.path.join(DATA_DIR, f'{sid}EC-Hypnogram.edf')
    return os.path.exists(psg) and os.path.exists(hyp)


def download_subject(sid):
    for fname in [f'{sid}E0-PSG.edf', f'{sid}EC-Hypnogram.edf']:
        dest = os.path.join(DATA_DIR, fname)
        if os.path.exists(dest):
            continue
        wfdb.dl_files(DATABASE, DATA_DIR, [f'sleep-cassette/{fname}'])
        sub = os.path.join(DATA_DIR, 'sleep-cassette', fname)
        if os.path.exists(sub):
            shutil.move(sub, dest)
    return is_complete(sid)


print("Sleep-EDF Expanded — Study 1 Download")
print(f"Target: N={N_TARGET} subjects")
print()

already = [s for s in STUDY1_SUBJECTS if is_complete(s)]
print(f"Already downloaded: {len(already)}/{N_TARGET}")

for sid in STUDY1_SUBJECTS:
    if len([s for s in STUDY1_SUBJECTS if is_complete(s)]) >= N_TARGET:
        break
    if is_complete(sid):
        continue
    print(f"  {sid}: downloading...", end=' ', flush=True)
    try:
        ok = download_subject(sid)
        if ok:
            mb = sum(os.path.getsize(os.path.join(DATA_DIR, f'{sid}{s}')) / 1e6
                     for s in ['E0-PSG.edf', 'EC-Hypnogram.edf']
                     if os.path.exists(os.path.join(DATA_DIR, f'{sid}{s}')))
            print(f"✓  ({mb:.0f} MB)")
        else:
            print("✗  incomplete")
    except Exception as e:
        print(f"✗  {e}")

complete = [s for s in STUDY1_SUBJECTS if is_complete(s)]
print(f"\n{'='*45}")
print(f"Ready: {len(complete)}/{N_TARGET} subjects")
for s in complete:
    print(f"  {s}")

if len(complete) >= 10:
    print("\n✓  Run next:  python 01_preprocess.py")
else:
    print(f"\n⚠  Only {len(complete)} subjects. PhysioNet requires a credentialed account.")
    print("   See: https://physionet.org/register/")
