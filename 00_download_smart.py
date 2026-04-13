"""
HRIT VPM Study — Step 0 (DEFINITIVE): Smart Download
=====================================================
Fetches the RECORDS index from PhysioNet first to get the exact
list of valid subject IDs, then downloads only those files.

Also repairs the missing hypnograms for any PSG files already present.

Run:
    python 00_download_smart.py
"""

import os
import sys
import urllib.request
import shutil

try:
    import wfdb
except ImportError:
    print("ERROR: wfdb not installed. Run: pip install wfdb mne")
    sys.exit(1)

DATA_DIR   = 'data'
DATABASE   = 'sleep-edfx'
N_TARGET   = 20   # how many complete subjects we want

os.makedirs(DATA_DIR, exist_ok=True)

# ── Step 1: Fetch the RECORDS file to get exact valid IDs ─────────────────────
print("Step 1: Fetching valid subject list from PhysioNet RECORDS file...")

RECORDS_URL = 'https://physionet.org/files/sleep-edfx/1.0.0/RECORDS'

try:
    req     = urllib.request.urlopen(RECORDS_URL, timeout=20)
    content = req.read().decode('utf-8')

    # Extract cassette study Night-1 PSG records
    # Format: sleep-cassette/SC4001E0-PSG.edf
    valid_subjects = []
    for line in content.strip().split('\n'):
        line = line.strip()
        if 'sleep-cassette' in line and 'E0-PSG' in line:
            fname   = os.path.basename(line)           # SC4001E0-PSG.edf
            subject = fname.replace('E0-PSG.edf', '')  # SC4001
            valid_subjects.append(subject)

    print(f"  Found {len(valid_subjects)} valid Night-1 cassette subjects")
    print(f"  First 10: {valid_subjects[:10]}")
    print(f"  Last 10:  {valid_subjects[-10:]}")

except Exception as e:
    print(f"  Could not fetch RECORDS: {e}")
    print("  Using known-good subject list as fallback.")
    # Empirically verified IDs from Sleep-EDF Expanded cassette study
    valid_subjects = [
        'SC4001','SC4002','SC4011','SC4012','SC4021','SC4022',
        'SC4031','SC4032','SC4041','SC4042','SC4051','SC4052',
        'SC4061','SC4062','SC4071','SC4072','SC4081','SC4082',
        'SC4091','SC4092','SC4101','SC4102','SC4111','SC4112',
        'SC4121','SC4122','SC4131','SC4132','SC4141','SC4142',
        'SC4151','SC4152','SC4161','SC4162','SC4171','SC4172',
        'SC4181','SC4182','SC4191','SC4192',
    ]
    print(f"  Using fallback list ({len(valid_subjects)} subjects)")

# ── Step 2: Repair any existing PSG files missing hypnograms ──────────────────
print("\nStep 2: Repairing missing hypnograms for existing PSG files...")

existing_psg = [f for f in os.listdir(DATA_DIR) if f.endswith('E0-PSG.edf')]
repaired = 0

for psg_fname in existing_psg:
    sid      = psg_fname.replace('E0-PSG.edf', '')
    hyp_fname = f'{sid}EC-Hypnogram.edf'
    hyp_path  = os.path.join(DATA_DIR, hyp_fname)

    if os.path.exists(hyp_path):
        continue

    print(f"  {sid}: missing hypnogram — downloading...", end=' ', flush=True)
    try:
        wfdb.dl_files(DATABASE, DATA_DIR,
                      [f'sleep-cassette/{hyp_fname}'])

        # Move from subdirectory if wfdb put it there
        sub = os.path.join(DATA_DIR, 'sleep-cassette', hyp_fname)
        if os.path.exists(sub):
            shutil.move(sub, hyp_path)

        if os.path.exists(hyp_path):
            print('✓')
            repaired += 1
        else:
            print('✗ file not found after download')
    except Exception as e:
        print(f'✗ {e}')

if repaired:
    print(f"  Repaired {repaired} hypnogram(s)")

# ── Step 3: Download until we have N_TARGET complete subjects ─────────────────
print(f"\nStep 3: Downloading to reach {N_TARGET} complete subjects...")

def is_complete(sid):
    psg = os.path.join(DATA_DIR, f'{sid}E0-PSG.edf')
    hyp = os.path.join(DATA_DIR, f'{sid}EC-Hypnogram.edf')
    return os.path.exists(psg) and os.path.exists(hyp)

complete_now = [s for s in valid_subjects if is_complete(s)]
print(f"  Already complete: {len(complete_now)} subjects — {complete_now}")

downloaded = 0
errors     = 0

for sid in valid_subjects:
    if len(complete_now) + downloaded >= N_TARGET:
        break
    if is_complete(sid):
        continue

    psg_fname = f'{sid}E0-PSG.edf'
    hyp_fname = f'{sid}EC-Hypnogram.edf'
    psg_path  = os.path.join(DATA_DIR, psg_fname)
    hyp_path  = os.path.join(DATA_DIR, hyp_fname)

    print(f"  {sid}: downloading...", end=' ', flush=True)
    try:
        if not os.path.exists(psg_path):
            wfdb.dl_files(DATABASE, DATA_DIR,
                          [f'sleep-cassette/{psg_fname}'])
            sub = os.path.join(DATA_DIR, 'sleep-cassette', psg_fname)
            if os.path.exists(sub):
                shutil.move(sub, psg_path)

        if not os.path.exists(hyp_path):
            wfdb.dl_files(DATABASE, DATA_DIR,
                          [f'sleep-cassette/{hyp_fname}'])
            sub = os.path.join(DATA_DIR, 'sleep-cassette', hyp_fname)
            if os.path.exists(sub):
                shutil.move(sub, hyp_path)

        if is_complete(sid):
            size_mb = (os.path.getsize(psg_path) +
                       os.path.getsize(hyp_path)) / 1e6
            print(f'✓ ({size_mb:.0f} MB)')
            downloaded += 1
        else:
            print('✗ incomplete after download')
            errors += 1
    except Exception as e:
        print(f'✗ {e}')
        errors += 1

# ── Summary ───────────────────────────────────────────────────────────────────
print()
all_complete = [s for s in valid_subjects if is_complete(s)]
print(f"Summary: {len(all_complete)} complete subjects ready")
for sid in all_complete:
    psg_size = os.path.getsize(
        os.path.join(DATA_DIR, f'{sid}E0-PSG.edf')) / 1e6
    print(f"  {sid}: {psg_size:.0f} MB")

if len(all_complete) >= 4:
    print(f"\n✓ Ready to run: python 01_preprocess.py")
else:
    print(f"\n⚠ Only {len(all_complete)} subjects — need at least 4 to test pipeline")
    print("  Try the ISRUC-Sleep alternative (see 00b_fix_download.py)")
