"""Epochs of sub-MS001 with |value| > 1e6 µV, in the BIDS file and the sourcedata .mat.

Usage: python 2026-10-22_invalid-epoch.py <nm000455 root> <PR #3 checkout>
"""

import csv
import hashlib
import subprocess
import sys
from pathlib import Path

import h5py
import mne
import numpy as np

STEM = "sub-MS001/ses-1/ieeg/sub-MS001_ses-1_task-auditorythreshold"
SRC = "sourcedata/icEEGandelectrodes/MS001/sorted_trials_sounds_restricted_CPonly_all_gray_matter.mat"
LIMIT = 1e6


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while chunk := f.read(1 << 24):
            h.update(chunk)
    return h.hexdigest()


def git(root: Path, *args: str) -> str:
    out = subprocess.run(["git", "-C", root, *args], capture_output=True, text=True)
    return out.stdout.strip() if out.returncode == 0 else "n/a"


def main(root: Path, pr: Path) -> None:
    print(f"mne {mne.__version__}, numpy {np.__version__}, h5py {h5py.__version__}")
    print(f"git HEAD of data root: {git(root, 'rev-parse', 'HEAD')}")
    print(f"git HEAD of PR checkout: {git(pr, 'rev-parse', 'HEAD')}")
    for rel in (f"{STEM}_ieeg.eeg", SRC):
        key = git(root, "cat-file", "-p", f"HEAD:{rel}").rsplit("--", 1)[-1]
        digest = sha256(root / rel)
        print(
            f"{Path(rel).name}: sha256 {digest}, = git-annex key: {key.startswith(digest)}"
        )

    raw = mne.io.read_raw_brainvision(
        root / f"{STEM}_ieeg.vhdr", preload=True, verbose="error"
    )
    # stored float32 value x resolution 1 µV
    x = raw.get_data(units="uV")
    with (pr / f"{STEM}_events.tsv").open(encoding="utf-8") as f:
        events = list(csv.DictReader(f, delimiter="\t"))
    n_ch, n_ep = len(raw.ch_names), len(events)
    length = raw.n_times // n_ep
    assert length * n_ep == raw.n_times, f"{raw.n_times} samples, {n_ep} epochs"
    print(
        f"\n{raw.n_times} samples = {n_ep} epochs x {length} samples, {n_ch} channels"
    )

    src = h5py.File(root / SRC, "r")
    rows = []
    for ev in events:
        a = int(ev["sample"])
        epoch = x[:, a : a + length]
        # h5py shape of a page: (samples, channels)
        page = src[ev["source_variable"]][int(ev["source_trial_index"]) - 1]
        finite = np.abs(np.nan_to_num(epoch, nan=0.0))
        rows.append(
            {
                **ev,
                "n_big": int((finite > LIMIT).sum()),
                "n_ch_big": int((finite > LIMIT).all(axis=1).sum()),
                "absmax": float(np.nanmax(np.abs(epoch))),
                "src_big": int((np.abs(np.nan_to_num(page)) > LIMIT).sum()),
                "equal_f32": np.array_equal(
                    np.float32(page.T), np.float32(epoch), equal_nan=True
                ),
            }
        )
    src.close()

    print(
        "\n| onset (s) | sample | source page | samples > 1e6 µV | channels with all"
        " samples > 1e6 µV | max abs (µV) | source page samples > 1e6 | status |"
        " status_description (PR head) |"
    )
    print("|---:|---:|---|---:|---:|---:|---:|---|---|")
    for r in rows:
        if r["n_big"] or r["status"] != "good":
            print(
                f"| {r['onset']} | {r['sample']} | {r['source_variable']} "
                f"{r['source_trial_index']} | {r['n_big']} of {n_ch * length} | "
                f"{r['n_ch_big']} of {n_ch} | {r['absmax']:.3g} | {r['src_big']} | "
                f"{r['status']} | {r['status_description'].replace('|', r'\|')} |"
            )
    rest = [r for r in rows if not r["n_big"]]
    print(
        f"\nOther epochs: {len(rest)} of {n_ep}, all status good:"
        f" {all(r['status'] == 'good' for r in rest)}; max |value|"
        f" {max(r['absmax'] for r in rest):.1f} µV;"
        f" source pages with |value| > 1e6: {sum(r['src_big'] > 0 for r in rest)}"
    )
    print(
        f"float32(source page) == BIDS epoch (NaN included): "
        f"{sum(r['equal_f32'] for r in rows)} of {n_ep} epochs"
    )


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
