# Receipt: Flag the sub-MS001 epoch with invalid values

**Claim**: one sub-MS001 epoch (BC_trials page 75, sample 296000, onset 148 s) has all 244000 samples above 1e6 µV (max 1.96e38), as in the source; every other epoch is within 1161 µV.

**PR**: nemarDatasets/nm000455#3, receipt on top of head c20aa748.

## Source
- NEMAR nm000455 `origin/main` 8d42762debf76d2b167d25ec01977ce75b1be7e2: `sub-MS001/ses-1/ieeg/sub-MS001_ses-1_task-auditorythreshold_ieeg.eeg` (+ `.vhdr`) and `sourcedata/icEEGandelectrodes/MS001/sorted_trials_sounds_restricted_CPonly_all_gray_matter.mat`; sha256 of both equals their git-annex keys (below). `events.tsv` read at the PR head.
- Upstream: Dryad doi:10.5061/dryad.h44j0zpww. Per `sourcedata/README_sourcedata.txt`, the sourcedata `.mat` is the Dryad file (SHA-256 verified against Dryad at download) with only the 116-byte MAT text header date changed.

## Method
- Reader: MNE 1.13.2 `read_raw_brainvision`, `get_data(units="uV")` (stored float32 x resolution 1 µV); source pages via h5py 3.16 (MAT v7.3, float64).
- Epochs: the 135 `events.tsv` rows, 4000 samples each from `sample`; each row's `source_variable` / `source_trial_index` gives its source page.
- Per epoch: samples with |value| > 1e6 µV (NaN not counted), channels with all samples above it, max |value| (NaN ignored); same count on the source page; float32(source page) == BIDS epoch, NaN included.

## Result
- sub-MS001_ses-1_task-auditorythreshold_ieeg.eeg: sha256 1ed022861655a6d89ebadf5a0b4e8785460d3013359ba2dc8444d6927c0fc0f3, = git-annex key: True
- sorted_trials_sounds_restricted_CPonly_all_gray_matter.mat: sha256 075cab214db6c5924ec38f1c1b963c3dbf7ce066411f21e587adafed574cf633, = git-annex key: True

| onset (s) | sample | source page | samples > 1e6 µV | channels with all samples > 1e6 µV | max abs (µV) | source page samples > 1e6 | status | status_description (PR head) |
|---:|---:|---|---:|---:|---:|---:|---|---|
| 148.0 | 296000 | BC_trials 75 | 244000 of 244000 | 61 of 61 | 1.96e+38 | 244000 | bad | invalid source values: \|x\| > 1e6 µV in all 61 channels (max 1.96e38) |

- Other epochs: 134 of 135, all status good: True; max |value| 1160.6 µV; source pages with |value| > 1e6: 0
- float32(source page) == BIDS epoch (NaN included): 135 of 135 epochs

## Coverage
1 of 1 sub-MS001 recording, 135 of 135 epochs, all 61 channels; every epoch compared with its source page.

## Reproduce
`uv run --no-project --python 3.12 --with mne==1.13.2 --with h5py python code/receipts/2026-10-22_invalid-epoch.py <nm000455 root with the two annexed files> <checkout of this branch>`
(Slurm job 13448588, 2026-10-22; root = NEMAR `origin/main` 8d42762.)
