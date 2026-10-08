"""d3-dryad-auditory-conscious (Dryad doi:10.5061/dryad.h44j0zpww, CC0-1.0) -> BIDS-iEEG DERIVATIVE dataset (Bruno 2026-10-08: epochs count).

Source: per subject folder <id>[_2]/sorted_trials_*.mat (MATLAB v7.3) with BC_trials (Perceived), WG_trials (Not Perceived),
blank_TN_trials (true-negative blanks): MATLAB channels x samples x trials float64, 2-s epochs (Dryad README) of grey-matter
electrodes; NaN = channel-epochs removed by the authors. L_/R_MontageMap.mat: [electrode #, 3 coordinates] per row.
Output per folder: sub-<id>/ses-<1|2>/ieeg/sub-<id>_ses-<k>_task-auditorythreshold_ieeg (BrainVision IEEE float32, µV):
all epochs concatenated in the order BC (MATLAB pages 1..n), WG, blank; each epoch = one events.tsv row (trial_type
perceived / not_perceived / blank) + a BrainVision 'New Segment' marker at every epoch boundary. RecordingType 'epoched'.
Float64 -> float32 cast (bitwise equal to float32 of the source, NaN kept); max |float32 - float64| reported. Originals in sourcedata/.
Usage: python f2_convert_aud.py <work_dir> <bids_root>   (env ONLY=comma list of folders for a smoke test; SKIP_SOURCEDATA=1)
"""
import glob
import json
import os
import re
import shutil
import sys

import h5py
import numpy as np
import scipy.io as sio

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b3w3_common import scrub_mat_header, wjson, wtsv, write_vhdr  # noqa: E402

W, OUT = sys.argv[1], sys.argv[2]
X = os.path.join(W, "extract")
DL = os.path.join(W, "sourcedata_dl")
SKIP_SRC = os.environ.get("SKIP_SOURCEDATA") == "1"
ONLY = [x for x in os.environ.get("ONLY", "").split(",") if x]
TASK = "auditorythreshold"
CONDS = [("BC_trials", "perceived"), ("WG_trials", "not_perceived"), ("blank_TN_trials", "blank")]
CH_COLS = ["name", "type", "units", "low_cutoff", "high_cutoff", "sampling_frequency", "reference", "group", "status",
           "status_description", "hemisphere", "description"]
SITE = {  # article section 2.4 (Christison-Lagay et al. 2025)
    "yale": dict(site="Yale New Haven Hospital", man="Natus", model="Natus Neurolink/Braintronics amplifier, NeuroWorks export", hp=0.1, lp=1330.0, stated_fs=4096),
    "yale1024": dict(site="Yale New Haven Hospital", man="Natus", model="Natus Neurolink/Braintronics amplifier, NeuroWorks export", hp=0.1, lp=400.0, stated_fs=1024),
    "northwell": dict(site="Northwell Health", man="Tucker-Davis Technologies", model="PZ5M amplifier", hp=0.4, lp=1500.0, stated_fs=3000),
    "sinai": dict(site="Mount Sinai Hospital", man="Natus", model="Natus NeuroWorks", hp=0.1, lp=670.0, stated_fs=2048),
}
rep = {"recordings": [], "notes": [], "sourcedata": [], "mat_strings": {}}


class _V5Var:
    """MAT v5 array presented like an h5py dataset of a v7.3 file (axes reversed: (trial, samp, ch))."""

    def __init__(self, a):
        self.a = np.asarray(a)
        self.shape = tuple(reversed(self.a.shape))

    def __getitem__(self, k):
        assert k == ()
        return self.a.T


class _V5File(dict):
    """Read-only stand-in for h5py.File for MAT v5 sorted_trials files (one release file, 390, is MAT v5)."""

    def __init__(self, p):
        super().__init__({k: _V5Var(v) for k, v in sio.loadmat(p).items() if not k.startswith("__")})

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def open_mat(p):
    if h5py.is_hdf5(p):
        return h5py.File(p, "r")
    rep["notes"].append(f"{os.path.relpath(p, X)}: MAT v5 file (read with scipy.io.loadmat)")
    return _V5File(p)


def mat_char_strings(p):
    out = []
    if not h5py.is_hdf5(p):
        for k, v in sio.loadmat(p).items():
            if k.startswith("__"):
                continue
            a = np.asarray(v)
            if a.dtype.kind in "US" or (a.dtype.kind == "O" and any(np.asarray(e).dtype.kind in "USO" for e in a.ravel())):
                out.append(f"{k}: {str(v)[:200]}")
        return out
    with h5py.File(p, "r") as f:
        def vis(n, o):
            if isinstance(o, h5py.Dataset):
                c = o.attrs.get("MATLAB_class", b"")
                c = c.decode() if isinstance(c, bytes) else str(c)
                if c == "char":
                    out.append("".join(chr(x) for x in np.asarray(o[()]).ravel(order="F")))
        f.visititems(vis)
    return out


def montage(p):
    if not os.path.exists(p):
        return None
    if not h5py.is_hdf5(p):  # MAT v5 (e.g. an empty map)
        m = {k: v for k, v in sio.loadmat(p).items() if not k.startswith("__")}
        rep["notes"].append(f"{os.path.relpath(p, X)}: MAT v5 file, variables {[(k, np.shape(v)) for k, v in m.items()]}")
        assert len(m) == 1, (p, list(m))
        a = np.asarray(next(iter(m.values())), dtype=float)
        return a if a.ndim == 2 and a.shape[1] == 4 and a.shape[0] > 0 else np.zeros((0, 4))
    with h5py.File(p, "r") as f:
        k = [x for x in f.keys() if not x.startswith("#")]
        assert len(k) == 1, (p, k)
        o = f[k[0]]
        if o.attrs.get("MATLAB_empty", 0) or o[()].ndim < 2:
            return np.zeros((0, 4))
        return np.asarray(o[()]).T  # MATLAB rows: [elec#, c1, c2, c3]


folders = sorted(glob.glob(os.path.join(X, "Part*of4icEEGandelectrodes", "*")), key=lambda d: os.path.basename(d))
folders = [d for d in folders if os.path.isdir(d)]
if ONLY:
    folders = [d for d in folders if os.path.basename(d) in ONLY]
parts = {}
for d in folders:
    fold = os.path.basename(d)
    m = re.fullmatch(r"([A-Za-z]*\d+)(?:_(\d))?", fold)
    assert m, fold
    sub, ses = m.group(1), (m.group(2) or "1")
    st = sorted(glob.glob(os.path.join(d, "sorted_trials*.mat")))
    assert len(st) == 1, (fold, st)
    src = st[0]
    S = f"sub-{sub}"
    D = os.path.join(OUT, S, f"ses-{ses}", "ieeg")
    os.makedirs(D, exist_ok=True)
    stem = f"{S}_ses-{ses}_task-{TASK}"
    rep["mat_strings"][f"{fold}/{os.path.basename(src)}"] = mat_char_strings(src)
    with open_mat(src) as f:
        shapes = {c: tuple(reversed(f[c].shape)) for c, _ in CONDS if c in f}  # MATLAB (ch, samp, trials)
        miss = [c for c, _ in CONDS if c not in f]
        if miss:
            rep["notes"].append(f"{fold}: variables missing {miss}")
        nch = {s[0] for s in shapes.values() if len(s) == 3}
        nsm = {s[1] for s in shapes.values() if len(s) == 3}
        assert len(nch) == 1 and len(nsm) == 1, (fold, shapes)
        nch, n = nch.pop(), nsm.pop()
        assert n % 2 == 0, (fold, n)
        fs = n / 2.0  # README: 'samples contain 2 seconds at the patient's sampling rate'
        # site from folder prefix and rate
        if sub.startswith("NS"):
            site = "northwell"
        elif sub.startswith("MS"):
            site = "sinai"
        else:
            site = "yale1024" if fs == 1024 else "yale"
        si = SITE[site]
        if abs(fs - si["stated_fs"]) > 0.5:
            rep["notes"].append(f"{fold}: samples/2 s = {fs} Hz differs from the article's rate for {si['site']} ({si['stated_fs']} Hz)")
        # montage -> channel names / coordinates
        L = montage(os.path.join(d, "L_MontageMap.mat"))
        R = montage(os.path.join(d, "R_MontageMap.mat"))
        coords = {}
        if L is not None or R is not None:
            for hemi, M in (("L", L), ("R", R)):
                if M is None:
                    continue
                for row in M:
                    e = int(round(row[0]))
                    assert e not in coords, (fold, e)
                    coords[e] = (hemi, row[1:4].tolist())
        if coords and sorted(coords) == list(range(1, nch + 1)):
            mon_ok = True
        else:
            mon_ok = False
            rep["notes"].append(f"{fold}: montage electrode numbers {sorted(coords)[:5]}..{sorted(coords)[-3:] if coords else []} "
                                f"(n={len(coords)}) do not cover channels 1..{nch}; coordinates kept only for listed numbers")
        names = [f"E{i:03d}" for i in range(1, nch + 1)]
        # write .eeg epoch by epoch; NaN count per (epoch, channel)
        eeg = os.path.join(D, stem + "_ieeg.eeg")
        ev, markers, nanall, k0 = [], [], np.zeros(nch, dtype=np.int64), 0
        maxerr, exact = 0.0, True
        maxrel, maxabs, nbig = 0.0, 0.0, 0  # relative float32 error; largest |value|; samples with |value| > 1e6 µV
        with open(eeg, "wb") as fh:
            for var, lab in CONDS:
                if var not in f or len(shapes[var]) != 3:
                    if var in f:
                        rep["notes"].append(f"{fold}: {var} shape {shapes[var]} (no epochs)")
                    continue
                full = np.asarray(f[var][()], dtype=np.float64)  # h5 (trial, samp, ch); chunks span all trials -> read whole variable
                for t in range(shapes[var][2]):
                    a = full[t]  # (samp, ch)
                    assert a.shape == (n, nch), (fold, var, t, a.shape)
                    a32 = a.astype(np.float32)
                    fin = np.isfinite(a)
                    if fin.any():
                        dif = np.abs(a32[fin].astype(np.float64) - a[fin])
                        av = np.abs(a[fin])
                        maxerr = max(maxerr, float(dif.max()))
                        maxabs = max(maxabs, float(av.max()))
                        nz = av > 0
                        if nz.any():
                            maxrel = max(maxrel, float((dif[nz] / av[nz]).max()))
                        nb_ = int((av > 1e6).sum())
                        nbig += nb_
                        if nb_:
                            rep["notes"].append(f"{fold}: {var} page {t + 1} ({lab}_{t + 1}) has {nb_} samples with |value| > 1e6 µV "
                                                 f"(max {float(av.max()):.3g}); kept as released")
                    a32.astype("<f4").tofile(fh)
                    chnan = np.isnan(a).all(axis=0)
                    part = np.isnan(a).any(axis=0) & ~chnan
                    nanall += chnan
                    rej = [names[i] for i in np.nonzero(chnan)[0]]
                    ev.append([round(k0 / fs, 6), 2.0, lab, k0, var, t + 1, len(rej), " ".join(rej) if rej else "n/a",
                               int(part.sum())])
                    if k0 > 0:
                        markers.append(("New Segment", "", k0, 1))
                    markers.append(("Comment", f"{lab}_{t + 1}", k0, n))
                    k0 += n
                del full
        ntot = k0
        # read back bit-exact
        mm = np.memmap(eeg, dtype="<f4", mode="r").reshape(-1, nch)
        assert mm.shape[0] == ntot
        k = 0
        for var, lab in CONDS:
            if var not in f or len(shapes[var]) != 3:
                continue
            full = np.asarray(f[var][()], dtype=np.float64)
            for t in range(shapes[var][2]):
                b = full[t].astype(np.float32)
                if not np.array_equal(mm[k:k + n], b, equal_nan=True):
                    exact = False
                k += n
            del full
        del mm
    write_vhdr(os.path.join(D, stem + "_ieeg"), nch, fs, names, ["µV"] * nch, markers=markers,
               comment=f"f2_convert_aud.py: {fold}/{os.path.basename(src)} BC/WG/blank_TN epochs concatenated, float64 -> float32")
    # MNE check
    import mne
    rr = mne.io.read_raw_brainvision(os.path.join(D, stem + "_ieeg.vhdr"), preload=False, verbose="error")
    A = rr.get_data(start=0, stop=min(ntot, 3 * n)) * 1e6
    B = np.fromfile(eeg, dtype="<f4", count=min(ntot, 3 * n) * nch).reshape(-1, nch).T.astype(np.float64)
    mne_ok = bool(rr.n_times == ntot and len(rr.ch_names) == nch and np.allclose(A, B, rtol=1e-6, atol=1e-9, equal_nan=True))
    nep = len(ev)
    # sidecars
    rows = []
    for i, nm in enumerate(names):
        hemi, xyz = coords.get(i + 1, ("n/a", None))
        allnan = nanall[i] == nep
        rows.append([nm, "EEG", "µV", si["hp"], si["lp"], fs, "referential: single reference electrode chosen by the clinical team (not identified in the release)",
                     "n/a", "bad" if allnan else "good",
                     "all epochs NaN in the release" if allnan else (f"NaN (removed by the authors) in {int(nanall[i])} of {nep} epochs" if nanall[i] else "n/a"),
                     hemi, f"grey-matter intracranial contact, electrode #{i + 1} of the release (grid, strip or depth not stated)"])
    wtsv(os.path.join(D, stem + "_channels.tsv"), CH_COLS, rows)
    wjson(os.path.join(D, stem + "_channels.json"), {
        "hemisphere": {"Description": "Hemisphere from the montage file the electrode is listed in (L_MontageMap / R_MontageMap); n/a if not listed"},
        "type": {"Description": "EEG = intracranial contact; the release does not say which contacts are subdural grid, strip or depth electrodes"},
        "low_cutoff": {"Description": "Hardware high-pass corner of the recording site (article section 2.4)", "Units": "Hz"},
        "high_cutoff": {"Description": "Hardware low-pass corner of the recording site (article section 2.4)", "Units": "Hz"}})
    wtsv(os.path.join(D, stem + "_events.tsv"),
         ["onset", "duration", "trial_type", "sample", "source_variable", "source_trial_index", "n_channels_rejected", "channels_rejected",
          "n_channels_partial_nan"], ev)
    wjson(os.path.join(D, stem + "_ieeg.json"), {
        "TaskName": TASK,
        "TaskDescription": "Auditory perceptual threshold task: a 75 ms target sound (waterdrop, whistle or 'laser') embedded in continuous white noise at the "
                           "participant's calibrated 50 % detection threshold (14 % blank trials), followed by two questions (heard? which sound?). Article section 2.2.",
        "Instructions": "Report whether a target sound was heard (yes/no buttons, counterbalanced), then identify it (guess if not heard).",
        "InstitutionName": si["site"],
        "Manufacturer": si["man"], "ManufacturersModelName": si["model"],
        "SamplingFrequency": fs, "PowerLineFrequency": 60,
        "HardwareFilters": {"Bandpass": {"Description": f"{si['hp']}-{si['lp']} Hz corner frequencies stated in the article for {si['site']}"
                                                         + (" (the 1024 Hz Yale patient)" if site == "yale1024" else "")}},
        "SoftwareFilters": {"AuthorsProcessing": {"Description": "Authors' epoching and artifact screening as released (sorted_trials file '" + os.path.basename(src)
                                                  + "'); see README. No filtering stated for the released epochs" + (
                                                      "; file name states 'no_TTL_artifact_subtract_30Hz_movmean' (authors' TTL artifact handling, not documented further)"
                                                      if "TTL" in src else "") + "."}},
        "iEEGReference": "Referential montage of the original acquisition: a single reference electrode selected by the clinical team (article section 2.4)",
        "iEEGGround": "n/a", "iEEGPlacementScheme": "Clinically determined subdural grids, strips and depth electrodes; only grey-matter contacts released",
        "ElectrodeManufacturer": "n/a",
        "ECOGChannelCount": 0, "SEEGChannelCount": 0, "EEGChannelCount": nch,
        "RecordingType": "epoched", "EpochLength": 2.0, "RecordingDuration": round(ntot / fs, 6),
        "SourceFile": f"{fold}/{os.path.basename(src)}",
    })
    # electrodes (coordinates exactly as released; frame: see coordsystem)
    if True:  # electrodes.tsv is required for iEEG; n/a rows when the release has no montage entry
        er = []
        for i, nm in enumerate(names):
            if i + 1 in coords:
                hemi, xyz = coords[i + 1]
                er.append([nm, *[round(float(v), 6) for v in xyz], "n/a", hemi])
            else:
                er.append([nm, "n/a", "n/a", "n/a", "n/a", "n/a"])
        wtsv(os.path.join(D, f"{S}_ses-{ses}_space-Other_electrodes.tsv"), ["name", "x", "y", "z", "size", "hemisphere"], er)
        wjson(os.path.join(D, f"{S}_ses-{ses}_space-Other_electrodes.json"), {
            "x": {"Description": "Montage map column 2 as released"}, "y": {"Description": "Montage map column 3 as released"},
            "z": {"Description": "Montage map column 4 as released"},
            "hemisphere": {"Description": "L_MontageMap (L) or R_MontageMap (R)"},
            "SpatialReference": "Authors' template image: Colin 22 brain template as used by BioImage Suite; the exact image and the "
                                "voxel-to-MNI transform were not released. Coordinates are copied unchanged from L_/R_MontageMap.mat"
                                + ("" if coords else "; this participant has no montage entries in the release, so all coordinates are n/a") + "."})
        wjson(os.path.join(D, f"{S}_ses-{ses}_space-Other_coordsystem.json"), {
            "iEEGCoordinateSystem": "Other", "iEEGCoordinateUnits": "n/a",
            "iEEGCoordinateSystemDescription": "Columns 2-4 of the authors' L_/R_MontageMap.mat, copied unchanged. The Dryad README calls them '3d MNI coordinates' and the "
                                               "article says BioImage Suite placed electrodes 'in MNI space (Colin 22 brain template)'. The values are all positive "
                                               "(about 10-200), i.e. image (voxel-index) coordinates in the template image rather than MNI millimetres, and the "
                                               "transform to MNI mm is not given. They are therefore NOT labelled MNI here; do not use them for atlas lookups "
                                               "without that transform."})
        if not mon_ok:
            rep["notes"].append(f"{fold}: electrodes.tsv has n/a rows for channels missing from the montage maps")
    nh = sum(1 for e in ev) * 2.0 / 3600
    parts.setdefault(sub, {"site": si["site"], "sessions": []})["sessions"].append(ses)
    rr_ = dict(folder=fold, sub=S, ses=ses, file=os.path.basename(src), n_channels=nch, sfreq=fs, n_epochs=nep,
               n_by_cond={lab: sum(1 for e in ev if e[2] == lab) for _, lab in CONDS}, hours=nh, exact_float32=exact, mne_ok=mne_ok,
               max_abs_err_vs_float64_uV=maxerr, max_rel_err_vs_float64=maxrel, max_abs_value_uV=maxabs, n_samples_abs_gt_1e6_uV=nbig,
               channels_all_nan=int((nanall == nep).sum()), nan_channel_epochs=int(nanall.sum()),
               montage_complete=mon_ok, roundtrip_ok=bool(exact and mne_ok))
    rep["recordings"].append(rr_)
    print("REC", json.dumps(rr_), flush=True)

# participants
wtsv(os.path.join(OUT, "participants.tsv"), ["participant_id", "age", "sex", "site", "sessions"],
     [[f"sub-{s}", "n/a", "n/a", v["site"], ",".join(sorted(v["sessions"]))] for s, v in sorted(parts.items())])
wjson(os.path.join(OUT, "participants.json"), {
    "age": {"Description": "Not released per participant (article: 38 enrolled, mean age 38.8, range 16-69 years)", "Units": "year"},
    "sex": {"Description": "Not released per participant (article: 16 of 38 enrolled female)"},
    "site": {"Description": "Recording site, from the subject code (numeric = Yale, NS = Northwell, MS = Mount Sinai) and the article"},
    "sessions": {"Description": "ses-2 = the release folder with suffix '_2' (repeat participant with a second implantation, article 2.1)"}})
wjson(os.path.join(OUT, f"task-{TASK}_events.json"), {
    "onset": {"Description": "Start of the 2-s epoch in the concatenated file", "Units": "s"},
    "duration": {"Description": "Epoch length (Dryad README: 2 s)", "Units": "s"},
    "trial_type": {"Description": "Trial category of the epoch", "Levels": {
        "perceived": "BC_trials: target presented, reported heard and correctly identified (article 2.3)",
        "not_perceived": "WG_trials: target presented, reported not heard and incorrectly identified",
        "blank": "blank_TN_trials: no target presented, reported not heard (true negative)"}},
    "sample": {"Description": "First sample of the epoch (0-based)"},
    "source_variable": {"Description": "MATLAB variable in the sorted_trials file"},
    "source_trial_index": {"Description": "Page (3rd dimension, 1-based) of the MATLAB array; trials are in order of appearance within their category (Dryad README)"},
    "n_channels_rejected": {"Description": "Channels whose samples are all NaN in this epoch (removed by the authors' artifact screening)"},
    "channels_rejected": {"Description": "Names of those channels"},
    "n_channels_partial_nan": {"Description": "Channels with some (not all) NaN samples in this epoch"}})

# sourcedata
S_ = os.path.join(OUT, "sourcedata")
os.makedirs(S_, exist_ok=True)
for fn in ("README.md", "ReadMe.txt"):
    shutil.copyfile(os.path.join(DL, fn), os.path.join(S_, fn))
if not SKIP_SRC:
    for d in folders:
        fold = os.path.basename(d)
        for p in sorted(glob.glob(os.path.join(d, "*"))):
            dst = os.path.join(S_, "icEEGandelectrodes", fold, os.path.basename(p))
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copyfile(p, dst)
            ch = None
            if p.endswith(".mat"):
                with open(dst, "r+b") as fh:
                    head = fh.read(128)
                    nh_, ch = scrub_mat_header(head)
                    fh.seek(0)
                    fh.write(nh_[:128])
            rep["sourcedata"].append({"file": os.path.relpath(dst, OUT), "header_date_day01": ch, "bytes": os.path.getsize(dst)})
        print("sourcedata", fold, flush=True)
# behaviour: file names and the subject_idnumber field carry participant initials -> renamed, field replaced by the study code
B = os.path.join(S_, "behavioraldata")
os.makedirs(B, exist_ok=True)
bmap = []
for p in sorted(glob.glob(os.path.join(X, "behavioraldata", "*.mat"))):
    bn = os.path.basename(p)
    mm_ = re.fullmatch(r"(NS\d+|\d+)([A-Za-z]*)(\d?)_(?:B\d+_)?(\d+|combined)\.mat", bn)
    assert mm_, bn
    code, sess2, tail = mm_.group(1), mm_.group(3), mm_.group(4)
    lab_ = f"{code}{'-' + sess2 if sess2 else ''}"
    new = f"{lab_}_{tail}.mat"
    m = sio.loadmat(p, squeeze_me=False)
    v = {k: val for k, val in m.items() if not k.startswith("__")}
    if "subject_idnumber" in v:
        v["subject_idnumber"] = np.array(lab_)
    dst = os.path.join(B, new)
    assert not os.path.exists(dst), (bn, new)
    sio.savemat(dst, v, do_compression=False, oned_as="column")
    with open(dst, "r+b") as fh:
        head = fh.read(128)
        nh_, _ = scrub_mat_header(head)
        fh.seek(0)
        fh.write(nh_[:128])
    bmap.append([new, ";".join(sorted(v.keys()))])
wtsv(os.path.join(B, "behavioraldata_files.tsv"), ["file", "variables"], bmap)
open(os.path.join(S_, "README_sourcedata.txt"), "w").write(
    "Original Dryad release doi:10.5061/dryad.h44j0zpww (CC0-1.0), downloaded 2026-10 via the Dryad API; every file SHA-256 verified against Dryad.\n"
    "- icEEGandelectrodes/<folder>/: the authors' MATLAB v7.3 files (sorted_trials_*.mat, L_/R_MontageMap.mat), identical to Dryad except the\n"
    "  116-byte MAT text header, whose file-creation date was set to day 01 (privacy rule). The four 'PartNof4' zips are merged.\n"
    "- behavioraldata/: the authors' behaviour files, re-saved (scipy.io.savemat, MAT v5) because the original file names and the\n"
    "  'subject_idnumber' field contain participant initials. New names: <study code>[-<digit that followed the initials>]_<suffix>.mat\n"
    "  (suffix as released: last session number or 'combined'); the subject_idnumber field now holds '<study code>[-<digit>]';\n"
    "  all other variables are unchanged (see behavioraldata_files.tsv).\n"
    "  Column meanings are in README.md (Dryad README). An optional infix of the form B<number> in one Northwell file name was dropped.\n"
    "- README.md / ReadMe.txt: the Dryad README files (unchanged).\n"
    "- parcellationmaps.zip (100 cortical parcellations used for statistics, no neural data) is not included; get it from Dryad.\n")

# code + metadata
os.makedirs(os.path.join(OUT, "code"), exist_ok=True)
here = os.path.dirname(os.path.abspath(__file__))
for x_ in ("f2_convert_aud.py", "b3w3_common.py"):
    shutil.copy(os.path.join(here, x_), os.path.join(OUT, "code", x_))
recs = rep["recordings"]
hrs = sum(r["hours"] for r in recs)
rt_all = all(r["roundtrip_ok"] for r in recs)
rep["summary"] = {"folders": len(recs), "participants": len(parts), "epochs": sum(r["n_epochs"] for r in recs), "hours": hrs, "rt_all": rt_all,
                  "max_abs_err_vs_float64_uV": max(r["max_abs_err_vs_float64_uV"] for r in recs),
                  "max_rel_err_vs_float64": max(r["max_rel_err_vs_float64"] for r in recs),
                  "big_value_recordings": {r["sub"] + "_ses-" + r["ses"]: [r["n_samples_abs_gt_1e6_uV"], r["max_abs_value_uV"]] for r in recs if r["n_samples_abs_gt_1e6_uV"]},
                  "no_montage": [r["sub"] + "_ses-" + r["ses"] for r in recs if not r["montage_complete"]], "rates": sorted({r["sfreq"] for r in recs})}
json.dump(rep, open(os.path.join(OUT, "code", "conversion_report.json"), "w"), indent=1, default=str)
AUTH = ["Kate L. Christison-Lagay", "Aya Khalaf", "Noah C. Freedman", "Christopher Micek", "Sharif I. Kronemer", "Mariana M. Gusso", "Lauren Kim",
        "Sarit Forman", "Julia Ding", "Mark Aksen", "Ahmad Abdel-Aty", "Hunki Kwon", "Noah Markowitz", "Erin Yeagle", "Elizabeth Espinal", "Jose Herrero",
        "Stephan Bickel", "James Young", "Ashesh Mehta", "Kun Wu", "Jason Gerrard", "Eyiyemisi Damisah", "Dennis Spencer", "Jiayin Qu", "Hal Blumenfeld"]
FUND = ["NIH R01 NS134655 (HB)", "NIH R01 DC007657 (Europe PMC grant list)", "Mark Loughridge and Michelle Williams Foundation (HB)",
        "Betsy and Jonathan Blattmachr Family (HB)", "CAPES 88887.147295/2017-00 and 88881.186875/2018-01 (MG)", "Fundação Araucária and CAPES 88887.185226/2018-00 (MG)"]
wjson(os.path.join(OUT, "dataset_description.json"), {
    "Name": "The neural activity of auditory conscious perception: intracranial EEG epochs (Christison-Lagay et al. 2025; derivative: authors' 2-s grey-matter epochs)",
    "BIDSVersion": "1.10.0", "DatasetType": "derivative", "License": "CC0-1.0", "Authors": AUTH,
    "Acknowledgements": "Data released on Dryad by Kate L. Christison-Lagay, Hal Blumenfeld et al. (Yale University). Repackaged to BIDS for NEMAR.",
    "HowToAcknowledge": "Cite Christison-Lagay et al. (2025) NeuroImage 308:121041, doi:10.1016/j.neuroimage.2025.121041, and the Dryad dataset doi:10.5061/dryad.h44j0zpww.",
    "Funding": FUND,
    "EthicsApprovals": ["Institutional review boards of Yale University, The Mount Sinai Hospital and Northwell Health (article section 2.1); written informed consent"],
    "ReferencesAndLinks": ["https://doi.org/10.5061/dryad.h44j0zpww", "https://doi.org/10.1016/j.neuroimage.2025.121041", "https://doi.org/10.1101/2023.01.12.523829"],
    "Keywords": ["conscious perception", "auditory perception", "perceptual threshold", "intracranial EEG", "epilepsy", "epochs", "derivative"],
    "SourceDatasets": [{"DOI": "doi:10.5061/dryad.h44j0zpww", "URL": "https://datadryad.org/dataset/doi:10.5061/dryad.h44j0zpww", "Version": "Dryad version 369244 (published 2025-06-09)"}],
    "GeneratedBy": [
        {"Name": "Authors' epoching and artifact screening (MATLAB 2019a)", "Description": "Continuous icEEG (not released) segmented around the target sound; epochs with "
         "clinical or subclinical epileptic events removed by visual inspection; four-stage artifact rejection (Welch PSD peak > 200 µV/Hz above 30 Hz; MSE < 200 µV^2; "
         "MSE to the mean trace > 3000 µV^2, at most 20 % per electrode; samples > 5 SD); trials sorted into Perceived / Not Perceived / blank; grey-matter electrodes only "
         "(article sections 2.3-2.5). Released as 2-s epochs (Dryad README); removed channel-epochs are NaN. Some file names state further steps, e.g. "
         "'sounds_restricted', 'CPonly', 'recalibrated', 'no_TTL_artifact_subtract_30Hz_movmean' (not documented further).",
         "CodeURL": "https://doi.org/10.5061/dryad.h44j0zpww"},
        {"Name": "f2_convert_aud.py", "Description": "Repackaging of BC_trials / WG_trials / blank_TN_trials (float64) to BrainVision IEEE float32 (bit-exact to float32 of "
         "the source, NaN kept), one file per participant and implantation with the epochs concatenated and marked in events.tsv; montage-map values as released "
         "(space-Other). Code in code/.", "Version": "2026-10-08"}],
})
s = rep["summary"]
big_txt = ("" if not s["big_value_recordings"] else
           " CAVEAT: some released epochs contain implausibly large values (|value| > 1e6 µV; recording: [number of samples, largest |value|]): "
           + "; ".join(f"{k}: [{v[0]}, {v[1]:.3g}]" for k, v in s["big_value_recordings"].items())
           + ". They are kept exactly as released (float32 rounding only); treat them as artifacts.")
notes_txt = "\n".join(f"  - {x}" for x in rep["notes"]) or "  - none"
open(os.path.join(OUT, "README"), "w").write(f"""# The neural activity of auditory conscious perception (Christison-Lagay et al. 2025): intracranial EEG epochs. DERIVATIVE

Intracranial EEG from adults with drug-resistant epilepsy (subdural grids, strips and depth electrodes; Yale New Haven Hospital,
Northwell Health, Mount Sinai) performing an auditory perceptual threshold task. A 75 ms target sound (waterdrop, whistle or
'laser') was embedded in white noise at each participant's 50 % detection threshold (14 % blank trials). Each trial asked whether a
sound was heard and which one. Perceived = heard and correctly identified; Not Perceived = not heard and wrongly identified.

## THIS IS A DERIVATIVE DATASET: the continuous recordings were not released
- The release holds only the authors' **epoched, artifact-screened** trials, sorted by category (sorted_trials_*.mat). The
  continuous icEEG recordings and the TTL trigger channels were not released by the authors.
- Each epoch is **2 s** long ("samples contain 2 seconds at the patient's sampling rate", Dryad README). The release does not state
  where the target sound falls within the 2-s epoch. The article describes 4-s epochs centred on the target with the outer 500 ms
  later pruned for the gamma analysis. Treat the epoch timing relative to the stimulus as unknown unless you verify it.
- **Grey-matter electrodes only** (article: 2,805 of 5,793 contacts), as chosen by the authors from co-registered MRI/CT.
- Authors' processing (article sections 2.3-2.5): trials sorted by behaviour; epochs with clinical or subclinical epileptic events removed
  by visual inspection; four-stage artifact rejection (Welch PSD peak > 200 µV/Hz above 30 Hz; mean-square < 200 µV²; MSE to the electrode's
  mean trace > 3000 µV², at most 20 % of epochs; any sample > 5 SD). Removed channel-epochs are **NaN** in the release and stay NaN
  here (events.tsv `channels_rejected`; channels.tsv `status_description`). Several source file names name further processing that is
  not documented further: `sounds_restricted`, `CPonly`, `recalibrated`, `redone`, `extended_norun6or8`, `_r_`, and for two
  participants (436, 451) **`no_TTL_artifact_subtract_30Hz_movmean`** (TTL artifact handling with a 30 Hz moving mean). The source
  file of each recording is listed in its `_ieeg.json` (`SourceFile`).
- Reference: the referential montage of the original acquisition (one reference electrode chosen by the clinical team; not identified).

## Contents
- {s['participants']} participants, {s['folders']} recordings (`ses-2` = second implantation, release folders with suffix `_2`),
  {s['epochs']} epochs, {s['hours']:.2f} h of epochs. Sampling rates {', '.join(f'{x:g}' for x in s['rates'])} Hz (= samples per epoch / 2 s).
- One file per participant and session: `task-{TASK}` with all epochs concatenated in the order Perceived (BC_trials), Not
  Perceived (WG_trials), blank (blank_TN_trials), each in the order of the MATLAB pages (order of appearance within the category).
  **The file is not a continuous recording**: every epoch is one `events.tsv` row (onset, 2 s, trial_type), and BrainVision
  'New Segment' markers mark every epoch boundary. The sidecar says `RecordingType: epoched`, `EpochLength: 2`.
- Single trials cannot be linked to the behaviour files from the release alone (artifact rejection and per-sound exclusions removed
  trials). The behaviour files are in sourcedata.
- Channels `E001`... = electrode # of the release (channel i of the arrays = electrode #i of the montage maps). Type `EEG` because
  the release does not say which contacts are grid, strip or depth. µV (unit implied by the article's µV thresholds; not stated in the files).
- Data format: BrainVision IEEE float32. The float64 source values are stored as float32 (largest relative difference
  {s['max_rel_err_vs_float64']:.3g}, i.e. float32 rounding).{big_txt} Every file was read back bit-identical to float32 of the source, NaN included, and read with MNE.
  The float64 originals are in sourcedata/.
- Electrodes: `space-Other` with columns 2-4 of the authors' L_/R_MontageMap.mat, unchanged. CAVEAT: the README calls them "3d MNI
  coordinates" and the article says BioImage Suite placed electrodes in "MNI space (Colin 22 brain template)". But the values are all
  positive (about 10-200), i.e. image coordinates in the template volume, not MNI millimetres. The transform is not given, so they are not
  labelled MNI here. `hemisphere` = which montage file lists the electrode. Recordings whose electrodes are not (all) listed in the
  montage maps have n/a coordinates for those electrodes: {', '.join(s['no_montage']) or 'none'}.
- Hardware filters and sites (article 2.4): Yale 4096 Hz, 0.1-1330 Hz (one patient 1024 Hz, 0.1-400 Hz); Northwell (NS*) TDT PZ5M 3 kHz,
  0.4-1500 Hz; Mount Sinai (MS*) 2048 Hz, 0.1-670 Hz.
- participants.tsv: age and sex n/a (article: 38 enrolled, 16 female, mean age 38.8, range 16-69; 31 analysed).
- Conversion notes:
{notes_txt}
- sourcedata/: all sorted_trials and montage files (MAT header dates set to day 01), behaviour files (renamed and the
  subject_idnumber field replaced because they contained initials), Dryad README files. Parcellation maps not included.

## Source
- Dryad: {', '.join(AUTH)}. "The neural activity of auditory conscious perception". doi:10.5061/dryad.h44j0zpww (CC0-1.0).
- Article: Christison-Lagay KL et al. (2025) NeuroImage 308:121041, doi:10.1016/j.neuroimage.2025.121041. Preprint doi:10.1101/2023.01.12.523829.
- Funding: {'; '.join(FUND)}.
- Ethics: IRBs of Yale University, The Mount Sinai Hospital and Northwell Health; written informed consent (article 2.1).

## Privacy
Subject labels are the authors' study codes. The MAT files hold no strings in the signal and montage files. The MAT header file-creation
dates were set to day 01. The behaviour file names and their subject_idnumber field contained participant initials. They were removed,
and the key status is recorded as retained.
""")
open(os.path.join(OUT, "CHANGES"), "w").write("1.0.0 2026-10-08\n  - Derivative repackaging of the sorted-trials epochs of Dryad doi:10.5061/dryad.h44j0zpww for NEMAR.\n")
open(os.path.join(OUT, ".bidsignore"), "w").write("code/\n")
print(json.dumps(s))
print("notes", rep["notes"])
print("mat_strings_nonempty", {k: v for k, v in rep["mat_strings"].items() if v})
print("CONVERT_DONE")
