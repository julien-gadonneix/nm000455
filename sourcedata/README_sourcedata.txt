Original Dryad release doi:10.5061/dryad.h44j0zpww (CC0-1.0), downloaded 2026-10 via the Dryad API; every file SHA-256 verified against Dryad.
- icEEGandelectrodes/<folder>/: the authors' MATLAB v7.3 files (sorted_trials_*.mat, L_/R_MontageMap.mat), identical to Dryad except the
  116-byte MAT text header, whose file-creation date was set to day 01 (privacy rule). The four 'PartNof4' zips are merged.
- behavioraldata/: the authors' behaviour files, re-saved (scipy.io.savemat, MAT v5) because the original file names and the
  'subject_idnumber' field contain participant initials. New names: <study code>[-<digit that followed the initials>]_<suffix>.mat
  (suffix as released: last session number or 'combined'); the subject_idnumber field now holds '<study code>[-<digit>]';
  all other variables are unchanged (see behavioraldata_files.tsv).
  Column meanings are in README.md (Dryad README). An optional infix of the form B<number> in one Northwell file name was dropped.
- README.md / ReadMe.txt: the Dryad README files (unchanged).
- parcellationmaps.zip (100 cortical parcellations used for statistics, no neural data) is not included; get it from Dryad.
