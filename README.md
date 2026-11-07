Code and data availability
All data are publicly available through the openFDA device event API. Two scripts reproduce the analysis, both using the Python standard library only:
1.	maude_stent_dislodgement.py: the full extraction and analysis pipeline (download, deduplication, classification, variable counts, proportional reporting ratios, and dislodgement share by device-year).
2.	maude_stent_dislodgement_MH.py: the same pipeline, unchanged, with the stratified Mantel-Haenszel analysis appended. It also prints a reconciliation against the manuscript counts, caches the raw download, and writes the strata counts used in Supplementary Table S2.
Because openFDA accrues reports retrospectively, later snapshots yield slightly larger counts than those reported here
