# Current Supplementary Data S1

`Supplementary_Data_S1_machine_readable_20261015.zip` is the machine-readable
companion to the revised manuscript and supplementary material, synchronized
on 2 November 2026. The original filename is retained. `archive_manifest.json`
identifies this revision by byte count, SHA-256, and member count.

The archive contains:

- Full-precision records and selection rules for Supplementary Tables S1-S11,
  including the simple weighted baselines and their recorded intervals.
- Mappings for manuscript Figs. 1-9 and Supplementary Figs. S1-S16, distinguishing
  current Word display images from archived plotting versions.
- Seven matched-window pairs: 14 complete uniform 5 x 5 x 3 windows, 1,050
  coarse-cell values, and 8,400 original fine-cell porosity/facies records.
- The three verified search-support tables (197,280 facies-valid windows),
  selected-window verification, complete pair-screening reconstruction,
  Table S11 exports, and editable SVG/TIFF Figure S16.

The manuscript Fig. 8/Table 4 pair is `dataset=z45_56, role=strict`, with centres
A=(52,17,47) and B=(42,7,48) in one-based processed coarse-grid coordinates.
The archive's protocol and scripts give the matching, selection, coordinate,
facies, and numerical definitions. This is post-processing of existing model
outputs; the original experiment results remain in
[`current_results/2026-10-12/`](../current_results/2026-10-12/README.md).

## Verify and extract

From the repository root, Python 3.10 or later can verify the archive using only
the standard library:

```text
python Paper-data/tools/verify_supplementary_data.py
```

Then extract it to a separate work directory, for example in PowerShell:

```powershell
Expand-Archive -LiteralPath Paper-data/supplementary_data_s1/Supplementary_Data_S1_machine_readable_20261015.zip -DestinationPath reviewer-output/Data_S1
```

From `reviewer-output/Data_S1/matched_fourth_order_windows/`, run:

```text
python scripts/verify_selected_windows.py --output-dir ../matched_window_checks
python scripts/export_table_s11.py --check
python scripts/make_figure_s16.py --output-dir ../figure_s16_rebuilt
```

These three commands use only the Python standard library. To reconstruct the
complete candidate search, install the matched-window `requirements.txt`
(NumPy and SciPy), then run:

```text
python scripts/rebuild_screening.py --output-dir ../matched_window_screening --save-candidates
```

The extracted README, `matched_fourth_order_windows/PROTOCOL.md`, and scripts
describe all input and output options. The archive retains historical source
paths as provenance; executable inputs resolve relative to the extracted
package or are supplied on the command line. Stanford VI-E source data retain
the attribution and terms described in [the raw-data README](../raw_data/README.txt).
