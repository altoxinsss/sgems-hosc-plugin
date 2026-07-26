Paper-data

This package contains the processed datasets and numerical outputs supporting the paper "An integrated SGeMS implementation for high-order spatial cumulant calculation and local geological attribute mapping".

Directory contents

raw_data/
  Provenance and download information for the Stanford VI-E source model. The raw third-party dataset remains available from SCRF at https://github.com/SCRFpublic/Stanford-VI-E.

processed_data/
  stanford_vi_e_z62_2d.txt: complete 75 x 100 coarse-grid slice at z = 62 (7,500 cells).
  stanford_vi_e_z10_z21_3d_full.txt: complete 75 x 100 x 12 coarse-grid subvolume for z = 10-21 (90,000 cells).
  stanford_vi_e_z10_z21_3d_high_bands_ge025_cc300.txt: high-porosity-band subset retaining values >= 0.25 in 26-neighbor connected components with at least 300 nodes.

template_results/
  Numerical output and metadata for Fig. 4.10, together with supplementary three-dimensional template outputs. Fig. 4.10 is the 2D lag section. The F411 files are retained as supplementary calculations and are not the spatial map shown as Fig. 4.11 in the current paper. Exact settings are in parameter_description.txt and the .meta.txt files.

moving_window_results/
  Numerical output and metadata for Fig. 4.11 and the Section 4.5 A/B/C comparison, together with the associated figures and revision_metrics.csv. Exact settings are in parameter_description.txt.

figures/from_Paper.docx/
  All 18 image resources embedded in the current Paper.docx, converted to PNG where necessary. figure_index.txt records the confirmed Fig. 4.10 and Fig. 4.11 mappings.

Data format

Each .txt dataset or result uses the SGeMS text convention: object name, property count, one property name per line, then tab-separated rows. Coordinates are one-based; x varies fastest, followed by y and z. Porosity values are stored with six decimal places. The processed datasets contain no missing values.

Reproducibility

The calculation code is publicly available at https://github.com/altoxinsss/sgems-hosc-plugin. Result files are supplied so that the reported figures and numerical comparisons can be checked without rebuilding the software.

License and attribution

The Stanford VI-E source model remains subject to its original repository terms. Please cite the Stanford VI-E source and the paper when reusing these derived datasets or results.
