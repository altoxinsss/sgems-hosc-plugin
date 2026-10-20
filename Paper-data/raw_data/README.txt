Stanford VI-E source-data provenance

The authoritative source is the Stanford Center for Earth Resources Forecasting (SCRF) repository:
https://github.com/SCRFpublic/Stanford-VI-E

The recorded preprocessing consumes Porosity/porosity.txt and Facies/facies.dat and checks the coarsened porosity against Reservoir Simulation (Eclipse)/poro.inc. The original grid has 150 x 200 x 200 cells. Contiguous 2 x 2 x 2 blocks yield a 75 x 100 x 100 coarse grid: porosity is averaged and channel facies are counted. Channel occupancy is channel_count / 8; the strict-majority channel label uses channel_count > 4; exact ties (channel_count = 4) are excluded by the label-valid mask.

The current formal run uses coarse-grid slabs z=15-26, z=45-56, and z=75-86. The separate Section 4.3 benchmark uses the recorded z=62 slice. The preprocessing manifest records the original input hashes, axis mapping (including the y/z reversals), coordinate fixtures, output hashes, and grid properties. Coordinates in the coarse grids are one-based, with x varying fastest.

See ../current_results/2026-10-12/sections_4_4_to_4_6/configs/preprocessing_config.json and ../current_results/2026-10-12/sections_4_4_to_4_6/preprocessing/preprocessing_manifest.json for the archived settings. The runnable preprocessing implementation is in standalone/experiments/stanford_vie_preprocessing/ at the repository root. Set source paths for the local Stanford VI-E files before regeneration.

Raw third-party files and the complete derived grids are not duplicated in this compact result package. Their source terms and attribution continue to apply.
