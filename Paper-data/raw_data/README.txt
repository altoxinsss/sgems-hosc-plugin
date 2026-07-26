Raw-data provenance

The source data are the Stanford VI-E porosity model distributed by the Stanford Center for Earth Resources Forecasting (SCRF):
https://github.com/SCRFpublic/Stanford-VI-E

The local source file used during processing was the Stanford VI-E file named porosity.dat.

The original fine grid contains 150 x 200 x 200 cells. A 2 x 2 x 2 linear block average produces the 75 x 100 x 100 coarse flow grid used in the paper. The public source repository is the authoritative location for the raw model; this Zenodo package stores the derived subsets and calculation outputs rather than duplicating the third-party raw file.

Processing and coordinate conventions are documented in processed_data/ and in the top-level README.txt. Coordinates are one-based, with x varying fastest, followed by y and z.
