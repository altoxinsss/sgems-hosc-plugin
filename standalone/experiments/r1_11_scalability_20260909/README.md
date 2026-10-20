# R1 11 standalone scalability benchmark

This experiment exercises the current `WithoutSGeMS` command-line programs. It does not call SGeMS and it does not modify the manuscript.

The generated INI files explicitly use the SGeMS-compatible defaults documented in `SGEMS_COMPATIBILITY.md`:

- complete moving windows only;
- complete data support inside each window;
- anchor presence required;
- output at the window centre;
- source origin and spacing fixed explicitly;
- no output inversion or percentile clipping;
- the selected source property is `value`.

Every run creates a new timestamped directory under `output/r1_11_scalability_20260909`. The directory contains all generated configurations, raw stdout and stderr logs, per-run timing and peak-memory records, Bell-versus-AVL error metrics, source and executable SHA-256 hashes, and summary tables.

Build the current standalone code first, then run:

```powershell
python run_benchmarks.py --build-dir ../../build/Release --data /path/to/exhaustive_image_z62.txt --repeats 5 --warmups 1 --timeout 1800
```

Use `--smoke --repeats 1 --warmups 0` for a short r=4/r=8 Bell--AVL compatibility check before the full matrix.

Run that command from this experiment directory. Use `../../build` with a
single configuration generator. The required `--data` argument identifies the
original Stanford VI-E z=62 grid; `--paper` is optional and records only a hash.
The common [reproduction guide](../README.md) describes the overall setup.

The full matrix covers:

- fourth-order template workload cases with 10 by 10, 15 by 15, 20 by 20, and 25 by 25 output sections;
- fourth-order moving-window cases with 5 by 5, 9 by 9, and 13 by 13 windows;
- fixed-workload template cumulant orders 4 through 8;
- direction counts 2 through 4 at fixed order 6;
- four order-6 occurrence-vector arrangements;

The order-extension cases retain the Stanford VI-E z=62 slice, three lag directions,
`NLAG=(6,20,20)`, `LAG=(3,3,3)`, `AZM=(180,90,0)`, a fixed first-direction lag,
and a 20 by 20 output section. Their occurrence vectors are `(1,1,1)`, `(2,1,1)`,
`(2,2,1)`, `(2,2,2)`, and `(3,2,2)` for orders 4 through 8. The order-4 20 by 20
case is measured once and reused by both the workload and order-extension summaries.

An experiment is marked `PASS` only when Bell and AVL produce identical coordinate support, all compared
values are finite, and the two routes satisfy

```text
abs(AVL - Bell) <= 1e-12 + 1e-10 * abs(Bell)
```

for every output value. A failed numerical gate remains in the evidence directory and is not silently discarded.

After completion, `RUN_REPORT.md` summarizes the fourth-order workload/window series, the
fixed-workload r=4--8 series, the sixth-order structure cases, peak working-set measurements,
and the Bell--AVL numerical gate. The CSV files remain the authoritative raw evidence.
