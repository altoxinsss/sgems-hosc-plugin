# Current paper experiments

[Data](Data) is the compact evidence package for Sections 4.2–4.5 and Supplementary Methods S5, Table S12 and Fig. S5, synchronized on 2026-11-25. It contains three input datasets, 138 SGeMS parameter files, 17 native result grids and 14 result/performance CSV files. It includes the 9,000-observation Dataset 3, Fig. 9 event diagnostic and Fig. 10 window sampling density. This guide is kept outside the data directory.

The main experiments use 98 parameter files and 15 grids. The supplementary comparison adds 40 order-one-to-four parameter files, two multi-attribute grids and two CSV tables. Fifth-order results are reused from the matching main experiments. `Parameters/Supp_weighted_attributes/Supp_samples.dat` is a copy of Dataset 3 with its archived weights appended, not a fourth dataset.

## Load the inputs and results

Use the 64-bit SGeMS host and preserve the supplied object names when loading inputs. Parameter files do not load their input objects automatically.

| Input | File under `Data/Inputs.prj` | Loading instructions |
|---|---|---|
| Dataset 1: 75 × 100 slice at model layer 62 | `Dataset1_z62.sgems` | Load the native SGeMS object directly. |
| Dataset 2: complete 75 × 100 × 12 reference, layers 61–72 | `Dataset2_reference_metres.dat` | Import as a GSLIB point set; columns 1–3 are X, Y and Z, and column 4 is porosity. |
| Dataset 3: 9,000 observations in 750 complete vertical columns | `Dataset3_samples_metres.dat` | Use the same point-set column mapping as Dataset 2. |

The point-set coordinates are in metres. Keep the filenames without extensions as the object names. For Dataset 1, the coordinates are model-grid indices. The `.prj` paths are project directories whose internal structure should be retained. Result objects can be loaded from [Data/Results.prj](Data/Results.prj) without rerunning the experiments.

For supplementary sample calculations, also import `Data/Parameters/Supp_weighted_attributes/Supp_samples.dat` as the object `Supp_samples`. Columns 1–3 are coordinates; import the remaining columns as properties, including porosity and the three observation-weight properties.

## Experiment parameters

Choose the corresponding plugin, use **Load Parameters**, and select the loaded source object and porosity property.

| Parameter directory under `Data/Parameters` | Files | Experiment |
|---|---:|---|
| `S42` | 2 | Bell/AVL numerical comparison for Fig. 4 |
| `S43` | 32 | Sixteen configurations, each evaluated with Bell and AVL |
| `S44` | 24 | Template moment/cumulant and moving-window local attributes |
| `S45_template` | 20 | Four template settings and five estimates per setting |
| `S45_window` | 20 | Four window settings and five estimates per setting |
| `Supp_weighted_attributes` | 40 | Supplementary order-one-to-four weighted attributes |

The fixed 500 cell-origin offsets are supplied in `Data/Parameters/origin_offsets_500.csv`.

## Use another run directory

The package matches the current manuscript data directory byte-for-byte, including its checksum list. Its parameters retain the recorded local output and cell-offset paths. From the repository root, generate path-adjusted copies for a new ASCII-only run directory:

```powershell
python tools/prepare_paper_parameters.py --output-root "D:\HOSC-paper-runs"
```

This creates a new run directory with path-adjusted parameters, the fixed cell offsets and the supplementary weight-bearing input. It changes only run-directory and origin-offset paths, and refuses to overwrite an existing directory. Scientific parameters, datasets and results are unchanged. Load the three inputs from `paper/Data/Inputs.prj` and, for the supplementary sample cases, the copied `Supp_samples.dat` before running the corresponding parameters.

Alternatively, change **Record directory** and the cell-origin offset file path in the plugin interface. The offset file requires an absolute path. Keep the scientific parameters unchanged when reproducing the reported experiments.

## Full-precision result tables

All files below are in [Data/CSV](Data/CSV). Each CSV starts with a header row; explicit configuration columns distinguish experiments within a file.

| File | Contents |
|---|---|
| `S42_S43_accuracy_pairs.csv` | Paired Bell/AVL values for Section 4.2 and the extended configurations |
| `S42_S43_accuracy_summary.csv` | Numerical agreement summaries |
| `Performance_runs.csv` | Individual timing and independent-process memory measurements |
| `Performance_summary.csv` | Timing means, sample standard deviations, Bell/AVL time ratios and available memory summaries |
| `S44_template.csv` | Fig. 6 template results and moment terms |
| `S44_windows.csv` | Six local attributes for window sizes 3, 5 and 9; `window_cells` identifies the side length |
| `Fig08_AB_values.csv` | All 50 porosity values from the two selected 5 × 5 windows |
| `S45_templates.csv` | Four template settings, with 100 positive horizontal lag pairs per setting, estimates and support diagnostics |
| `Dataset3_point_weights.csv` | Three observation-weight vectors for the 9,000 Dataset 3 source records |
| `S45_windows.csv` | Estimates and support diagnostics at all valid window centers for four settings |
| `Fig09_event_contributions.csv` | Lag-domain dense-region event contributions and error improvements |
| `Fig10_sampling_density.csv` | Window sampling density at all 21,000 centers of the displayed configuration |
| `Supp_weighted_attributes.csv` | Full-precision local attributes for the two supplementary window scales |
| `Supp_weighted_attributes_metrics.csv` | Overall and layer-specific errors, including comparisons on the centers shared by both scales |

In result-property names, R denotes the complete-model reference, U uniform weighting, C cell declustering, P polygonal declustering and K global kriging declustering. The `source_record_1based` column in the point-weight table identifies data rows after the header of `Dataset3_samples_metres.dat`, starting at one. Source-point weights are normalized over the input records; template and window weights are normalized over the support of each calculation.

The Fig. 9 diagnostic uses horizontal lag coordinates. Its `excess_contribution_pp` column stores 100(Fs − Fr), where Fs and Fr are the mean dense-region point fractions in sample and reference events. Fig. 10 sampling density is the sample count within a 41 × 41 × 3-cell window divided by 5,043 reference nodes and by the domain sampling fraction 0.1. Window grids retain metre coordinates; the manuscript converts the displayed horizontal axes to kilometres.

## Precision and performance records

Native SGeMS result properties are stored as float32. Use the full-precision CSV files for numerical agreement checks and error calculations. Supplementary fifth-order results reuse the matching main-experiment fields; `kappa5_n` and `kappa5_ESS` retain their support diagnostics. The supplementary comparison reports errors for all valid centers at each scale and for the shared 21,000-center set.

Runtime is in seconds and memory is in MiB. In `Performance_runs.csv`, `measurement` distinguishes timing from memory records; the same `repeat` number across these groups does not identify a single execution. Timing summaries use five `formal` runs per case and evaluator, excluding `warmup` and `extra` records. Memory summaries use the largest working-set increase among five independent-process formal runs, calculated as the run peak minus the pre-execution working set. Higher-order configurations without memory measurements have empty memory fields.

## Verify the package

`Data/SHA256SUMS.txt` lists the SHA-256 hash of each data file, excluding the checksum list itself. From the repository root, check the package and documentation links with:

```powershell
python tools/check_repository.py
```

The data directory contains inputs, parameters, result grids, CSV evidence and the checksum list. Processing scripts and historical experiments are not part of this package.

## Original data source

The source reservoir model is [Stanford VI-E](https://github.com/SCRFpublic/Stanford-VI-E), provided by the Stanford Center for Earth Resources Forecasting. Retain the source attribution and applicable upstream terms when reusing derived inputs.
