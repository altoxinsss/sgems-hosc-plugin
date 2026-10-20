# Quick start: install, calculate, inspect

## 1. Requirements

- Windows with PowerShell 5.1 or newer.
- A working **32-bit SGeMS 2.0** installation using **Qt 3.3.5 and VC7.1**.
  The application folder must contain `sgems.exe` and its runtime DLLs.
- Write access to the SGeMS `plugins/Geostat` folder and to the example-output
  directory.

Obtain `sgems_2.0_installer.exe` from the
[official SGeMS 2.0 Windows installer directory](https://sourceforge.net/projects/sgems/files/sgems/v2.0%20Windows%20installer/).
The installer validates the selected host before copying plugin files. These
prebuilt DLLs do not require Visual Studio, a compiler, CMake, or Python.

## 2. Install the plugin package

Download and extract the repository ZIP, or clone the repository. Close SGeMS,
then double-click `Install.cmd`. Enter the directory containing `sgems.exe`
when prompted. The directory is the SGeMS application root, not its
`plugins/Geostat` subdirectory.

Alternatively, open PowerShell in the repository directory:

```powershell
.\Install.cmd "C:\SGeMS"
```

For a preflight without installing:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\tools\install_plugin.ps1 -SGeMSRoot "C:\SGeMS" -CheckOnly
```

The matched files installed under `plugins/Geostat` are:

| DLL | Interface | Algorithm shown in SGeMS |
|---|---|---|
| `HOSC2008.dll` | `HOSC_TEMPLATE.ui` | `HOSC_TEMPLATE` |
| `HOSCWindow.dll` | `HOSC_WINDOW.ui` | `HOSC_WINDOW` |

Existing files are backed up before replacement. Keep the backup path printed
by the installer if a later rollback is needed. Open SGeMS after installation;
both algorithms should be available in the algorithm selector.

## 3. Run the supplied examples

From the repository directory:

```powershell
.\RunExamples.cmd "C:\SGeMS"
```

For an explicit output location:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\tools\run_examples.ps1 -SGeMSRoot "C:\SGeMS" -OutputDirectory "C:\HOSC_examples"
```

The runner uses a separate SGeMS test environment and the supplied synthetic
Cartesian grid. It calculates window and exact-template examples, exports
results, and checks them against reference values. The console prints the
output directory and validation result. Inspect the generated validation
manifest and result files before using a new installation with research data.
Each execution creates a separate `run_<timestamp>_<id>` directory beneath the
output directory:

| Output | Contents |
|---|---|
| `RESULT.json` | Overall validation result and execution metadata |
| `reported_issues.csv` | Expected field diagnostics from deliberately invalid parameter cases |
| `template_result_checks.csv` | Numerical and support checks for the exact template |
| `input.sgems` | Source Cartesian grid used by the examples |
| `valid_template.sgems` | Exported template result grid |
| `valid_window.sgems` | Exported moving-window result grid |
| `template_3d.sgems` | Result from the shipped exact-template parameter file |
| `window_uniform.sgems` | Result from the shipped uniform fourth-order parameter file |
| `window_coordinate_gaussian.sgems` | Result from the shipped spatially weighted sixth-order parameter file |
| `window_soft_high_porosity.sgems` | Result from the shipped value-weighted sixth-order parameter file |
| `gui_example_result_checks.csv` | Checks for all 54 output nodes from the four shipped parameter files |
| `executed_parameters/` | Exact copies of the four parameter files used in this run |
| `audit/<run>/window_results.csv` | Double-precision window results and validity/support diagnostics |

The runner also checks that invalid inputs produce specific diagnostics and
that corrected inputs succeed in the same process. Four recorded field
diagnostics are expected in `reported_issues.csv`; the overall success signal
is `"status": "passed"` in `RESULT.json` and `PASS` in the console.

The four shipped parameter files are also executed without changing their
parameter values. Their file hashes, selected algorithms, weight modes, and
maximum numerical discrepancies are recorded in `RESULT.json`.

Load the exported `.sgems` grids in SGeMS to inspect their properties. The
optional `-WorkDirectory` parameter selects a separate working directory for
the isolated host; use an ASCII-only writable path when specifying it.

The examples are small software checks. They do not launch the full reservoir
experiments. Existing interactive projects are not used as example inputs.

## 4. Repeat a moving-window calculation in the interface

Load `tests/fixtures/manual/selftest_3d_cartesian.sgems` through SGeMS's object
loading command. This native file defines a **Cartesian grid named `test`**,
with dimensions `5 × 5 × 3`, unit cell sizes, origin `(0,0,0)`, and one property,
`phi`. It avoids the grid-type choices required by a GSLIB import.

Parameter files in `examples/` can be loaded with **Load Parameters** in the
SGeMS algorithm panel after loading this grid:

| Parameter file | Algorithm and configuration |
|---|---|
| `window_uniform.par` | `HOSC_WINDOW`, fourth-order uniform calculation |
| `window_coordinate_gaussian.par` | `HOSC_WINDOW`, sixth-order coordinate-Gaussian calculation |
| `window_soft_high_porosity.par` | `HOSC_WINDOW`, sixth-order soft-high-porosity calculation |
| `template_3d.par` | `HOSC_TEMPLATE`, exact three-dimensional template |

Choose the corresponding algorithm before loading its parameters, confirm that
the source is `test / phi`, and use a distinct target property for each run.
If audit output is enabled, select an existing writable audit directory before
running. The settings below also allow the baseline to be entered manually.

Select `HOSC_WINDOW` and set:

| Page | Field | Value |
|---|---|---|
| General | Calculate Target | `Cumulant` |
| General | Calculate Algorithm | `Bell` |
| General | Source Grid / Source Property | `test / phi` |
| General | Target Property | `window_phi_r4` |
| General | Target Origin Coordinates Ox/Oy/Oz | `0 / 0 / 0` |
| Sliding window | Moment/Cumulant Order | `4` |
| Sliding window | Window Size X/Y/Z | `3 / 3 / 3` |
| Sliding window | Move Step X/Y/Z | `1 / 1 / 1` |
| Sliding window | Invert primary statistic | `false` |
| Sliding window | Percentile clipping | `false` |
| Weighted Moment | Enable normalized weighting | `Disable` |
| Weighted Moment | Weight Mode | `uniform` |
| Weighted Moment | Normalize weights (required) | `true` |
| Weighted Moment | Legacy centering (unsupported) | `false` |
| Weighted Moment | Center cumulant inputs for order >= 2 | `true` |
| Weighted Moment | Write standardized cumulant lambda | `true` |
| Weighted Moment | Lambda variance epsilon | `1e-12` |
| Weighted Moment | Write support and ESS diagnostics | `true` |
| Weighted Moment | Coordinate Gaussian sigma | `1` |
| Weighted Moment | Minimum weight | `0` |
| Weighted Moment | High-porosity threshold | `0.25` |
| Weighted Moment | Soft transition width | `0.02` |
| Weighted Moment | Soft weight strength | `0.8` |
| Weighted Moment | Coordinate center mode | `window_center` |
| Weighted Moment | Write double audit sidecar | `true` |
| Weighted Moment | Existing audit base directory | An existing writable directory, e.g. `C:/HOSC_examples/audit` |

Create the audit directory before running. Weight-shape parameters do not affect
this unweighted baseline. The output grid has dimensions `3 × 3 × 1`, with
nine window centers from `(1,1,1)` to `(3,3,1)` in X-fastest order.

| Quantity | Expected value at every center |
|---|---|
| Raw fourth cumulant | `-6.2306722005208369e-5` |
| Variance, `kappa2` | `0.009114583333333334` |
| Standardized fourth cumulant, `lambda` | `-0.75` |
| Support / ESS | `27 / 27` |
| Numeric, display, and lambda validity | `1` |

The main property is the raw cumulant. Additional properties use suffixes
`_kappa2`, `_lambda`, `_lambda_valid`, `_support`, `_ess`, `_numeric_valid`,
and `_display_valid`. The audit run directory contains `resolved_config.txt`
and `window_results.csv`; the CSV preserves double-precision results, whereas
SGeMS display properties use single precision.

Repeat with Algorithm=`AVL` and a new Target Property. The values, support, and
centers should agree. To test weighting, enable normalized weighting and select
one of the three weight modes. With `uniform`, the same baseline is recovered.

## 5. Repeat an exact-template calculation

Use the same `test / phi` grid. Select `HOSC_TEMPLATE` and set:

| Page/group | Field | Value |
|---|---|---|
| General | Calculate Target / Algorithm | `Cumulant / AVL` |
| General | Source Grid / Source Property | `test / phi` |
| General | Target Property | `template_phi_r4` |
| General | Target Origin Coordinates Ox/Oy/Oz | `0 / 0 / 0` |
| Template shape | Number of occurrences | `1,1,1` |
| Template shape | Number of lags | `3,3,3` |
| Template shape | Lag separation | `1,1,1` |
| Template shape | Lag tolerance, azimuth, azimuth tolerance, azimuth bandwidth, dip, dip tolerance, dip bandwidth | `0,0,0` in each of the seven fields |
| Target Crosssections | Direction X / Y / Z | `1 / 2 / 3` |
| Target Crosssections | Position | `1` |
| 3D Slice | Enable 3D Slice | `Disable` |
| Weighted Moment | Enable Weighted Cumulant | `Disable` |
| Exact Template | Sampling mode | `exhaustive_integer_offsets` |
| Exact Template | Integer direction vectors | `1,0,0;0,1,0;0,0,1` |

The origin contributes one occurrence, so `1,1,1` produces total cumulant
order four. The three integer vectors select X, Y, and Z. Each lag count of
three enumerates indices `0,1,2`, with a separation of one grid cell.

The result is a `3 × 3 × 3` **lag-combination grid**, not a map of the original
reservoir nodes. For lag indices `(i,j,k)`, support is
`(5-i)(5-j)(3-k)`: 75 at `(0,0,0)` and 9 at `(2,2,2)`. All 27 cells have
complete-template support. Repeat with Bell to compare evaluators.

Three-dimensional grids require the exact mode. The exact mode requires
weighted calculation to be disabled, including when the selected weight mode
is `uniform`. Angle-search controls do not define the exact template; the
integer vectors and lag steps do.

## 6. Use a new dataset

For a GSLIB property file, choose **Cartesian grid**, provide its dimensions,
cell sizes and origin, and set the appropriate no-data value. Rows must follow
X-fastest, then Y, then Z order. A property named `x`, `y`, or `z` is not a
substitute for grid geometry. Point-set inputs are unsupported.

The optional `selftest_2d.gslib` and `selftest_3d.gslib` fixtures contain six
properties for the more detailed [manual tests](MANUAL_TESTING.md). Import
their dimensions as `5 × 5 × 1` and `5 × 5 × 3`, respectively; use unit cell
sizes, origin zero, and no-data value `-999`.

Choose a window that fits inside the source grid. Two-dimensional window
calculations use size Z=`1` and move step Z=`0`; three-dimensional calculations
use a positive Z step. For templates, choose lag spans that leave complete
translations inside the source grid. A lag with no support is uninformed.

Save the SGeMS project and parameter settings after a successful calculation.
For quantitative window comparisons, retain the audit files as well as the
project. See [the numerical contract](NUMERICAL_CONTRACT.md) for weights,
missing data, standardized cumulants, and validity masks.
