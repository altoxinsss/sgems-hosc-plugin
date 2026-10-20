# Troubleshooting

Start with the installer preflight and the supplied examples. These separate
host compatibility from dataset or parameter problems.

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\tools\install_plugin.ps1 -SGeMSRoot "C:\SGeMS" -CheckOnly
.\RunExamples.cmd "C:\SGeMS"
```

## Installation and startup

| Symptom | Check and action |
|---|---|
| `sgems.exe` not found | Pass the folder containing the executable, not `plugins/Geostat`. Extract a complete SGeMS distribution before installing. |
| Host compatibility check fails | Use the Windows x86 SGeMS 2.0 distribution with Qt 3.3.5 and VC7.1 runtime dependencies. A 64-bit Windows operating system is acceptable; a 64-bit SGeMS executable is a different binary target. |
| Missing `MSVCR71.dll`, `MSVCP71.dll`, or a SGeMS/Qt DLL | Restore the dependencies from the matching SGeMS distribution. Keep those DLLs with their host application. |
| Access denied or a plugin cannot be replaced | Save and close all SGeMS sessions. Verify write access to its installation folder; an installation under a protected system folder may require administrator privileges. |
| PowerShell blocks a downloaded script | Run the supplied `.cmd` launcher. It invokes its local script with a process-scoped execution-policy setting. Check any organizational policy or file blocking if execution is still denied. |
| Algorithm missing after installation | Restart SGeMS and check that each DLL and matching `.ui` is in `plugins/Geostat`. Rerun preflight to verify the installed host and source package. |
| New DLL but old interface, or vice versa | Reinstall the complete package; the DLLs and interfaces are a matched set. |
| An example path is not usable by the old host | Give `tools/run_examples.ps1` an existing ASCII-only `-WorkDirectory`, such as `C:\HOSC_work`. Repository and output directories may contain Unicode; the isolated host command file uses the ASCII working directory. |
| Example log includes `'import site' failed` | This message can come from SGeMS's embedded Python startup in the isolated copy. These HOSC examples use native commands. Assess the run using its exit status, complete numerical checks, and `RESULT.json`; a missing or failed result is not a pass. |

The plugin package contains HOSC binaries and interfaces. SGeMS, Qt, and the
Microsoft runtime are supplied separately by their respective distributions.
Building with a newer Visual C++ version does not make a DLL interchangeable
with this host's VC7.1 binary interface.

To restore replaced files, close SGeMS and pass the backup manifest printed by
the installer:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\tools\uninstall_plugin.ps1 -SGeMSRoot "C:\SGeMS" -BackupManifest "C:\SGeMS\plugin_backups\<installation-backup>\backup_manifest.json"
```

Replace `<installation-backup>` with the actual backup directory name. Keep
SGeMS closed until the restore completes.

## Parameter validation

The HOSC interface displays a validation summary, marks invalid fields, and
selects the relevant page. SGeMS may also display its standard parameter-error
dialog; close that dialog to read the detailed HOSC summary. Correct all
reported fields, then run again.

| Message or symptom | Meaning and correction |
|---|---|
| Source must be a Cartesian grid | Reimport the dataset as a regular Cartesian grid with explicit dimensions and geometry. Point sets are unsupported. |
| `legacy_search` is not supported for 3D | Select `exhaustive_integer_offsets` for a three-dimensional source. |
| Exact mode supports unweighted Cumulant only | Select `Cumulant` and set **Enable Weighted Cumulant = Disable**. A selected `uniform` mode does not override an enabled weighting switch. |
| Integer direction vectors invalid | Supply one nonzero integer vector for each direction. A 3D example is `1,0,0;0,1,0;0,0,1`; a 2D vector has two components. |
| Template lists have inconsistent lengths | Occurrences, lag counts, lag steps, and all angle/tolerance lists need the same direction count. |
| Template order exceeds eight | Total order equals one origin occurrence plus the sum of the specified direction occurrences. Reduce the occurrences. |
| Window step or size invalid | A window must fit the source. X/Y steps must be positive; use Z step zero for a 2D source and a positive Z step for a 3D source. |
| Weight normalization required | Set `Normalize weights = true` for enabled weighting. |
| Coordinate center mode invalid | Use `window_center` for coordinate-Gaussian window weighting. |
| Audit directory invalid | Create an existing writable directory, then select it, or disable the audit sidecar. |
| Moment with standardized cumulant enabled | Disable standardized cumulant output when calculating moments. |

Exact templates require several compatible settings together. A complete
working set is provided in [the quick start](QUICKSTART.md#5-repeat-an-exact-template-calculation).

## Results and display

| Observation | Interpretation and check |
|---|---|
| Window output is smaller than the source | Only full windows are evaluated. Its nodes represent window centers, with spacing determined by the move step. |
| Template output looks unlike the source map | Its axes represent lag combinations. The support property counts complete translations for each combination. |
| A result cell is blank | Inspect support and validity. Incomplete windows, unsupported lags, or non-finite results are uninformed rather than valid zeros. |
| Lambda is blank but the raw cumulant exists | Variance may be at or below the configured epsilon. Inspect `_kappa2` and `_lambda_valid`; standardized cumulants are not defined for zero variance. |
| A cumulant is negative | Higher-order cumulants can be negative. Check its order, support, and expected value rather than treating its sign as a failure. |
| Most colors look identical | Check the numerical range and the selected color limits. The synthetic linear fixture deliberately produces the same fourth-order window cumulant at all nine centers. |
| Display values differ slightly from the audit CSV | SGeMS stores display properties as floats; the audit preserves doubles. Use the numerical tolerances in [manual testing](MANUAL_TESTING.md). |
| Clipping appears not to change the main property | Clipping is written to a separate `<output>_clipped` property. The raw primary statistic is retained. |
| A nonuniform calculation matches uniform | Confirm that weighting is enabled and the intended mode is selected. Soft strength zero recovers uniform weighting. |

Use a new output name when comparing runs. Preserve parameter files with
results so that changes in order, weights, centering, or output switches remain
traceable.

## Information to retain when reporting an issue

Record the repository revision or package manifest, Windows and SGeMS versions,
algorithm name, saved parameters, source dimensions and property name, and
exact reproduction steps. Include the validation message or crash details,
the example validation manifest if available, and window audit files
(`resolved_config.txt`, `window_results.csv`). A small synthetic input that
reproduces the issue is preferable when the original data cannot be shared.
