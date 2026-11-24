# Quick start

The supplied binaries are for the verified SGeMS x64 Beta / VC2010 / Qt 4.7.3 host. A 32-bit SGeMS installation cannot load them.

## Installation

Close the application, then run from this repository:

```powershell
.\Install.cmd "D:\SGeMS-x64-Beta"
```

To check compatibility without copying anything:

```powershell
.\tools\install_x64.ps1 -TargetRoot "D:\SGeMS-x64-Beta" -VerifyOnly
```

Previous plugin files are backed up under the target host's `plugin_backups` directory. The installation record identifies the four copied files and their prior hashes. To restore an older installation, close SGeMS and copy those backed-up files back to `plugins/Geostat`.

## First calculation

1. Load `examples/Inputs.prj/grid.sgems` with **Load Object**. Its object name is `grid`; the statistical property is `value`.
2. Select **Utilities → HOSC_TEMPLATE**.
3. Load `examples/Parameters/grid_template_uniform.par`.
4. Select a writable English-path **Record directory** on **Records**, then run.
5. Inspect the new `HOSC_TEMPLATE_grid_template_uniform` result object.

For a moving-window example, load `grid_window_uniform.par` in **HOSC_WINDOW**. For point-set examples, first load `points.sgems` and select the corresponding parameters. Each input/parameter pairing is listed in [examples](../examples/README.md).

The two external-weight examples select the existing `weight` property. The template weighting checkbox must remain enabled. Parameters do not import inputs; load the correct object first.

An empty record directory uses a temporary location. Select a persistent location to retain full-precision outputs and the automatic per-run performance CSV. No command-line calculation program is involved.
