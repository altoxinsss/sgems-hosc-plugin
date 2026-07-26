# WithSGeMS

Source code for the SGeMS high-order spatial cumulant plugins.

## Contents

- `MoveWindow/HOSCWindow`: moving-window moment and cumulant plugin.
- `Template/HOSC2008`: template-based moment and cumulant plugin.
- `Template/GsTL`, `Template/GsTLAppli`, `Template/boost`, `Template/Eigen`, and `Template/tools`: source dependencies included by the projects.

## Build requirements

The projects target the legacy 32-bit SGeMS 2.0 plugin interface and use Visual C++ 2008 project files. A compatible SGeMS 2.0 development tree and its Qt3 libraries are required for compilation.

The project files contain local include and library paths that may need to be adjusted after cloning. The generated DLLs should be copied to the SGeMS `plugins/Geostat` directory.

This repository contains source files only. Build outputs, local IDE settings, binary libraries, generated executables, logs, and test data are intentionally excluded.
