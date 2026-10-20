# Current installation and build instructions

Install the tested Windows binaries using [Quick start](QUICKSTART.md).
No compiler is required to install or run the plugins.

The stock SGeMS 2.0 host used by this release uses the VC7.1 C++ runtime.
Use [Building from source](BUILD_FROM_SOURCE.md) for a matching VC7.1 build.
The legacy `.vcproj` files are retained as historical project metadata;
they are not the release build entry point. A VC2008/VC9 DLL is not a
compatible replacement for the distributed VC7.1 DLLs.
