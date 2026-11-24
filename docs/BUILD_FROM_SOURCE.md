# Build the 64-bit SGeMS plugins

The host boundary must use the **VC2010 x64 ABI**, Release runtime and Qt 4.7.3. A modern compiler is suitable for the optional development reference, but is not a drop-in replacement for this host interface.

The external SDK is not included in this repository. Its required layout and provenance are described in [SDK sources](SDK.md).

```powershell
python tools/build_x64.py --sdk-root "D:\Development\DLL\DLL\WithSGeMS64\sdk" --check-only
python tools/build_x64.py --sdk-root "D:\Development\DLL\DLL\WithSGeMS64\sdk"
```

Alternatively set `HOSC_SDK_ROOT` and omit `--sdk-root`. The build reads `src/native`, `src/shared` and the two small plugin entry points. It does not generate numerical sources from a command-line implementation. Output and logs go to `build/x64`; installed files and `bin/win64` are not overwritten by a build.

The script uses `/EHsc /GR /MD /O2 /DNDEBUG /D_ITERATOR_DEBUG_LEVEL=0` and links against the target host's import libraries. `tests/native_core_sha256.csv` binds the imported numerical core to the previously verified source. The build manifest records the files actually compiled.

To compile the Qt4 interface and per-run measurement checks too:

```powershell
python tools/build_x64.py --sdk-root "D:\Development\DLL\DLL\WithSGeMS64\sdk" --with-tests
```

See [tests](../tests/README.md) for running the outputs in an isolated host. Do not copy newly built DLLs into a running SGeMS process.

For optional independent numerical comparisons, the small CMake project is under [dev/reference](../dev/reference/README.md). It is not needed to build or use the plugins.
