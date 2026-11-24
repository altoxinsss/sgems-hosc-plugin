# Optional numerical reference for developers

This C++17 project is retained only for regression comparisons with the native SGeMS plugins. It is not the recommended user workflow, is not required by the plugin build, and supplies no current paper performance measurements.

The old experiment drivers, output directories and user examples have been removed from this component. Numerical source, dependency notices and focused tests are preserved.

```powershell
cmake -S dev/reference -B build/reference -A x64 -DPython3_EXECUTABLE=<python.exe>
cmake --build build/reference --config Release --parallel
ctest --test-dir build/reference -C Release --output-on-failure
```

Python 3 is required for the integration checks. The build produces `template_cli` and `movewindow_cli` for tests that replay exported plugin inputs. Native SGeMS installation and normal calculations do not invoke them. See [the host tests](../../tests/README.md).
