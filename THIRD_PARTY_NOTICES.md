# Third-party dependencies and data

The current plugin package targets SGeMS x64 Beta, Qt 4.7.3 and the Visual C++ 2010 x64 runtime. SGeMS/Qt/Boost/GsTL headers, host import libraries, compiler tools and host runtime binaries are external development dependencies; they are not bundled in this repository. Sources and compatibility details are in [SDK sources](docs/SDK.md).

The four files in `bin/win64` are the two HOSC plugin DLLs and their interface forms, not an SGeMS runtime. Follow the applicable upstream licenses when obtaining or redistributing dependencies. Source-level notices and the existing upstream notices in `licenses/` are retained. No new repository-wide license is assigned to third-party code or data.

This product includes software developed by The Stanford Center for Reservoir Forecasting, Stanford University.

Stanford VI-E is provided by the [SCRF repository](https://github.com/SCRFpublic/Stanford-VI-E). The source model and derived paper data remain subject to applicable upstream terms and attribution. The synthetic introductory fixtures in `examples/` are separate from that reservoir model.

`dev/reference/vendor` holds the existing numerical reference sources used by developer checks, not an archived SDK or a runtime dependency of the plugins.
