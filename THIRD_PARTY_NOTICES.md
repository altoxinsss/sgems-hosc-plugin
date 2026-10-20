# Third-party dependencies and data

The SGeMS plugin projects require a compatible external SGeMS SDK and Qt3
development environment. The SDK supplies GsTL, GsTLAppli, Boost, and its
transitive headers; the build also requires the import libraries listed in
[the build guide](docs/BUILD_FROM_SOURCE.md). The host application, SDK,
Qt libraries, Microsoft compiler and Microsoft runtime binaries are not
distributed in this repository. The two HOSC plugin DLLs built from the
included plugin source are distributed under `bin/win32`.

Obtain these dependencies from their respective providers and follow their
licenses. GsTL and SGeMS headers carry their own notices; Boost uses the Boost
Software License; Qt3 has GPL/QPL/commercial licensing alternatives. Coin,
SoQt, Qwt, SimVoleon, and compiler/runtime components retain their own terms.
Selecting a compatible build environment does not change these licenses.

The retained plugin and standalone source files preserve their existing
file-level notices. `standalone/vendor` contains the numerical implementation
required by the current standalone build; it is not an archived SDK tree.

Stanford VI-E is provided by the [SCRF repository](https://github.com/SCRFpublic/Stanford-VI-E).
The source model and derived reservoir results remain subject to applicable
upstream terms and attribution. The small synthetic software-test fixtures
are identified separately from that model.

The upstream SGeMS license notices are retained in `licenses/`.
This product includes software developed by The Stanford Center for
Reservoir Forecasting, Stanford University.

No new repository-wide license is assigned to third-party components or data.
