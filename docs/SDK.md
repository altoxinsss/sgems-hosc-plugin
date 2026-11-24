# SDK sources and compatibility

The local verified kit is at `D:\Development\DLL\DLL\WithSGeMS64\sdk`. It combines compatible headers with import libraries generated from the target host. It is not an official complete Beta SDK and does not rebuild SGeMS itself.

| Component | Directory below the SDK root | Source |
|---|---|---|
| VC2010 x64 compiler and libraries | `winsdk71_extracted/vc_stdamd64/Program Files(64)/Microsoft Visual Studio 10.0/VC` | Microsoft Windows SDK 7.1 |
| VC2010 headers | `winsdk71_extracted/vc_stdx86/Program Files/Microsoft Visual Studio 10.0/VC/include` | Same SDK |
| Windows headers and x64 system libraries | `winsdk71_extracted/WinSDKBuild_amd64/Program Files/Microsoft SDKs/Windows/v7.1` | Same SDK |
| SGeMS headers | `SGeMS-5d87afac331e12b5cdf6c1c1b2f0e6e20e9485d7` | [2011 source mirror](https://github.com/fnavarrov/SGeMS/tree/5d87afac331e12b5cdf6c1c1b2f0e6e20e9485d7) |
| Qt headers | `qt-4.7.3/include` | [Qt v4.7.3](https://github.com/qt/qt/tree/v4.7.3) |
| Boost headers | `boost_1_47_0` | [Boost 1.47.0 archive](https://archives.boost.io/release/1.47.0/source/) |
| GsTL headers | `gstl/GsTL` | Compatible existing SGeMS developer source, with its original license |
| Host import libraries | `host_imports` | Generated from the verified host's actual exports |

Microsoft's SDK information is available at https://www.microsoft.com/en-us/download/details.aspx?id=8442. Download/extraction provenance remains with the external kit; compiler and SDK binaries are not redistributed here.

`tools/prepare_host_imports.py --sdk-root <kit> --host <host>` regenerates import libraries from that host, recording DLL hashes. This does not prove compatibility with another SGeMS release: matching headers, ABI and real-host tests remain necessary.

`tools/prepare_qt_headers.py --sdk-root <kit>` creates forwarding headers from the Qt source tree. A source-only Qt checkout also needs the compatible Qt build configuration; the existing verified kit already has it. Normal plugin builds do not need to regenerate the headers or convert the supplied Qt4 forms.

The package's `bin/win64/host_libraries.csv` records the supported host DLL fingerprints. No host DLL is replaced during installation.
