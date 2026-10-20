#ifndef WITHOUT_SGEMS_GSTL_MESSAGES_COMPAT_H
#define WITHOUT_SGEMS_GSTL_MESSAGES_COMPAT_H

#include <iostream>

// The algorithm cores retain a few diagnostic calls from their SGeMS origin.
// This minimal compatibility layer keeps the standalone CLI independent from
// SGeMS while preserving those diagnostics on stderr.
namespace gstlIO {
struct End {};
static const End end = End();
}

class WithoutSgemsDiagnosticStream {
public:
    template <typename T>
    WithoutSgemsDiagnosticStream& operator<<(const T& value) {
        std::cerr << value;
        return *this;
    }

    WithoutSgemsDiagnosticStream& operator<<(const gstlIO::End&) {
        std::cerr << std::endl;
        return *this;
    }
};

static WithoutSgemsDiagnosticStream GsTLcerr;

#endif
