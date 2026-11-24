#ifndef HOSC_NATIVE_CORE_H
#define HOSC_NATIVE_CORE_H
#include <string>
namespace hosc_native {
// Synchronous, in-process computation. No executable or helper DLL is loaded.
void run_config(bool is_template, const std::string &directory, const std::string &filename);
} // namespace hosc_native
#endif
