#ifndef HOSC_WINDOW_PLUGIN_H
#define HOSC_WINDOW_PLUGIN_H
#include "../../shared/EnginePlugin.h"
class __declspec(dllexport) HOSC_WINDOW : public hosc_sgems::EnginePlugin {
public:
    HOSC_WINDOW(): hosc_sgems::EnginePlugin(false) {}
    virtual ~HOSC_WINDOW() {}
    static Named_interface* create_new_interface(std::string&) { return new HOSC_WINDOW; }
};
#endif
