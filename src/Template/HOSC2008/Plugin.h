#ifndef HOSC_TEMPLATE_PLUGIN_H
#define HOSC_TEMPLATE_PLUGIN_H
#include "../../shared/EnginePlugin.h"
class __declspec(dllexport) HOSC_TEMPLATE : public hosc_sgems::EnginePlugin {
public:
    HOSC_TEMPLATE(): hosc_sgems::EnginePlugin(true) {}
    virtual ~HOSC_TEMPLATE() {}
    static Named_interface* create_new_interface(std::string&) { return new HOSC_TEMPLATE; }
};
#endif
