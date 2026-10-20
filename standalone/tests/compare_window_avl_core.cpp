#include "CumulantCal.h"
#include "Initialization.h"

#include <cstdlib>
#include <iomanip>
#include <iostream>

namespace {

void clear_roots() {
    if (RootAVL2 == NULL) {
        return;
    }
    for (int i = 1; i <= 1; ++i) {
        if (RootAVL2[i] != NULL) {
            FreeAVLTree(&RootAVL2[i]);
            RootAVL2[i] = NULL;
        }
    }
    std::free(RootAVL2);
    RootAVL2 = NULL;
}

void clear_moments() {
    std::free(G_Moment);
    G_Moment = NULL;
    if (G_MomentIndex != NULL) {
        if (G_SeriesNum > 0) {
            Free2DVec_char(&G_MomentIndex, G_SeriesNum);
        }
        else {
            std::free(G_MomentIndex);
            G_MomentIndex = NULL;
        }
    }
    std::free(G_MomentHash);
    G_MomentHash = NULL;
    G_MomentHashSize = 0;
    G_MomentHashReady = false;
    G_MomentFull = false;
    G_SeriesNum = 0;
    G_CurrMomentPos = 0;
}

}

int main() {
    const double values[] = {1.0, 2.0, 4.0, 7.0};
    const int order = 4;

    G_Template = static_cast<Template*>(std::calloc(4, sizeof(Template)));
    G_TemplateNum = 4;
    for (int i = 0; i < G_TemplateNum; ++i) {
        G_Template[i].Value = static_cast<double*>(std::calloc(2, sizeof(double)));
        G_Template[i].Value[1] = values[i];
    }
    G_Template->index = static_cast<int*>(std::calloc(2, sizeof(int)));
    G_Template->index[0] = order;
    G_Template->index[1] = order;

    clear_moments();
    clear_roots();
    if (InitRoot2(1) == -1) {
        return 2;
    }
    const double result = BetterRecCumulant3(G_Template->index, 1);
    RefreshTree(1);
    clear_roots();
    clear_moments();

    for (int i = 0; i < G_TemplateNum; ++i) {
        std::free(G_Template[i].Value);
    }
    std::free(G_Template->index);
    std::free(G_Template);
    G_Template = NULL;
    G_TemplateNum = 0;

    std::cout << std::setprecision(17) << "avl_cumulant4=" << result << "\n";
    return 0;
}
