#ifndef __FUNCTIONCALL_H__
#define __FUNCTIONCALL_H__
#include "CumulantCal.h"
#include "WeightedMomentCore.h"

#ifdef __cplusplus
extern "C" {
#endif


int RecCalculation(const char* KeyName, const char* GridName);


int RecCalculation4Dll(const char* KeyName, double** grid);


double** RecCalculation4DLL2(Keys* Key, double**& grid, int* gridSize, int* crossection, int crossectionCount);


double** Moment4Dll(Keys* Key, double**& grid, int* gridSize, int* crossection, int crossectionCount);

double** RecCalculation4DLL3(Keys* Key, double**& grid, int* gridSize, int* crossection, int crossectionCount);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
double** WeightedMoment4Dll(
    Keys* Key,
    double**& grid,
    int* gridSize,
    int* crossection,
    int crossectionCount,
    const weighted_moment_core::WeightedMomentOptions& options,
    bool momentTarget);

double** WeightedRecCalculation4DLL2(
    Keys* Key,
    double**& grid,
    int* gridSize,
    int* crossection,
    int crossectionCount,
    const weighted_moment_core::WeightedMomentOptions& options,
    bool momentTarget);
#endif

#endif
