#ifndef __FUNCTIONCALL_H__
#define __FUNCTIONCALL_H__
#include "CumulantCal.h"

#ifdef __cplusplus
extern "C" {
#endif

//       (    )
int RecCalculation(const char* KeyName, const char* GridName);

// dll  
int RecCalculation4Dll(const char* KeyName, double** grid);

// dll         
double** RecCalculation4DLL2(Keys* Key, double**& grid, int* gridSize, int* crossection, int crossectionCount);

//calculate moment
double** Moment4Dll(Keys* Key, double**& grid, int* gridSize, int* crossection, int crossectionCount);

double** RecCalculation4DLL3(Keys* Key, double**& grid, int* gridSize, int* crossection, int crossectionCount);

void ResetMomentWorkspaceForDLL(void);
void ClearAVLWorkspaceForDLL(int order);
int InitAVLWorkspaceForDLL(int order);
void ClearOneEvaluationWorkspaceForDLL(int order);

#ifdef __cplusplus
}
#endif

#endif
