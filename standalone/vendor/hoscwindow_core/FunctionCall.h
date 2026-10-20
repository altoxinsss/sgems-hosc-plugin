#ifndef __FUNCTIONCALL_H__
#define __FUNCTIONCALL_H__
#include "CumulantCal.h"

#ifdef __cplusplus
extern "C" {
#endif

//改进的递归算法(循环实现)
int RecCalculation(const char* KeyName, const char* GridName);

//为dll实现
int RecCalculation4Dll(const char* KeyName, double** grid);

//为dll实现（第二个版本）
double** RecCalculation4DLL2(Keys* Key, double**& grid, int* gridSize, int* crossection, int crossectionCount);

//calculate moment
double** Moment4Dll(Keys* Key, double**& grid, int* gridSize, int* crossection, int crossectionCount);

double** RecCalculation4DLL3(Keys* Key, double**& grid, int* gridSize, int* crossection, int crossectionCount);

#ifdef __cplusplus
}
#endif

#endif