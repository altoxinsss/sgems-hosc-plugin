#ifndef __CUMULANTCALL_H__
#define __CUMULANTCALL_H__

#include "SearchNode.h"
#include "AVL.h"

#define MaxCumulantNum 10000
#define MaxFloor 2
#define DEFAULT_FST_FLOOR_NUM 300

extern int CurrCumulantNum;
extern int AllCumulantNum;
extern int FoundCumulantNum;
extern CumulantIndex** Floor2;


double BetterRecCumulant3(int* index, int Order);


double LowOrderCumulant(int* index, int Order);




void SetLowOrderCumulantMomentCacheMode(int enabled);


int PrecomputeMomentsForIndex(const int* index, int Order);


void SetFactorSwitchesForRun(int avlReuse, int iterativeSchedule, int hashLookup, int momentPrecompute);


void ResetFactorSwitchesForRun(void);


int InitialCoeff(CumulantIndex* CurrStruct);


double TrueCumulant(int* index, int Order);


double BellCumulantFromMultiplicityIndex(const int* index, int Order);


double BellRawMomentFromMultiplicityIndex(const int* index, int Order);

void backtrack(int* arr, int n, int k, int index, Partition* current, int used_blocks);
void generate_k_partitions(int* arr, int n, int k);
void print_partitions(int k);

#endif
