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

// Iterative cumulant evaluation used by the optimized pipeline.
double BetterRecCumulant3(int* index, int Order);

// Direct low-order cumulant evaluation.
double LowOrderCumulant(int* index, int Order);

// Use the precomputed moment cache for low-order cumulants. This is needed by
// weighted AVL runs because first-order weighted moments may be non-zero.
void SetLowOrderCumulantMomentCacheMode(int enabled);

// Precompute all moments required by the current target index.
int PrecomputeMomentsForIndex(const int* index, int Order);

// Explicitly configure A/I/H/M factor switches for one run.
void SetFactorSwitchesForRun(int avlReuse, int iterativeSchedule, int hashLookup, int momentPrecompute);

// Reset factor switches so the next call reloads from environment.
void ResetFactorSwitchesForRun(void);

// Initialize coefficient and child-cumulant terms for one node.
int InitialCoeff(CumulantIndex* CurrStruct);

// True cumulant evaluation for high-order validation path.
double TrueCumulant(int* index, int Order);

// Bell cumulant evaluation on the multiplicity-index representation used by the template plugin.
double BellCumulantFromMultiplicityIndex(const int* index, int Order);

// Bell raw moment for diagnostics and validation.
double BellRawMomentFromMultiplicityIndex(const int* index, int Order);

void backtrack(int* arr, int n, int k, int index, Partition* current, int used_blocks);
void generate_k_partitions(int* arr, int n, int k);
void print_partitions(int k);

#endif
