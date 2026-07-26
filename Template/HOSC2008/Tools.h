


#ifndef __TOOLS_H__
#define __TOOLS_H__
#include "Initialization.h"

#define myEmpty -999



#ifdef __cplusplus
extern "C" {
#endif


double CumulantOrder1(Grid* pt_Grids);


double CumulantOrder2(Template* Templates, int TemplateNum);


double CumulantOrder3(Template* Templates, int TemplateNum);


int Expect(double* Variable, int Size, double* Result);


double Combinatorial(int a, int b);


double Moment(int** index, int Order);


double Moment2(int** index, int Order);


double Moment3(int** index, int Order);


int** Catalan(int n, int* returnSize);


double** GetAllMoments_test();
double** GetAllMoments();


int IndexRank(int* index, int Order);


void Merge(int arr[], int left, int mid, int right);


void MergeSort(int arr[], int n);




int GetZeroMean(double* Variable, int Size);


int SaveResult(double* CumulantResult, int** LAGS, int LagNum, int DirectionNum, int start, int end);

double** SaveResult2Grid(double* CumulantResult, int** LAGS, int LagNum, int DirectionNum, int start, int end);


int RemoveZero(int** index, int* Order);

int Init3DVec(double**** VEC, int size1, int Size2, int Size3);

int Free3DVec(double**** VEC, int Size1, int Size2);

int Free3DVec_Int(int**** VEC, int Size1, int Size2);

int Free2DVec(double*** VEC, int Size);

int Free2DVec_int(int*** VEC, int Size);

int FreeTemplate(Template** pt_Template);

int FreeCumulantIndex(CumulantIndex** pt_Cumulant);

int Free2DVec_char(char*** VEC, int Size);


double Moment4(int** index, int* pOrder);

double Moment5(int* index, int Order);


int GetBell(int Order);

int Factorial(int n);

bool myStrcmp(int* index, int indexLen, char* str2);





#ifdef __cplusplus
}
#endif


#endif