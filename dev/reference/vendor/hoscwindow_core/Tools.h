/*
此处存放功能性函数
*/
#ifndef __TOOLS_H__
#define __TOOLS_H__
#include "Initialization.h"

#define myEmpty -999

//------------------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif

//一阶累积量直接调用
double CumulantOrder1(Grid* pt_Grids);

//二阶累积量直接调用
double CumulantOrder2(Template* Templates, int TemplateNum);

//三阶累积量直接调用
double CumulantOrder3(Template* Templates, int TemplateNum);

//获取数组的均值
int Expect(double* Variable, int Size, double* Result);

//组合数运算函数
double Combinatorial(int a, int b);

//矩计算函数(优化前)
double Moment(int** index, int Order);

//矩计算函数(优化后)
double Moment2(int** index, int Order);

//矩计算函数(在外部一次性计算出所有矩，该函数用于找到指定的矩并直接返回值)
double Moment3(int** index, int Order);

//一次性求得所有不重复的索引序列（利用Catalan数）
int** Catalan(int n, int* returnSize);

//一次性求得所有矩
double** GetAllMoments_test();
double** GetAllMoments();

//求出指定索引序列在全部Catlan数序列中的位置，从而定位到指定索引下的矩
int IndexRank(int* index, int Order);

// 合并两个已排序的子数组
void Merge(int arr[], int left, int mid, int right);

// 迭代式归并排序，不使用递归
void MergeSort(int arr[], int n);

//------------------------------------------------------------------------------------

//将通过空间模板构造函数生成的随机变量，修改为0均值随机变量
int GetZeroMean(double* Variable, int Size);

//保存累积量结果
int SaveResult(double* CumulantResult, int** LAGS, int LagNum, int DirectionNum, int start, int end);
//for dll
double** SaveResult2Grid(double* CumulantResult, int** LAGS, int LagNum, int DirectionNum, int start, int end);

//累积量索引重新排序，并去零
int RemoveZero(int** index, int* Order);

int Init3DVec(double**** VEC, int size1, int Size2, int Size3);

int Free3DVec(double**** VEC, int Size1, int Size2);

int Free3DVec_Int(int**** VEC, int Size1, int Size2);

int Free2DVec(double*** VEC, int Size);

int Free2DVec_int(int*** VEC, int Size);

int FreeTemplate(Template** pt_Template);

int FreeCumulantIndex(CumulantIndex** pt_Cumulant);

int Free2DVec_char(char*** VEC, int Size);

//记录了矩阶数的变化
double Moment4(int** index, int* pOrder);

double Moment5(int* index, int Order);

//在指定阶数下索引可被划分的数量
int GetBell(int Order);

int Factorial(int n);

bool myStrcmp(int* index, int indexLen, char* str2);





#ifdef __cplusplus
}
#endif


#endif