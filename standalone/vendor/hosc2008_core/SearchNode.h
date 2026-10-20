/*
在GRID中根据KEYS文件指定的条件
构造空间模板
*/
#ifndef __SEARCHNODE_H__
#define __SEARCHNODE_H__
#include "Tools.h"

#ifdef __cplusplus
extern "C" {
#endif

//替换模板
//int ChangeTemplate(const char* filename, int Col);

//给出指定索引下的全部模板生成结果
int SaveTemplateResult(Template* Templates, int Order, int TemplateNum);

//给出所有可能的滞后距组合
int InitialLagMat(int*** LagMat, int* TotalDirectionsNum);

//Get max or min number's index, get Max number's index when MAX_MIN == true 
int GetIndex(int* Num, int Size, bool MAX_MIN);

//释放指定索引下生成的全部模板
int FreeTemplate(Template** pt_Template);

//根据指定索引生成该索引下的二维图像的所有模板，通过Templates传出
int TemplateConstruct_2D(int* pint32_LAG, int Order, Template** Templates, int* TemplateNums, int range1, int range2);

//根据指定索引生成该索引下的三维图像所有模板，通过Templates传出
int TemplateConstruct_3D(int* pint32_LAG, int Order, Template** Templates, int* TemplateNums, int range1, int range2);

//根据指定索引生成该索引下的二维图像的所有模板，通过Templates传出
int TemplateConstruct_2D_old(int* pint32_LAG, int Order, Template** Templates, int* TemplateNums, int range1, int range2);

//根据指定索引生成该索引下的三维图像所有模板，通过Templates传出
int TemplateConstruct_3D_old(int* pint32_LAG, int Order, Template** Templates, int* TemplateNums, int range1, int range2);

//在确定累积量矩的情况下，一次性求出所有矩
//double* GetAllMoment();

//在3D容忍区内搜索满足条件的邻接点
//计算两点的差向量
Point Subtract(Point a, Point b);

//求向量的点积
double Dot(Point a, Point b);

//判断目标点是否在容忍区内
bool PointInRectangle(Point p, Point v0, Point v1, Point v3);

//在容忍区内搜寻点
int SearchInTolerentZone(double** TolerentZone, int* CoorX, int* CoorY, int CoorZ);

#ifdef __cplusplus
}
#endif

#endif