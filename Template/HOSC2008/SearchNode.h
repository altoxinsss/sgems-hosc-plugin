



#ifndef __SEARCHNODE_H__
#define __SEARCHNODE_H__
#include "Tools.h"

#ifdef __cplusplus
extern "C" {
#endif





int SaveTemplateResult(Template* Templates, int Order, int TemplateNum);


int InitialLagMat(int*** LagMat, int* TotalDirectionsNum);


int GetIndex(int* Num, int Size, bool MAX_MIN);


int FreeTemplate(Template** pt_Template);


int TemplateConstruct_2D(int* pint32_LAG, int Order, Template** Templates, int* TemplateNums, int range1, int range2);


int TemplateConstruct_3D(int* pint32_LAG, int Order, Template** Templates, int* TemplateNums, int range1, int range2);


int TemplateConstruct_2D_old(int* pint32_LAG, int Order, Template** Templates, int* TemplateNums, int range1, int range2);


int TemplateConstruct_3D_old(int* pint32_LAG, int Order, Template** Templates, int* TemplateNums, int range1, int range2);






Point Subtract(Point a, Point b);


double Dot(Point a, Point b);


bool PointInRectangle(Point p, Point v0, Point v1, Point v3);


int SearchInTolerentZone(double** TolerentZone, int* CoorX, int* CoorY, int CoorZ);

#ifdef __cplusplus
}
#endif

#endif