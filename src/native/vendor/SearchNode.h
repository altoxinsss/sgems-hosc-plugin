/*
 GRID   KEYS       
      
*/
#ifndef __SEARCHNODE_H__
#define __SEARCHNODE_H__
#include "Tools.h"

#ifdef __cplusplus
extern "C" {
#endif

//    
//int ChangeTemplate(const char* filename, int Col);

//                
int SaveTemplateResult(Template* Templates, int Order, int TemplateNum);

//            
int InitialLagMat(int*** LagMat, int* TotalDirectionsNum);

//Get max or min number's index, get Max number's index when MAX_MIN == true 
int GetIndex(int* Num, int Size, bool MAX_MIN);

//              
int FreeTemplate(Template** pt_Template);

//                         Templates  
int TemplateConstruct_2D(int* pint32_LAG, int Order, Template** Templates, int* TemplateNums, int range1, int range2);

//                        Templates  
int TemplateConstruct_3D(int* pint32_LAG, int Order, Template** Templates, int* TemplateNums, int range1, int range2);

// Optional raw-sample trace for the Cartesian weighted adapter. Null callbacks
// preserve the legacy selection and centering path. Coordinates are zero-based
// dense-grid indices; values are supplied before the legacy centering pass.
typedef void (*TemplateRawTrace)(int variables, const int* x, const int* y,
    const int* z, const double* values, void* context);
void SetTemplateRawTrace(TemplateRawTrace callback, void* context,
    const unsigned char* informed_mask, int nx, int ny);

//                         Templates  
int TemplateConstruct_2D_old(int* pint32_LAG, int Order, Template** Templates, int* TemplateNums, int range1, int range2);

//                        Templates  
int TemplateConstruct_3D_old(int* pint32_LAG, int Order, Template** Templates, int* TemplateNums, int range1, int range2);

//                    
//double* GetAllMoment();

// 3D              
//        
Point Subtract(Point a, Point b);

//      
double Dot(Point a, Point b);

//            
bool PointInRectangle(Point p, Point v0, Point v1, Point v3);

//        
int SearchInTolerentZone(double** TolerentZone, int* CoorX, int* CoorY, int CoorZ);

#ifdef __cplusplus
}
#endif

#endif
