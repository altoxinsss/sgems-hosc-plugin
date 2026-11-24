/*
           ?
           KEYS GRID    ?
*/

#ifndef __Initialization_h__
#define __Initialization_h__
#include "stdio.h"
#include "stdlib.h"
#include "assert.h"
#include "string.h"
//#include <stdbool.h>
#include <math.h>
#include <fstream>
#include <GsTLAppli/utils/gstl_messages.h>


#define PI 3.14159265

//    
extern double** ValueMap_2D;
extern double*** ValueMap_3D;
extern int TotalNodeNum;
extern int G_TemplateNum;
//       ?
extern double** G_Moments;
//         ?
extern int** G_IndexSeries;
//        
extern int G_SeriesNum;

//                                 ?
extern double* G_Moment;
//                               ?
extern int G_CurrMomentPos;
//                    ?
extern char** G_MomentIndex;
//        ?
extern bool G_MomentFull;
//      
extern int* G_MomentHash;
extern int G_MomentHashSize;
extern bool G_MomentHashReady;

//                   ?
extern int MaxX;
extern int MaxY;
extern int MaxZ;

//      
#define MAX_PARTITIONS 10000   //       
#define MAX_BLOCK_SIZE 32      //         ?
#define MAX_BLOCKS     10      //     ?

//         k        MAX_BLOCK_SIZE    ?
typedef struct {
    int blocks[MAX_BLOCKS][MAX_BLOCK_SIZE];
    int block_sizes[MAX_BLOCKS]; //           
} Partition;

//       ?
extern Partition g_all_partitions[MAX_PARTITIONS];
extern int g_total_partitions;




typedef struct 
{
    int DIM;
    int REG_IRREG;
    int NX, NY, NZ, DIMX, DIMY, DIMZ;
    int ORDER;
    int NDIR;
	int* NOD;
    int* NLAG;
    int* LAG;
    int* LTOL;  //LAG TOLERANCE IN EACH DIRECTION
    int* AZM;   //AZIMUTH
    //int* ELEVATION;
    int* ATOL;  //AZIMUTH TOLERANCE IN EACH DIRECTION(AZIMUTH)
    int* BANDWH;    //WIDTH OF TOLERANCE'S ZONE
    int* DIP;   //NO USE
    int* DTOL;  //NO USE
    int* BANDWD;    //HEIGHT OF TOLERANCE'S ZONES

} Keys;



typedef struct 
{
    /* data */
    int x;
    int y;
    int z;
    double val;
    int index;
    bool Flag_Used;
}Nodes;


typedef struct 
{
    int ND;
    int NVAR;
    Nodes* pt_Nodes;
} Grid;

extern Keys* pt_Keys;
extern Grid* pt_Grids;

typedef struct
{
    double* Value;
    int* Lag;
    double Mean;
    int* index;
} Template;

extern Template* G_Template;


//          value         ?
//value               ?
//1.  InitialCumulant                           0
//2.        InitialCumulant                               value
typedef struct CumulantIndex
{
    double* moment;
    int* MomentOrder;
    double* combine;
    double value;
    int Order;
    bool Calculated;
    bool exist;
    struct CumulantIndex* nxtCumulant;
    int FormulaNum;
    struct CumulantIndex* lstCumulant;
    int* index;
}CumulantIndex;

//        
//       
typedef struct {
    double x;
    double y;
} Point;


#ifdef __cplusplus
extern "C" {
#endif

//   KEYS    ?
int GetKeys(FILE* FileKeys);

//   GRID    ?
int GetGrid(FILE* FileGrid);

//   KEYS GRID          GetKeys  GetGrid
int InitialStruct();

//  KEYS GRID           
int Destruction();

//                             ?
int GetMaxElem(Nodes* pt_Nodes, int NodeNum);

//for dll
int InitialStruct2(const char* filename1, const char* filename2);

//    dll  ,.KEYS        ?GRID  SGeMS  Object  
int InitialStruct3(const char* fileKeys, double** Grid);

//    dll  
int GetGrid2(double** Grid);
int InitialStruct4(Keys* Key, double** grid);

//Test
int TestPTGrid();

void mylog(const char* msg, int num);


#ifdef __cplusplus
}
#endif


#endif
