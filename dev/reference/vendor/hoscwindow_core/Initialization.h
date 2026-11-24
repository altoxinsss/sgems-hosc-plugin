/*
声明程序需要的结构�?
根据输入文件名称初始化KEYS及GRID结构�?
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

//全局变量
extern double** ValueMap_2D;
extern double*** ValueMap_3D;
extern int TotalNodeNum;
extern int G_TemplateNum;
//存储所有的�?
extern double** G_Moments;
//存储所有矩的索�?
extern int** G_IndexSeries;
//存储所有矩的总数
extern int G_SeriesNum;

//暂且使用一维数组存储矩，在确保没有冗余计算的前提下使算法简单一�?
extern double* G_Moment;
//在未储存完所有矩之前，该值一直递增，记录当前存储到哪一个矩�?
extern int G_CurrMomentPos;
//在存储矩的值的同时，还需要存储矩的索�?
extern char** G_MomentIndex;
//判断矩是否存�?
extern bool G_MomentFull;
//矩索引哈希表
extern int* G_MomentHash;
extern int G_MomentHashSize;
extern bool G_MomentHashReady;

//记录原始图像在三个坐标轴上的最大维�?
extern int MaxX;
extern int MaxY;
extern int MaxZ;

//用于处理划分
#define MAX_PARTITIONS 10000   //最多的划分数量
#define MAX_BLOCK_SIZE 32      //每块最大元素个�?
#define MAX_BLOCKS     10      //最多块�?

//保存一个划分结果，k个子集，每个最多MAX_BLOCK_SIZE个元�?
typedef struct {
    int blocks[MAX_BLOCKS][MAX_BLOCK_SIZE];
    int block_sizes[MAX_BLOCKS]; // 每个子集的实际元素数
} Partition;

//所有划分结�?
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


//该结构体的核心是计算value的值，即累积量�?
//value值只能通过以下两种方式计算�?
//1. 在InitialCumulant函数过程中发现特殊索引，将下一层的对应累积量直接赋值为0
//2. 完成所有结构体InitialCumulant后，将结构体中的矩、组合数、下一层累积量相乘并相加得到当前层的value
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

//搜寻容忍区内的点
// 定义点结构体
typedef struct {
    double x;
    double y;
} Point;


#ifdef __cplusplus
extern "C" {
#endif

//初始化KEYS结构�?
int GetKeys(FILE* FileKeys);

//初始化GRID结构�?
int GetGrid(FILE* FileGrid);

//初始化KEYS，GRID结构体，内部分别调用GetKeys以及GetGrid
int InitialStruct();

//释放KEYS，GRID结构体以及部分全局变量
int Destruction();

//获取数组中最大的元素，用来确定每个方向上最大的元素是什�?
int GetMaxElem(Nodes* pt_Nodes, int NodeNum);

//for dll
int InitialStruct2(const char* filename1, const char* filename2);

//为第二个dll设计,.KEYS通过文件读取�?GRID通过SGeMS里的Object读取
int InitialStruct3(const char* fileKeys, double** Grid);

//为第二个dll设计
int GetGrid2(double** Grid);
int InitialStruct4(Keys* Key, double** grid);

//Test
int TestPTGrid();

void mylog(const char* msg, double num);


#ifdef __cplusplus
}
#endif


#endif
