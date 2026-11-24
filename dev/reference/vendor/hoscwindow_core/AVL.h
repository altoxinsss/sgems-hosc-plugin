#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <stdbool.h>
#include "Tools.h"

//AVL成功查找的次数
extern int TreeFounded;

//AVL树结构体
typedef struct Node {
    double value;
    char* index;
    //double index;
    struct Node *left;
    struct Node *right;
    int height;
    bool existed;
} Node;

//栈元素结构体，用于释放AVL树
typedef struct
{
    Node* node;
    bool visited;
} AVLStackElement;

//栈结构体，用于释放AVL树
typedef struct
{
    AVLStackElement* elements;
    int capacity;
    int top;
}AVLStack;

//初始化AVL栈
int InitAVLStack(AVLStack* avlStack, int initCapacity);

//压栈（自动扩容）
int PushAVLStack(AVLStack* avlStack, AVLStackElement avlElement);
//出栈
AVLStackElement* PopAVLStack(AVLStack* avlStack);

//释放AVL树栈内存
int FreeAVLStack(AVLStack* avlStack);

//释放AVL树内存
int FreeAVLTree(Node** node);

//全局变量记得更新(单AVL树，不对阶数进行限制)
extern Node* RootAVL;

//多AVL树，不同的树对应不同阶数
extern Node** RootAVL2;

//初始化根节点(单AVL树)
int InitRoot(int* index, int value, int Order);

//初始化根节点(多AVL树)，一次初始化会初始化Order个AVL树
int InitRoot2(int Order);

//AVL树的更新，重新创造一个Order阶的AVL树
int RefreshTree(int Order);

// 获取节点的高度
int height(Node *node);

// 计算平衡因子
int balanceFactor(Node *node);

// 右旋转操作
Node *rightRotate(Node *y);

// 左旋转操作
Node *leftRotate(Node *x);

// 插入节点到树中
Node *insert(Node* node, int* indexI, int value, int Order, bool existed);


//将累积量索引转换成十进制，方便存储和进行索引，暂时使用int，可能出现越界情况
double Index2Dex(int* indexI, int Order);

//将累积量滞后距按索引次序转换成字符串，方便存储和进行索引
char* Index2String(int* indexI, int Order);

//将累积量滞后距按索引次序转换成字符串，并在后面加上索引本身
char* Index2String2(int* index, int Order);

//在树中搜索指定索引(单AVL树)
Node* SearchIndex(int* indexI, int Order);

//在树中搜索指定索引(多AVL树)
Node* SearchIndex2(int* indexI, int Order);

//遍历查看树
int TraverseTree(Node *root);

int max(int a, int b);



