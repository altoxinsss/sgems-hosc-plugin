#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Tools.h"


extern int TreeFounded;


typedef struct Node {
    double value;
    char* index;

    struct Node *left;
    struct Node *right;
    int height;
    bool existed;
} Node;


typedef struct
{
    Node* node;
    bool visited;
} AVLStackElement;


typedef struct
{
    AVLStackElement* elements;
    int capacity;
    int top;
}AVLStack;


int InitAVLStack(AVLStack* avlStack, int initCapacity);


int PushAVLStack(AVLStack* avlStack, AVLStackElement avlElement);

AVLStackElement* PopAVLStack(AVLStack* avlStack);


int FreeAVLStack(AVLStack* avlStack);


int FreeAVLTree(Node** node);


extern Node* RootAVL;


extern Node** RootAVL2;


int InitRoot(int* index, int value, int Order);


int InitRoot2(int Order);


int RefreshTree(int Order);


int height(Node *node);


int balanceFactor(Node *node);


Node *rightRotate(Node *y);


Node *leftRotate(Node *x);


Node *insert(Node* node, int* indexI, int value, int Order, bool existed);



double Index2Dex(int* indexI, int Order);


char* Index2String(int* indexI, int Order);


char* Index2String2(int* index, int Order);


Node* SearchIndex(int* indexI, int Order);


Node* SearchIndex2(int* indexI, int Order);


int TraverseTree(Node *root);

int max(int a, int b);



