#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <stdbool.h>
#include "Tools.h"

//AVL       
extern int TreeFounded;

//AVL    
typedef struct Node {
    double value;
    char* index;
    //double index;
    struct Node *left;
    struct Node *right;
    int height;
    bool existed;
} Node;

//           AVL 
typedef struct
{
    Node* node;
    bool visited;
} AVLStackElement;

//         AVL 
typedef struct
{
    AVLStackElement* elements;
    int capacity;
    int top;
}AVLStack;

//   AVL 
int InitAVLStack(AVLStack* avlStack, int initCapacity);

//        
int PushAVLStack(AVLStack* avlStack, AVLStackElement avlElement);
//  
AVLStackElement* PopAVLStack(AVLStack* avlStack);

//  AVL    
int FreeAVLStack(AVLStack* avlStack);

//  AVL   
int FreeAVLTree(Node** node);

//        ( AVL          )
extern Node* RootAVL;

// AVL            
extern Node** RootAVL2;

//      ( AVL )
int InitRoot(int* index, int value, int Order);

//      ( AVL )          Order AVL 
int InitRoot2(int Order);

//AVL           Order  AVL 
int RefreshTree(int Order);

//        
int height(Node *node);

//       
int balanceFactor(Node *node);

//      
Node *rightRotate(Node *y);

//      
Node *leftRotate(Node *x);

//        
Node *insert(Node* node, int* indexI, int value, int Order, bool existed);


//                           int         
double Index2Dex(int* indexI, int Order);

//                            
char* Index2String(int* indexI, int Order);

//                             
char* Index2String2(int* index, int Order);

//         ( AVL )
Node* SearchIndex(int* indexI, int Order);

//         ( AVL )
Node* SearchIndex2(int* indexI, int Order);

//     
int TraverseTree(Node *root);

int max(int a, int b);



