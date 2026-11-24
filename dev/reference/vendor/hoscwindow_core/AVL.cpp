#include "AVL.h"
#include "ExplainCounters.h"

Node* RootAVL;
Node** RootAVL2;
int TreeFounded;

int max(int a, int b)
{
    return (a > b) ? a : b;
}

int InitRoot(int* indexI, int value, int Order)
{
    if(RootAVL != NULL)
    {
        printf("Root Of AVL is not NULL\n");
        return 0;
    }
    if(indexI == NULL)
    {
        fprintf(stderr, "Error in function InitRoot: Parameter is NULL\n");
        return -1;
    }
    //double index = Index2Dex(indexI, Order);
    char* index = Index2String(indexI, Order);
    RootAVL = (Node*)malloc(sizeof(Node));
    if(RootAVL == NULL)
    {
        perror("Error in function InitRoot: Malloc return NULL\n");
        return -1;
    }
    RootAVL->height = 1;
    RootAVL->left = NULL;
    RootAVL->right = NULL;
    RootAVL->index = index;
    RootAVL->value = value;
    return 1;
}

int InitRoot2(int Order)
{
    if(RootAVL2 != NULL)
    {
        printf("Root Of AVL is not NULL\n");
        return 0;
    }
    
    RootAVL2 = (Node**)calloc(Order+1, sizeof(Node*));
    if(RootAVL2 == NULL)
    {
        perror("Error in function InitRoot2: Malloc return NULL\n");
        return -1;
    }
    for(int i = 1; i <= Order; i++)
    {
        RootAVL2[i] = (Node*)malloc(sizeof(Node));
        if(RootAVL2[i] == NULL)
        {
            perror("Error in function InitRoot2: Failed to Malloc\n");
            return -1;
        }
        RootAVL2[i]->height = 1;
        RootAVL2[i]->index = (char*)malloc(2);  // <-- 先分配内存
        if(RootAVL2[i]->index == NULL)
        {
            perror("Error in function InitRoot2: malloc for index failed\n");
            return -1;
        }
        strcpy(RootAVL2[i]->index, "0\0");
        //RootAVL2[i]->index = 0;
        RootAVL2[i]->left = NULL;
        RootAVL2[i]->right = NULL;
        RootAVL2[i]->value = 0;
    }
    return 1;
}

int RefreshTree(int OOrder)
{
    for(int Order = 1; Order <= OOrder; Order++)
    {
        FreeAVLTree(&RootAVL2[Order]);
        RootAVL2[Order] = (Node*)malloc(sizeof(Node));
        if(RootAVL2[Order] == NULL)
        {
            perror("Error in function RefreshTree: Malloc return NULL\n");
            return -1;
        }
        RootAVL2[Order]->height = 1;
        RootAVL2[Order]->index = (char*)malloc(2);
        if(RootAVL2[Order]->index == NULL)
        {
            perror("Error in function RefreshTree: malloc for index failed\n");
            free(RootAVL2[Order]);
            RootAVL2[Order] = NULL;
            return -1;
        }
        strcpy(RootAVL2[Order]->index, "0\0");
        RootAVL2[Order]->left = NULL;
        RootAVL2[Order]->right = NULL;
        RootAVL2[Order]->value = 0;
        RootAVL2[Order]->existed = false;
    }
    return 0;
}

// 获取节点的高度
int height(Node* node) 
{
    if (node == NULL) {
        return 0;
    }
    return node->height;
}

// 计算平衡因子
int balanceFactor(Node* node) 
{
    if (node == NULL) {
        return 0;
    }
    return height(node->left) - height(node->right);
}

// 右旋转操作
Node*rightRotate(Node* y) 
{
    Node*x = y->left;
    Node*T2 = x->right;

    x->right = y;
    y->left = T2;

    // 更新高度
    y->height = 1 + (height(y->left) > height(y->right) ? height(y->left) : height(y->right));
    x->height = 1 + (height(x->left) > height(x->right) ? height(x->left) : height(x->right));

    return x; // 新的根节点
}

// 左旋转操作
Node*leftRotate(Node* x) 
{
    Node*y = x->right;
    Node*T2 = y->left;

    y->left = x;
    x->right = T2;

    // 更新高度
    x->height = 1 + (height(x->left) > height(x->right) ? height(x->left) : height(x->right));
    y->height = 1 + (height(y->left) > height(y->right) ? height(y->left) : height(y->right));

    return y;
}


Node* insert(Node* root, int* indexI, int value, int Order, bool existed) {
    //printf("%d\n", Order);
    // 计算新节点的 index 值
    //double index = Index2Dex(indexI, Order);
    char* index = Index2String(indexI, Order);
    if(index == NULL)
    {
        return root;
    }
    // 如果树为空，则创建新节点并返回
    if (root == NULL) {
        Node* newNode = (Node*)malloc(sizeof(Node));
        newNode->value = value;
        newNode->index = index;
        newNode->height = 1;
        newNode->existed = false;
        newNode->left = newNode->right = NULL;
        g_hosc_counters.avl_insert_calls++;
        return newNode;
    }

    // 定义一个数组作为栈来记录从根到插入位置的路径
    Node* path[100];
    int top = 0;
    Node* curr = root;
    Node* parent = NULL;

    // 查找插入位置（查找路径）
    while (curr != NULL) {
        path[top++] = curr;  // 保存当前节点
        parent = curr;
        //if (index < curr->index)
        if(strncmp(index, curr->index, 128) < 0)
            curr = curr->left;
        //else if (index > curr->index)
        else if(strncmp(index, curr->index, 128) > 0)
            curr = curr->right;
        else
        {
            // 如果已存在相同的 index，直接返回原树
            free(index);
            return root;
        }
    }

    // 创建新节点
    Node* newNode = (Node*)malloc(sizeof(Node));
    //newNode->value = value;
    newNode->value = -999;
    newNode->index = index;
    newNode->height = 1;
    newNode->existed = false;
    newNode->left = newNode->right = NULL;
    g_hosc_counters.avl_insert_calls++;

    // 将新节点挂到父节点的相应分支上
    //if (index < parent->index)
    if(strncmp(index, parent->index, 128) < 0)
        parent->left = newNode;
    else
        parent->right = newNode;

    // 从插入点向上回溯，更新高度并检查是否失衡
    for (int i = top - 1; i >= 0; i--) {
        Node* currNode = path[i];
        // 更新当前节点高度
        currNode->height = 1 + max(height(currNode->left), height(currNode->right));
        int balance = balanceFactor(currNode);

        // 根据失衡情况执行旋转

        // 左左情况
        //if (balance > 1 && index < currNode->left->index) {
        if(balance > 1 && strncmp(index, currNode->left->index, 128) <0) {
            Node* rotated = rightRotate(currNode);
            // 如果当前节点不是根节点，则更新其父节点的指针
            if (i == 0)
                root = rotated;
            else {
                Node* parentNode = path[i - 1];
                if (parentNode->left == currNode)
                    parentNode->left = rotated;
                else
                    parentNode->right = rotated;
            }
            // 更新当前栈中记录的节点
            path[i] = rotated;
        }
        // 右右情况
        //else if (balance < -1 && index > currNode->right->index) {
        else if(balance < -1 && strncmp(index, currNode->right->index, 128) > 0) {
            Node* rotated = leftRotate(currNode);
            if (i == 0)
                root = rotated;
            else {
                Node* parentNode = path[i - 1];
                if (parentNode->left == currNode)
                    parentNode->left = rotated;
                else
                    parentNode->right = rotated;
            }
            path[i] = rotated;
        }
        // 左右情况
        //else if (balance > 1 && index > currNode->left->index) {
        else if(balance > 1 && strncmp(index, currNode->left->index, 128) > 0) {
            currNode->left = leftRotate(currNode->left);
            Node* rotated = rightRotate(currNode);
            if (i == 0)
                root = rotated;
            else {
                Node* parentNode = path[i - 1];
                if (parentNode->left == currNode)
                    parentNode->left = rotated;
                else
                    parentNode->right = rotated;
            }
            path[i] = rotated;
        }
        // 右左情况
        //else if (balance < -1 && index < currNode->right->index) {
        else if(balance < -1 && strncmp(index, currNode->right->index, 128) < 0) {
            currNode->right = rightRotate(currNode->right);
            Node* rotated = leftRotate(currNode);
            if (i == 0)
                root = rotated;
            else {
                Node* parentNode = path[i - 1];
                if (parentNode->left == currNode)
                    parentNode->left = rotated;
                else
                    parentNode->right = rotated;
            }
            path[i] = rotated;
        }
        // 如果平衡，则不需要旋转，继续向上回溯更新
    }

    return root;
}

#if 0
double Index2Dex(int* index, int Order)
{
    //传进来的index必须是已经排好序且没有0的
    if(index == NULL)
    {
        return -999;
    }
    if(Order == 0)
    {
        return 0;
    }
    double result = 0;
    int CurrNum;
    int digit = 1;
  
    for(int i = Order; i >= 1; i--)
    {
        CurrNum = index[i];
        //CurrNum = G_Template->Lag[index[i]];
        //if(CurrNum == 0) CurrNum = index[i] + pt_Keys->NX;
        if(result == 0)
            result += CurrNum;
        
        // else if(CurrNum == 0)
        // {
        //     digit++;
        //     continue;
        // }
        
        else
            result += CurrNum * pow(10, digit);
        while(true)
        {
            if((int)(result / pow(10, digit)) == 0)
                break;
            digit++;
        }
    }
    
    return result;
}
#endif

double Index2Dex(int* index, int Order)
{
    if (index == NULL) return -999;
    if (Order == 0) return 0;

    double result = 0;
    for (int i = 1; i <= Order; ++i) {
        int d = index[i];
        if (d < 0 || d > 9) {
            // 如果 index 中可能 >=10，需要改用更大进制或字符串
            // 这里暂以简单行为处理
        }
        result = result * 10.0 + (double)d;
    }
    return result;
}

#if 0
char* Index2String(int* index, int Order)
{
    if(index == NULL)
    {
        return NULL;
    }
    if(Order == 0)
    {
        return strdup("0");
    }
    char* buffer = (char*)malloc(128*sizeof(char));
    int offset = 0;
    for (int i = 0; i < Order; i++) {
        offset += snprintf(buffer + offset, 128 - offset, "%d ", G_Template->Lag[index[i]]);
    }
    buffer[offset - 1] = '\0'; // 去除末尾空格
    return buffer;
}
#endif

char* Index2String(int* index, int Order) 
{
    if (index == NULL) return NULL;
    if (Order <= 0) return strdup("0");
    char* buf = (char*)malloc(Order + 1);
    if (!buf) return NULL;
    for (int i = 1; i <= Order; ++i) {
        int d = index[i];
        if (d < 0 || d > 9) { free(buf); return NULL; } // 或者处理错误
        buf[i-1] = '0' + d;
    }
    buf[Order] = '\0';
    return buf;
}


char* Index2String2(int* index, int Order) {
    if (index == NULL) return NULL;
    if (Order == 0) return strdup("0");  // 返回常量字符串需谨慎（不可修改）

    // 动态计算缓冲区大小（更安全）
    const int max_int_len = 12;  // 假设每个整数最多占12字符（含空格）
    const int buffer_size = 2 * Order * max_int_len + 1;  // Lag和索引各占一半
    char* buffer = (char*)malloc(buffer_size);
    if (!buffer) return NULL;  // 内存分配失败处理

    int offset = 0;
    
    // 1. 拼接Lag数组元素
    for (int i = 0; i < Order; i++) {
        offset += _snprintf(buffer + offset, buffer_size - offset, 
                          "%d ", G_Template->Lag[index[i]]);
        if (offset >= buffer_size) {  // 溢出检查
            free(buffer);
            return NULL;
        }
    }

    // 2. 拼接索引数组元素
    for (int i = 0; i < Order; i++) {
        offset += _snprintf(buffer + offset, buffer_size - offset, 
                          "%d ", index[i]);
        if (offset >= buffer_size) {  // 溢出检查
            free(buffer);
            return NULL;
        }
    }

    // 处理末尾空格
    if (offset > 0) buffer[offset - 1] = '\0';  // 替换最后一个空格[3](@ref)
    return buffer;
}


Node* SearchIndex(int* indexI, int Order)
{
    if(indexI == NULL || RootAVL == NULL)
    {
        if(indexI == NULL && RootAVL != NULL)
            fprintf(stderr, "Error in function SearchIndex 286: IndexI is NULL\n");
        else if(RootAVL == NULL && indexI != NULL)
            fprintf(stderr, "Error in function SearchIndex 288: RootAVL is NULL\n");
        else
            fprintf(stderr, "Error in function SearchIndex 290: Parameter is NULL\n");
        return NULL;
    }
    //double index = Index2Dex(indexI, Order);
    char* index = Index2String(indexI, Order);
    Node* CurrNode = RootAVL;
    while(true)
    {
        if(index == CurrNode->index)
        //if(strncmp(index, CurrNode->index, 128)==0)
        {
            //free(index);
            return CurrNode;
        }
        else if(index < CurrNode->index)
        //else if(strncmp(index, CurrNode->index, 128)<0)
        {
            if(CurrNode->left != NULL)
            {
                CurrNode = CurrNode->left;
                continue;
            }
            else
                return NULL;
        }
        else if(index > CurrNode->index)
        //else if(strncmp(index, CurrNode->index, 128)>0)
        {
            if(CurrNode->right != NULL)
            {
                CurrNode = CurrNode->right;
                continue;
            }
            else
                return NULL;
        }
    }
}

Node* SearchIndex2(int* indexI, int Order)
{
    g_hosc_counters.avl_search_calls++;
    if(indexI == NULL || RootAVL2 == NULL)
    {
        if(indexI == NULL && RootAVL2 != NULL)
            fprintf(stderr, "Error in function SearchIndex 286: IndexI is NULL\n");
        else if(RootAVL2 == NULL && indexI != NULL)
            fprintf(stderr, "Error in function SearchIndex 288: RootAVL is NULL\n");
        else
            fprintf(stderr, "Error in function SearchIndex 290: Parameter is NULL\n");
        return NULL;
    }
    
    char* index = Index2String(indexI, Order);
    if(index == NULL)
    {
        return NULL;
    }
    //double test = Index2Dex(indexI, Order);
    //char* index = Index2String2(indexI, Order);
    Node* CurrNode = RootAVL2[Order];
    if(CurrNode == NULL || CurrNode->index == NULL)
    {
        free(index);
        return NULL;
    }
    while(true)
    {
        //if(index == CurrNode->index)
        if(strncmp(index, CurrNode->index, 128) == 0)
        {
            free(index);
            TreeFounded++;
            g_hosc_counters.avl_search_hits++;
            return CurrNode;
        }
        //else if(index < CurrNode->index)
        else if(strncmp(index, CurrNode->index, 128) < 0)
        {
            if(CurrNode->left != NULL)
            {
                CurrNode = CurrNode->left;
                continue;
            }
            else
            {
                //printf("Can't find index:%s\n", index);
                free(index);
                return NULL;
            }
        }
        //else if(index > CurrNode->index)
        else if(strncmp(index, CurrNode->index, 128) > 0)
        {
            if(CurrNode->right != NULL)
            {
                CurrNode = CurrNode->right;
                continue;
            }
            else
            {
                //printf("Can't find index:%s\n", index);
                free(index);
                return NULL;
            }
        }
    }
}

int InitAVLStack(AVLStack* avlStack, int initCapacity)
{
    avlStack->elements = (AVLStackElement*)calloc(initCapacity, sizeof(AVLStackElement));
    if(avlStack == NULL)
    {
        perror("Error in function InitAVLStack: Calloc return NULL\n");
        return -1;
    }
    avlStack->capacity = initCapacity;
    avlStack->top = -1;
    return 1;
}

int PushAVLStack(AVLStack* avlStack, AVLStackElement avlElement)
{
    if(avlStack->top == avlStack->capacity - 1)
    {
        avlStack->capacity *= 2;
        avlStack->elements = (AVLStackElement*)realloc(avlStack->elements, avlStack->capacity * sizeof(AVLStackElement));
        if(avlStack->elements == NULL)
        {
            perror("Error in function PushAVLStack: realloc return NULL\n");
            return -1;
        }
    }
    if (avlElement.visited)
        free(avlElement.node->index);
    avlStack->elements[++(avlStack->top)] = avlElement;
    return 1;
}

AVLStackElement* PopAVLStack(AVLStack* avlStack)
{
    if(avlStack->top == -1)
    {
        fprintf(stderr, "Error in function: PopAVLStack: Stack has no element\n");
        return NULL;
    }
    return &avlStack->elements[avlStack->top--];
}

int FreeAVLStack(AVLStack* avlStack)
{
    //free(avlStack->elements->node->index);
    free(avlStack->elements);
    avlStack->capacity = 0;
    avlStack->top = -1;
    return 1;
}

int FreeAVLTree(Node** node1)
{
    Node* node = *node1;
    if(node == NULL)
    {
        return 0;
    }
    AVLStack avlStack;
    InitAVLStack(&avlStack, 10);

    AVLStackElement root = {node, false};
    PushAVLStack(&avlStack, root);

    while(avlStack.top != -1)
    {
        AVLStackElement avlElement = *(PopAVLStack(&avlStack));
        if(avlElement.node == NULL)
            continue;
        if(avlElement.visited)
            free(avlElement.node);
        else
        {
            avlElement.visited = true;
            PushAVLStack(&avlStack, avlElement);
            if(avlElement.node->right != NULL)
            {
                AVLStackElement right = {avlElement.node->right, false};
                PushAVLStack(&avlStack, right);
            }
            if(avlElement.node->left != NULL)
            {
                AVLStackElement left = {avlElement.node->left, false};
                PushAVLStack(&avlStack, left);
            }
        }
    }
    FreeAVLStack(&avlStack);
    *node1 = NULL;
    return 1;
}


int TraverseTree(Node *root) {
    // 假设树中节点数量不超过100个
    Node *stack[100];
    int top = -1;  // 栈顶索引，初始为空
    Node *current = root;
    
    while (current != NULL || top != -1) {
        // 先将所有左子树节点入栈
        while (current != NULL) {
            stack[++top] = current;
            current = current->left;
        }
        
        // 弹出栈顶节点并访问
        current = stack[top--];
        printf("Index is:%f value is:%f\n", current->index, current->value);
        
        // 转向右子树
        current = current->right;
    }
    return 0;
}



