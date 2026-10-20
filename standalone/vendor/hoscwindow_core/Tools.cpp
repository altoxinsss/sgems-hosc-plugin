#include "Tools.h"
#include "ExplainCounters.h"
//New approach of calculating N-Order cumulant

int** G_IndexSeries;
int G_SeriesNum;
HOSCExplainCounters g_hosc_counters;

void HOSCResetCounters(void)
{
    memset(&g_hosc_counters, 0, sizeof(g_hosc_counters));
}

void HOSCDumpCounters(FILE* out)
{
    FILE* fp = out == NULL ? stdout : out;
    fprintf(fp,
            "[HOSC_COUNTERS] template_attempts=%llu template_success=%llu template_fail=%llu "
            "cumulant_eval_attempts=%llu cumulant_eval_success=%llu cumulant_eval_fail=%llu "
            "moment_precompute_calls=%llu moment_precompute_fail=%llu moment_eval_calls=%llu "
            "original_recursive_calls=%llu original_recursive_node_visits=%llu "
            "original_recursive_cache_queries=%llu original_recursive_cache_hits=%llu "
            "optimized_iterative_calls=%llu avl_search_calls=%llu avl_search_hits=%llu avl_insert_calls=%llu\n",
            g_hosc_counters.template_attempts,
            g_hosc_counters.template_success,
            g_hosc_counters.template_fail,
            g_hosc_counters.cumulant_eval_attempts,
            g_hosc_counters.cumulant_eval_success,
            g_hosc_counters.cumulant_eval_fail,
            g_hosc_counters.moment_precompute_calls,
            g_hosc_counters.moment_precompute_fail,
            g_hosc_counters.moment_eval_calls,
            g_hosc_counters.original_recursive_calls,
            g_hosc_counters.original_recursive_node_visits,
            g_hosc_counters.original_recursive_cache_queries,
            g_hosc_counters.original_recursive_cache_hits,
            g_hosc_counters.optimized_iterative_calls,
            g_hosc_counters.avl_search_calls,
            g_hosc_counters.avl_search_hits,
            g_hosc_counters.avl_insert_calls);
}

double CumulantOrder1(Grid* pt_Grids)
{
    if(pt_Grids == NULL)
    {
        fprintf(stderr, "Error in function CumulantOrder1: Grid is NULL\n");
        return -1;
    }
    double Sum = 0;
    double Mean = 0;
    for(int i = 0; i < pt_Grids->ND; i++)
    {
        Sum += pt_Grids->pt_Nodes->val;
    }
    Mean = Sum / pt_Grids->ND;
    if(Mean != 0)
    {
        Sum = 0;
        for(int i = 0; i < pt_Grids->ND; i++)
        {
            Sum = Sum + pt_Grids->pt_Nodes[i].val - Mean;
        }
        Mean = Sum / pt_Grids->ND;
    }
    return Mean;
}

double CumulantOrder2(Template* Templates, int Size)
{
    if(Templates == NULL)
    {
        fprintf(stderr, "Error in function CumulantOrder2: Templates is NULL\n");
        return -1;
    }
    if(Size == 0)
    {
        fprintf(stderr, "Error in function CumulantOrder2: Size is NULL\n");
        return -1;
    }
    double* Temp1 = (double*)calloc(Size, sizeof(double));
    double* Temp2 = (double*)calloc(Size, sizeof(double));
    if(Temp1 == NULL || Temp2 == NULL)
    {
        perror("Error in function CumulantOrder2: Calloc return NULL\n");
        return -1;
    }
    for(int i = 0; i < Size; i++)
    {
        Temp1[i] = Templates[i].Value[0];
        Temp2[i] = Templates[i].Value[1];
    }
    GetZeroMean(Temp1, Size);
    GetZeroMean(Temp2, Size);
    for(int i = 0; i < Size; i++)
    {
        Temp1[i] = Temp1[i] * Temp2[i];
    }
    double Result;
    if(Expect(Temp1, Size, &Result) == -1)
    {
        fprintf(stderr, "Error in function CumulantOrder2: Expect return -1\n");
        return -1;
    }
    free(Temp1);
    free(Temp2);
    if(FreeTemplate(&Templates) == -1)
    {
        fprintf(stderr, "Failed to free template\n");
        return myEmpty;
    }
    Temp1 = NULL;
    Temp2 = NULL;
    //Templates = NULL;
    return Result;
}

double CumulantOrder3(Template* Templates, int Size)
{
    if(Templates == NULL)
    {
        fprintf(stderr, "Error in function CumulantOrder3: Templates is NULL\n");
        return -1;
    }
    if(Size == 0)
    {
        fprintf(stderr, "Error in function CumulantOrder3: Size is NULL\n");
        return -1;
    }
    double result;
    double* Temp1 = (double*)calloc(Size, sizeof(double));
    double* Temp2 = (double*)calloc(Size, sizeof(double));
    double* Temp3 = (double*)calloc(Size, sizeof(double));
    for(int i = 0; i < Size; i++)
    {
        Temp1[i] = Templates[i].Value[0];
        Temp2[i] = Templates[i].Value[1];
        Temp3[i] = Templates[i].Value[2];
    }
    GetZeroMean(Temp1, Size);
    GetZeroMean(Temp2, Size);
    GetZeroMean(Temp3, Size);
    for(int i = 0; i < Size; i++)
    {
        Temp1[i] = Temp1[i] * Temp2[i] * Temp3[i];
    }
    if(Expect(Temp1, Size, &result) == -1)
    {
        fprintf(stderr, "Error in function CumulantOrder3: Expect return -1\n");
        return -1;
    }
    free(Temp1);
    free(Temp2);
    free(Temp3);
    FreeTemplate(&Templates);
    Temp1 = NULL;
    Temp2 = NULL;
    Temp3 = NULL;
    //Templates = NULL;
    return result;
}

int Expect(double* Variable, int Size, double* Result)
{
    if(Variable == NULL)
    {
        fprintf(stderr, "Error in function Expect: Varible is NULL\n");
        return -1;
    }
    if(Size == 0)
    {
        fprintf(stderr, "Error in function Expect: Size is 0\n");
        return -1;
    }
    double result;
    double Sum = 0;
    for(int i = 0; i < Size; i++)
    {
        Sum+=Variable[i];
    }
    result = Sum / Size;
    *Result = result;
    if(result == 0)
        return 1;
    return 0;
}

int GetZeroMean(double* Variable, int Size)
{
    if(Variable == NULL)
    {
        fprintf(stderr, "Error in function GetZeroMean: Varible is NULL\n");
        return -1;
    }
    if(Size == 0)
    {
        fprintf(stderr, "Error in function GetZeroMean: Size is 0\n");
        return -1;
    }
    double result;
    if(Expect(Variable, Size, &result) == 1)
        return 0;
    printf("%lf", result);
    for(int i =0; i < Size; i++)
    {
        Variable[i]-=result;
    }
    return 1;
}


int SaveResult(double* CumulantResult, int** LAGS, int LagNum, int DirectionNum, int start, int end)
{
    if(CumulantResult == NULL)
    {
        fprintf(stderr, "Error in function SaveResult: CumulantResult is NULL\n");
        return -1;
    }
    if(LagNum == 0 || DirectionNum == 0)
    {
        fprintf(stderr, "Error in function SaveResult: No direction input\n");
    }
    FILE* file = fopen("result.csv", "w");
    if(file == NULL)
    {
        perror("Error opening file");
        return -1;
    }
    //fprintf(file, "Directions: Dir1, Dir2, Dir3 ......\n");
    fprintf(file, "Index,   ");
    for(int i = 0; i < DirectionNum; i++)
    {
        fprintf(file, "h%d,   ", i+1);
    }
    fprintf(file, "Value\n");
    for(int i = start; i <= end; i++)
    {
        fprintf(file, "%d,   ", i);
        for(int j = 0; j < DirectionNum; j++)
        {
            fprintf(file, "%d,   ", LAGS[i][j]);
        }
       fprintf(file, "%f\n", CumulantResult[i]);
    }
    fclose(file);
    return 0;
}

int SaveResult2(double* CumulantResult, int** LAGS, int LagNum, int DirectionNum, int start, int end, const char* fileName)
{
    if(CumulantResult == NULL)
    {
        fprintf(stderr, "Error in function SaveResult: CumulantResult is NULL\n");
        return -1;
    }
    if(LagNum == 0 || DirectionNum == 0)
    {
        fprintf(stderr, "Error in function SaveResult: No direction input\n");
    }
    //char fileName[100];
    //snprintf(fileName, sizeof(fileName), "stanford%d.txt", DirectionNum+1);
    char NewFileName[100];
    sprintf(NewFileName, "%s%s%s", "TempData/", fileName, ".txt");
    FILE* file = fopen(NewFileName, "w");
    if(file == NULL)
    {
        printf("%s\n", NewFileName);
        perror("Error opening file");
        return -1;
    }
    //fprintf(file, "Directions: Dir1, Dir2, Dir3 ......\n");
    fprintf(file, "Order%d\n", DirectionNum+1);
    fprintf(file, "%d\n", DirectionNum+2);
    fprintf(file, "index\n");
    for(int i = 0; i < DirectionNum; i++)
    {
        fprintf(file, "h%d\n", i+1);
    }
    fprintf(file, "Value\n");
    for(int i = start; i <= end; i++)
    {
        fprintf(file, "%d    ", i);
        for(int j = 0; j < DirectionNum; j++)
        {
            fprintf(file, "%d    ", LAGS[i][j]);
        }
       fprintf(file, "%.15lf\n", CumulantResult[i]);
    }
    fclose(file);
    return 0;
}

double** SaveResult2Grid(double* CumulantResult, int** LAGS, int LagNum, int DirectionNum, int start, int end)
{
    if(end < start)
    {
        return NULL;
    }
    if(CumulantResult == NULL)
    {
        fprintf(stderr, "Error in function SaveResult: CumulantResult is NULL\n");
        return NULL;
    }
    if(LagNum == 0 || DirectionNum == 0)
    {
        fprintf(stderr, "Error in function SaveResult: No direction input\n");
    }

    int rowsNum = end - start + 1;
    double** targetGrid = (double**)calloc(rowsNum, sizeof(double*));
    if(targetGrid == NULL)
    {
        return NULL;
    }

    for(int i = 0; i < rowsNum; i++)
    {
        targetGrid[i] = (double*)calloc(DirectionNum + 1, sizeof(double));
        if(targetGrid[i] == NULL)
        {
            return NULL;
        }
        for(int j = 0; j < DirectionNum; j++)
        {
            targetGrid[i][j] = LAGS[i + start][j];
        }
        targetGrid[i][DirectionNum] = CumulantResult[i];
    }
    return targetGrid;
}

double Combinatorial(int a, int b)
{
    if(a == 0 && b == 0)
        return 1;
    if(a == 0 && b != 0)
        return myEmpty;
    if(a > 0 && b == 0)
        return 1;
    if(a < 0 || b < 0)
        return 0;
    if(a < b)
        return 0;
    if(a == b)
        return 1;
    double FactoriA = (double)a;
    double FactoriB = (double)b;
    b--;
    for( ;b >= 1; --b)
    {
        FactoriA *= (--a);
        FactoriB *= b;
    }
    return FactoriA/FactoriB;
}

double Moment(int** index, int Order)
{
    g_hosc_counters.moment_eval_calls++;
    if(G_Template == NULL || index == NULL)
    {
        fprintf(stderr, "Error in function Moment: No Elements\n");
        return myEmpty;
    }
    //m00...000 = 1(Peter J. Smith)
    if(Order == 0)
    {
        return 1;
    }

    double** TempElems = (double**)calloc(Order, sizeof(double*));
    if(TempElems == NULL)
    {
        perror("Error in function Moment: Calloc return NULL\n");
        return myEmpty;
    }
    int numOfNoZero = 0;
    int OrderOfNoZero = Order;
    for(int i = 0; i < Order; i++)
    {
        if(numOfNoZero == Order)
            break;
        TempElems[numOfNoZero] = (double*)calloc(G_TemplateNum, sizeof(double));
        if(TempElems[numOfNoZero] == NULL)
        {
            perror("Error in function Moment: Calloc return NULL\n");
            for(int j = 0; j < numOfNoZero; j++)
            {
                free(TempElems[j]);
                TempElems[j] = NULL;
            }
            TempElems[numOfNoZero] = NULL;
            free(TempElems);
            TempElems = NULL;
            return myEmpty;
        }
        //每一个TempElems都是一个新的随机变�?
        if(index[numOfNoZero][0] == index[numOfNoZero][1])
        {
            i--;
            OrderOfNoZero--;
        }
        else
        {
            for(int j = 0; j < G_TemplateNum; j++)
            {
                TempElems[i][j] = G_Template[j].Value[index[numOfNoZero][0]-index[numOfNoZero][1]];
            }
        }
        numOfNoZero++;
    }
    //printf("%d\n", OrderOfNoZero);
    if(OrderOfNoZero == 0)
    {
        Free2DVec(&TempElems, Order);
        return 1;
    }
    //零均值随机变量，一阶矩�?
    if(OrderOfNoZero == 1)
    {
        Free2DVec(&TempElems, Order);
        return 0;
    }

    for(int i = 0; i < G_TemplateNum; i++)
    {
        double TempProduct = 1;
        for(int j = 0; j < OrderOfNoZero; j++)
        {
            TempProduct *= TempElems[j][i];
        }
        TempElems[0][i] = TempProduct;
    }
    double TempMean = 0;
    for(int i = 0; i < G_TemplateNum; i++)
    {
        TempMean += TempElems[0][i];
    }
    TempMean /= G_TemplateNum;
    Free2DVec(&TempElems, Order);
    return TempMean;
    
    //return 1;
}

double Moment2(int** index, int Order)
{
    g_hosc_counters.moment_eval_calls++;
    if(G_Template == NULL || index == NULL)
    {
        fprintf(stderr, "Error in function Moment2: No Elements\n");
        return myEmpty;
    }
    //m00...000 = 1(Peter J. Smith)
    if(Order == 0)
    {
        return 1;
    }
    else if(Order == 1)
    {
        return 0;
    }
    int NoZeroOrder = Order;
    for(int i = Order - 1; i >= 0; i--)
    {
        if(index[i][1] == index[i][0])
            NoZeroOrder--;
        else
            break;
    }
    if(NoZeroOrder == 0)
        return 1;
    else if(NoZeroOrder == 1)
        return 0;
    int* TempIndex = (int*)calloc(NoZeroOrder, sizeof(int));
    int TempOrder = NoZeroOrder;
    int TempI = 0;
    for(int i = 0; i < TempOrder; i++)
    {
        TempIndex[TempI] = index[i][0] - index[i][1];
        if(TempIndex[TempI] == 0)
            NoZeroOrder--;
        else
            TempI++;
    }
    double Accumulate = 0;
    double Multiply = 1;
    for(int i = 0; i < G_TemplateNum; i++)
    {
        for(int j = 0; j < NoZeroOrder; j++)
        {
            Multiply *= G_Template[i].Value[TempIndex[j]];
        }
        Accumulate += Multiply;
        Multiply = 1;
    }
    double result = Accumulate / G_TemplateNum;
    free(TempIndex);
    return result;
}

bool myStrcmp(int* index, int indexLen, char* str2)
{
    char* str1 = (char*)calloc(indexLen+1, sizeof(char));
    if(str1 == NULL)
    {
        perror("Initial: calloc return NULL\n");
        return myEmpty;
    }
    bool Skipped = false;
    int zeroNum = 0;
    for(int i = indexLen; i >= 1; i--)
    {
        if(index[i] == 0 && Skipped == false)
        {
            zeroNum++;
            continue;
        }
        if(index[i] != 0 && Skipped == false)
            Skipped = true;
        str1[i-1] = '0' + index[i];
    }
    str1[indexLen-zeroNum] = '\0';
    if(strcmp(str1, str2) == 0)
    {
        free(str1);
        return true;
    }
    free(str1);
    return false;
}

double Moment5(int* index, int Order)
{
    if(G_Template == NULL || index == NULL || G_MomentIndex == NULL)
    {
        fprintf(stderr, "Error in function Moment5: No Elements || G_Moments == NULL\n");
        return myEmpty;
    }
    
    if(Order == 0)
    {
        return 1;
    }

    if(index[0] == 1)
    {
        return 0;
    }
    

    double Accumulate = 0;
    double Multiply = 1;
    for(int i = 0; i < G_TemplateNum; i++)
    {
        for(int j = 1; j <= Order; j++)
        {
            Multiply *= pow(G_Template[i].Value[j], (double)index[j]);
        }
        Accumulate += Multiply;
        Multiply = 1;
    }
    double result = Accumulate / G_TemplateNum;

    G_MomentIndex[G_CurrMomentPos] = (char*)calloc(Order+1, sizeof(char));
    if(G_MomentIndex[G_CurrMomentPos] == NULL)
    {
        perror("Initial: calloc return NULL\n");
        return myEmpty;
    }
    bool Skipped = false;
    int zeroNum = 0;
    for(int i = Order; i >= 1; i--)
    {
        if(index[i] == 0 && Skipped == false)
        {
            zeroNum++;
            continue;
        }
        if(index[i] != 0 && Skipped == false)
            Skipped = true;
         G_MomentIndex[G_CurrMomentPos][i-1] = '0' + index[i];
    }
    G_MomentIndex[G_CurrMomentPos][Order-zeroNum] = '\0';
    //G_CurrMomentPos++;

    return result;
}

int RemoveZero(int** index, int* Order)
{
    if(index == NULL)
    {
        perror("Error in function RemoveZero: index is NULL\n");
        return myEmpty;
    }
    int ZeroNum = 0;
    int* TempIndex = *index;
    int TempOrder = *Order;
    //TempOrder--;
    //MergeSort(TempIndex, TempOrder);
    for(int i = TempOrder; i >=1; i--)
    {
        if(TempIndex[i] == 0)
            ZeroNum++;
        else
            break;
    }
    if(ZeroNum == 0)
    {
        int TempSum = 0;
        for(int i = TempOrder; i >= 1; i--)
        {
            TempSum += TempIndex[i];
        }
        if(TempIndex[0] != TempSum)
        {
            TempIndex[0] = TempSum;
            *index = TempIndex;
        }
        return 0;
    }
    TempOrder -= ZeroNum;
    if(TempOrder == 0)
    {
        //所有index全为0
        *Order = 0;
        TempIndex[0] = 0;
        return -1;
    }
    int* FinalIndex = (int*)calloc(TempOrder+1, sizeof(int));
    if(FinalIndex == NULL)
    {
        perror("Error in function TureCumulantLoop: Calloc return NULL\n");
        return myEmpty;
    }
    int jointOrder = 0;
    for(int i = 1; i <= TempOrder; i++)
    {
        FinalIndex[i] = TempIndex[i];
        jointOrder += TempIndex[i];
    }
    FinalIndex[0] = jointOrder;
    free(TempIndex);
    *index = FinalIndex;
    *Order = TempOrder;

    return 0;
}

int Init3DVec(double**** VEC, int Size1, int Size2, int Size3)
{
    double*** DestVec = (double***)calloc(Size1, sizeof(double**));
    if(DestVec == NULL)
    {
        perror("Error in function Init3DVec: calloc return NULL\n");
        return -1;
    }
    for(int i = 0; i < Size1; i++)
    {
        DestVec[i] = (double**)calloc(Size2, sizeof(double*));
        if(DestVec[i] == NULL)
        {
            perror("Error in function Init3DVec: calloc return NULL\n");
            return -1;
        }
        for(int j = 0; j < Size2; j++)
        {
            DestVec[i][j] = (double*)calloc(Size3, sizeof(double));
            if(DestVec[i][j] == NULL)
            {
                perror("Error in function Init3DVec: calloc return NULL\n");
                return -1;
            }
        }
    }
    *VEC = DestVec;
    DestVec = NULL;
    return 0;
}

int Free3DVec(double**** VEC, int Size1, int Size2)
{
    if(VEC == NULL)
    {
        fprintf(stderr, "Error in function Free3DVec: VEC is NULL\n");
        return -1;
    }
    double*** DestVec = *VEC;
    for(int i = 0; i < Size1; i++)
    {
        if(DestVec[i] == NULL)
        {
            fprintf(stderr, "Error in functoin Free3DVec: Can't free because DestVec[i] is NULL\n");
            return -1;
        }
        for(int j = 0; j < Size2; j++)
        {
            if(DestVec[i][j] == NULL)
            {
                fprintf(stderr, "Error in functoin Free3DVec: Can't free because DestVec[i][j] is NULL\n");
                return -1;
            }
            free(DestVec[i][j]);
            DestVec[i][j] = NULL;
        }
        free(DestVec[i]);
        DestVec[i] = NULL;
    }
    free(DestVec);
    return 0;
}

int Free3DVec_Int(int**** VEC, int Size1, int Size2)
{
    if(VEC == NULL)
    {
        fprintf(stderr, "Error in function Free3DVec: VEC is NULL\n");
        return -1;
    }
    int*** DestVec = *VEC;
    for(int i = 0; i < Size1; i++)
    {
        if(DestVec[i] == NULL)
        {
            fprintf(stderr, "Error in functoin Free3DVec: Can't free because DestVec[i] is NULL\n");
            return -1;
        }
        for(int j = 0; j < Size2; j++)
        {
            if(DestVec[i][j] == NULL)
            {
                fprintf(stderr, "Error in functoin Free3DVec: Can't free because DestVec[i][j] is NULL\n");
                return -1;
            }
            free(DestVec[i][j]);
        }
        free(DestVec[i]);
    }
    free(DestVec);
    return 0;
}

int Free2DVec(double*** VEC, int Size)
{
    if(VEC == NULL)
    {
        fprintf(stderr, "Error in function Free2DVec: VEC is NULL\n");
        return -1;
    }
    double** VECTemp = *VEC;
    
    //共有size个模�?
    for(int i = 0; i < Size; i++)
    {
        if(VECTemp[i] == NULL)
        {
            fprintf(stderr, "Error in function Free2DVec: Can't free because VecTemp[i] is NULL\n");
            return -1;
        }
        free(VECTemp[i]);
        VECTemp[i] = NULL;
    }
    free(VECTemp);
    VECTemp = NULL;
    return 0;
}

int Free2DVec_int(int*** VEC, int Size)
{
    if(VEC == NULL)
    {
        fprintf(stderr, "Error in function Free2DVec: VEC is NULL\n");
        return -1;
    }
    int** VECTemp = *VEC;
    
    //共有size个模�?
    for(int i = 0; i < Size; i++)
    {
        if(VECTemp[i] == NULL)
        {
            fprintf(stderr, "Error in function Free2DVec: Can't free because VecTemp[i] is NULL\n");
            return -1;
        }
        free(VECTemp[i]);
        VECTemp[i] = NULL;
    }
    free(VECTemp);
    VECTemp = NULL;
    return 0;
}

int Free2DVec_char(char*** VEC, int Size)
{
    if(VEC == NULL)
    {
        fprintf(stderr, "Error in function Free2DVec: VEC is NULL\n");
        return -1;
    }
    char** VECTemp = *VEC;
    if(VECTemp == NULL)
    {
        return 0;
    }
    
    //共有size个模�?
    for(int i = 0; i < Size; i++)
    {
        if(VECTemp[i] == NULL)
        {
            continue;
        }
        free(VECTemp[i]);
        VECTemp[i] = NULL;
    }
    free(VECTemp);
    *VEC = NULL;
    return 0;
}

int FreeTemplate(Template** pt_Template)
{
    Template* templateData = *pt_Template;
    if(*pt_Template == NULL)
    {
        //fprintf(stderr, "Error in function FreeTemplate: Template is NULL\n");
        return -1;
    }
    for(int i = 0; i < TotalNodeNum; i++)
    {
        if(templateData[i].Value != NULL)
        {
            free(templateData[i].Value);
            templateData[i].Value = NULL;
        }
        if(templateData[i].Lag != NULL)
        {
            free(templateData[i].Lag);
            templateData[i].Lag = NULL;
        }
    }
    free(templateData);
    //*pt_Template = NULL;
    return 0;
}


int FreeCumulantIndex(CumulantIndex** pt_Cumulant)
{
    CumulantIndex* TempCumulant = *pt_Cumulant;
    if(pt_Cumulant == NULL)
    {
        fprintf(stderr, "Error in function FreeCumulantIndex: Parameter is NULL\n");
        return myEmpty;
    }
    if(TempCumulant->moment != NULL)
    {
        free(TempCumulant->moment);
        TempCumulant->moment = NULL;
    }
    if(TempCumulant->MomentOrder != NULL)
    {
        free(TempCumulant->MomentOrder);
        TempCumulant->MomentOrder = NULL;
    }
    if(TempCumulant->combine != NULL)
    {
        free(TempCumulant->combine);
        TempCumulant->combine = NULL;
    }
    if(TempCumulant->index != NULL)
    {
        free(TempCumulant->index);
        TempCumulant->index = NULL;
    }
    return 0;
}

int** Catalan(int n, int* returnSize)  
{
    int cap = 128;
    int **res = (int**)malloc(sizeof(int*) * cap);
    *returnSize = 0;

    // 当前序列缓冲区，1-based 下标（更容易处理�?
    int *a = (int*)malloc(sizeof(int) * (n + 1));

    // 初始化为�?0（从 0 开始）
    for (int i = 1; i <= n; i++) 
        a[i] = 0;

    while (1) 
    {
        // 扩容
        if (*returnSize >= cap) 
        {
            cap <<= 1;
            res = (int**)realloc(res, sizeof(int*) * cap);
        }

        // 存储当前序列（转 0-based 存入�?
        int *copy = (int*)malloc(sizeof(int) * n);
        for (int i = 1; i <= n; i++) 
            copy[i-1] = a[i];
        res[(*returnSize)++] = copy;

        // 找到最右侧 still‐incrementable 的位�?k
        int k;
        for (k = n; k >= 1; k--) 
        {
            int limit = (k < n ? k : n);  // 每位最大值规�?
            if (a[k] < limit) 
            {
                a[k]++;  // 当前位增�?
                // 保持非递减：后续都设为 a[k]
                for (int j = k + 1; j <= n; j++) 
                    a[j] = a[k];
                break;
            }
        }

        if (k == 0)  // 所有可能枚举完
            break;
    }
    free(a);
    return res;
}


double** GetAllMoments_test()
{
    double** Moments = (double**)calloc(pt_Keys->ORDER+1, sizeof(double*));
    if(Moments == NULL)
    {
        perror("Error in function Catalan: calloc return NULL\n");
        return NULL;
    }
    for(int i = 2; i <= pt_Keys->ORDER; i++)
    {
        int SeriesNum;
        int** IndexSeries = Catalan(i, &SeriesNum);
        #if 0
        for(int test = 0; test < SeriesNum; test++)
        {
            for(int test2 = 0; test2 < i; test2++)
            {
                printf("%d ", IndexSeries[test][test2]);
            }
            printf("\n");
        }
        #endif
        Moments[i] = (double*)calloc(SeriesNum, sizeof(double));
        if(Moments[i] == NULL)
        {
            perror("Error in function GetAllMoments\n");
            return NULL;
        }
        for(int j = 0; j < SeriesNum; j++)
        {
            double Sum = 0;
            for(int k = 0; k < G_TemplateNum; k++)
            {
                double Mutiply = 1;
                for(int m = 0; m < i; m++)
                {
                    Mutiply *= G_Template[k].Value[IndexSeries[j][m]];
                }
                Sum += Mutiply;
            }
            Sum /= G_TemplateNum;
            Moments[i][j] = Sum;

            #if 0
            for(int k = 0; k < i; k++)
            {
                printf("%d ", IndexSeries[j][k]);
            }
            printf(" = %f\n", Sum);
            #endif
        }
        for(int j = 0; j < SeriesNum; j++)
        {
            free(IndexSeries[j]);
        }
        free(IndexSeries);
    }
    return Moments;
}

double** GetAllMoments()
{
    double** Moments = (double**)calloc(pt_Keys->ORDER+1, sizeof(double*));
    if(Moments == NULL)
    {
        perror("Error in function GetAllMoments: Calloc return NULL\n");
        return NULL;
    }
    Moments[0] = (double*)calloc(pt_Keys->ORDER+1, sizeof(double));
    if(Moments[0] == NULL)
    {
        perror("Error in function GetAllMoments: Calloc return NULL\n");
        return NULL;
    }
    int SeriesNum;
    int** IndexSeries = Catalan(pt_Keys->ORDER, &SeriesNum);
    Moments[1] = (double*)calloc(SeriesNum, sizeof(double));
    if(Moments[1] == NULL)
    {
        perror("Error in function GetAllMoments: Calloc return NULL\n");
        return NULL;
    }
    //计算每一阶矩有多少个索引以及记录IndexSeries中每一组索引是多少阶的
    int SkipNum = 0;
    for(int i = 0; i < SeriesNum; i++)
    {
        int TempOrder = pt_Keys->ORDER;
        if(RemoveZero(&IndexSeries[i], &TempOrder) == -1 || TempOrder == 1)
        {
            Moments[1][i] = TempOrder;
            SkipNum++;
            continue;
        }
        Moments[0][TempOrder]++;
        Moments[1][i] = TempOrder;
    }
    //初始化记录矩的二维矩�?
    for(int i = 2; i <= pt_Keys->ORDER; i++)
    {
        Moments[i] = (double*)calloc(Moments[0][i], sizeof(double));
    }
    int* MomentIndex = (int*)calloc(pt_Keys->ORDER+1, sizeof(int));
    //计算所有矩
    for(int i = SkipNum; i < SeriesNum; i++)
    {
        //double Mutiply = 1;
        double Sum = 0;
        for(int j = 0; j < G_TemplateNum; j++)
        {
            double Mutiply = 1;
            for(int k = 0; k < Moments[1][i]; k++)
            {
                Mutiply *= G_Template[j].Value[IndexSeries[i][k]];
                //printf("%lf\n", G_Template[j].Value[IndexSeries[i][k]]);
            }
            Sum += Mutiply;
        }
        Sum /= G_TemplateNum;
        Moments[(int)Moments[1][i]][MomentIndex[(int)Moments[1][i]]] = Sum;
        MomentIndex[(int)Moments[1][i]]++;
    }
    //Free2DVec_int(&IndexSeries, SeriesNum);
    G_IndexSeries = IndexSeries;
    G_SeriesNum = SeriesNum;
    free(MomentIndex);
    //Test Moments
    #if 0
    for(int i = 2; i < pt_Keys->ORDER+1; i++)
    {
        printf("Curr Order: %d\n", i);
        for(int j = 0; j < Moments[0][i]; j++)
        {
            printf("%.15lf ", Moments[i][j]);
        }
        printf("\n");
    }
    #endif
    return Moments;
}

void Merge(int arr[], int left, int mid, int right) 
{
    int n1 = mid - left + 1;
    int n2 = right - mid;
    int *L = (int *)malloc(n1 * sizeof(int));
    int *R = (int *)malloc(n2 * sizeof(int));
    if (L == NULL || R == NULL) 
    {
        perror("Error in function merge: malloc return NULL\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < n1; i++)
        L[i] = arr[left + i];
    for (int j = 0; j < n2; j++)
        R[j] = arr[mid + 1 + j];
    
    int i = 0, j = 0, k = left;
    while (i < n1 && j < n2) 
    {
        if (L[i] <= R[j])
            arr[k++] = L[i++];
        else
            arr[k++] = R[j++];
    }
    while (i < n1)
        arr[k++] = L[i++];
    while (j < n2)
        arr[k++] = R[j++];
    
    free(L);
    free(R);
}

void MergeSort(int arr[], int n) 
{
    for (int curr_size = 1; curr_size < n; curr_size *= 2) 
    {
        for (int left_start = 0; left_start < n - 1; left_start += 2 * curr_size) 
        {
            int mid = left_start + curr_size - 1;
            if(mid >= n) break;  // 若左子数组不存在则退�?
            int right_end = left_start + 2 * curr_size - 1;
            if (right_end >= n)
                right_end = n - 1;
            Merge(arr, left_start, mid, right_end);
        }
    }
}

#if 0
int IndexRank(int* index, int Order) 
{
    MergeSort(index, Order);
    if (index[0] != 0) 
    {
        fprintf(stderr, "Error in function IndexRank: First index must be 0\n");
        return 1;
    }
    for (int i = 1; i < Order; i++) 
    {
        if (index[i] < index[i - 1] || index[i] > i + 1) 
        {
            fprintf(stderr, "Error in function IndexRank: The %dth number %d in seq is illegle(should be in range of [%d, %d])\n", i + 1, index[i], index[i - 1], i);
            printf("Order is: %d\n", Order);
            printf("Index is: ");
            for(int j = 0; j < Order; j++)
            {
                printf("%d ", index[j]);
            }
            printf("\n");
            return -1;
        }
    }

    int dp_rows = Order + 2;
    int dp_cols = Order + 2;
    long long** dp = (long long**)malloc(dp_rows * sizeof(long long*));
    if (dp == NULL) 
    {
        perror("Error in function IndexRank: malloc return NULL");
        return -1;
    }
    for (int i = 0; i < dp_rows; i++) 
    {
        dp[i] = (long long*)malloc(dp_cols * sizeof(long long));
        if (dp[i] == NULL) 
        {
            printf("内存分配失败。\n");
            for (int j = 0; j < i; j++) free(dp[j]);
            free(dp);
            return -1;
        }
    }

    // 边界：dp[n+1][prev] = 1
    for (int prev = 0; prev < dp_cols; prev++) 
    {
        dp[Order + 1][prev] = 1;
    }

    // 自底向上计算 dp[i][prev]
    for (int i = Order; i >= 2; i--) 
    {
        for (int prev = 0; prev < dp_cols; prev++) 
        {
            long long sum = 0;
            for (int d = prev; d <= i; d++) {
                sum += dp[i + 1][d];
            }
            dp[i][prev] = sum;
        }
    }

    // 计算字典序排�?
    long long rank = 0;
    int currentPrev = index[0];
    for (int pos = 2; pos <= Order; pos++) 
    {
        int a_val = index[pos - 1];
        for (int d = currentPrev; d < a_val; d++) 
        {
            if (d <= pos) 
            {
                rank += dp[pos + 1][d];
            }
        }
        currentPrev = a_val;
    }
    rank++;

    //printf("%lld\n", rank);

    for (int i = 0; i < dp_rows; i++) 
    {
        free(dp[i]);
    }
    free(dp);
    return rank-1;
}
#endif

int IndexRank(int *seq, int n) 
{
    if(G_Moments[n] == NULL)
    {
        fprintf(stderr, "Error in function IndexRank: G_Moments is NULL\n");
        return myEmpty;
    }

    int idx = -1;

    for (int i = 0; i < G_SeriesNum; i++) 
    {
        int match = 1;
        if(G_Moments[1][i] != n)
            continue;
        for (int j = 0; j < n; j++) {
            if (G_IndexSeries[i][j] != seq[j]) {
                match = 0;
                break;
            }
        }
        if (match) {
            idx = i;
            break;
        }
    }

    int TempOrder = n - 1;
    int* TempSeries = (int*)calloc(TempOrder, sizeof(int));
    if(TempSeries == NULL)
    {
        perror("Error in function IndexRank: Calloc return NULL\n");
        return myEmpty;
    }
    for(int i = TempOrder-1 ; i >= 0; i--)
    {
        if(i == TempOrder-1)
        {
            TempSeries[i] = pt_Keys->ORDER;
        }
        else
        {
           TempSeries[i] = pt_Keys->ORDER - TempOrder + i + 1;
        }
    }
    int idx2 = -1;
    for (int i = 0; i < G_SeriesNum; i++) 
    {
        if(G_Moments[1][i] != TempOrder)
            continue;
        int match = 1;
        for (int j = 0; j < TempOrder; j++) {
            if (G_IndexSeries[i][j] != TempSeries[j]) {
                match = 0;
                break;
            }
        }
        if (match) {
            idx2 = i;
            break;
        }
    }
    if(idx2 == -1)
    {
        fprintf(stderr, "Something wrong in function IndexRank when Calculated idx2\n");
        return myEmpty;
    }
    idx = idx - idx2 - 1;
    free(TempSeries);
    return idx;
}

double Moment3(int** index, int Order)
{
    g_hosc_counters.moment_eval_calls++;
    if(G_Template == NULL || index == NULL)
    {
        fprintf(stderr, "Error in function Moment2: No Elements\n");
        return myEmpty;
    }

    //m00...000 = 1(Peter J. Smith)
    if(Order == 0)
    {
        return 1;
    }
    else if(Order == 1)
    {
        return 0;
    }

    int* TempIndex = (int*)calloc(Order, sizeof(int));
    if(TempIndex == NULL)
    {
        perror("Error in function Moment3: calloc return NULL\n");
        return -1;
    }
    for(int i = 0; i < Order; i++)
    {
        if(i == Order - 1)
        {
             TempIndex[i] = index[i][0] + 1 - index[i][1];
        }
        else
            TempIndex[i] = index[i][0] - index[i][1];
    }
    //TempIndex[Order-1]++;
    //对TempIndex进行排序、去零处�?
    if(RemoveZero(&TempIndex, &Order) == -1)
    {
        //m00...000 = 1(Peter J. Smith)
        return 1;
    }
    if(Order == 1)
    {
        //0均�?
        return 0;
    }
    int MomentIndex = IndexRank(TempIndex, Order);
    double result = G_Moments[Order][MomentIndex];
    free(TempIndex);
    return result;
}

double Moment4(int** index, int* pOrder)
{
    int Order= *pOrder;
    if(G_Template == NULL || index == NULL)
    {
        fprintf(stderr, "Error in function Moment2: No Elements\n");
        return myEmpty;
    }

    //m00...000 = 1(Peter J. Smith)
    if(Order == 0)
    {
        return 1;
    }
    else if(Order == 1)
    {
        return 0;
    }

    int* TempIndex = (int*)calloc(Order, sizeof(int));
    if(TempIndex == NULL)
    {
        perror("Error in function Moment3: calloc return NULL\n");
        return -1;
    }
    for(int i = 0; i < Order; i++)
    {
        if(i == Order - 1)
        {
             TempIndex[i] = index[i][0] + 1 - index[i][1];
        }
        else
            TempIndex[i] = index[i][0] - index[i][1];
    }
    //TempIndex[Order-1]++;
    //对TempIndex进行排序、去零处�?
    if(RemoveZero(&TempIndex, &Order) == -1)
    {
        //m00...000 = 1(Peter J. Smith)
        return 1;
    }
    if(Order == 1)
    {
        //0均�?
        return 0;
    }
    int MomentIndex = IndexRank(TempIndex, Order);
    double result = G_Moments[Order][MomentIndex];
    free(TempIndex);
    *pOrder = Order;
    return result;
}

int GetBell(int Order)
{
    if(Order == 1)
        return 0;
    int PartitionNum = 0;
    for(int i = 0; i < Order; i++)
    {
        PartitionNum += (int)Combinatorial(Order-1, i);
    }
    return PartitionNum;
}

int Factorial(int n) {
    if (n < 0) 
    {
        printf("错误：阶乘只定义在非负整数上。\n");
        return 0;
    }
    int result = 1;
    for (int i = 2; i <= n; i++) 
    {
        result *= i;
    }
    return result;
}


int CheckIndex(int** LagMat, int LagNum)
{
    for(int i = 0; i < LagNum; i++)
    {
        for(int j = 0; j < pt_Keys->NDIR; j++)
        {
            printf("%d ", LagMat[i][j]);
        }
        printf("\n");
    }
    return 0;
}
