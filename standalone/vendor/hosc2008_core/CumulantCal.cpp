#include "CumulantCal.h"

int CurrCumulantNum;

int AllCumulantNum;
int FoundCumulantNum;

char** G_MomentIndex;



CumulantIndex** Floor2;


Partition g_all_partitions[MAX_PARTITIONS];
int g_total_partitions = 0;

static double g_max_moment_err = 0.0;
static char g_max_moment_index[256];
#ifndef MOMENT_COMPARE_ENABLE
#define MOMENT_COMPARE_ENABLE 0
#endif

static bool g_factor_switch_ready = false;
static bool g_factor_avl_reuse = true;
static bool g_factor_iterative_schedule = true;
static bool g_factor_hash_lookup = true;
static bool g_factor_moment_precompute = true;
static bool g_low_order_uses_moment_cache = false;

static int ReadBinarySwitch(const char* name, int defaultValue)
{
    const char* v = getenv(name);
    if(v == NULL || v[0] == '\0')
    {
        return defaultValue;
    }
    if(strcmp(v, "0") == 0 || strcmp(v, "false") == 0 || strcmp(v, "FALSE") == 0 || strcmp(v, "off") == 0 || strcmp(v, "OFF") == 0)
    {
        return 0;
    }
    if(strcmp(v, "1") == 0 || strcmp(v, "true") == 0 || strcmp(v, "TRUE") == 0 || strcmp(v, "on") == 0 || strcmp(v, "ON") == 0)
    {
        return 1;
    }
    return defaultValue;
}

static void InitFactorSwitches(void)
{
    if(g_factor_switch_ready == true)
    {
        return;
    }
    g_factor_avl_reuse = ReadBinarySwitch("HOSC_FACTOR_A", 1) == 1;
    g_factor_iterative_schedule = ReadBinarySwitch("HOSC_FACTOR_I", 1) == 1;
    g_factor_hash_lookup = ReadBinarySwitch("HOSC_FACTOR_H", 1) == 1;
    g_factor_moment_precompute = ReadBinarySwitch("HOSC_FACTOR_M", 1) == 1;
    g_factor_switch_ready = true;
}

void SetFactorSwitchesForRun(int avlReuse, int iterativeSchedule, int hashLookup, int momentPrecompute)
{
    g_factor_avl_reuse = avlReuse != 0;
    g_factor_iterative_schedule = iterativeSchedule != 0;
    g_factor_hash_lookup = hashLookup != 0;
    g_factor_moment_precompute = momentPrecompute != 0;
    g_factor_switch_ready = true;
}

void ResetFactorSwitchesForRun(void)
{
    g_factor_switch_ready = false;
}

void SetLowOrderCumulantMomentCacheMode(int enabled)
{
    g_low_order_uses_moment_cache = enabled != 0;
}

static double Moment5ValueOnly(const int* index, int Order)
{
    if(G_Template == NULL || index == NULL)
    {
        fprintf(stderr, "Error in function Moment5ValueOnly: No Elements\n");
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
    return Accumulate / G_TemplateNum;
}

static double BellMomentValueOnly(const int* index, int Order)
{
    if(G_Template == NULL || index == NULL)
    {
        fprintf(stderr, "Error in function BellMomentValueOnly: No Elements\n");
        return myEmpty;
    }
    if(Order == 0)
    {
        return 1;
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
    return Accumulate / G_TemplateNum;
}

double BellRawMomentFromMultiplicityIndex(const int* index, int Order)
{
    return BellMomentValueOnly(index, Order);
}

static void LogMomentCompare(const int* index, int Order, double direct, double lookup)
{
    double err = fabs(direct - lookup);
    if(err <= g_max_moment_err)
    {
        return;
    }
    g_max_moment_err = err;

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
        g_max_moment_index[i-1] = '0' + index[i];
    }
    g_max_moment_index[Order - zeroNum] = '\0';

    FILE* f = fopen("moment_compare.log", "a");
    if(f != NULL)
    {
        fprintf(f, "max_err=%.15g order=%d index=%s direct=%.15g lookup=%.15g\n",
            g_max_moment_err, Order, g_max_moment_index, direct, lookup);
        fclose(f);
    }
}

static unsigned long HashMomentKey(const char* s)
{
    unsigned long h = 2166136261UL;
    while(*s)
    {
        h ^= (unsigned char)(*s);
        h *= 16777619UL;
        s++;
    }
    return h;
}

static int BuildMomentKey(const int* index, int Order, char* out, int outSize)
{
    if(outSize < Order + 1)
    {
        return -1;
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
        out[i-1] = '0' + index[i];
    }
    out[Order - zeroNum] = '\0';
    return 0;
}

static int BuildMomentHash(void)
{
    if(G_MomentIndex == NULL || G_SeriesNum <= 0)
    {
        return -1;
    }
    if(G_MomentHash != NULL)
    {
        free(G_MomentHash);
        G_MomentHash = NULL;
        G_MomentHashSize = 0;
        G_MomentHashReady = false;
    }

    int size = 1;
    while(size < G_SeriesNum * 2)
    {
        size <<= 1;
    }
    G_MomentHash = (int*)malloc(sizeof(int) * size);
    if(G_MomentHash == NULL)
    {
        perror("Error in function BuildMomentHash: malloc return NULL\n");
        return -1;
    }
    for(int i = 0; i < size; i++)
    {
        G_MomentHash[i] = -1;
    }

    int mask = size - 1;
    for(int i = 0; i < G_SeriesNum; i++)
    {
        if(G_MomentIndex[i] == NULL || G_MomentIndex[i][0] == '\0')
        {
            continue;
        }
        unsigned long h = HashMomentKey(G_MomentIndex[i]);
        int pos = (int)(h & (unsigned long)mask);
        while(true)
        {
            if(G_MomentHash[pos] == -1)
            {
                G_MomentHash[pos] = i;
                break;
            }
            if(strcmp(G_MomentIndex[G_MomentHash[pos]], G_MomentIndex[i]) == 0)
            {
                break;
            }
            pos = (pos + 1) & mask;
        }
    }

    G_MomentHashSize = size;
    G_MomentHashReady = true;
    return 0;
}

static int MomentLookupIndex(const int* index, int Order)
{
    InitFactorSwitches();
    if(G_MomentIndex == NULL || G_SeriesNum <= 0)
    {
        return -1;
    }

    if(g_factor_hash_lookup == true && G_MomentHashReady == true && G_MomentHash != NULL && G_MomentHashSize > 0)
    {
        char* key = (char*)calloc(Order + 1, sizeof(char));
        if(key != NULL && BuildMomentKey(index, Order, key, Order + 1) == 0)
        {
            int mask = G_MomentHashSize - 1;
            unsigned long h = HashMomentKey(key);
            int pos = (int)(h & (unsigned long)mask);
            int start = pos;
            while(true)
            {
                int idx = G_MomentHash[pos];
                if(idx == -1)
                {
                    break;
                }
                if(G_MomentIndex[idx] != NULL && strcmp(G_MomentIndex[idx], key) == 0)
                {
                    return idx;
                }
                pos = (pos + 1) & mask;
                if(pos == start)
                {
                    break;
                }
            }
        }
        if(key != NULL)
        {
            free(key);
        }
    }

    for(int i = 0; i < G_SeriesNum; i++)
    {
        if(G_MomentIndex[i] == NULL)
        {
            continue;
        }
        if(myStrcmp((int*)index, Order, G_MomentIndex[i]) == true)
        {
            return i;
        }
    }
    return -1;
}

static double RootMomentOrDirect(const int* index, int Order)
{
    if(g_factor_moment_precompute == true && G_Moment != NULL && G_SeriesNum > 0)
    {
        return G_Moment[0];
    }
    return Moment5ValueOnly(index, Order);
}

static double MomentFromCacheOrDirect(const int* index, int Order)
{
    if(index == NULL)
    {
        return myEmpty;
    }
    if(Order == 0 || index[0] == 0)
    {
        return 1.0;
    }
    int indexPos = MomentLookupIndex(index, Order);
    if(indexPos >= 0 && G_Moment != NULL)
    {
        return G_Moment[indexPos];
    }
    return Moment5ValueOnly(index, Order);
}

static double LowOrderCumulantFromMoments(const int* index, int Order)
{
    if(index == NULL)
    {
        return myEmpty;
    }
    if(Order == 0 || index[0] == 0)
    {
        return 0.0;
    }

    const int jointOrder = index[0];
    if(jointOrder == 1)
    {
        return MomentFromCacheOrDirect(index, Order);
    }

    if(jointOrder == 2)
    {
        int firstLabel = -1;
        int secondLabel = -1;
        for(int variable = 1; variable <= Order; variable++)
        {
            for(int repeat = 0; repeat < index[variable]; repeat++)
            {
                if(firstLabel < 0)
                {
                    firstLabel = variable;
                }
                else
                {
                    secondLabel = variable;
                }
            }
        }
        if(firstLabel < 1 || secondLabel < 1)
        {
            return myEmpty;
        }

        double m12 = MomentFromCacheOrDirect(index, Order);
        int* idx1 = (int*)calloc(Order + 1, sizeof(int));
        int* idx2 = (int*)calloc(Order + 1, sizeof(int));
        if(idx1 == NULL || idx2 == NULL)
        {
            free(idx1);
            free(idx2);
            return myEmpty;
        }
        idx1[0] = 1;
        idx2[0] = 1;
        idx1[firstLabel] = 1;
        idx2[secondLabel] = 1;
        double m1 = MomentFromCacheOrDirect(idx1, Order);
        double m2 = MomentFromCacheOrDirect(idx2, Order);
        free(idx1);
        free(idx2);
        if(m12 == myEmpty || m1 == myEmpty || m2 == myEmpty)
        {
            return myEmpty;
        }
        return m12 - m1 * m2;
    }

    if(jointOrder == 3)
    {
        int labels[3] = {0, 0, 0};
        int cursor = 0;
        for(int variable = 1; variable <= Order; variable++)
        {
            for(int repeat = 0; repeat < index[variable]; repeat++)
            {
                if(cursor >= 3)
                {
                    return myEmpty;
                }
                labels[cursor++] = variable;
            }
        }
        if(cursor != 3)
        {
            return myEmpty;
        }

        int* idx1 = (int*)calloc(Order + 1, sizeof(int));
        int* idx2 = (int*)calloc(Order + 1, sizeof(int));
        int* idx3 = (int*)calloc(Order + 1, sizeof(int));
        int* idx12 = (int*)calloc(Order + 1, sizeof(int));
        int* idx13 = (int*)calloc(Order + 1, sizeof(int));
        int* idx23 = (int*)calloc(Order + 1, sizeof(int));
        if(idx1 == NULL || idx2 == NULL || idx3 == NULL || idx12 == NULL || idx13 == NULL || idx23 == NULL)
        {
            free(idx1);
            free(idx2);
            free(idx3);
            free(idx12);
            free(idx13);
            free(idx23);
            return myEmpty;
        }

        idx1[0] = 1;
        idx2[0] = 1;
        idx3[0] = 1;
        idx12[0] = 2;
        idx13[0] = 2;
        idx23[0] = 2;
        idx1[labels[0]]++;
        idx2[labels[1]]++;
        idx3[labels[2]]++;
        idx12[labels[0]]++;
        idx12[labels[1]]++;
        idx13[labels[0]]++;
        idx13[labels[2]]++;
        idx23[labels[1]]++;
        idx23[labels[2]]++;

        double m123 = MomentFromCacheOrDirect(index, Order);
        double m12 = MomentFromCacheOrDirect(idx12, Order);
        double m13 = MomentFromCacheOrDirect(idx13, Order);
        double m23 = MomentFromCacheOrDirect(idx23, Order);
        double m1 = MomentFromCacheOrDirect(idx1, Order);
        double m2 = MomentFromCacheOrDirect(idx2, Order);
        double m3 = MomentFromCacheOrDirect(idx3, Order);

        free(idx1);
        free(idx2);
        free(idx3);
        free(idx12);
        free(idx13);
        free(idx23);

        if(m123 == myEmpty || m12 == myEmpty || m13 == myEmpty || m23 == myEmpty ||
           m1 == myEmpty || m2 == myEmpty || m3 == myEmpty)
        {
            return myEmpty;
        }
        return m123 - m12 * m3 - m13 * m2 - m23 * m1 + 2.0 * m1 * m2 * m3;
    }

    return myEmpty;
}

static bool IsLegacyBellRootIndex(const int* index, int Order)
{
    if(index == NULL || Order <= 0)
    {
        return false;
    }
    for(int i = 0; i < Order; i++)
    {
        if(index[i] != i + 1)
        {
            return false;
        }
    }
    return true;
}

static double LegacyBellBlockMoment(const int* labels, int blockSize, int fullOrder)
{
    int* multiplicity = (int*)calloc(fullOrder + 1, sizeof(int));
    if(multiplicity == NULL)
    {
        perror("Error in function LegacyBellBlockMoment: calloc return NULL\n");
        return myEmpty;
    }

    multiplicity[0] = blockSize;
    for(int i = 0; i < blockSize; i++)
    {
        int label = labels[i];
        if(label < 1 || label > fullOrder)
        {
            fprintf(stderr, "Error in function LegacyBellBlockMoment: label %d out of range [1, %d]\n", label, fullOrder);
            free(multiplicity);
            return myEmpty;
        }
        multiplicity[label]++;
    }

    double result = BellMomentValueOnly(multiplicity, fullOrder);
    free(multiplicity);
    return result;
}

static void BellBacktrackAllPartitions(int* arr, int n, int k, int index, Partition* current, int used_blocks)
{
    if(index == n)
    {
        if(used_blocks == k)
        {
            if(g_total_partitions < MAX_PARTITIONS)
            {
                g_all_partitions[g_total_partitions++] = *current;
            }
        }
        return;
    }

    for(int i = 0; i < used_blocks; i++)
    {
        int size = current->block_sizes[i];
        current->blocks[i][size] = arr[index];
        current->block_sizes[i]++;
        BellBacktrackAllPartitions(arr, n, k, index + 1, current, used_blocks);
        current->block_sizes[i]--;
    }

    if(used_blocks < k)
    {
        current->blocks[used_blocks][0] = arr[index];
        current->block_sizes[used_blocks] = 1;
        BellBacktrackAllPartitions(arr, n, k, index + 1, current, used_blocks + 1);
        current->block_sizes[used_blocks] = 0;
    }
}

static void GenerateBellAllPartitions(int* arr, int n, int k)
{
    g_total_partitions = 0;
    Partition init;
    memset(&init, 0, sizeof(init));
    BellBacktrackAllPartitions(arr, n, k, 0, &init, 0);
}

static double LegacyBellCumulantDirect(const int* index, int Order)
{
    double total = 0;
    for(int blockCount = 1; blockCount <= Order; blockCount++)
    {
        GenerateBellAllPartitions((int*)index, Order, blockCount);
        double partitionSum = 0;
        for(int partitionIndex = 0; partitionIndex < g_total_partitions; partitionIndex++)
        {
            double product = 1.0;
            for(int blockIndex = 0; blockIndex < blockCount; blockIndex++)
            {
                double blockMoment = LegacyBellBlockMoment(
                    g_all_partitions[partitionIndex].blocks[blockIndex],
                    g_all_partitions[partitionIndex].block_sizes[blockIndex],
                    Order);
                if(blockMoment == myEmpty)
                {
                    return myEmpty;
                }
                product *= blockMoment;
            }
            partitionSum += product;
        }
        partitionSum *= Factorial(blockCount - 1);
        if(blockCount % 2 != 0)
        {
            total += partitionSum;
        }
        else
        {
            total -= partitionSum;
        }
    }
    return total;
}

double BellCumulantFromMultiplicityIndex(const int* index, int Order)
{
    if(index == NULL)
    {
        fprintf(stderr, "Error in function BellCumulantFromMultiplicityIndex: index is NULL\n");
        return myEmpty;
    }
    if(Order <= 0)
    {
        return 0;
    }

    int jointOrder = index[0];
    if(jointOrder <= 0)
    {
        return 0;
    }

    int* labels = (int*)calloc(jointOrder, sizeof(int));
    if(labels == NULL)
    {
        perror("Error in function BellCumulantFromMultiplicityIndex: calloc return NULL\n");
        return myEmpty;
    }

    int cursor = 0;
    for(int variable = 1; variable <= Order; variable++)
    {
        if(index[variable] < 0)
        {
            fprintf(stderr, "Error in function BellCumulantFromMultiplicityIndex: negative multiplicity at %d\n", variable);
            free(labels);
            return myEmpty;
        }
        for(int repeat = 0; repeat < index[variable]; repeat++)
        {
            if(cursor >= jointOrder)
            {
                fprintf(stderr, "Error in function BellCumulantFromMultiplicityIndex: multiplicity sum exceeds joint order\n");
                free(labels);
                return myEmpty;
            }
            labels[cursor++] = variable;
        }
    }

    if(cursor != jointOrder)
    {
        fprintf(stderr, "Error in function BellCumulantFromMultiplicityIndex: multiplicity sum %d does not match joint order %d\n", cursor, jointOrder);
        free(labels);
        return myEmpty;
    }

    double total = 0;
    for(int blockCount = 1; blockCount <= jointOrder; blockCount++)
    {
        GenerateBellAllPartitions(labels, jointOrder, blockCount);
        double partitionSum = 0;
        for(int partitionIndex = 0; partitionIndex < g_total_partitions; partitionIndex++)
        {
            double product = 1.0;
            for(int blockIndex = 0; blockIndex < blockCount; blockIndex++)
            {
                double blockMoment = LegacyBellBlockMoment(
                    g_all_partitions[partitionIndex].blocks[blockIndex],
                    g_all_partitions[partitionIndex].block_sizes[blockIndex],
                    Order);
                if(blockMoment == myEmpty)
                {
                    free(labels);
                    return myEmpty;
                }
                product *= blockMoment;
            }
            partitionSum += product;
        }
        partitionSum *= Factorial(blockCount - 1);
        if(blockCount % 2 != 0)
        {
            total += partitionSum;
        }
        else
        {
            total -= partitionSum;
        }
    }

    free(labels);
    return total;
}




double BetterRecCumulant3(int* index, int Order)
{
    InitFactorSwitches();
    if(Order == 0 || index == NULL)
    {
        fprintf(stderr, "Error in function RecCumulant3: Wrong Parameters\n");
        return myEmpty;
    }

    #if 0
    if(Order >= 5)
    {
        double result = TrueCumulant(index, Order);
        memset(g_all_partitions, 0, sizeof(Partition) * MAX_PARTITIONS);
        g_total_partitions = 0;
        return result;
    }
    #endif

    //鍒濆鍖栧叏灞€鍙橀噺
    CurrCumulantNum = 0;
    //Floor2 = (CumulantIndex**)calloc(FstFloorNum, sizeof(CumulantIndex*));

    CumulantIndex* FstCumulant = (CumulantIndex*)calloc(1, sizeof(CumulantIndex));
    if(FstCumulant == NULL)
    {
        perror("Error in function RecCumulant3: calloc return NULL\n");
        return myEmpty;
    }
    FstCumulant->index = (int*)calloc(Order + 1, sizeof(int));
    if(FstCumulant->index == NULL)
    {
        free(FstCumulant);
        perror("Error in function RecCumulant3: calloc return NULL\n");
        return myEmpty;
    }
    memcpy(FstCumulant->index, index, sizeof(int) * (Order + 1));
    FstCumulant->Order = Order;

    //璁板綍绗竴灞傚眰鍏辨湁澶氬皯涓瓙绱Н閲?
    int NuminFstFloor = InitialCoeff(FstCumulant);
    //if(FstFloorNum < NuminFstFloor)
        //Floor2 = (CumulantIndex**)realloc(Floor2, (NuminFstFloor + 1) * sizeof(CumulantIndex*));
    if(NuminFstFloor == myEmpty)
    {
        fprintf(stderr, "Error in function RecCumulant3: Failed to Initial Cumulant\n");
        return myEmpty;
    }
    if(FstCumulant->Calculated == true)
    {
        return FstCumulant->value;
    }
    if(FstCumulant->nxtCumulant == NULL)
    {   
        printf("Value is %f\n", FstCumulant->value);
        fprintf(stderr, "Error in function RecCumulant3: this cumulant has no nxtCumulant\n");
        return myEmpty;
    }

    //鍒嗘壒杩涜瀛愮疮绉噺鐨勮绠楋紝璇ユ壒娆＄殑褰撳墠浣嶇疆锛屾瘡璁＄畻涓€涓瓙绱Н閲忓氨鍔犱竴
    int FstFloorStart = 0;
    //鍒嗘壒杩涜瀛愮疮绉噺鐨勮绠楋紝璇ユ壒娆＄殑鏈€鍚庝竴涓綅缃紝 FstFloorNum闄愬埗浜嗙涓€灞傛爤鐨勫閲?
    int FirstFloorBatch = g_factor_iterative_schedule == true ? DEFAULT_FST_FLOOR_NUM : NuminFstFloor;
    if(FirstFloorBatch <= 0)
    {
        FirstFloorBatch = NuminFstFloor;
    }
    int FstFloorEnd = FstFloorStart + FirstFloorBatch;
    //璁板綍姣忎竴鎵圭殑瀛愮疮绉噺鐩稿姞鐨勭粨鏋?
    double GroupResult = 0;
    //璁板綍鏄惁鏉ュ埌鏈€鍚庝竴涓壒娆?
    bool ReachEnd = false;
    //鎬绘爤
    CumulantIndex*** Floor;
    //璁板綍褰撳墠鎵€鍦ㄦ€绘爤鐨勫眰鏁?
    int CurrFloor;
    //璁板綍褰撳墠灞傜殑绱Н閲忔€绘暟
    int FloorCumNum[MaxFloor];
    //int TempTime = 0;
    while(true)
    {
        //TempTime++;
        //printf("%d\n", TempTime);
        if(FstFloorEnd >= NuminFstFloor)
        {
            ReachEnd = true;
            FstFloorEnd = NuminFstFloor;
        }  
        
        //ResultIndex = 0;

        //鍒濆鍖栨€绘爤
        Floor = (CumulantIndex***)calloc(MaxFloor, sizeof(CumulantIndex**));
        if(Floor == NULL)
        {
            perror("Error in function RecCumulant: Calloc return NULL\n");
            return myEmpty;
        }

        for(int i = 0; i < MaxFloor; i++)
        {
            Floor[i] = (CumulantIndex**)calloc(FstFloorEnd - FstFloorStart, sizeof(CumulantIndex*));
            if(Floor[i] == NULL)
            {
                perror("Error in function RecCumulant: Calloc return NULL\n");
                return myEmpty;
            }
        }

        CurrFloor = 0;
        //璁板綍褰撳墠澶勭悊鐨勭疮绉噺
        CumulantIndex* CurrCumulant = FstCumulant;
        //璁板綍褰撳墠澶勭悊绱Н閲忕殑绱㈠紩
        int CurrIndex = 0;
        //瀛樺偍褰撳墠灞傜殑鎵€鏈夊瓙绱Н閲?
        //鍒ゆ柇褰撳墠灞傛墍鏈夌疮绉噺鏄惁宸茬粡琚绠?
        bool AllTrue = true;
        int NuminCurrFloor = FstFloorEnd - FstFloorStart;
        for(int i = 0; i < NuminCurrFloor; i++)
        //for(int i = NuminCurrFloor - 1; i >= 0; i--)
        {
            Floor[CurrFloor][i] = &CurrCumulant->nxtCumulant[i + FstFloorStart];
            
            if(Floor[CurrFloor][i]->Calculated == true)
                continue;
            if(Floor[CurrFloor][i]->index == NULL)
            {
                printf("here %d %d", CurrFloor, i);
                getchar();
            }
            Node* CheckExist = SearchIndex2(Floor[CurrFloor][i]->index, Floor[CurrFloor][i]->Order);
            if(CheckExist == NULL)
            {
                //鐩存帴鍚戞爲涓彃鍏ョ储寮?
                RootAVL2[CurrCumulant->nxtCumulant[i + FstFloorStart].Order] = insert(RootAVL2[CurrCumulant->nxtCumulant[i + FstFloorStart].Order], CurrCumulant->nxtCumulant[i + FstFloorStart].index, CurrCumulant->nxtCumulant[i + FstFloorStart].value, CurrCumulant->nxtCumulant[i + FstFloorStart].Order, CurrCumulant->nxtCumulant[i + FstFloorStart].Calculated);
                //fgetc(stdin);
            }
            else if(g_factor_avl_reuse == true)
            {
                Floor[CurrFloor][i]->exist = true;
                //fgetc(stdin);
            }
            if(CurrCumulant->nxtCumulant[i + FstFloorStart].Calculated == false)
            {
                AllTrue = false;
            }
        }
        
        //璁板綍褰撳墠灞傦紙绗竴灞傦級鍏辨湁澶氬皯涓疮绉噺
        FloorCumNum[CurrFloor] = NuminCurrFloor;

        if(AllTrue == true)
        {
            //褰撳墠灞傛墍鏈夊瓙绱Н閲忛兘宸茬粡琚绠?
            //姝ゅ鍙互鐩存帴璁＄畻鍑烘渶缁堢粨鏋?
            //涓€鑸笉鍙兘鍑虹幇璇ユ儏鍐?
            printf("Warning in function RecCumulant: May have something wrong\n");
            if(ReachEnd == true)
            {
                FstCumulant->value = 0;
                FstCumulant->Calculated = true;
                goto FreeStack;
            }
            else
            {
                FstFloorStart = FstFloorEnd;
                FstFloorEnd += FirstFloorBatch;
                for(int i = 0; i < FloorCumNum[CurrFloor]; i++)
                {
                    if(FreeCumulantIndex(&Floor[CurrFloor][i]) == myEmpty)
                    {
                        fprintf(stderr, "Error in function: BetterRecCumulant3: Cant free Cumulant Index\n");
                        return myEmpty;
                    }
                                
                }
                for(int i = 0; i < MaxFloor; i++)
                {
                    free(Floor[i]);
                }
                free(Floor);
                continue;
            }
            //return 0;
        }
        //璁板綍褰撳墠绱Н閲忔湁澶氬皯涓瓙绱Н閲忥紙鍗曞眰锛?
        int NuminFloor[MaxCumulantNum];

      
        //鐖剁疮绉噺鐨勭粨鏋勪綋鍦ㄤ笂涓€灞傜殑绱㈠紩
        int cumulantIndex = 0;
        int nxtStart = 0;
        CurrFloor++;
        //NuminCurrFloor = 0;
        int SubCumulantNum = 0;
        //鍒濆鍖栫涓€灞傛爤锛屽嵆瀵圭涓€灞傚瓙绱Н閲忚皟鐢╥nit
        while(true)
        {
            //鏉ュ埌寰幆鍑哄彛锛屾洿鏂癗uminCurrFloor
            if(cumulantIndex >= NuminCurrFloor)
            {
                //鏇存柊绗簩灞傛爤涓厓绱犵殑鏁伴噺
                FloorCumNum[CurrFloor] = SubCumulantNum;
                //printf("Current Cumulant number: %d\n", NuminCurrFloor);
                break;
            }
            //鍒ゆ柇涓婁竴灞傜殑绱Н閲忔槸鍚﹀凡缁忚璁＄畻鍑虹粨鏋滐紙鍦ㄦ湭璋冪敤InitialCoeff涔嬪墠灏辩畻鍑虹粨鏋滐紝璇存槑绱Н閲忓€间负0锛夛紝鑻ユ槸鍒欒烦杩?
            if(Floor[CurrFloor - 1][cumulantIndex]->Calculated != true && Floor[CurrFloor - 1][cumulantIndex]->exist == false)
                CurrCumulant = Floor[CurrFloor - 1][cumulantIndex];
            else
            {
                //NuminFloor鐢ㄤ簬璁板綍褰撳墠鐖剁疮绉噺鏈夊灏戜釜瀛愮疮绉噺锛屾渶鍚庡皢鍏剁疮鍔犲緱鍒颁笅涓€灞傛爤鐨勭疮绉噺鎬绘暟
                //瀵逛簬宸茬粡涓簍rue鐨勭埗绱Н閲忚€岃█锛屽叾瀛愮疮绉噺鏁伴噺涓?
                NuminFloor[cumulantIndex] = 0;
                //RootAVL2[Floor[CurrFloor - 1][cumulantIndex]->Order] = insert(RootAVL2[Floor[CurrFloor - 1][cumulantIndex]->Order], Floor[CurrFloor - 1][cumulantIndex]->index, Floor[CurrFloor - 1][cumulantIndex]->value, Floor[CurrFloor - 1][cumulantIndex]->Order);
                cumulantIndex++;
                continue;
            }

            NuminFloor[cumulantIndex] = InitialCoeff(CurrCumulant);

            if(NuminFloor[cumulantIndex] == myEmpty)
            {
                printf("here %d %d", CurrFloor-1, cumulantIndex);
                getchar();
                fprintf(stderr, "Error in function RecCumulant3: Failed to Initial Cumulant\n");
                return myEmpty;
            }
            //鍒ゆ柇涓婁竴灞傜殑绱Н閲忔槸鍚﹀凡缁忚璁＄畻鍑虹粨鏋滐紝鑻ユ槸鍒欒烦杩? 杩欓噷閽堝灏忎簬鍥涢樁锛屽ぇ浜庨浂闃剁殑绱Н閲?
            if(CurrCumulant->Calculated == true)
            {
                NuminFloor[cumulantIndex] = 0;
                //RootAVL2[CurrCumulant->Order] = insert(RootAVL2[CurrCumulant->Order], CurrCumulant->index, CurrCumulant->value, CurrCumulant->Order);
                cumulantIndex++;
                //printf("%d\n", cumulantIndex);
                Node* TempNode = SearchIndex2(CurrCumulant->index, CurrCumulant->Order);
                if(TempNode == NULL)
                {
                    fprintf(stderr, "Error in function BetterRecCumulant3: Cant find index\n");
                    return myEmpty;
                }
                TempNode->value = CurrCumulant->value;
                TempNode->existed = true;
                continue;
            }
            
            

            SubCumulantNum += NuminFloor[cumulantIndex];
            //if(SubCumulantNum > FstFloorNum)
            if(SubCumulantNum > FstFloorEnd - FstFloorStart)
            {
                //鍔ㄦ€佸垎閰嶆爤
                Floor[CurrFloor] = (CumulantIndex**)realloc(Floor[CurrFloor], SubCumulantNum * sizeof(CumulantIndex*));
                if(Floor[CurrFloor] == NULL)
                {
                    perror("Error in function BetterRecCumulat3: realloc return NULL\n");
                    return myEmpty;
                }
            }
            for(int i = nxtStart; i < NuminFloor[cumulantIndex] + nxtStart; i++)
            //for(int i = NuminFloor[cumulantIndex] + nxtStart - 1; i >= nxtStart; i--)
            {
                //printf("nxtstart is %d\n", nxtStart);
                Floor[CurrFloor][i] = &CurrCumulant->nxtCumulant[i - nxtStart];
                //Floor[CurrFloor][i]->lstCumulant = CurrCumulant;

                if(i + 1 == NuminFloor[cumulantIndex] + nxtStart)
                {
                    //鏇存柊涓嬩竴涓埗绱Н閲忕殑绗竴涓瓙绱閲忓湪涓嬩竴灞傛爤鐨勮捣濮嬩綅缃?
                    nxtStart = i + 1;
                    break;
                }
            }
            cumulantIndex++;
        }
        
        //璁板綍褰撳墠鎵规鏄惁璁＄畻瀹屾瘯
        bool TempFinish = false;

        //瀛樺偍褰撳墠绱姞涓殑涓€椤圭粨鏋?
        //鑷虫宸插畬鎴愭鍚戝惊鐜?
        //涓嬮潰寮€濮嬩粠椤跺悜涓嬪紑濮嬮€嗗悜寰幆
        int testI = 0;
        while(true)
        {
            //涓婂眰绱Н閲忕储寮?
            CurrIndex = 0;
            //涓嬪眰绯绘暟绱㈠紩
            int CurrCumulantIndex = 0;
            if(CurrFloor == 0 && GroupResult != 0)
            {
                CurrCumulantIndex = FstFloorStart;
            }
            double CurrValue = 0;
            //鐢ㄤ簬閲婃斁绗簩灞傚瓙绱Н閲?
            int CurrPos = 0;
            int LstPos = 0;
            
            int TempNum = 0;

            while(true)
            {
                if(Floor[CurrFloor][CurrIndex] == NULL)
                {
                    CurrFloor--;
                    if(GroupResult != 0)
                        CurrCumulantIndex = FstFloorStart;
                }
                CurrCumulant = Floor[CurrFloor][CurrIndex];
                //閲婃斁绗竴灞傜殑涓嬩竴灞傚瓙绱Н閲?
                if(CurrCumulant->nxtCumulant != NULL)
                {
                    for(; CurrPos < NuminFloor[CurrIndex] + LstPos; CurrPos++)
                    {
                        if(FreeCumulantIndex(&Floor[CurrFloor+1][CurrPos]) == myEmpty)
                        {
                            fprintf(stderr, "Error in function RecCumulant3: Cant free Stack\n");
                            return myEmpty;
                        }
                    }
                    free(CurrCumulant->nxtCumulant);
                    LstPos = CurrPos;
                }
                if(CurrCumulant == NULL)
                {
                    return myEmpty;
                }
                else if(CurrCumulant->lstCumulant == NULL)
                {
                    return myEmpty;
                }
                else if(CurrCumulant->lstCumulant->moment == NULL)
                {
                    return myEmpty;
                }
                else if(CurrCumulant->lstCumulant->combine == NULL)
                {
                    return myEmpty;
                }

                if(CurrCumulant->Calculated == false)
                {
                    Node* CheckExist = SearchIndex2(CurrCumulant->index, CurrCumulant->Order);
                    if(CheckExist == NULL)
                    {
                        printf("Curr Order is %d\n", CurrCumulant->Order);
                        printf("Curr index is %d\n", CurrCumulant->index[Order-1]);
                        fprintf(stderr, "Error in function RecCumulant3:1208 Cant find index in AVLTree\n");
                        return myEmpty;
                    }
                    
                    CurrCumulant->value = CheckExist->value;
                    if(CheckExist->existed == false)
                    {
                        fprintf(stderr, "Error in function BetterRecCumulant3: Value has not been recorded in tree\n");
                        return myEmpty;
                    }
                    else
                        CurrCumulant->Calculated = true;
                }
                

                
                CurrValue += CurrCumulant->lstCumulant->moment[CurrCumulantIndex] * CurrCumulant->lstCumulant->combine[CurrCumulantIndex] * CurrCumulant->value;
                testI++;
                
                #if 0
                printf("CurrFloor: %d ", CurrFloor);
                printf("Lst Index: ");
                for(int testJ = 1; testJ <= CurrCumulant->lstCumulant->Order; testJ++)
                {
                    printf("%d ", CurrCumulant->lstCumulant->index[testJ]);
                }
                printf(" ");
                printf("Curr cumulant Index: ");
                for(int testJ = 1; testJ <= CurrCumulant->Order; testJ++)
                {
                    printf("%d ", CurrCumulant->index[testJ]);
                }
                printf(" ");
                printf("%d: Moment:%f  combine:%f CurrCumulant:%f CurrValue:%f\n", testI, CurrCumulant->lstCumulant->moment[CurrCumulantIndex], CurrCumulant->lstCumulant->combine[CurrCumulantIndex], CurrCumulant->value, CurrValue);
                #endif

                CurrIndex++;
                CurrCumulantIndex++;
                
                //鏈眰宸茶绠楀畬姣?
                if((CurrIndex == FloorCumNum[CurrFloor] && CurrFloor != 0) || (CurrIndex == FstFloorEnd - FstFloorStart && CurrFloor == 0))
                {
                    CurrCumulant->lstCumulant->Calculated = true;

                    //int MomentIndex = IndexRank(CurrCumulant->lstCumulant->index, CurrCumulant->lstCumulant->Order);
                    //double except = G_Moments[CurrCumulant->lstCumulant->Order][MomentIndex];
                    #if 0
                    int* tempIndex = NULL;
                    if(CurrFloor != 0)
                    {
                        if(CurrCumulant->lstCumulant->lstCumulant == NULL)
                        {
                            fprintf(stderr, "Error in function BetterRecCumulant3: there is no lst-lst cumulant\n");
                            return myEmpty;
                        }
                        int lst_lst_order = CurrCumulant->lstCumulant->lstCumulant->Order;
                        int lst_order = CurrCumulant->lstCumulant->Order;
                        tempIndex = (int*)calloc(lst_lst_order+1, sizeof(int));
                        if(tempIndex == NULL)
                        {
                            perror("Error in function BetterRecCumulant3 : Calloc return NULL\n");
                            return myEmpty;
                        }
                        printf("Lst Cumulant moment index: ");
                        for(int i = 1; i <=lst_lst_order; i++)
                        {
                            if(i > lst_order)
                                tempIndex[i] = CurrCumulant->lstCumulant->lstCumulant->index[i];
                            else
                            {
                                tempIndex[i] = CurrCumulant->lstCumulant->lstCumulant->index[i] - CurrCumulant->lstCumulant->index[i];
                            }
                            printf("%d ", tempIndex[i]);
                        }
                    }
                    
                    int* myTempIndex = CurrCumulant->lstCumulant->index;
                    int myTempOrder = CurrCumulant->lstCumulant->Order;
                    if(tempIndex != NULL)
                    {
                        myTempIndex = tempIndex;
                        myTempOrder = CurrCumulant->lstCumulant->lstCumulant->Order;
                    }
                    //鏂扮殑鏌ユ壘绱㈠紩鐨勫€煎搴旂殑鐭?
                    bool fouded = false;
                    int indexPos;
                    for(indexPos = 0; indexPos < G_SeriesNum; indexPos++)
                    {
                        
                        if(myStrcmp(myTempIndex, myTempOrder, G_MomentIndex[indexPos]) == true)
                        {
                            fouded = true;
                            break;
                        }
                    }
                    if(fouded == false)
                    {
                        fprintf(stderr, "Failed to find moment index\n");
                        return myEmpty;
                    }
                    double except = G_Moment[indexPos];
                    printf("value is: %f\n", except);

                    CurrCumulant->lstCumulant->value = except - CurrValue;

                    if(tempIndex != NULL)
                    {
                        free(tempIndex);
                        tempIndex = NULL;
                    }
                    #endif


                    //鏂扮殑鏌ユ壘绱㈠紩鐨勫€煎搴旂殑鐭?
                    int* myTempIndex = CurrCumulant->lstCumulant->index;
                    int myTempOrder = CurrCumulant->lstCumulant->Order;
                    int indexPos = MomentLookupIndex(myTempIndex, myTempOrder);
                    double except = myEmpty;
                    if(indexPos < 0)
                    {
                        except = Moment5ValueOnly(myTempIndex, myTempOrder);
                    }
                    else
                    {
                        except = G_Moment[indexPos];
                    }
                    if(except == myEmpty)
                    {
                        fprintf(stderr, "Failed to calculate moment 1\n");
                        return myEmpty;
                    }
                    //printf("value is: %f\n", except);

                    CurrCumulant->lstCumulant->value = except - CurrValue;
                    CurrCumulant->lstCumulant->Calculated = true;



                    if(CurrFloor != 0)
                    {
                    Node* TempNode = SearchIndex2(CurrCumulant->lstCumulant->index, CurrCumulant->lstCumulant->Order);
                if(TempNode == NULL)
                {
                    fprintf(stderr, "Error in function BetterRecCumulant3: Cant find index\n");
                    return myEmpty;
                }
                TempNode->value = CurrCumulant->lstCumulant->value;
                TempNode->existed = true;
            }





                    //鏇存柊AVL鏍戜腑鐨勫€?
                    if(CurrFloor != 0)
                    {
                        Node* RefreshNode = SearchIndex2(CurrCumulant->lstCumulant->index, CurrCumulant->lstCumulant->Order);
                        if(RefreshNode == NULL)
                        {
                            fprintf(stderr, "Error in function BetterRecCumulant3:1238: Failed to refresh value in AVL\n");
                            return myEmpty;
                        }
                        RefreshNode->value = CurrCumulant->lstCumulant->value;
                        RefreshNode->existed = true;
                        CurrFloor--;
                        break;
                    }
                    else
                    {
                        //褰撳墠灞傛墍鏈夌疮绉噺閮藉搴斾簬涓嬩竴灞傜殑涓€涓疮绉噺锛岃鏄庡凡缁忚绠楀埌浜嗗簳灞?
                        if(ReachEnd == true && GroupResult == 0)
                        {
                            FstCumulant->Calculated = true;
                            
                            //int MomentIndex = IndexRank(FstCumulant->index, FstCumulant->Order);
                            //double except = G_Moments[FstCumulant->Order][MomentIndex];
                            
                            //鏂扮殑鏌ユ壘绱㈠紩鐨勫€煎搴旂殑鐭?
                            double except1 = RootMomentOrDirect(FstCumulant->index, FstCumulant->Order);
                            
                            FstCumulant->value = except1 - CurrValue;
                            
                            goto FreeStack;
                        }
                        //瀹屾垚鎵€鏈夋壒娆＄疮璁￠噺鐨勮绠?
                        else if(ReachEnd == true && GroupResult != 0)
                        {
                            GroupResult += CurrValue;
                            FstCumulant->Calculated = true;
                            
                            //int MomentIndex = IndexRank(FstCumulant->index, FstCumulant->Order);
                            //double except = G_Moments[FstCumulant->Order][MomentIndex];
                            
                            //鏂扮殑鏌ユ壘绱㈠紩鐨勫€煎搴旂殑鐭?
                            double except2 = RootMomentOrDirect(FstCumulant->index, FstCumulant->Order);
                            
                            FstCumulant->value = except2 - GroupResult;

                            goto FreeStack;
                        }
                        //浠呭畬鎴愪竴鎵圭疮绉噺鐨勮绠?
                        else
                        {
                            GroupResult += CurrValue;
                            TempFinish = true;
                            FstFloorStart = FstFloorEnd;
                            FstFloorEnd += FirstFloorBatch;
                            for(int i = 0; i < FloorCumNum[CurrFloor]; i++)
                            {
                                if(FreeCumulantIndex(&Floor[CurrFloor][i]) == myEmpty)
                                {
                                    fprintf(stderr, "Error in function:1279: BetterRecCumulant3: Cant free Cumulant Index\n");
                                    return myEmpty;
                                }
                                
                            }
                            for(int i = 0; i < MaxFloor; i++)
                            {
                                free(Floor[i]);
                            }
                            free(Floor);
                            break;
                        }
                    }
                }

                //if(TempFinish == true)
                //{
                    //break;
                //}

                //涓婂眰鐨勪竴涓疮绉噺宸茶绠楀畬姣?
                if(Floor[CurrFloor][CurrIndex]->lstCumulant != CurrCumulant->lstCumulant)
                {
                    //閫氳繃涓€绯诲垪瀛愮疮绉噺瀹屾垚涓婁竴灞傜殑鐨勪竴涓疮绉噺璁＄畻
                    //杩涘叆鍚屼竴灞備笅涓€涓疮绉噺鐨勫瓙绱Н閲?
                    CurrCumulantIndex = 0;
                    CurrCumulant->lstCumulant->Calculated = true;
                    
                    //int MomentIndex = IndexRank(CurrCumulant->lstCumulant->index, CurrCumulant->lstCumulant->Order);
                    //double except = G_Moments[CurrCumulant->lstCumulant->Order][MomentIndex];
                    
                    #if 0
                    //鏇挎崲涓虹煩鐨勭储寮?
                    int* tempIndex = NULL;
                    if(CurrFloor != 0)
                    {
                        if(CurrCumulant->lstCumulant->lstCumulant == NULL)
                        {
                            fprintf(stderr, "Error in function BetterRecCumulant3: there is no lst-lst cumulant\n");
                            return myEmpty;
                        }
                        int lst_lst_order = CurrCumulant->lstCumulant->lstCumulant->Order;
                        int lst_order = CurrCumulant->lstCumulant->Order;
                        tempIndex = (int*)calloc(lst_lst_order+1, sizeof(int));
                        if(tempIndex == NULL)
                        {
                            perror("Error in function BetterRecCumulant3 : Calloc return NULL\n");
                            return myEmpty;
                        }
                        printf("Lst Cumulant moment index: ");
                        for(int i = 1; i <=lst_lst_order; i++)
                        {
                            if(i > lst_order)
                                tempIndex[i] = CurrCumulant->lstCumulant->lstCumulant->index[i];
                            else
                            {
                                tempIndex[i] = CurrCumulant->lstCumulant->lstCumulant->index[i] - CurrCumulant->lstCumulant->index[i];
                            }
                            printf("%d ", tempIndex[i]);
                        }
                        //printf("\n");
                    }
                    
                    int* myTempIndex = CurrCumulant->lstCumulant->index;
                    int myTempOrder = CurrCumulant->lstCumulant->Order;
                    if(tempIndex != NULL)
                    {
                        myTempIndex = tempIndex;
                        myTempOrder = CurrCumulant->lstCumulant->lstCumulant->Order;
                    }

                    //鏂扮殑鏌ユ壘绱㈠紩鐨勫€煎搴旂殑鐭?
                    bool fouded = false;
                    int indexPos;
                    for(indexPos = 0; indexPos < G_SeriesNum; indexPos++)
                    {
                        if(myStrcmp(myTempIndex, myTempOrder, G_MomentIndex[indexPos]) == true)
                        {
                            fouded = true;
                            break;
                        }
                    }
                    if(fouded == false)
                    {
                        fprintf(stderr, "Failed to find moment index\n");
                        return myEmpty;
                    }
                    double except = G_Moment[indexPos];
                    printf("value is: %f\n", except);
                    
                    CurrCumulant->lstCumulant->value = except - CurrValue;

                    if(tempIndex != NULL)
                    {
                        free(tempIndex);
                        tempIndex = NULL;
                    }
                    #endif



                    //鏂扮殑鏌ユ壘绱㈠紩鐨勫€煎搴旂殑鐭?
                    int* myTempIndex = CurrCumulant->lstCumulant->index;
                    int myTempOrder = CurrCumulant->lstCumulant->Order;
                    int indexPos = MomentLookupIndex(myTempIndex, myTempOrder);
                    double except = myEmpty;
                    if(indexPos < 0)
                    {
                        except = Moment5ValueOnly(myTempIndex, myTempOrder);
                    }
                    else
                    {
                        except = G_Moment[indexPos];
                    }
                    if(except == myEmpty)
                    {
                        fprintf(stderr, "Failed to calculate moment 2\n");
                        return myEmpty;
                    }
                    //printf("value is: %f\n", except);
                    
                    CurrCumulant->lstCumulant->value = except - CurrValue;
                    CurrCumulant->lstCumulant->Calculated = true;




                    



                    //鏇存柊AVL鏍戜腑鐨勫€?
                    Node* RefreshNode = SearchIndex2(CurrCumulant->lstCumulant->index, CurrCumulant->lstCumulant->Order);
                    if(RefreshNode == NULL)
                    {
                        fprintf(stderr, "Error in function BetterRecCumulant3:1301 Failed to refresh value in AVL\n");
                        return myEmpty;
                    }
                    else
                    {
                        RefreshNode->value = CurrCumulant->lstCumulant->value;
                        RefreshNode->existed = true;
                    }
                    TempNum++;
                    CurrValue = 0;
                }
            }
            
            if(TempFinish == true)
            {
                break;
            }
            
        }
    }
    FreeStack:
    for(int i = 0; i < FloorCumNum[CurrFloor]; i++)
    {
        if(FreeCumulantIndex(&Floor[CurrFloor][i]) == myEmpty)
        {
            fprintf(stderr, "Error in function: BetterRecCumulant3: Cant free Cumulant Index\n");
            return myEmpty;
        }
        
    }
    
    for(int i = 0; i < MaxFloor; i++)
    {
        free(Floor[i]);
    }
    free(Floor);

    #if 0
    for(int i = 0; i < CurrCumulantNum; i++)
    {
        free(Floor2[i]);
    }
    free(Floor2);
    #endif
    
    //free(G_Template->index);
    //G_Moments = NULL;
    //G_Template->index = NULL;
    FreeCumulantIndex(&FstCumulant);
    free(FstCumulant->nxtCumulant);
    if(FstCumulant->Calculated == true)
    {
        double value = FstCumulant->value;
        free(FstCumulant);
        //TraverseTree(RootAVL2[4]);
        return value;
    }
    return myEmpty;
}

double LowOrderCumulant(int* index, int Order)
{
    RemoveZero(&index, &Order);
    //if(Order == 0 || index[0] == 1)
    if(Order == 0)
    {
        return 0;
    }
    if(g_low_order_uses_moment_cache == true)
    {
        return LowOrderCumulantFromMoments(index, Order);
    }
    for(int i = 0; i < G_TemplateNum; i++)
    {
        double TempProduct = 1;
        for(int j = 1; j <= Order; j++)
        {
            TempProduct *= pow(G_Template[i].Value[j], (double) index[j]);
        }
        G_Template[i].Value[0] = TempProduct;
    }
    double TempMean = 0;
    for(int i = 0; i < G_TemplateNum; i++)
    {
        TempMean += G_Template[i].Value[0];
        
        //瀹夊叏璧疯灏咷_Templte[i].value[0]閲嶆柊璁剧疆涓?
        G_Template[i].Value[0] = 0;
    }
    TempMean /= G_TemplateNum;
    return TempMean;
}

int PrecomputeMomentsForIndex(const int* index, int Order)
{
    if(index == NULL || Order <= 0)
    {
        fprintf(stderr, "Error in function PrecomputeMomentsForIndex: Wrong parameters\n");
        return myEmpty;
    }
    if(G_MomentFull == true)
    {
        return 0;
    }
    int* limits = (int*)calloc(Order + 1, sizeof(int));
    double** powerTable = (double**)calloc(Order + 1, sizeof(double*));
    if(limits == NULL || powerTable == NULL)
    {
        perror("Error in function PrecomputeMomentsForIndex: calloc return NULL\n");
        free(limits);
        free(powerTable);
        return myEmpty;
    }

    int formulaNum = 1;
    for(int j = 1; j <= Order; j++)
    {
        limits[j] = index[j];
        if(j == Order)
        {
            limits[j] = index[j] - 1;
        }
        formulaNum *= (limits[j] + 1);
    }
    formulaNum--;

    int totalNum = formulaNum + 1;
    G_Moment = (double*)calloc(totalNum, sizeof(double));
    G_MomentIndex = (char**)calloc(totalNum, sizeof(char*));
    if(G_Moment == NULL || G_MomentIndex == NULL)
    {
        perror("Error in function PrecomputeMomentsForIndex: calloc return NULL\n");
        free(limits);
        free(powerTable);
        return myEmpty;
    }
    G_SeriesNum = totalNum;
    G_CurrMomentPos = 0;

    for(int j = 1; j <= Order; j++)
    {
        int maxExp = index[j];
        int tableSize = (maxExp + 1) * G_TemplateNum;
        powerTable[j] = (double*)calloc(tableSize, sizeof(double));
        if(powerTable[j] == NULL)
        {
            perror("Error in function PrecomputeMomentsForIndex: calloc return NULL\n");
            for(int k = 1; k < j; k++)
            {
                free(powerTable[k]);
            }
            free(powerTable);
            free(limits);
            return myEmpty;
        }
        for(int i = 0; i < G_TemplateNum; i++)
        {
            powerTable[j][i] = 1.0;
        }
        for(int e = 1; e <= maxExp; e++)
        {
            int currOff = e * G_TemplateNum;
            int prevOff = (e - 1) * G_TemplateNum;
            for(int i = 0; i < G_TemplateNum; i++)
            {
                powerTable[j][currOff + i] = powerTable[j][prevOff + i] * G_Template[i].Value[j];
            }
        }
    }

    {
        double sum = 0;
        for(int i = 0; i < G_TemplateNum; i++)
        {
            double prod = 1.0;
            for(int j = 1; j <= Order; j++)
            {
                prod *= powerTable[j][index[j] * G_TemplateNum + i];
            }
            sum += prod;
        }
        G_Moment[0] = sum / G_TemplateNum;
        G_MomentIndex[0] = (char*)calloc(Order + 1, sizeof(char));
        if(G_MomentIndex[0] == NULL || BuildMomentKey(index, Order, G_MomentIndex[0], Order + 1) != 0)
        {
            fprintf(stderr, "Error in function PrecomputeMomentsForIndex: Failed to build root index\n");
            for(int j = 1; j <= Order; j++) free(powerTable[j]);
            free(powerTable);
            free(limits);
            return myEmpty;
        }
    }
    G_CurrMomentPos = 1;

    if(formulaNum > 0)
    {
        int radix = 1;
        int* base = (int*)calloc(Order + 1, sizeof(int));
        if(base == NULL)
        {
            perror("Error in function PrecomputeMomentsForIndex: calloc return NULL\n");
            for(int j = 1; j <= Order; j++) free(powerTable[j]);
            free(powerTable);
            free(limits);
            return myEmpty;
        }
        for(int j = 1; j <= Order; j++)
        {
            base[j] = limits[j] + 1;
            radix *= base[j];
        }

        int allocFailed = 0;
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) reduction(max:allocFailed)
        #endif
        for(int code = 1; code < radix; code++)
        {
            int* localExp = (int*)calloc(Order + 1, sizeof(int));
            if(localExp == NULL)
            {
                allocFailed = 1;
                continue;
            }
            localExp[0] = 0;
            int rem = code;
            for(int j = Order; j >= 1; j--)
            {
                localExp[j] = rem % base[j];
                rem /= base[j];
            }

            double sum = 0;
            for(int i = 0; i < G_TemplateNum; i++)
            {
                double prod = 1.0;
                for(int j = 1; j <= Order; j++)
                {
                    prod *= powerTable[j][localExp[j] * G_TemplateNum + i];
                }
                sum += prod;
            }

            G_Moment[code] = sum / G_TemplateNum;
            G_MomentIndex[code] = (char*)calloc(Order + 1, sizeof(char));
            if(G_MomentIndex[code] != NULL)
            {
                BuildMomentKey(localExp, Order, G_MomentIndex[code], Order + 1);
            }
            else
            {
                allocFailed = 1;
            }
            free(localExp);
        }
        if(allocFailed != 0)
        {
            fprintf(stderr, "Error in function PrecomputeMomentsForIndex: calloc return NULL\n");
            free(base);
            for(int j = 1; j <= Order; j++) free(powerTable[j]);
            free(powerTable);
            free(limits);
            return myEmpty;
        }
        G_CurrMomentPos = radix;
        free(base);
    }

    for(int j = 1; j <= Order; j++)
    {
        free(powerTable[j]);
    }
    free(powerTable);
    free(limits);

    G_MomentFull = true;
    InitFactorSwitches();
    if(g_factor_hash_lookup == true)
    {
        BuildMomentHash();
    }
    return 0;
}

int InitialCoeff(CumulantIndex* CurrStruct)
{
    InitFactorSwitches();
    if(CurrStruct == NULL)
    {
        fprintf(stderr, "Error in function InitialCoeff: Parameter is NULL\n");
        return myEmpty;
    }
    //澶栭儴搴斿皢CurrStruct鐨刬ndex鍜孫rder澶勭悊瀹屾瘯
    if(CurrStruct->index == NULL)
    {
        fprintf(stderr, "Error in function InitialCoeff: Wrong parameter\n");
        return myEmpty;
    }

    //鍏堝浼犺繘鏉ョ殑绱Н閲忕储寮曡繘琛屾帓搴忋€佸幓闆跺鐞嗐€傝嫢绱Н閲忛樁鏁颁负0锛屽垯杩斿洖0
    if(RemoveZero(&CurrStruct->index, &CurrStruct->Order) == -1 || CurrStruct->Order == 0)
    {
        CurrStruct->value = 0;
        CurrStruct->Calculated = true;
        return 0;
    }

    //灏嗗垵濮嬬疮绉噺璁＄畻鏍囩缃负false
    CurrStruct->Calculated = false;
    
    int* index = CurrStruct->index;
    int Order = CurrStruct->Order;

    if(index[0] <= 3)
    {
            CurrStruct->value = LowOrderCumulant(index, Order);
            CurrStruct->Calculated = true;
            return 0;
    }

    //鍒濆鍖栫储寮曚俊鎭?
    int* indexJ = (int*)calloc(Order+1, sizeof(int));
    if(indexJ == NULL)
    {
        perror("Error in function InitialCoeff 379: calloc return NULL\n");
        system("PAUSE");
        return myEmpty;
    }
    indexJ[Order] = 1;

    int** TempIndex = (int**)calloc(Order+1, sizeof(int*));
    if(TempIndex == NULL)
    {
        perror("Error in function InitialCoeff 388: calloc return NULL\n");
        system("PAUSE");
        return myEmpty;
    }
    //printf("indexi: ");
    for(int i = 0; i < Order+1; i++)
    {
        TempIndex[i] = (int*)calloc(2, sizeof(int));
        if(TempIndex[i] == NULL)
        {
            perror("Error in function InitialCoeff 397: calloc return NULL\n");
            system("PAUSE");
            for(int j = 0; j < i; j++)
            {
                free(TempIndex[j]);
                TempIndex[j] = NULL;
            }
            free(TempIndex);
            TempIndex = NULL;
            return myEmpty;
        }
        if(i == Order)
        {
            TempIndex[i][0] = index[i] - 1;
            indexJ[i] = 0;
            for(int j = i; j >= 1; j--)
            {
                if(TempIndex[j][0] == 0)
                    continue;
                else
                {
                    indexJ[j] = 1;
                    break;
                }
            }
        }
        else
            TempIndex[i][0] = index[i];
        //printf("%d ", TempIndex[i][0]);
    }
    //printf("\n");

    //璁＄畻绱姞娆℃暟
    int formulaNum = 1;
    for(int i = 1; i <= Order; i++)
    {
        //formulaNum *= index[i] + 1 - indexJ[i];
        formulaNum *= TempIndex[i][0] + 1;
    }
    //formulaNum *= index[Order - 1] - indexJ[Order - 1] + 1;
    formulaNum--;
    
    if(g_factor_moment_precompute == true && G_Moment == NULL)
    {
        G_Moment = (double*)calloc(formulaNum+1+formulaNum, sizeof(double));
        G_MomentIndex = (char**)calloc(formulaNum+1+formulaNum, sizeof(char*));
        G_SeriesNum = formulaNum+1+formulaNum;
        if(G_Moment == NULL || G_MomentIndex == NULL)
        {
            perror("InitialCoeff : calloc return NULL\n");
            return myEmpty;
        }
        //G_Moment[G_CurrMomentPos] = Moment5(index, Order);
        //G_CurrMomentPos++;
    }
    if(g_factor_moment_precompute == true && G_MomentFull == false)
    {
        G_Moment[G_CurrMomentPos] = Moment5(index, Order);
        G_CurrMomentPos++;
    }

    //鍒濆鍖栧綋鍓嶇疮绉噺鐨勭粨鏋勪綋
    CurrStruct->nxtCumulant = (CumulantIndex*)calloc(formulaNum, sizeof(CumulantIndex));
    

    if(CurrStruct->nxtCumulant == NULL)
    {
        perror("Error in function InitialCoeff 427: calloc return NULL\n");
        system("PAUSE");
        return myEmpty;
    }
    CurrStruct->moment = (double*)calloc(formulaNum, sizeof(double));
    if(CurrStruct->moment == NULL)
    {
        free(CurrStruct->nxtCumulant);
        perror("Error in function InitialCoeff 435: calloc return NULL\n");
        system("PAUSE");
        return myEmpty;
    }

    //remove this
    CurrStruct->MomentOrder = (int*)calloc(formulaNum, sizeof(int));
    if(CurrStruct->MomentOrder == NULL)
    {
        free(CurrStruct->nxtCumulant);
        perror("Error in function InitialCoeff 435: calloc return NULL\n");
        system("PAUSE");
        return myEmpty;
    }

    CurrStruct->combine = (double*)calloc(formulaNum, sizeof(double));
    if(CurrStruct->combine == NULL)
    {
        free(CurrStruct->nxtCumulant);
        free(CurrStruct->moment);
        free(CurrStruct->MomentOrder);
        perror("Error in function InitialCoeff: calloc return NULL\n");
        system("PAUSE");
        return myEmpty;
    }


    //绱姞杩囩▼鍋滄鏍囪
    bool AccumulateStop = false;
    //绱姞杩囩▼涓储寮曪紝璁板綍绱姞娆℃暟
    //int AccumulateTime = 0;
    int AccumulateTime = formulaNum-1;

    while(AccumulateStop == false)
    {
        CurrStruct->nxtCumulant[AccumulateTime].lstCumulant = CurrStruct;
        CurrStruct->combine[AccumulateTime] = 1;
        //printf("indexJ: ");
        for(int i = 1; i <= Order; i++)
        {
            TempIndex[i][1] = indexJ[i];
            CurrStruct->combine[AccumulateTime] *= Combinatorial(TempIndex[i][0], TempIndex[i][1]);
            //printf("%d ", TempIndex[i][1]);
        }
        //printf("\n");

        //绗竴娆¤皟鐢↖nitial鍚庯紝G_CurrMomentPos == formulaNum+1
        if(g_factor_moment_precompute == true && G_MomentFull == false)
        {   
            //printf("%d\n", CurrStruct->lstCumulant->nxtCumulant[0].Order);
            CurrStruct->moment[AccumulateTime] = Moment5(indexJ, Order);
            //printf("%d\n", CurrStruct->lstCumulant->nxtCumulant[0].Order);
            //printf("ok");
            G_Moment[G_CurrMomentPos] = CurrStruct->moment[AccumulateTime];
            G_CurrMomentPos++;
        }
        else
        {
            //鏌ユ壘绱㈠紩鐨勫€煎搴旂殑鐭?
            int indexPos = MomentLookupIndex(indexJ, Order);
            if(indexPos < 0)
            {
                CurrStruct->moment[AccumulateTime] = Moment5ValueOnly(indexJ, Order);
                if(CurrStruct->moment[AccumulateTime] == myEmpty)
                {
                    fprintf(stderr, "Failed to calculate moment in InitialCoeff\n");
                    return myEmpty;
                }
            }
            else
            {
                #if MOMENT_COMPARE_ENABLE
                {
                    double direct = Moment5ValueOnly(indexJ, Order);
                    double lookup = G_Moment[indexPos];
                    if(direct != myEmpty)
                    {
                        LogMomentCompare(indexJ, Order, direct, lookup);
                    }
                }
                #endif
                CurrStruct->moment[AccumulateTime] = G_Moment[indexPos];
            }
        }

        //CurrStruct->moment[AccumulateTime] = Moment4(TempIndex, &TempOrder);
        //remove this
        //CurrStruct->MomentOrder[AccumulateTime] = TempOrder; 

        if(CurrStruct->moment[AccumulateTime] == myEmpty || CurrStruct->combine[AccumulateTime] == myEmpty)
        {
            free(CurrStruct->combine);
            free(CurrStruct->nxtCumulant);
            free(CurrStruct->moment);
            fprintf(stderr, "Failed to calculate moment or combine\n");
            return myEmpty;
        }

        CurrStruct->nxtCumulant[AccumulateTime].Calculated = false;
        CurrStruct->nxtCumulant[AccumulateTime].value = 0;
    carryon:
        CurrStruct->nxtCumulant[AccumulateTime].Order = Order;  //姝ゅOrder鍜宨ndex鍦ㄤ笅闈㈣繘琛屼簡鏇存柊
        CurrStruct->nxtCumulant[AccumulateTime].index = (int*)calloc(Order+1, sizeof(int));
        CurrStruct->nxtCumulant[AccumulateTime].exist =false;
        //CurrStruct->nxtCumulant[AccumulateTime].index = indexJ;
        if(CurrStruct->nxtCumulant[AccumulateTime].index == NULL)
        {
            perror("Error in function InitialCoeff 491: Calloc return NULL\n");
            system("PAUSE");
            return myEmpty;
        }
        for(int i = 1; i < Order+1; i++)
        {
            int CurrIndex = index[i] - indexJ[i];
            //CurrStruct->nxtCumulant[AccumulateTime].index[i] = indexJ[i];
            CurrStruct->nxtCumulant[AccumulateTime].index[i] = CurrIndex;
        }
        //printf("\n");
        //瀵瑰瓙绱Н閲忕殑index杩涜鍘婚浂澶勭悊
        if(RemoveZero(&CurrStruct->nxtCumulant[AccumulateTime].index, &CurrStruct->nxtCumulant[AccumulateTime].Order) == -1)
        {
            CurrStruct->nxtCumulant[AccumulateTime].Calculated = true;
            CurrStruct->nxtCumulant[AccumulateTime].value = 0;
            //return 0;
        }
        //printf("%d\n", CurrStruct->nxtCumulant[AccumulateTime].Order);
        //瀹屾垚涓€娆¤凯浠ｅ悗杩涜鍗曟鍔犳硶锛岃€冭檻杩涗綅
        //AccumulateTime++;
        AccumulateTime--;
        indexJ[Order]++;
        int AccumulateIndex = Order;
        if(indexJ[AccumulateIndex] > index[AccumulateIndex] - 1)
        {
            indexJ[AccumulateIndex] = 0;
            indexJ[--AccumulateIndex]++;
        }
        while(true)
        {
            if(AccumulateIndex < 1)
            {
                AccumulateStop = true;
                break;
            }
            if(indexJ[AccumulateIndex] > index[AccumulateIndex])
            {
                indexJ[AccumulateIndex] = 0;
                if(AccumulateIndex == 1)
                {
                    AccumulateStop = true;
                    break;
                }
                indexJ[--AccumulateIndex]++;
            }
            else
            {
                break;
            }
        }
    }
    Free2DVec_int(&TempIndex, Order+1);
    free(indexJ);
    if(g_factor_moment_precompute == true && G_MomentFull == false)
    {
        G_MomentFull = true;
        InitFactorSwitches();
        if(g_factor_hash_lookup == true && G_MomentHashReady == false)
        {
            BuildMomentHash();
        }
    }
    return formulaNum;
}

double TrueCumulant(int* index, int Order)
{
    if(index == NULL)
    {
        fprintf(stderr, "Index is NULL\n");
        return myEmpty;
    }
    if(Order <= 0)
    {
        return 0;
    }
    if(IsLegacyBellRootIndex(index, Order) == true)
    {
        return LegacyBellCumulantDirect(index, Order);
    }
    if(Order == 1)
    {
        int MomentIndex = IndexRank(index, Order);
        return G_Moments[Order][MomentIndex];
    }
    //int MomentIndex = IndexRank(index, Order);
    //double except = G_Moments[Order][MomentIndex];
    //鎵€鏈夊垝鍒嗛」鐨勬€诲拰
    double CurrValue1 = 0;
    for(int i = 1; i <= Order; i++)
    {
        generate_k_partitions(index, Order, i);
        //鍗曚釜鍒掑垎椤圭殑鍜?
        double CurrValue2 = 0;
        for(int j = 0; j < g_total_partitions; j++)
        {
            int** TempIndex = (int**)calloc(i, sizeof(int*));
            if(TempIndex == NULL)
            {
                perror("Error in function TrueCumulant: calloc return NULL\n");
                return myEmpty;
            }
            //姣忎竴椤圭殑涔樼Н
            double TempCum = 1;
            for(int k = 0; k < i; k++)
            {
                TempIndex[k] = (int*)calloc(g_all_partitions[j].block_sizes[k], sizeof(int));
                if(TempIndex[k] == 0)
                {
                    perror("Error in function TrueCumulant: calloc return NULL\n");
                    return myEmpty;
                }
                for(int l = 0; l < g_all_partitions[j].block_sizes[k]; l++)
                {
                    TempIndex[k][l] = g_all_partitions[j].blocks[k][l];
                }
                int TempMomentIndex = IndexRank(TempIndex[k], g_all_partitions[j].block_sizes[k]);
                TempCum *= G_Moments[g_all_partitions[j].block_sizes[k]][TempMomentIndex];
            }
            CurrValue2 += TempCum; 
            Free2DVec_int(&TempIndex, i);
        }
        CurrValue2 *= Factorial(i-1);
        if(i % 2 != 0)
            CurrValue1 += CurrValue2;
        else
            CurrValue1 -= CurrValue2;
    }
    return CurrValue1;
}

void backtrack(int* arr, int n, int k, int index, Partition* current, int used_blocks) 
{
    if (index == n) 
    {
        if (used_blocks == k) 
        {
            // 妫€鏌ユ槸鍚﹀瓨鍦ㄦ煇涓?block 鍙湁涓€涓厓绱?
            int has_singleton = 0;
            for (int i = 0; i < k; ++i) 
            {
                if (current->block_sizes[i] == 1) 
                {
                    has_singleton = 1;
                    break;
                }
            }

            if (!has_singleton)
            {
                g_all_partitions[g_total_partitions++] = *current;
            }
        }
        return;
    }

    for (int i = 0; i < used_blocks; ++i) 
    {
        int size = current->block_sizes[i];
        current->blocks[i][size] = arr[index];
        current->block_sizes[i]++;
        backtrack(arr, n, k, index + 1, current, used_blocks);
        current->block_sizes[i]--;  // 鍥炴函
    }

    if (used_blocks < k) 
    {
        current->blocks[used_blocks][0] = arr[index];
        current->block_sizes[used_blocks] = 1;
        backtrack(arr, n, k, index + 1, current, used_blocks + 1);
        current->block_sizes[used_blocks] = 0;  // 鍥炴函
    }
}

void generate_k_partitions(int* arr, int n, int k) 
{
    g_total_partitions = 0;
    Partition init;
    memset(&init, 0, sizeof(init));
    backtrack(arr, n, k, 0, &init, 0);
}

void print_partitions(int k) 
{
    for (int p = 0; p < g_total_partitions; ++p) 
    {
        printf("Partition %d:\n", p + 1);
        for (int b = 0; b < k; ++b) {
            printf("  Block %d: ", b + 1);
            for (int i = 0; i < g_all_partitions[p].block_sizes[b]; ++i) 
            {
                printf("%d ", g_all_partitions[p].blocks[b][i]);
            }
            printf("\n");
        }
        printf("\n");
    }
}


