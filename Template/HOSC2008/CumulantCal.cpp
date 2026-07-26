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


    CurrCumulantNum = 0;


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


    int NuminFstFloor = InitialCoeff(FstCumulant);


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


    int FstFloorStart = 0;

    int FirstFloorBatch = g_factor_iterative_schedule == true ? DEFAULT_FST_FLOOR_NUM : NuminFstFloor;
    if(FirstFloorBatch <= 0)
    {
        FirstFloorBatch = NuminFstFloor;
    }
    int FstFloorEnd = FstFloorStart + FirstFloorBatch;

    double GroupResult = 0;

    bool ReachEnd = false;

    CumulantIndex*** Floor;

    int CurrFloor;

    int FloorCumNum[MaxFloor];

    while(true)
    {


        if(FstFloorEnd >= NuminFstFloor)
        {
            ReachEnd = true;
            FstFloorEnd = NuminFstFloor;
        }




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

        CumulantIndex* CurrCumulant = FstCumulant;

        int CurrIndex = 0;


        bool AllTrue = true;
        int NuminCurrFloor = FstFloorEnd - FstFloorStart;
        for(int i = 0; i < NuminCurrFloor; i++)

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

                RootAVL2[CurrCumulant->nxtCumulant[i + FstFloorStart].Order] = insert(RootAVL2[CurrCumulant->nxtCumulant[i + FstFloorStart].Order], CurrCumulant->nxtCumulant[i + FstFloorStart].index, CurrCumulant->nxtCumulant[i + FstFloorStart].value, CurrCumulant->nxtCumulant[i + FstFloorStart].Order, CurrCumulant->nxtCumulant[i + FstFloorStart].Calculated);

            }
            else if(g_factor_avl_reuse == true)
            {
                Floor[CurrFloor][i]->exist = true;

            }
            if(CurrCumulant->nxtCumulant[i + FstFloorStart].Calculated == false)
            {
                AllTrue = false;
            }
        }


        FloorCumNum[CurrFloor] = NuminCurrFloor;

        if(AllTrue == true)
        {



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

        }

        int NuminFloor[MaxCumulantNum];



        int cumulantIndex = 0;
        int nxtStart = 0;
        CurrFloor++;

        int SubCumulantNum = 0;

        while(true)
        {

            if(cumulantIndex >= NuminCurrFloor)
            {

                FloorCumNum[CurrFloor] = SubCumulantNum;

                break;
            }

            if(Floor[CurrFloor - 1][cumulantIndex]->Calculated != true && Floor[CurrFloor - 1][cumulantIndex]->exist == false)
                CurrCumulant = Floor[CurrFloor - 1][cumulantIndex];
            else
            {


                NuminFloor[cumulantIndex] = 0;

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

            if(CurrCumulant->Calculated == true)
            {
                NuminFloor[cumulantIndex] = 0;

                cumulantIndex++;

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

            if(SubCumulantNum > FstFloorEnd - FstFloorStart)
            {

                Floor[CurrFloor] = (CumulantIndex**)realloc(Floor[CurrFloor], SubCumulantNum * sizeof(CumulantIndex*));
                if(Floor[CurrFloor] == NULL)
                {
                    perror("Error in function BetterRecCumulat3: realloc return NULL\n");
                    return myEmpty;
                }
            }
            for(int i = nxtStart; i < NuminFloor[cumulantIndex] + nxtStart; i++)

            {

                Floor[CurrFloor][i] = &CurrCumulant->nxtCumulant[i - nxtStart];


                if(i + 1 == NuminFloor[cumulantIndex] + nxtStart)
                {

                    nxtStart = i + 1;
                    break;
                }
            }
            cumulantIndex++;
        }


        bool TempFinish = false;




        int testI = 0;
        while(true)
        {

            CurrIndex = 0;

            int CurrCumulantIndex = 0;
            if(CurrFloor == 0 && GroupResult != 0)
            {
                CurrCumulantIndex = FstFloorStart;
            }
            double CurrValue = 0;

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


                if((CurrIndex == FloorCumNum[CurrFloor] && CurrFloor != 0) || (CurrIndex == FstFloorEnd - FstFloorStart && CurrFloor == 0))
                {
                    CurrCumulant->lstCumulant->Calculated = true;



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

                        if(ReachEnd == true && GroupResult == 0)
                        {
                            FstCumulant->Calculated = true;





                            double except1 = RootMomentOrDirect(FstCumulant->index, FstCumulant->Order);

                            FstCumulant->value = except1 - CurrValue;

                            goto FreeStack;
                        }

                        else if(ReachEnd == true && GroupResult != 0)
                        {
                            GroupResult += CurrValue;
                            FstCumulant->Calculated = true;





                            double except2 = RootMomentOrDirect(FstCumulant->index, FstCumulant->Order);

                            FstCumulant->value = except2 - GroupResult;

                            goto FreeStack;
                        }

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







                if(Floor[CurrFloor][CurrIndex]->lstCumulant != CurrCumulant->lstCumulant)
                {


                    CurrCumulantIndex = 0;
                    CurrCumulant->lstCumulant->Calculated = true;




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


                    CurrCumulant->lstCumulant->value = except - CurrValue;
                    CurrCumulant->lstCumulant->Calculated = true;









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




    FreeCumulantIndex(&FstCumulant);
    free(FstCumulant->nxtCumulant);
    if(FstCumulant->Calculated == true)
    {
        double value = FstCumulant->value;
        free(FstCumulant);

        return value;
    }
    return myEmpty;
}

double LowOrderCumulant(int* index, int Order)
{
    RemoveZero(&index, &Order);

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

    if(CurrStruct->index == NULL)
    {
        fprintf(stderr, "Error in function InitialCoeff: Wrong parameter\n");
        return myEmpty;
    }


    if(RemoveZero(&CurrStruct->index, &CurrStruct->Order) == -1 || CurrStruct->Order == 0)
    {
        CurrStruct->value = 0;
        CurrStruct->Calculated = true;
        return 0;
    }


    CurrStruct->Calculated = false;

    int* index = CurrStruct->index;
    int Order = CurrStruct->Order;

    if(index[0] <= 3)
    {
            CurrStruct->value = LowOrderCumulant(index, Order);
            CurrStruct->Calculated = true;
            return 0;
    }


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

    }



    int formulaNum = 1;
    for(int i = 1; i <= Order; i++)
    {

        formulaNum *= TempIndex[i][0] + 1;
    }

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


    }
    if(g_factor_moment_precompute == true && G_MomentFull == false)
    {
        G_Moment[G_CurrMomentPos] = Moment5(index, Order);
        G_CurrMomentPos++;
    }


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



    bool AccumulateStop = false;


    int AccumulateTime = formulaNum-1;

    while(AccumulateStop == false)
    {
        CurrStruct->nxtCumulant[AccumulateTime].lstCumulant = CurrStruct;
        CurrStruct->combine[AccumulateTime] = 1;

        for(int i = 1; i <= Order; i++)
        {
            TempIndex[i][1] = indexJ[i];
            CurrStruct->combine[AccumulateTime] *= Combinatorial(TempIndex[i][0], TempIndex[i][1]);

        }



        if(g_factor_moment_precompute == true && G_MomentFull == false)
        {

            CurrStruct->moment[AccumulateTime] = Moment5(indexJ, Order);


            G_Moment[G_CurrMomentPos] = CurrStruct->moment[AccumulateTime];
            G_CurrMomentPos++;
        }
        else
        {

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
        CurrStruct->nxtCumulant[AccumulateTime].Order = Order;
        CurrStruct->nxtCumulant[AccumulateTime].index = (int*)calloc(Order+1, sizeof(int));
        CurrStruct->nxtCumulant[AccumulateTime].exist =false;

        if(CurrStruct->nxtCumulant[AccumulateTime].index == NULL)
        {
            perror("Error in function InitialCoeff 491: Calloc return NULL\n");
            system("PAUSE");
            return myEmpty;
        }
        for(int i = 1; i < Order+1; i++)
        {
            int CurrIndex = index[i] - indexJ[i];

            CurrStruct->nxtCumulant[AccumulateTime].index[i] = CurrIndex;
        }


        if(RemoveZero(&CurrStruct->nxtCumulant[AccumulateTime].index, &CurrStruct->nxtCumulant[AccumulateTime].Order) == -1)
        {
            CurrStruct->nxtCumulant[AccumulateTime].Calculated = true;
            CurrStruct->nxtCumulant[AccumulateTime].value = 0;

        }



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



    double CurrValue1 = 0;
    for(int i = 1; i <= Order; i++)
    {
        generate_k_partitions(index, Order, i);

        double CurrValue2 = 0;
        for(int j = 0; j < g_total_partitions; j++)
        {
            int** TempIndex = (int**)calloc(i, sizeof(int*));
            if(TempIndex == NULL)
            {
                perror("Error in function TrueCumulant: calloc return NULL\n");
                return myEmpty;
            }

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
        current->block_sizes[i]--;
    }

    if (used_blocks < k)
    {
        current->blocks[used_blocks][0] = arr[index];
        current->block_sizes[used_blocks] = 1;
        backtrack(arr, n, k, index + 1, current, used_blocks + 1);
        current->block_sizes[used_blocks] = 0;
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


