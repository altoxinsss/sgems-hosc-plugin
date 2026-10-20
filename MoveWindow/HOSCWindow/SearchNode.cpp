#include "SearchNode.h"
#include <float.h>

int G_TemplateNum;
Template* G_Template;
double** G_Moments;
double* G_Moment;
int G_CurrMomentPos;
bool G_MomentFull;
int* G_MomentHash;
int G_MomentHashSize;
bool G_MomentHashReady;

int ChangeTemplate(const char* filename, int Col)
{
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("无法打开文件");
        return 0;
    }

    char line[1024];
    int test1 = -1;
    int test2 = 1;
    while (fgets(line, sizeof(line), file)) {
        char* token = strtok(line, ",\n");
        if(test1 == -1)
        {
            test1++;
            G_TemplateNum = atoi(token);
            continue;
        }
        int count = 0;
        while (token && count < Col) {
            double value = atof(token);
            G_Template[test1].Value[test2] = value;
            test2++;
            token = strtok(NULL, ",\n");
            count++;
        }
        test2 = 0;
        test1++;

        if (count != Col) {
            fprintf(stderr, "[警告] 此行只有 %d 列，期望 %d 列。\n", count, Col);
        }
    }

    fclose(file);
    return 0;
}

int GetIndex(int* Num, int Size, bool MAX_MIN)
{
    int MaxNum = 0;
    int Index = 0;
    if(MAX_MIN == true)
    for(int i = 0; i < Size; i++)
    {
        if(Num[i] > MaxNum)
        {
            MaxNum = Num[i];
            Index = i;
        }
    }
    else
    for(int i = 0; i < Size; i++)
    {
        if(Num[i] < MaxNum)
        {
            MaxNum = Num[i];
            Index = i;
        }
    }
    return Index;
}

int InitialLagMat(int*** LagMat, int* TotalDirectionsNum)
{
    int Directions = 1;
    for(int i = 0; i < pt_Keys->NDIR; i++)
    {
        Directions = Directions * pt_Keys->NLAG[i];
    }

    int** LagMats = (int**)calloc(Directions, sizeof(int*));
    if(LagMats == NULL)
    {
        perror("Error in function LagMats: calloc return NULL\n");
        return -1;
    }
    for(int i = 0; i < Directions; i++)
    {
        LagMats[i] = (int*)calloc(pt_Keys->NDIR, sizeof(int));
        if(LagMats[i] == NULL)
        {
            perror("Error in function LagMats: calloc return NULL\n");
            return -1;
        }
    }
    int CurrPos = 0;
    for(int i = 0; i < Directions; i++)
    {

        if(i == 0)
            continue;
        CurrPos = pt_Keys->NDIR - 1;
        LagMats[i][CurrPos]+=pt_Keys->LAG[CurrPos];
        if(LagMats[i][CurrPos] >= pt_Keys->LAG[CurrPos] * pt_Keys->NLAG[CurrPos])
        {
            LagMats[i][CurrPos] = 0;
            while(true)
            {
                CurrPos--;
                LagMats[i][CurrPos]+=pt_Keys->LAG[CurrPos];
                if(LagMats[i][CurrPos] < pt_Keys->LAG[CurrPos] * pt_Keys->NLAG[CurrPos])
                {
                    break;
                }
                else if(CurrPos == 0)
                {

                    *LagMat = LagMats;
                    *TotalDirectionsNum = Directions;
                    return 0;
                }
                else
                {
                    LagMats[i][CurrPos] = 0;
                }
            }

        }

        if(i + 1< Directions)
            for(int j = 0; j < pt_Keys->NDIR; j++)
            {
                LagMats[i+1][j] = LagMats[i][j];
            }

    }
    *LagMat = LagMats;
    *TotalDirectionsNum = Directions;
    return 0;
}

int SaveTemplateResult(Template* Templates, int Order, int TemplateNum)
{
    if(Templates == NULL || TemplateNum == 0 || Order == 0)
    {
        fprintf(stderr, "Error in function SaveTemplateResult: Templates is NULL or Order is 0 or TemplateNum is 0\n");
        return -1;
    }
    FILE* file = fopen("Template.txt", "w");
    if(file == NULL)
    {
        perror("Error opening file");
        return -1;
    }

    fprintf(file, "%d\n", Order);

    for(int i = 0; i < TemplateNum; i++)
    {
        if(Templates[i].Value == NULL)
        {
            fprintf(stderr, "Error in function SaveTemplateResult: Templates[%d].Value is NULL\n", i);
            fclose(file);
            return -1;
        }
        for(int j = 1; j <= Order; j++)
        {
            if(j == Order)
            {
                fprintf(file, "%f\n", Templates[i].Value[j]);

            }
            else
            {
                fprintf(file, "%f    ", Templates[i].Value[j]);

            }
        }
    }
    fclose(file);
    return 0;
}


int TemplateConstruct_2D(int* pint32_LAG, int Order, Template** Templates, int* TemplateNums, int range1, int range2)
{
    if(pt_Keys == NULL)
    {
        fprintf(stderr, "Error in function TemplateConstruct_2D: pt_Keys is NULL\n");
        return -1;
    }
    if(pt_Grids == NULL)
    {
        fprintf(stderr, "Error in function TemplateConstruct_2D: pt_Grids is NULL\n");
        return -1;
    }
    if(ValueMap_2D == NULL)
    {
        fprintf(stderr, "Error in function TemplateConstruct_2D: ValueMap_2D is NULL\n");
        return -1;
    }

    double** ValueMap;
        ValueMap = (double**)calloc(pt_Keys->NX, sizeof(double*));
        if(ValueMap == NULL)
        {
            perror("Error in function TemplateConstruct_2D: calloc return NULL\n");
            return -1;
        }
        for(int i = 0; i < pt_Keys->NX; i++)
        {
            ValueMap[i] = (double*)calloc(pt_Keys->NY, sizeof(double));
            if(ValueMap[i] == NULL)
            {
                perror("Error in function TemplateConstruct_2D: calloc return NULL\n");
                free(ValueMap);
                ValueMap = NULL;
                return -1;
            }
        }

        int TemplateNum = 0;
        int Templated = 0;
        int TemplateSize = range2;

        Template* templates = (Template*)calloc(TemplateSize, sizeof(Template));
        if(templates == NULL)
        {
            perror("Error in function TemplateConstruct_2D: calloc return NULL\n");
            if(Free2DVec(&ValueMap, pt_Keys->NX) == -1)
            {
                fprintf(stderr, "Failed to free 2DVec\n");
            }
            ValueMap = NULL;
            return -1;
        }

        for(int i = range1; i < range2; i++)
        {
            templates[i].Lag = (int*)calloc(Order + 1, sizeof(int));
            templates[i].Value = (double*)calloc(Order + 1, sizeof(double));
            if(templates[i].Lag == NULL || templates[i].Value == NULL)
            {
                perror("Error in function TemplateConstruct_2D: calloc return NULL\n");
                Free2DVec(&ValueMap, pt_Keys->NX);
                free(templates);
                templates = NULL;
                ValueMap = NULL;
                return -1;
            }

        }


        double*** TolerentZone;
        if(Init3DVec(&TolerentZone, Order - 1, 4, 2) == -1)
        {
            fprintf(stderr, "Failed to free 3DVec\n");
        }

        int* TaggedX = (int*)calloc(Order-1, sizeof(int));
        int* TaggedY = (int*)calloc(Order-1, sizeof(int));
        if(TaggedX == NULL || TaggedY == NULL)
        {
            perror("Error in function TemplateConstruct_2D: calloc return NULL\n");
            FreeTemplate(&templates);
            Free2DVec(&ValueMap, pt_Keys->NX);
            Free3DVec(&TolerentZone, Order-1, 4);
            templates = NULL;
            ValueMap = NULL;
            TolerentZone = NULL;
            return -1;
        }

        for(int i = 0; i < TemplateSize; i++)
        {

            if(ValueMap[pt_Grids->pt_Nodes[i].x][pt_Grids->pt_Nodes[i].y] < 0)
            {
                continue;
            }

            Templated = 0;

            for(int j = 0; j < Order; j++)
            {
                if(j == Order - 1 && Templated == Order - 1)
                {
                    ValueMap[pt_Grids->pt_Nodes[i].x][pt_Grids->pt_Nodes[i].y] = -i;
                    for(int k = 0; k < Order - 1; k++)
                    {
                        ValueMap[TaggedX[k]][TaggedY[k]] = -i;
                    }

                    templates[TemplateNum].Value[1] = pt_Grids->pt_Nodes[i].val;
                    templates[TemplateNum].Lag[1] = 0;

                    TemplateNum++;

                    break;
                }
                else if(j == Order - 1 && Templated < Order - 1)
                {
                    break;
                }
                int LAG = pint32_LAG[j];
                templates[TemplateNum].Lag[j+2] = LAG;
                int LAG_t = pt_Keys->LTOL[j];
                double ATOL = pt_Keys->ATOL[j] * PI / 180;
                double AZM = pt_Keys->AZM[j] * PI / 180;


                int lag = LAG - LAG_t;

                int MiddleZone_x = pt_Grids->pt_Nodes[i].x + LAG * cos(AZM);
                int MiddleZone_y = pt_Grids->pt_Nodes[i].y + LAG * sin(AZM);


                TolerentZone[j][0][0] = pt_Grids->pt_Nodes[i].x + lag / cos(ATOL) * cos(AZM + ATOL);
                TolerentZone[j][0][1] = pt_Grids->pt_Nodes[i].y + lag / cos(ATOL) * sin(AZM + ATOL);
                TolerentZone[j][1][0] = pt_Grids->pt_Nodes[i].x + lag / cos(ATOL) * cos(AZM - ATOL);
                TolerentZone[j][1][1] = pt_Grids->pt_Nodes[i].y + lag / cos(ATOL) * sin(AZM - ATOL);
                TolerentZone[j][2][0] = 2 * MiddleZone_x - TolerentZone[j][0][0];
                TolerentZone[j][2][1] = 2 * MiddleZone_y - TolerentZone[j][0][1];
                TolerentZone[j][3][0] = 2 * MiddleZone_x - TolerentZone[j][1][0];
                TolerentZone[j][3][1] = 2 * MiddleZone_y - TolerentZone[j][1][1];


                if(TolerentZone[j][0][0] >= 0 && TolerentZone[j][0][1] >= 0 && TolerentZone[j][0][0] <= pt_Keys->DIMX && TolerentZone[j][0][1] <= pt_Keys->DIMY)
                {
                    if(ValueMap[(int)TolerentZone[j][0][0]][(int)TolerentZone[j][0][1]] >= 0)
                    {
                        if(ValueMap_2D[(int)TolerentZone[j][0][0]][(int)TolerentZone[j][0][1]] != -1)
                        {


                            templates[TemplateNum].Value[j+2] = ValueMap_2D[(int)TolerentZone[j][0][0]][(int)TolerentZone[j][0][1]];

                            TaggedX[j] = (int)TolerentZone[j][0][0];
                            TaggedY[j] = (int)TolerentZone[j][0][1];
                            Templated++;
                            continue;
                        }
                        else
                            break;
                    }
                    else
                        break;
                }
                else if(TolerentZone[j][1][0] >= 0 && TolerentZone[j][1][1] >= 0 && TolerentZone[j][1][0] <= pt_Keys->DIMX && TolerentZone[j][1][1] <= pt_Keys->DIMY)
                {
                    if(ValueMap[(int)TolerentZone[j][1][0]][(int)TolerentZone[j][1][1]] >= 0)
                    {
                        if(ValueMap_2D[(int)TolerentZone[j][1][0]][(int)TolerentZone[j][1][1]] != -1)
                        {


                            templates[TemplateNum].Value[j+2] = ValueMap_2D[(int)TolerentZone[j][1][0]][(int)TolerentZone[j][1][1]];

                            TaggedX[j] = (int)TolerentZone[j][1][0];
                            TaggedY[j] = (int)TolerentZone[j][1][1];
                            Templated++;
                            continue;
                        }
                        else
                            break;
                    }
                    else
                        break;
                }

                else if(TolerentZone[j][2][0] >= 0 && TolerentZone[j][2][1] >= 0 && TolerentZone[j][2][0] <= pt_Keys->DIMX && TolerentZone[j][2][1] <= pt_Keys->DIMY)
                {
                    if(ValueMap[(int)TolerentZone[j][2][0]][(int)TolerentZone[j][2][1]] >= 0)
                    {
                        if(ValueMap_2D[(int)TolerentZone[j][2][0]][(int)TolerentZone[j][2][1]] != -1)
                        {



                            templates[TemplateNum].Value[j+2] = ValueMap_2D[(int)TolerentZone[j][2][0]][(int)TolerentZone[j][2][1]];

                            TaggedX[j] = (int)TolerentZone[j][2][0];
                            TaggedY[j] = (int)TolerentZone[j][2][1];
                            Templated++;
                            continue;
                        }
                        else
                            break;
                    }
                    else
                        break;
                }
                else if(TolerentZone[j][3][0] >= 0 && TolerentZone[j][3][1] >= 0 && TolerentZone[j][3][0] <= pt_Keys->DIMX && TolerentZone[j][3][1] <= pt_Keys->DIMY)
                {
                    if(ValueMap[(int)TolerentZone[j][3][0]][(int)TolerentZone[j][3][1]] >= 0)
                    {
                        if(ValueMap_2D[(int)TolerentZone[j][3][0]][(int)TolerentZone[j][3][1]] != -1)
                        {



                            templates[TemplateNum].Value[j+2] = ValueMap_2D[(int)TolerentZone[j][3][0]][(int)TolerentZone[j][3][1]];

                            TaggedX[j] = (int)TolerentZone[j][3][0];
                            TaggedY[j] = (int)TolerentZone[j][3][1];
                            Templated++;
                            continue;
                        }
                        else
                            break;
                    }
                    else
                        break;
                }

                else
                {
                    break;
                }
            }
        }


        Free2DVec(&ValueMap, pt_Keys->NX);
        ValueMap = NULL;

        free(TaggedX);
        TaggedX = NULL;
        free(TaggedY);
        TaggedY = NULL;
        free(ValueMap);
        Free3DVec(&TolerentZone, Order-1, 4);
        TolerentZone = NULL;
        *TemplateNums += TemplateNum;
        G_TemplateNum = TemplateNum;
        if(TemplateNum == 0)
        {
            FreeTemplate(&templates);
            templates = NULL;
            return 1;
        }


        double TempMean;
        for(int i = 1; i < Order + 1; i++)
        {
            TempMean = 0;
            for(int j = 0; j < TemplateNum; j++)
            {
                TempMean += templates[j].Value[i];
            }

            if(TempMean != 0)
            {
                for(int j = 0; j < TemplateNum; j++)
                {
                    templates[j].Value[i] -= TempMean / TemplateNum;
                }
            }
        }


        templates->index = (int*)calloc(Order+1, sizeof(int));
        if(templates == NULL)
        {
            perror("Error in function TemplateConstruct_2D: calloc return NULL\n");
            return -1;
        }

        int tempOrder = 0;
        for(int i = 2; i <= Order; i++)
        {

            templates->index[i] = pt_Keys->NOD[i-2];
            tempOrder += templates->index[i];
        }

        templates->index[1] = 1;
        templates->index[0] = ++tempOrder;



        *Templates = templates;
        G_Template = templates;


        G_Moments = NULL;
        G_Moment = NULL;
        G_CurrMomentPos = 0;
        G_MomentFull = false;


    return 0;
}

int TemplateConstruct_3D(int* pint32_LAG, int Order, Template** Templates, int* TemplateNums, int range1, int range2)
{
    if(pt_Keys == NULL)
    {
        fprintf(stderr, "Error in function TemplateConstruct_2D: pt_Keys is NULL\n");
        return -1;
    }
    if(pt_Grids == NULL)
    {
        fprintf(stderr, "Error in function TemplateConstruct_2D: pt_Grids is NULL\n");
        return -1;
    }
    if(ValueMap_3D == NULL)
    {
        fprintf(stderr, "Error in function TemplateConstruct_2D: ValueMap_2D is NULL\n");
        return -1;
    }

    double*** ValueMap;
        ValueMap = (double***)calloc(MaxX, sizeof(double**));
        if(ValueMap == NULL)
        {
            perror("Error in function TemplateConstruct_3D: calloc return NULL\n");
            return -1;
        }
        for(int i = 0; i < MaxX; i++)
        {
            ValueMap[i] = (double**)malloc(MaxY * sizeof(double*));
            if(ValueMap[i] == NULL)
            {
                fprintf(stderr, "Error in function TemplateConstruct_3D: malloc return -1\n");
                free(ValueMap);
                ValueMap = NULL;
                return -1;
            }
            for(int j = 0; j < MaxY; j++)
            {
                ValueMap[i][j] = (double*)malloc(MaxZ * sizeof(double));
                if(ValueMap[i][j] == NULL)
                {
                    fprintf(stderr, "Error in function TemplateConstruct_3D: malloc return -1\n");
                    free(ValueMap);
                    ValueMap = NULL;
                    return -1;
                }
                for(int k = 0; k < MaxZ; k++)
                {
                    ValueMap[i][j][k] = ValueMap_3D[i][j][k];
                }
            }
        }

        int TemplateNum = 0;
        int Templated = 0;
        int TemplateSize = range2;

        Template* templates = (Template*)calloc(TemplateSize, sizeof(Template));
        if(templates == NULL)
        {
            perror("Error in function TemplateConstruct_2D: calloc return NULL\n");
            if(Free3DVec(&ValueMap, MaxX, MaxY) == -1)
            {
                fprintf(stderr, "Failed to free 2DVec\n");
            }
            ValueMap = NULL;
            return -1;
        }

        for(int i = range1; i < range2; i++)
        {
            templates[i].Lag = (int*)calloc(Order + 1, sizeof(int));
            templates[i].Value = (double*)calloc(Order + 1, sizeof(double));
            if(templates[i].Lag == NULL || templates[i].Value == NULL)
            {
                perror("Error in function TemplateConstruct_2D: calloc return NULL\n");
                Free3DVec(&ValueMap, MaxX, MaxY);
                free(templates);
                templates = NULL;
                ValueMap = NULL;
                return -1;
            }

        }


        double*** TolerentZone;
        if(pt_Keys->REG_IRREG == 0)
        {
            if(Init3DVec(&TolerentZone, Order - 1, 4, 3) == -1)
            {
                fprintf(stderr, "Failed to free 3DVec\n");
            }
        }
        int* TaggedX = (int*)calloc(Order-1, sizeof(int));
        int* TaggedY = (int*)calloc(Order-1, sizeof(int));
        int* TaggedZ = (int*)calloc(Order-1, sizeof(int));
        if(TaggedX == NULL || TaggedY == NULL || TaggedZ == NULL)
        {
            perror("Error in function TemplateConstruct_2D: calloc return NULL\n");
            FreeTemplate(&templates);
            Free3DVec(&ValueMap, MaxX, MaxY);
            if(pt_Keys->REG_IRREG == 0)
                Free3DVec(&TolerentZone, Order-1, 4);
            templates = NULL;
            ValueMap = NULL;
            TolerentZone = NULL;
            return -1;
        }

        for(int i = 0; i < TemplateSize; i++)
        {

            if(ValueMap[pt_Grids->pt_Nodes[i].x][pt_Grids->pt_Nodes[i].y][pt_Grids->pt_Nodes[i].z] < 0)
            {
                continue;
            }

            Templated = 0;

            for(int j = 0; j < Order; j++)
            {
                if(j == Order - 1 && Templated == Order - 1)
                {
                    ValueMap[pt_Grids->pt_Nodes[i].x][pt_Grids->pt_Nodes[i].y][pt_Grids->pt_Nodes[i].z] = -i;
                    for(int k = 0; k < Order - 1; k++)
                    {
                        ValueMap[TaggedX[k]][TaggedY[k]][TaggedZ[k]] = -i;
                    }

                    templates[TemplateNum].Value[1] = pt_Grids->pt_Nodes[i].val;
                    templates[TemplateNum].Lag[1] = 0;
                    TemplateNum++;

                    break;
                }
                else if(j == Order - 1 && Templated < Order - 1)
                {
                    break;
                }
                if(pt_Keys->REG_IRREG == 0)
                {
                    int LAG = pint32_LAG[j];
                    templates[TemplateNum].Lag[j+2] = LAG;
                    int LAG_t = pt_Keys->LTOL[j];

                    double ATOL = pt_Keys->ATOL[j] * PI / 180;
                    double AZM = pt_Keys->AZM[j] * PI / 180;

                    double ELEVATION = pt_Keys->DIP[j] * PI / 180;





                    int lag = LAG - LAG_t;
                    int TZoneHeight = 2 * LAG_t;
                    int MiddleZone_x = pt_Grids->pt_Nodes[i].x + LAG * cos(AZM);
                    int MiddleZone_y = pt_Grids->pt_Nodes[i].y + LAG * sin(AZM);


                    TolerentZone[j][0][0] = pt_Grids->pt_Nodes[i].x + lag / cos(ATOL) * cos(AZM + ATOL);
                    TolerentZone[j][0][1] = pt_Grids->pt_Nodes[i].y + lag / cos(ATOL) * sin(AZM + ATOL);
                    TolerentZone[j][0][2] = pt_Grids->pt_Nodes[i].z + LAG * sin(ELEVATION);
                    TolerentZone[j][1][0] = pt_Grids->pt_Nodes[i].x + lag / cos(ATOL) * cos(AZM - ATOL);
                    TolerentZone[j][1][1] = pt_Grids->pt_Nodes[i].y + lag / cos(ATOL) * sin(AZM - ATOL);
                    TolerentZone[j][1][2] = pt_Grids->pt_Nodes[i].z + LAG * sin(ELEVATION);
                    TolerentZone[j][2][0] = 2 * MiddleZone_x - TolerentZone[j][0][0];
                    TolerentZone[j][2][1] = 2 * MiddleZone_y - TolerentZone[j][0][1];
                    TolerentZone[j][2][2] = pt_Grids->pt_Nodes[i].z + LAG * sin(ELEVATION);
                    TolerentZone[j][3][0] = 2 * MiddleZone_x - TolerentZone[j][1][0];
                    TolerentZone[j][3][1] = 2 * MiddleZone_y - TolerentZone[j][1][1];
                    TolerentZone[j][3][2] = pt_Grids->pt_Nodes[i].z + LAG * sin(ELEVATION);

                    #if 0

                    int NodeX = -1;
                    int NodeY = -1;
                    int NodeZ = pt_Grids->pt_Nodes[i].z + LAG * sin(ELEVATION);

                    bool FoundZ = false;
                    for(int k = 0; k < pt_Grids->ND; k++)
                    {
                        if(abs(pt_Grids->pt_Nodes[i].z - NodeZ) <= pt_Keys->BANDWD[j])
                        {
                            NodeZ = pt_Grids->pt_Nodes[i].z;
                            FoundZ = true;
                            break;
                        }
                        else
                            continue;
                    }
                    if(FoundZ == false)
                        break;
                    if(SearchInTolerentZone(TolerentZone[j], &NodeX, &NodeY, NodeZ) == 0)
                    {
                        templates[TemplateNum].Value[j+2] = ValueMap_3D[NodeX][NodeY][NodeZ];
                        TaggedX[j] = NodeX;
                        TaggedY[j] = NodeY;
                        TaggedZ[j] = NodeZ;
                        Templated++;

                        continue;
                    }
                    #endif




                    if(TolerentZone[j][0][0] >= 0 && TolerentZone[j][0][1] >= 0 && TolerentZone[j][0][2] >= 0 && TolerentZone[j][0][0] <= pt_Keys->DIMX && TolerentZone[j][0][1] <= pt_Keys->DIMY && TolerentZone[j][0][2] <= pt_Keys->DIMZ + pt_Grids->pt_Nodes[0].z)
                    {

                        if(ValueMap[(int)TolerentZone[j][0][0]][(int)TolerentZone[j][0][1]][(int)TolerentZone[j][0][2]] >= 0)
                        if(ValueMap_3D[(int)TolerentZone[j][0][0]][(int)TolerentZone[j][0][1]][(int)TolerentZone[j][0][2]] != -1)
                        {


                            templates[TemplateNum].Value[j+2] = ValueMap_3D[(int)TolerentZone[j][0][0]][(int)TolerentZone[j][0][1]][(int)TolerentZone[j][0][2]];

                            TaggedX[j] = (int)TolerentZone[j][0][0];
                            TaggedY[j] = (int)TolerentZone[j][0][1];
                            TaggedZ[j] = (int)TolerentZone[j][0][2];
                            Templated++;
                            continue;
                        }
                    }
                    else if(TolerentZone[j][1][0] >= 0 && TolerentZone[j][1][1] >= 0 && TolerentZone[j][1][2] >= 0 && TolerentZone[j][1][0] <= pt_Keys->DIMX && TolerentZone[j][1][1] <= pt_Keys->DIMY && TolerentZone[j][1][2] <= pt_Keys->DIMZ + pt_Grids->pt_Nodes[0].z)
                    {
                        if(ValueMap[(int)TolerentZone[j][1][0]][(int)TolerentZone[j][1][1]][(int)TolerentZone[j][1][2]] >= 0)
                        if(ValueMap_3D[(int)TolerentZone[j][1][0]][(int)TolerentZone[j][1][1]][(int)TolerentZone[j][1][2]] != -1)
                        {


                            templates[TemplateNum].Value[j+2] = ValueMap_3D[(int)TolerentZone[j][1][0]][(int)TolerentZone[j][1][1]][(int)TolerentZone[j][1][2]];

                            TaggedX[j] = (int)TolerentZone[j][1][0];
                            TaggedY[j] = (int)TolerentZone[j][1][1];
                            TaggedZ[j] = (int)TolerentZone[j][1][2];
                            Templated++;
                            continue;
                        }
                    }
                    else if(TolerentZone[j][2][0] >= 0 && TolerentZone[j][2][1] >= 0 && TolerentZone[j][2][2] >= 0 && TolerentZone[j][2][0] <= pt_Keys->DIMX && TolerentZone[j][2][1] <= pt_Keys->DIMY && TolerentZone[j][2][2] <= pt_Keys->DIMZ + pt_Grids->pt_Nodes[0].z)
                    {
                        if(ValueMap[(int)TolerentZone[j][2][0]][(int)TolerentZone[j][2][1]][(int)TolerentZone[j][2][2]] >= 0)
                        if(ValueMap_3D[(int)TolerentZone[j][2][0]][(int)TolerentZone[j][2][1]][(int)TolerentZone[j][2][2]] != -1)
                        {



                            templates[TemplateNum].Value[j+2] = ValueMap_3D[(int)TolerentZone[j][2][0]][(int)TolerentZone[j][2][1]][(int)TolerentZone[j][2][2]];

                            TaggedX[j] = (int)TolerentZone[j][2][0];
                            TaggedY[j] = (int)TolerentZone[j][2][1];
                            TaggedZ[j] = (int)TolerentZone[j][2][2];
                            Templated++;
                            continue;
                        }
                    }
                    else if(TolerentZone[j][3][0] >= 0 && TolerentZone[j][3][1] >= 0 && TolerentZone[j][3][2] >= 0 && TolerentZone[j][3][0] <= pt_Keys->DIMX && TolerentZone[j][3][1] <= pt_Keys->DIMY && TolerentZone[j][3][2] <= pt_Keys->DIMZ + pt_Grids->pt_Nodes[0].z)
                    {
                        if(ValueMap[(int)TolerentZone[j][3][0]][(int)TolerentZone[j][3][1]][(int)TolerentZone[j][3][2]] >= 0)
                        if(ValueMap_3D[(int)TolerentZone[j][3][0]][(int)TolerentZone[j][3][1]][(int)TolerentZone[j][3][2]] != -1)
                        {



                            templates[TemplateNum].Value[j+2] = ValueMap_3D[(int)TolerentZone[j][3][0]][(int)TolerentZone[j][3][1]][(int)TolerentZone[j][3][2]];

                            TaggedX[j] = (int)TolerentZone[j][3][0];
                            TaggedY[j] = (int)TolerentZone[j][3][1];
                            TaggedZ[j] = (int)TolerentZone[j][3][2];
                            Templated++;
                            continue;
                        }
                    }

                    else
                    {
                        break;
                    }
                }
                else
                {
                    double TempX = 0;
                    double TempY = 0;
                    double TempZ = 0;
                    double AZM = pt_Keys->AZM[j] * PI / 180;
                    double ELEVATION = pt_Keys->DIP[j] * PI / 180;
                    int LAG = pint32_LAG[j];
                    TempX = pt_Grids->pt_Nodes[i].x + LAG * cos(AZM);
                    TempY = pt_Grids->pt_Nodes[i].y + LAG * sin(AZM);
                    TempZ = pt_Grids->pt_Nodes[i].z + LAG * sin(ELEVATION);
                    if(TempX >= 0 && TempY >= 0 && TempZ >= 0 && TempX <= pt_Keys->DIMX && TempY <= pt_Keys->DIMY && TempZ <= pt_Keys->DIMZ + pt_Grids->pt_Nodes[0].z)
                    {

                        if(ValueMap[(int)TempX][(int)TempY][(int)TempZ] >= 0)
                        if(ValueMap_3D[(int)TempX][(int)TempY][(int)TempZ] != -1)
                        {


                            templates[TemplateNum].Value[j+2] = ValueMap_3D[(int)TempX][(int)TempY][(int)TempZ];

                            TaggedX[j] = (int)TempX;
                            TaggedY[j] = (int)TempY;
                            TaggedZ[j] = (int)TempZ;
                            Templated++;
                            continue;
                        }
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }


        Free3DVec(&ValueMap, MaxX, MaxY);
        ValueMap = NULL;

        free(TaggedX);
        TaggedX = NULL;
        free(TaggedY);
        TaggedY = NULL;
        free(ValueMap);
        free(TaggedZ);
        TaggedZ = NULL;
        if(pt_Keys->REG_IRREG == 0)
            Free3DVec(&TolerentZone, Order-1, 4);
        TolerentZone = NULL;
        *TemplateNums += TemplateNum;
        G_TemplateNum = TemplateNum;
        if(TemplateNum == 0)
        {
            FreeTemplate(&templates);
            templates = NULL;
            return 1;
        }


        double TempMean;
        for(int i = 1; i < Order + 1; i++)
        {
            TempMean = 0;
            for(int j = 0; j < TemplateNum; j++)
            {
                TempMean += templates[j].Value[i];
            }

            if(TempMean != 0)
            {
                for(int j = 0; j < TemplateNum; j++)
                {
                    templates[j].Value[i] -= TempMean / TemplateNum;
                }
            }
        }


        templates->index = (int*)calloc(Order+1, sizeof(int));
        if(templates == NULL)
        {
            perror("Error in function TemplateConstruct_2D: calloc return NULL\n");
            return -1;
        }

        int tempOrder = 0;
        for(int i = 2; i <= Order; i++)
        {

            templates->index[i] = pt_Keys->NOD[i-2];
            tempOrder += templates->index[i];
        }

        templates->index[1] = 1;
        templates->index[0] = ++tempOrder;



        *Templates = templates;
        G_Template = templates;

        G_Moments = NULL;
        G_Moment = NULL;
         G_CurrMomentPos = 0;
		 G_MomentFull = false;


    return 0;
}




Point Subtract(Point a, Point b)
{
    Point result = {a.x - b.x, a.y - b.y};
    return result;
}

double Dot(Point a, Point b)
{
    return a.x * b.x + a.y * b.y;
}

bool PointInRectangle(Point p, Point v0, Point v1, Point v3)
{
    Point vec = Subtract(p, v0);
    Point edge1 = Subtract(v1, v0);
    Point edge2 = Subtract(v3, v0);

    double dot00 = Dot(edge1, edge1);
    double dot01 = Dot(edge1, edge2);
    double dot11 = Dot(edge2, edge2);
    double dot0v = Dot(edge1, vec);
    double dot1v = Dot(edge2, vec);

    double denom = (double)(dot00 * dot11 - dot01 * dot01);
    if (denom == 0)
    {
        return false;
    }
    double u = (dot11 * dot0v - dot01 * dot1v) / denom;
    double v_coef = (dot00 * dot1v - dot01 * dot0v) / denom;

    return (u >= 0.0 && u <= 1.0 && v_coef >= 0.0 && v_coef <= 1.0);
}

int SearchInTolerentZone(double** vertices, int* CoorX, int* CoorY, int CoorZ)
{
    if(vertices == NULL)
    {
        fprintf(stderr, "Error in function SearchInTolerentZone: TolerentZone is NULL\n");
        return -1;
    }
    Point v[4];
    for (int i = 0; i < 4; i++)
    {
        v[i].x = vertices[i][0];
        v[i].y = vertices[i][1];
    }

    int minX = v[0].x, maxX = v[0].x;
    int minY = v[0].y, maxY = v[0].y;
    for (int i = 1; i < 4; i++)
    {
        if (v[i].x < minX) minX = v[i].x;
        if (v[i].x > maxX) maxX = v[i].x;
        if (v[i].y < minY) minY = v[i].y;
        if (v[i].y > maxY) maxY = v[i].y;
    }
    for (int x = minX; x <= maxX; x++)
    {
        for (int y = minY; y <= maxY; y++)
        {
            Point p = {x, y};
            if (PointInRectangle(p, v[0], v[1], v[3]))
            {
                if(x >= 0 && y >=0 && x < MaxX && y < MaxY)
                    if (ValueMap_3D[x][y][CoorZ] >= 0)
                    {
                        *CoorX = x;
                        *CoorY = y;
                        return 0;
                    }
            }
        }
    }
    return -1;
}


int TemplateConstruct_2D_old(int* pint32_LAG, int Order, Template** Templates, int* TemplateNums, int range1, int range2)
{
    if(pt_Keys == NULL)
    {
        fprintf(stderr, "Error in function TemplateConstruct_2D: pt_Keys is NULL\n");
        return -1;
    }
    if(pt_Grids == NULL)
    {
        fprintf(stderr, "Error in function TemplateConstruct_2D: pt_Grids is NULL\n");
        return -1;
    }
    if(ValueMap_2D == NULL)
    {
        fprintf(stderr, "Error in function TemplateConstruct_2D: ValueMap_2D is NULL\n");
        return -1;
    }

    double** ValueMap;
        ValueMap = (double**)calloc(MaxX, sizeof(double*));
        if(ValueMap == NULL)
        {
            perror("Error in function TemplateConstruct_2D: calloc return NULL\n");
            return -1;
        }
        for(int i = 0; i < MaxX; i++)
        {
            ValueMap[i] = (double*)calloc(MaxY, sizeof(double));
            if(ValueMap[i] == NULL)
            {
                perror("Error in function TemplateConstruct_2D: calloc return NULL\n");
                free(ValueMap);
                ValueMap = NULL;
                return -1;
            }
        }

        int TemplateNum = 0;
        int Templated = 0;
        int TemplateSize = range2;

        Template* templates = (Template*)calloc(TemplateSize, sizeof(Template));
        if(templates == NULL)
        {
            perror("Error in function TemplateConstruct_2D: calloc return NULL\n");
            if(Free2DVec(&ValueMap, pt_Keys->NX) == -1)
            {
                fprintf(stderr, "Failed to free 2DVec\n");
            }
            ValueMap = NULL;
            return -1;
        }

        for(int i = range1; i < range2; i++)
        {
            templates[i].Lag = (int*)calloc(Order + 1, sizeof(int));
            templates[i].Value = (double*)calloc(Order + 1, sizeof(double));
            if(templates[i].Lag == NULL || templates[i].Value == NULL)
            {
                perror("Error in function TemplateConstruct_2D: calloc return NULL\n");
                Free2DVec(&ValueMap, pt_Keys->NX);
                free(templates);
                templates = NULL;
                ValueMap = NULL;
                return -1;
            }

        }


        double*** TolerentZone;
        if(Init3DVec(&TolerentZone, Order - 1, 4, 2) == -1)
        {
            fprintf(stderr, "Failed to free 3DVec\n");
        }

        int* TaggedX = (int*)calloc(Order-1, sizeof(int));
        int* TaggedY = (int*)calloc(Order-1, sizeof(int));
        if(TaggedX == NULL || TaggedY == NULL)
        {
            perror("Error in function TemplateConstruct_2D: calloc return NULL\n");
            FreeTemplate(&templates);
            Free2DVec(&ValueMap, pt_Keys->NX);
            Free3DVec(&TolerentZone, Order-1, 4);
            templates = NULL;
            ValueMap = NULL;
            TolerentZone = NULL;
            return -1;
        }

        for(int i = 0; i < TemplateSize; i++)
        {

            if(ValueMap[pt_Grids->pt_Nodes[i].x][pt_Grids->pt_Nodes[i].y] < 0)
            {
                continue;
            }

            Templated = 0;

            for(int j = 0; j < Order; j++)
            {
                if(j == Order - 1 && Templated == Order - 1)
                {
                    ValueMap[pt_Grids->pt_Nodes[i].x][pt_Grids->pt_Nodes[i].y] = -i;
                    for(int k = 0; k < Order - 1; k++)
                    {
                        ValueMap[TaggedX[k]][TaggedY[k]] = -i;
                    }

                    templates[TemplateNum].Value[1] = pt_Grids->pt_Nodes[i].val;
                    templates[TemplateNum].Lag[1] = 0;

                    TemplateNum++;

                    break;
                }
                else if(j == Order - 1 && Templated < Order - 1)
                {
                    break;
                }
                int LAG = pint32_LAG[j];
                templates[TemplateNum].Lag[j+2] = LAG;
                int LAG_t = pt_Keys->LTOL[j];
                double ATOL = pt_Keys->ATOL[j] * PI / 180;
                double AZM = pt_Keys->AZM[j] * PI / 180;


                int lag = LAG - LAG_t;

                int MiddleZone_x = pt_Grids->pt_Nodes[i].x + LAG * cos(AZM);
                int MiddleZone_y = pt_Grids->pt_Nodes[i].y + LAG * sin(AZM);


                TolerentZone[j][0][0] = pt_Grids->pt_Nodes[i].x + lag / cos(ATOL) * cos(AZM + ATOL);
                TolerentZone[j][0][1] = pt_Grids->pt_Nodes[i].y + lag / cos(ATOL) * sin(AZM + ATOL);
                TolerentZone[j][1][0] = pt_Grids->pt_Nodes[i].x + lag / cos(ATOL) * cos(AZM - ATOL);
                TolerentZone[j][1][1] = pt_Grids->pt_Nodes[i].y + lag / cos(ATOL) * sin(AZM - ATOL);
                TolerentZone[j][2][0] = 2 * MiddleZone_x - TolerentZone[j][0][0];
                TolerentZone[j][2][1] = 2 * MiddleZone_y - TolerentZone[j][0][1];
                TolerentZone[j][3][0] = 2 * MiddleZone_x - TolerentZone[j][1][0];
                TolerentZone[j][3][1] = 2 * MiddleZone_y - TolerentZone[j][1][1];


                if(TolerentZone[j][0][0] >= 0 && TolerentZone[j][0][1] >= 0 && TolerentZone[j][0][0] <= pt_Keys->DIMX && TolerentZone[j][0][1] <= pt_Keys->DIMY)
                {
                    if(ValueMap[(int)TolerentZone[j][0][0]][(int)TolerentZone[j][0][1]] >= 0)
                    {
                        if(ValueMap_2D[(int)TolerentZone[j][0][0]][(int)TolerentZone[j][0][1]] != -1)
                        {


                            templates[TemplateNum].Value[j+2] = ValueMap_2D[(int)TolerentZone[j][0][0]][(int)TolerentZone[j][0][1]];

                            TaggedX[j] = (int)TolerentZone[j][0][0];
                            TaggedY[j] = (int)TolerentZone[j][0][1];
                            Templated++;
                            continue;
                        }
                        else
                            break;
                    }
                    else
                        break;
                }
                else if(TolerentZone[j][1][0] >= 0 && TolerentZone[j][1][1] >= 0 && TolerentZone[j][1][0] <= pt_Keys->DIMX && TolerentZone[j][1][1] <= pt_Keys->DIMY)
                {
                    if(ValueMap[(int)TolerentZone[j][1][0]][(int)TolerentZone[j][1][1]] >= 0)
                    {
                        if(ValueMap_2D[(int)TolerentZone[j][1][0]][(int)TolerentZone[j][1][1]] != -1)
                        {


                            templates[TemplateNum].Value[j+2] = ValueMap_2D[(int)TolerentZone[j][1][0]][(int)TolerentZone[j][1][1]];

                            TaggedX[j] = (int)TolerentZone[j][1][0];
                            TaggedY[j] = (int)TolerentZone[j][1][1];
                            Templated++;
                            continue;
                        }
                        else
                            break;
                    }
                    else
                        break;
                }

                else if(TolerentZone[j][2][0] >= 0 && TolerentZone[j][2][1] >= 0 && TolerentZone[j][2][0] <= pt_Keys->DIMX && TolerentZone[j][2][1] <= pt_Keys->DIMY)
                {
                    if(ValueMap[(int)TolerentZone[j][2][0]][(int)TolerentZone[j][2][1]] >= 0)
                    {
                        if(ValueMap_2D[(int)TolerentZone[j][2][0]][(int)TolerentZone[j][2][1]] != -1)
                        {



                            templates[TemplateNum].Value[j+2] = ValueMap_2D[(int)TolerentZone[j][2][0]][(int)TolerentZone[j][2][1]];

                            TaggedX[j] = (int)TolerentZone[j][2][0];
                            TaggedY[j] = (int)TolerentZone[j][2][1];
                            Templated++;
                            continue;
                        }
                        else
                            break;
                    }
                    else
                        break;
                }
                else if(TolerentZone[j][3][0] >= 0 && TolerentZone[j][3][1] >= 0 && TolerentZone[j][3][0] <= pt_Keys->DIMX && TolerentZone[j][3][1] <= pt_Keys->DIMY)
                {
                    if(ValueMap[(int)TolerentZone[j][3][0]][(int)TolerentZone[j][3][1]] >= 0)
                    {
                        if(ValueMap_2D[(int)TolerentZone[j][3][0]][(int)TolerentZone[j][3][1]] != -1)
                        {



                            templates[TemplateNum].Value[j+2] = ValueMap_2D[(int)TolerentZone[j][3][0]][(int)TolerentZone[j][3][1]];

                            TaggedX[j] = (int)TolerentZone[j][3][0];
                            TaggedY[j] = (int)TolerentZone[j][3][1];
                            Templated++;
                            continue;
                        }
                        else
                            break;
                    }
                    else
                        break;
                }

                else
                {
                    break;
                }
            }
        }


        Free2DVec(&ValueMap, pt_Keys->NX);
        ValueMap = NULL;
        free(TaggedX);
        TaggedX = NULL;
        free(TaggedY);
        TaggedY = NULL;
        free(ValueMap);
        Free3DVec(&TolerentZone, Order-1, 4);
        TolerentZone = NULL;
        *TemplateNums += TemplateNum;
        G_TemplateNum = TemplateNum;
        if(TemplateNum == 0)
        {
            FreeTemplate(&templates);
            templates = NULL;
            return 1;
        }


        double TempMean;
        for(int i = 0; i < Order + 1; i++)
        {
            TempMean = 0;
            for(int j = 0; j < TemplateNum; j++)
            {
                TempMean += templates[j].Value[i];
            }

            if(TempMean != 0)
            {
                for(int j = 0; j < TemplateNum; j++)
                {
                    templates[j].Value[i] -= TempMean / TemplateNum;
                }
            }
        }


        templates->index = (int*)calloc(Order, sizeof(int));
        if(templates == NULL)
        {
            perror("Error in function TemplateConstruct_2D: calloc return NULL\n");
            return -1;
        }

        for(int i = 1; i <= Order; i++)
        {
            templates->index[i-1] = i;
        }


        *Templates = templates;
        G_Template = templates;

        G_Moments = GetAllMoments();

    return 0;
}

int TemplateConstruct_3D_old(int* pint32_LAG, int Order, Template** Templates, int* TemplateNums, int range1, int range2)
{
    if(pt_Keys == NULL)
    {
        fprintf(stderr, "Error in function TemplateConstruct_2D: pt_Keys is NULL\n");
        return -1;
    }
    if(pt_Grids == NULL)
    {
        fprintf(stderr, "Error in function TemplateConstruct_2D: pt_Grids is NULL\n");
        return -1;
    }
    if(ValueMap_3D == NULL)
    {
        fprintf(stderr, "Error in function TemplateConstruct_2D: ValueMap_2D is NULL\n");
        return -1;
    }

    double*** ValueMap;
        ValueMap = (double***)calloc(MaxX, sizeof(double**));
        if(ValueMap == NULL)
        {
            perror("Error in function TemplateConstruct_3D: calloc return NULL\n");
            return -1;
        }
        for(int i = 0; i < MaxX; i++)
        {
            ValueMap[i] = (double**)malloc(MaxY * sizeof(double*));
            if(ValueMap[i] == NULL)
            {
                fprintf(stderr, "Error in function TemplateConstruct_3D: malloc return -1\n");
                free(ValueMap);
                ValueMap = NULL;
                return -1;
            }
            for(int j = 0; j < MaxY; j++)
            {
                ValueMap[i][j] = (double*)malloc(MaxZ * sizeof(double));
                if(ValueMap[i][j] == NULL)
                {
                    fprintf(stderr, "Error in function TemplateConstruct_3D: malloc return -1\n");
                    free(ValueMap);
                    ValueMap = NULL;
                    return -1;
                }
                for(int k = 0; k < MaxZ; k++)
                {
                    ValueMap[i][j][k] = ValueMap_3D[i][j][k];
                }
            }
        }

        int TemplateNum = 0;
        int Templated = 0;
        int TemplateSize = range2;

        Template* templates = (Template*)calloc(TemplateSize, sizeof(Template));
        if(templates == NULL)
        {
            perror("Error in function TemplateConstruct_2D: calloc return NULL\n");
            if(Free3DVec(&ValueMap, MaxX, MaxY) == -1)
            {
                fprintf(stderr, "Failed to free 2DVec\n");
            }
            ValueMap = NULL;
            return -1;
        }

        for(int i = range1; i < range2; i++)
        {
            templates[i].Lag = (int*)calloc(Order + 1, sizeof(int));
            templates[i].Value = (double*)calloc(Order + 1, sizeof(double));
            if(templates[i].Lag == NULL || templates[i].Value == NULL)
            {
                perror("Error in function TemplateConstruct_2D: calloc return NULL\n");
                Free3DVec(&ValueMap, MaxX, MaxY);
                free(templates);
                templates = NULL;
                ValueMap = NULL;
                return -1;
            }

        }


        double*** TolerentZone;
        if(pt_Keys->REG_IRREG == 0)
        {
            if(Init3DVec(&TolerentZone, Order - 1, 4, 3) == -1)
            {
                fprintf(stderr, "Failed to free 3DVec\n");
            }
        }
        int* TaggedX = (int*)calloc(Order-1, sizeof(int));
        int* TaggedY = (int*)calloc(Order-1, sizeof(int));
        int* TaggedZ = (int*)calloc(Order-1, sizeof(int));
        if(TaggedX == NULL || TaggedY == NULL || TaggedZ == NULL)
        {
            perror("Error in function TemplateConstruct_2D: calloc return NULL\n");
            FreeTemplate(&templates);
            Free3DVec(&ValueMap, MaxX, MaxY);
            if(pt_Keys->REG_IRREG == 0)
                Free3DVec(&TolerentZone, Order-1, 4);
            templates = NULL;
            ValueMap = NULL;
            TolerentZone = NULL;
            return -1;
        }

        for(int i = 0; i < TemplateSize; i++)
        {

            if(ValueMap[pt_Grids->pt_Nodes[i].x][pt_Grids->pt_Nodes[i].y][pt_Grids->pt_Nodes[i].z] < 0)
            {
                continue;
            }

            Templated = 0;

            for(int j = 0; j < Order; j++)
            {
                if(j == Order - 1 && Templated == Order - 1)
                {
                    ValueMap[pt_Grids->pt_Nodes[i].x][pt_Grids->pt_Nodes[i].y][pt_Grids->pt_Nodes[i].z] = -i;
                    for(int k = 0; k < Order - 1; k++)
                    {
                        ValueMap[TaggedX[k]][TaggedY[k]][TaggedZ[k]] = -i;
                    }

                    templates[TemplateNum].Value[1] = pt_Grids->pt_Nodes[i].val;
                    templates[TemplateNum].Lag[1] = 0;

                    TemplateNum++;

                    break;
                }
                else if(j == Order - 1 && Templated < Order - 1)
                {
                    break;
                }
                if(pt_Keys->REG_IRREG == 0)
                {
                    int LAG = pint32_LAG[j];
                    templates[TemplateNum].Lag[j+2] = LAG;
                    int LAG_t = pt_Keys->LTOL[j];

                    double ATOL = pt_Keys->ATOL[j] * PI / 180;
                    double AZM = pt_Keys->AZM[j] * PI / 180;

                    double ELEVATION = pt_Keys->DIP[j] * PI / 180;





                    int lag = LAG - LAG_t;
                    int TZoneHeight = 2 * LAG_t;
                    int MiddleZone_x = pt_Grids->pt_Nodes[i].x + LAG * cos(AZM);
                    int MiddleZone_y = pt_Grids->pt_Nodes[i].y + LAG * sin(AZM);


                    TolerentZone[j][0][0] = pt_Grids->pt_Nodes[i].x + lag / cos(ATOL) * cos(AZM + ATOL);
                    TolerentZone[j][0][1] = pt_Grids->pt_Nodes[i].y + lag / cos(ATOL) * sin(AZM + ATOL);
                    TolerentZone[j][0][2] = pt_Grids->pt_Nodes[i].z + LAG * sin(ELEVATION);
                    TolerentZone[j][1][0] = pt_Grids->pt_Nodes[i].x + lag / cos(ATOL) * cos(AZM - ATOL);
                    TolerentZone[j][1][1] = pt_Grids->pt_Nodes[i].y + lag / cos(ATOL) * sin(AZM - ATOL);
                    TolerentZone[j][1][2] = pt_Grids->pt_Nodes[i].z + LAG * sin(ELEVATION);
                    TolerentZone[j][2][0] = 2 * MiddleZone_x - TolerentZone[j][0][0];
                    TolerentZone[j][2][1] = 2 * MiddleZone_y - TolerentZone[j][0][1];
                    TolerentZone[j][2][2] = pt_Grids->pt_Nodes[i].z + LAG * sin(ELEVATION);
                    TolerentZone[j][3][0] = 2 * MiddleZone_x - TolerentZone[j][1][0];
                    TolerentZone[j][3][1] = 2 * MiddleZone_y - TolerentZone[j][1][1];
                    TolerentZone[j][3][2] = pt_Grids->pt_Nodes[i].z + LAG * sin(ELEVATION);

                    #if 0

                    int NodeX = -1;
                    int NodeY = -1;
                    int NodeZ = pt_Grids->pt_Nodes[i].z + LAG * sin(ELEVATION);

                    bool FoundZ = false;
                    for(int k = 0; k < pt_Grids->ND; k++)
                    {
                        if(abs(pt_Grids->pt_Nodes[i].z - NodeZ) <= pt_Keys->BANDWD[j])
                        {
                            NodeZ = pt_Grids->pt_Nodes[i].z;
                            FoundZ = true;
                            break;
                        }
                        else
                            continue;
                    }
                    if(FoundZ == false)
                        break;
                    if(SearchInTolerentZone(TolerentZone[j], &NodeX, &NodeY, NodeZ) == 0)
                    {
                        templates[TemplateNum].Value[j+2] = ValueMap_3D[NodeX][NodeY][NodeZ];
                        TaggedX[j] = NodeX;
                        TaggedY[j] = NodeY;
                        TaggedZ[j] = NodeZ;
                        Templated++;

                        continue;
                    }
                    #endif




                    if(TolerentZone[j][0][0] >= 0 && TolerentZone[j][0][1] >= 0 && TolerentZone[j][0][2] >= 0 && TolerentZone[j][0][0] <= pt_Keys->DIMX && TolerentZone[j][0][1] <= pt_Keys->DIMY && TolerentZone[j][0][2] <= pt_Keys->DIMZ + pt_Grids->pt_Nodes[0].z)
                    {

                        if(ValueMap[(int)TolerentZone[j][0][0]][(int)TolerentZone[j][0][1]][(int)TolerentZone[j][0][2]] >= 0)
                        if(ValueMap_3D[(int)TolerentZone[j][0][0]][(int)TolerentZone[j][0][1]][(int)TolerentZone[j][0][2]] != -1)
                        {


                            templates[TemplateNum].Value[j+2] = ValueMap_3D[(int)TolerentZone[j][0][0]][(int)TolerentZone[j][0][1]][(int)TolerentZone[j][0][2]];

                            TaggedX[j] = (int)TolerentZone[j][0][0];
                            TaggedY[j] = (int)TolerentZone[j][0][1];
                            TaggedZ[j] = (int)TolerentZone[j][0][2];
                            Templated++;
                            continue;
                        }
                    }
                    else if(TolerentZone[j][1][0] >= 0 && TolerentZone[j][1][1] >= 0 && TolerentZone[j][1][2] >= 0 && TolerentZone[j][1][0] <= pt_Keys->DIMX && TolerentZone[j][1][1] <= pt_Keys->DIMY && TolerentZone[j][1][2] <= pt_Keys->DIMZ + pt_Grids->pt_Nodes[0].z)
                    {
                        if(ValueMap[(int)TolerentZone[j][1][0]][(int)TolerentZone[j][1][1]][(int)TolerentZone[j][1][2]] >= 0)
                        if(ValueMap_3D[(int)TolerentZone[j][1][0]][(int)TolerentZone[j][1][1]][(int)TolerentZone[j][1][2]] != -1)
                        {


                            templates[TemplateNum].Value[j+2] = ValueMap_3D[(int)TolerentZone[j][1][0]][(int)TolerentZone[j][1][1]][(int)TolerentZone[j][1][2]];

                            TaggedX[j] = (int)TolerentZone[j][1][0];
                            TaggedY[j] = (int)TolerentZone[j][1][1];
                            TaggedZ[j] = (int)TolerentZone[j][1][2];
                            Templated++;
                            continue;
                        }
                    }
                    else if(TolerentZone[j][2][0] >= 0 && TolerentZone[j][2][1] >= 0 && TolerentZone[j][2][2] >= 0 && TolerentZone[j][2][0] <= pt_Keys->DIMX && TolerentZone[j][2][1] <= pt_Keys->DIMY && TolerentZone[j][2][2] <= pt_Keys->DIMZ + pt_Grids->pt_Nodes[0].z)
                    {
                        if(ValueMap[(int)TolerentZone[j][2][0]][(int)TolerentZone[j][2][1]][(int)TolerentZone[j][2][2]] >= 0)
                        if(ValueMap_3D[(int)TolerentZone[j][2][0]][(int)TolerentZone[j][2][1]][(int)TolerentZone[j][2][2]] != -1)
                        {



                            templates[TemplateNum].Value[j+2] = ValueMap_3D[(int)TolerentZone[j][2][0]][(int)TolerentZone[j][2][1]][(int)TolerentZone[j][2][2]];

                            TaggedX[j] = (int)TolerentZone[j][2][0];
                            TaggedY[j] = (int)TolerentZone[j][2][1];
                            TaggedZ[j] = (int)TolerentZone[j][2][2];
                            Templated++;
                            continue;
                        }
                    }
                    else if(TolerentZone[j][3][0] >= 0 && TolerentZone[j][3][1] >= 0 && TolerentZone[j][3][2] >= 0 && TolerentZone[j][3][0] <= pt_Keys->DIMX && TolerentZone[j][3][1] <= pt_Keys->DIMY && TolerentZone[j][3][2] <= pt_Keys->DIMZ + pt_Grids->pt_Nodes[0].z)
                    {
                        if(ValueMap[(int)TolerentZone[j][3][0]][(int)TolerentZone[j][3][1]][(int)TolerentZone[j][3][2]] >= 0)
                        if(ValueMap_3D[(int)TolerentZone[j][3][0]][(int)TolerentZone[j][3][1]][(int)TolerentZone[j][3][2]] != -1)
                        {



                            templates[TemplateNum].Value[j+2] = ValueMap_3D[(int)TolerentZone[j][3][0]][(int)TolerentZone[j][3][1]][(int)TolerentZone[j][3][2]];

                            TaggedX[j] = (int)TolerentZone[j][3][0];
                            TaggedY[j] = (int)TolerentZone[j][3][1];
                            TaggedZ[j] = (int)TolerentZone[j][3][2];
                            Templated++;
                            continue;
                        }
                    }

                    else
                    {
                        break;
                    }
                }
                else
                {
                    double TempX = 0;
                    double TempY = 0;
                    double TempZ = 0;
                    double AZM = pt_Keys->AZM[j] * PI / 180;
                    double ELEVATION = pt_Keys->DIP[j] * PI / 180;
                    int LAG = pint32_LAG[j];
                    TempX = pt_Grids->pt_Nodes[i].x + LAG * cos(AZM);
                    TempY = pt_Grids->pt_Nodes[i].y + LAG * sin(AZM);
                    TempZ = pt_Grids->pt_Nodes[i].z + LAG * sin(ELEVATION);
                    if(TempX >= 0 && TempY >= 0 && TempZ >= 0 && TempX <= pt_Keys->DIMX && TempY <= pt_Keys->DIMY && TempZ <= pt_Keys->DIMZ + pt_Grids->pt_Nodes[0].z)
                    {

                        if(ValueMap[(int)TempX][(int)TempY][(int)TempZ] >= 0)
                        if(ValueMap_3D[(int)TempX][(int)TempY][(int)TempZ] != -1)
                        {


                            templates[TemplateNum].Value[j+2] = ValueMap_3D[(int)TempX][(int)TempY][(int)TempZ];

                            TaggedX[j] = (int)TempX;
                            TaggedY[j] = (int)TempY;
                            TaggedZ[j] = (int)TempZ;
                            Templated++;
                            continue;
                        }
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }


        Free3DVec(&ValueMap, MaxX, MaxY);
        ValueMap = NULL;
        free(TaggedX);
        TaggedX = NULL;
        free(TaggedY);
        TaggedY = NULL;
        free(ValueMap);
        free(TaggedZ);
        TaggedZ = NULL;
        if(pt_Keys->REG_IRREG == 0)
            Free3DVec(&TolerentZone, Order-1, 4);
        TolerentZone = NULL;
        *TemplateNums += TemplateNum;
        G_TemplateNum = TemplateNum;
        if(TemplateNum == 0)
        {
            FreeTemplate(&templates);
            templates = NULL;
            return 1;
        }


        double TempMean;
        for(int i = 0; i < Order + 1; i++)
        {
            TempMean = 0;
            for(int j = 0; j < TemplateNum; j++)
            {
                TempMean += templates[j].Value[i];
            }

            if(TempMean != 0)
            {
                for(int j = 0; j < TemplateNum; j++)
                {
                    templates[j].Value[i] -= TempMean / TemplateNum;
                }
            }
        }


        templates->index = (int*)calloc(Order, sizeof(int));
        if(templates == NULL)
        {
            perror("Error in function TemplateConstruct_2D: calloc return NULL\n");
            return -1;
        }
        for(int i = 1; i <= Order; i++)
        {
            templates->index[i-1] = i;
        }
        *Templates = templates;
        G_Template = templates;
        G_Moments = GetAllMoments();
    return 0;
}
