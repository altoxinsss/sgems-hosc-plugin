#include "Initialization.h"
#include "DebugLogger.h"


#define MaxStringLen 300
#define MaxNumLen 30



#define HOSCMaxStringLen 300
#define HOSCMaxNumLen 30

double** ValueMap_2D;
double*** ValueMap_3D;
Keys* pt_Keys;
Grid* pt_Grids;

int MaxX;
int MaxY;
int MaxZ;


int TotalNodeNum;

int InitialStruct()
{
    char* FileKeysName;
    char* FileGridName;
    printf("GIVE THE NAME OF INPUT DATA FILE (WITHOUT EXTENSION):");
    FileKeysName = (char*)calloc(MaxStringLen, sizeof(char));
    FileGridName = (char*)calloc(MaxStringLen, sizeof(char));
    assert(FileKeysName != NULL && FileGridName != NULL);
    scanf("%s", FileKeysName);
    strcpy(FileGridName, FileKeysName);

    assert(strncat(FileKeysName, ".KEYS", MaxStringLen-strlen(".KEYS")+1) != NULL);
    assert(strncat(FileGridName, ".GRID", MaxStringLen-strlen(".GRID")+1) != NULL);

    FILE* FileKeys = fopen(FileKeysName, "r");
    FILE* FileGrid = fopen(FileGridName, "r");

    assert(FileKeys != NULL && FileGrid != NULL);

    if(GetGrid(FileGrid) == -1)
    {
        fprintf(stderr, "Faild to get Grid\n");
        return -1;
    }

    if(GetKeys(FileKeys) == -1)
    {
        fprintf(stderr, "Faild to get keys\n");
        return -1;
    }

    if(FileKeys == NULL || FileGrid == NULL || FileKeysName == NULL || FileGridName == NULL)
        return -1;
    fclose(FileKeys);
    fclose(FileGrid);
    free(FileKeysName);
    free(FileGridName);

    return 0;
}

int GetKeys(FILE* FileKeys)
{

    if(FileKeys == NULL || pt_Keys == NULL || pt_Grids == NULL)
	{
        return -1;
	}



    char* pc_TempStr = (char*)calloc(MaxStringLen, sizeof(char));
    char* pc_TempParam = (char*)calloc(MaxNumLen, sizeof(char));
    if(pc_TempParam == NULL)
        return -1;


    int* p_TempKeys = (int* )pt_Keys;


    for(int i = 0; i < 10; i++)
    {
        if(i < 3 || i > 7)
            fgets(pc_TempStr, MaxStringLen, FileKeys);

        for(int j = 0; j < MaxNumLen; j++)
        {
            pc_TempParam[j] = fgetc(FileKeys);
            if(pc_TempParam[j] == ' ' || pc_TempParam[j] == '\0')
            {
                pc_TempParam[j] = '\0';
                break;
            }
        }
        *p_TempKeys = atoi(pc_TempParam);
        memset(pc_TempParam, 0, MaxNumLen*sizeof(char));
       p_TempKeys++;
    }

    if(pt_Keys->NDIR < 0 || pt_Keys->NDIR !=pt_Keys->ORDER - 1)
    {
        printf("NDIR IS TOO SMALL: CHECK PARAMETERS\n");
        return -1;
    }



    pt_Keys->NLAG = (int* )calloc(pt_Keys->NDIR, sizeof(int));
    pt_Keys->LAG = (int* )calloc(pt_Keys->NDIR, sizeof(int));
    pt_Keys->LTOL = (int* )calloc(pt_Keys->NDIR, sizeof(int));
    pt_Keys->AZM = (int* )calloc(pt_Keys->NDIR, sizeof(int));


    pt_Keys->ATOL = (int* )calloc(pt_Keys->NDIR, sizeof(int));
    pt_Keys->BANDWH = (int* )calloc(pt_Keys->NDIR, sizeof(int));
    pt_Keys->DIP = (int* )calloc(pt_Keys->NDIR, sizeof(int));
    pt_Keys->DTOL = (int* )calloc(pt_Keys->NDIR, sizeof(int));
    pt_Keys->BANDWD = (int* )calloc(pt_Keys->NDIR, sizeof(int));
    bool TES = false;







    for(int k = 0; k< 3; k++)
    {

            fgets(pc_TempStr, MaxStringLen, FileKeys);
        for(int i = 0; i < pt_Keys->NDIR; i++)
        {
            memset(pc_TempParam, 0, MaxNumLen*sizeof(char));
            for(int j = 0; j < MaxNumLen; j++)
            {
                pc_TempParam[j] = fgetc(FileKeys);
                if(pc_TempParam[j] == ' ' || pc_TempParam[j] == '\n')
                {
                    pc_TempParam[j] = '\0';



                    break;
                }
            }
            if(k == 0)
            {
                pt_Keys->NLAG[i] = atoi(pc_TempParam);
                if(pt_Keys->NLAG[i] > 1)
                    TES = true;
            }
            else if(k == 1)
            {
                pt_Keys->LAG[i] = atoi(pc_TempParam);
                if(pt_Keys->LAG[i] > 0)
                    TES = true;
            }
            else if(k == 2)
                pt_Keys->LTOL[i] = atoi(pc_TempParam);
        }
        if(TES == false && k == 0 && pt_Keys->ORDER != 1)
        {
            printf("XNLAG, YNLOG AND ZNLOG ARE TOO SMALL: CHECK PARAMETERS\n");
            return -1;
        }
        else if(TES == false && k == 1 && pt_Keys->ORDER != 1)
        {
            printf("XLAG AND YLOG AND ZLOG IS TOO SMALL: CHECK PARAMETER FILE\n");
            return -1;
        }
    }

    fgets(pc_TempStr, MaxStringLen, FileKeys);

    for(int i = 0; i < pt_Keys->NDIR; i++)
    {

        for(int k = 0; k < 6; k++)
        {
            memset(pc_TempParam, 0, MaxNumLen*sizeof(char));
            for(int j = 0; j < MaxNumLen; j++)
            {
                pc_TempParam[j] = fgetc(FileKeys);


                if(pc_TempParam[j] == ' ' || pc_TempParam[j] == '\n')
                {
                    pc_TempParam[j] = '\0';


                    break;
                }
            }
            #if 0
            if(k == 0)
                pt_Keys->AZM[i] = atoi(pc_TempParam);
            else if(k == 2)
            {
                pt_Keys->ELEVATION[i] = atoi(pc_TempParam);
                if(pt_Keys->DIM == 2 && pt_Keys->ELEVATION[i] != 0)
                {
                    fprintf(stderr, "In 2D space, Z is not 0\n");
                    return -1;
                }
            }
            else if(k == 1)
                pt_Keys->ATOL[i] = atoi(pc_TempParam);
            else if(k == 3)
                pt_Keys->BANDWH[i] = atoi(pc_TempParam);
            else if(k == 4)
                pt_Keys->DIP[i] = atoi(pc_TempParam);
            else if(k == 5)
                pt_Keys->DTOL[i] = atoi(pc_TempParam);
            else if(k == 6)
                pt_Keys->BANDWD[i] = atoi(pc_TempParam);
            #endif
            if(k == 0)
                pt_Keys->AZM[i] = atoi(pc_TempParam);
            else if(k == 1)
                pt_Keys->ATOL[i] = atoi(pc_TempParam);
            else if(k == 2)
                pt_Keys->BANDWH[i] = atoi(pc_TempParam);
            else if(k == 3)
                pt_Keys->DIP[i] = atoi(pc_TempParam);
            else if(k == 4)
                pt_Keys->DTOL[i] = atoi(pc_TempParam);
            else if(k == 5)
                pt_Keys->BANDWD[i] = atoi(pc_TempParam);
        }
    }

    if(pt_Keys->DIM == 2)
    {
        ValueMap_2D = (double**)calloc(pt_Keys->NX, sizeof(double*));
        if(ValueMap_2D == NULL)
            return -1;
        for(int i = 0; i < pt_Keys->NX; i++)
        {
            ValueMap_2D[i] = (double*)malloc(pt_Keys->NY * sizeof(double));
            if(ValueMap_2D[i] == NULL)
                return -1;
            memset(ValueMap_2D[i], -1, pt_Keys->NY * sizeof(double));
        }
        for(int i = 0; i < pt_Grids->ND; i++)
        {
            ValueMap_2D[pt_Grids->pt_Nodes[i].x][pt_Grids->pt_Nodes[i].y] = pt_Grids->pt_Nodes[i].val;
        }
    }
    else
    {
        MaxX = pt_Keys->NX;
        MaxY = pt_Keys->NY;
        MaxZ = pt_Keys->NZ;
        if(pt_Grids->pt_Nodes == NULL)
        {
            fprintf(stderr, "Error in function GetKeys: ");
            return -1;
        }
        GetMaxElem(pt_Grids->pt_Nodes, pt_Grids->ND);
        ValueMap_3D = (double***)calloc(MaxX, sizeof(double**));
        if(ValueMap_3D == NULL)
            return -1;
        for(int i = 0; i < MaxX; i++)
        {
            ValueMap_3D[i] = (double**)calloc(MaxY, sizeof(double*));
            if(ValueMap_3D[i] == NULL)
                return -1;
            for(int j = 0; j < MaxY; j++)
            {
                ValueMap_3D[i][j] = (double*)malloc(MaxZ * sizeof(double));
                if(ValueMap_3D[i][j] == NULL)
                    return -1;
				for (int k = 0; k < MaxZ; ++k)
				{
					ValueMap_3D[i][j][k] = -1.0;
				}
            }
        }
        for(int i = 0; i < pt_Grids->ND; i++)
        {
            ValueMap_3D[pt_Grids->pt_Nodes[i].x][pt_Grids->pt_Nodes[i].y][pt_Grids->pt_Nodes[i].z] = pt_Grids->pt_Nodes[i].val;
        }
    }

    free(pc_TempParam);
    free(pc_TempStr);
    return 0;
}

int GetGrid(FILE* FileGrid)
{
    if(FileGrid == NULL || pt_Grids == NULL)
	{
        return -1;
	}
    int TRANS = 0;
    char* pc_TempStr = (char*)calloc(MaxStringLen, sizeof(char));
    char* pc_TempParam = (char*)calloc(MaxNumLen, sizeof(char));
    fgets(pc_TempStr, MaxStringLen, FileGrid);

    for(int k = 0; k < 2; k++)
    {
        memset(pc_TempParam, 0, MaxNumLen*sizeof(char));
        while(true)
        {
            pc_TempParam[0] = fgetc(FileGrid);
            if(pc_TempParam[0] != ' ')
                break;
            if(pc_TempParam[0] == '\0')
            {
                printf("Didnt find index info");

				FILE* test1;
				test1 = fopen("t1", "w");

                return -1;
            }
        }
        for(int i = 1; i < MaxNumLen; i++)
        {
            pc_TempParam[i] = fgetc(FileGrid);
            if(pc_TempParam[i] == ' ')
            {
                pc_TempParam[i] = '\0';
                break;
            }
        }
        if(k == 0)
            pt_Grids->ND = atoi(pc_TempParam);
        if(k == 1)
            pt_Grids->NVAR = atoi(pc_TempParam);
    }
    fgets(pc_TempStr, MaxStringLen, FileGrid);
    fgets(pc_TempStr, MaxStringLen, FileGrid);
    TotalNodeNum = pt_Grids->ND;


    pt_Grids->pt_Nodes = (Nodes*)calloc(pt_Grids->ND, sizeof(Nodes));
    for(int i = 0; i < pt_Grids->ND; i++)
    {
        for(int j = 0; j < pt_Grids->NVAR + 1; j++)
        {
            memset(pc_TempParam, 0, MaxNumLen*sizeof(char));
            while(true)
            {
                pc_TempParam[0] = fgetc(FileGrid);
                if(pc_TempParam[0] != ' ')
                    break;
                if(pc_TempParam[0] == '\0')
                {

					FILE* test2;
					test2 = fopen("t2", "w");

                    return -1;
                }
            }
            if(j != 4)
                for(int k = 1; k < MaxNumLen; k++)
                {
                    pc_TempParam[k] = fgetc(FileGrid);
                    if(pc_TempParam[k] == ' ' || '\0')
                    {
                        pc_TempParam[k] = '\0';
                        break;
                    }
                }
            else if(j == 4)
            {
                fscanf(FileGrid, "%lf", &pt_Grids->pt_Nodes[i].val);
                fgets(pc_TempStr, MaxStringLen, FileGrid);
            }
            if(j == 0)
                pt_Grids->pt_Nodes[i].index = atoi(pc_TempParam);
            else if(j == 1)
            {
                pt_Grids->pt_Nodes[i].x = atoi(pc_TempParam);
                if(pt_Grids->pt_Nodes[i].x == -1)
                    return -1;
                if(i == 0 && pt_Grids->pt_Nodes[i].x == 1)
                {
                    pt_Grids->pt_Nodes[i].x--;
                    TRANS = 1;
                }
                else if(TRANS == 1)
                {
                    pt_Grids->pt_Nodes[i].x--;
                }

            }
            else if(j == 2)
            {
                pt_Grids->pt_Nodes[i].y = atoi(pc_TempParam);
                if(pt_Grids->pt_Nodes[i].y == -1)
                    return -1;
                if(TRANS == 1)
                    pt_Grids->pt_Nodes[i].y--;
            }
            else if(j == 3)
            {
                pt_Grids->pt_Nodes[i].z = atoi(pc_TempParam);
                if(TRANS == 1)
                    pt_Grids->pt_Nodes[i].z--;
            }
            else if(j == 4)
            {
                pt_Grids->pt_Nodes[i].val = atof(pc_TempParam);
                #if 0
                if(pt_Keys->DIM == 2)
                {
                    if(ValueMap_2D == NULL)
                        return -1;
                    ValueMap_2D[pt_Grids->pt_Nodes[i].x][pt_Grids->pt_Nodes[i].y] = pt_Grids->pt_Nodes[i].val;
                }
                else if(pt_Keys->DIM == 3)
                {
                    if(ValueMap_3D == NULL)
                        return -1;
                    ValueMap_3D[pt_Grids->pt_Nodes[i].x][pt_Grids->pt_Nodes[i].y][pt_Grids->pt_Nodes[i].z] = pt_Grids->pt_Nodes[i].val;
                }
                #endif
            }
            pt_Grids->pt_Nodes[i].Flag_Used = false;
        }

    }

    free(pc_TempParam);
    free(pc_TempStr);

    return 0;
}

int GetMaxElem(Nodes* pt_Nodes, int NodeNum)
{
    if(pt_Nodes == NULL)
    {
        fprintf(stderr, "Error in function GetMaxElem: Parameter is NULL\n");
        return -1;
    }
    int maxX = pt_Nodes[0].x;
    int maxY = pt_Nodes[0].y;
    int maxZ = pt_Nodes[0].z;
    for(int i = 1; i < NodeNum; i++)
    {
        maxX = (pt_Nodes[i].x > maxX) ? pt_Nodes[i].x : maxX;
        maxY = (pt_Nodes[i].y > maxY) ? pt_Nodes[i].y : maxY;
        maxZ = (pt_Nodes[i].z > maxZ) ? pt_Nodes[i].z : maxZ;
    }
    MaxX = (maxX > MaxX) ? maxX : MaxX;
    MaxY = (maxY > MaxY) ? maxY : MaxY;
    MaxZ = (maxZ > MaxZ) ? maxZ : MaxZ;
    MaxX = MaxX + 1;
    MaxY = MaxY + 1;
    MaxZ = MaxZ + 1;
    return 0;
}

int Destruction()
{
    if(pt_Keys == NULL || pt_Grids == NULL)
    {
        fprintf(stderr, "Error in function Destruction: Parameter is NULL\n");
        return -1;
    }
    if(pt_Keys->NLAG != NULL)
    {
        free(pt_Keys->NLAG);
        pt_Keys->NLAG = NULL;
    }
    if(pt_Keys->LAG != NULL)
    {
        free(pt_Keys->LAG);
        pt_Keys->LAG = NULL;
    }
    if(pt_Keys->LTOL != NULL)
    {
        free(pt_Keys->LTOL);
        pt_Keys->LTOL = NULL;
    }
    if(pt_Keys->AZM != NULL)
    {
        free(pt_Keys->AZM);
        pt_Keys->AZM = NULL;
    }
    #if 0
    if(pt_Keys->ELEVATION != NULL)
    {
        free(pt_Keys->ELEVATION);
        pt_Keys->ELEVATION = NULL;
    }
    #endif
    if(pt_Keys->ATOL != NULL)
    {
        free(pt_Keys->ATOL);
        pt_Keys->ATOL = NULL;
    }
    if(pt_Keys->BANDWH != NULL)
    {
        free(pt_Keys->BANDWH);
        pt_Keys->BANDWH = NULL;
    }
    if(pt_Keys->DIP != NULL)
    {
        free(pt_Keys->DIP);
        pt_Keys->DIP = NULL;
    }
    if(pt_Keys->DTOL != NULL)
    {
        free(pt_Keys->DTOL);
        pt_Keys->DTOL = NULL;
    }
    if(pt_Keys->BANDWD != NULL)
    {
        free(pt_Keys->BANDWD);
        pt_Keys->BANDWD = NULL;
    }


    if(pt_Grids->pt_Nodes != NULL)
    {
        free(pt_Grids->pt_Nodes);
        pt_Grids->pt_Nodes = NULL;
    }
    free(pt_Grids);
    pt_Grids = NULL;
    return 0;
}



int InitialStruct2(const char* FileKeysName, const char* FileGridName)
{
    pt_Keys = (Keys*)malloc(sizeof(Keys));
    pt_Grids = (Grid*)malloc(sizeof(Grid));

    FILE* FileKeys = fopen(FileKeysName, "r");
    FILE* FileGrid = fopen(FileGridName, "r");





    if (GetGrid(FileGrid) == -1)
    {
        fprintf(stderr, "Faild to get Grid\n");
        return -1;
    }

    if (GetKeys(FileKeys) == -1)
    {
        fprintf(stderr, "Faild to get keys\n");
        return -1;
    }

    if (FileKeys == NULL || FileGrid == NULL)
        return -1;
    fclose(FileKeys);
    fclose(FileGrid);

    return 0;
}

int InitialStruct3(const char* FileKeysName, double** grid)
{
	pt_Keys = (Keys*)malloc(sizeof(Keys));
    pt_Grids = (Grid*)malloc(sizeof(Grid));

	if(GetGrid2(grid) == -1)
	{
		return -1;
	}


	FILE* FileKeys = fopen(FileKeysName, "r");
	if (GetKeys(FileKeys) == -1)
    {
        return -1;
    }

	fclose(FileKeys);

	return 0;
}


int GetGrid2(double** grid)
{
	hosc_debug_log("GetGrid2.begin", "grid=%p", grid);



	if(pt_Grids == NULL)
	{
		hosc_debug_log("GetGrid2.fail", "pt_Grids is NULL");
		return -1;
	}
	if(grid == NULL)
	{

		hosc_debug_log("GetGrid2.fail", "grid is NULL");
		return -1;
	}


	pt_Grids->ND = (int)grid[0][0];

	pt_Grids->NVAR = (int)grid[0][1];
	TotalNodeNum = pt_Grids->ND;
	hosc_debug_log("GetGrid2.header", "ND=%d NVAR=%d", pt_Grids->ND, pt_Grids->NVAR);
	pt_Grids->pt_Nodes = (Nodes*)calloc(pt_Grids->ND, sizeof(Nodes));

	if(pt_Grids->pt_Nodes == NULL)
	{
		hosc_debug_log("GetGrid2.fail", "pt_Nodes allocation failed");
		return -1;
	}
	double tempMean = 0;
	for(int i = 1; i <= pt_Grids->ND; i++)
	{
		if(i == 1 || i == pt_Grids->ND)
		{
			hosc_debug_log("GetGrid2.row", "i=%d raw=(%.6f,%.6f,%.6f,%.6f)", i, grid[i][0], grid[i][1], grid[i][2], grid[i][3]);
		}
		pt_Grids->pt_Nodes[i-1].x = (int)grid[i][0];
		pt_Grids->pt_Nodes[i-1].y = (int)grid[i][1];
		pt_Grids->pt_Nodes[i-1].z = (int)grid[i][2];
		pt_Grids->pt_Nodes[i-1].val = grid[i][3];
		tempMean += pt_Grids->pt_Nodes[i-1].val;
		pt_Grids->pt_Nodes[i-1].Flag_Used = false;
	}
	tempMean /= pt_Grids->ND;
	hosc_debug_log("GetGrid2.mean", "mean=%.17g coordinate_basis=zero_based", tempMean);
	if(tempMean != 0)
	{
		for(int i = 0; i < pt_Grids->ND; i++)
		{
			pt_Grids->pt_Nodes[i].val -= tempMean;

		}
	}
	#if 0
	else
	{
		for(int i = 0; i <= pt_Grids->ND; i++)
		{

			delete[] grid[i];
		}
	}
	#endif
	for(int i = 0; i <= pt_Grids->ND; i++)
	{
		delete[] grid[i];
	}
	delete[] grid;
	hosc_debug_log("GetGrid2.end", "done");
	return 0;
}

int TestPTGrid()
{
	FILE* testGrid = fopen("testGrid.txt", "w");
	for(int i = 0; i < pt_Grids->ND; i++)
	{
		fprintf(testGrid, "%d, %d, %d, %lf\n", pt_Grids->pt_Nodes[i].x, pt_Grids->pt_Nodes[i].y, pt_Grids->pt_Nodes[i].z, pt_Grids->pt_Nodes[i].val);
	}
	fclose(testGrid);
	return 0;
}

void mylog(const char* msg, int num) {
    std::ofstream f("debug_log.txt", std::ios::app);
    f << msg << num<<std::endl;
    f.flush();
}


int InitialStruct4(Keys* Key, double** grid)
{
	hosc_debug_log("InitialStruct4.begin", "Key=%p grid=%p", Key, grid);
	if(Key == NULL || grid == NULL)
	{
		hosc_debug_log("InitialStruct4.fail", "Key or grid is NULL");
		return -1;
	}
	if(Key == NULL)
	{
		return -1;
	}
	pt_Keys = Key;
    pt_Grids = (Grid*)malloc(sizeof(Grid));
	if(pt_Grids == NULL)
	{
		hosc_debug_log("InitialStruct4.fail", "pt_Grids allocation failed");
		return -1;
	}

	if(GetGrid2(grid) == -1)
	{
		hosc_debug_log("InitialStruct4.fail", "GetGrid2 failed");
		return -1;
	}
	hosc_debug_log(
		"InitialStruct4.after_grid",
		"dim=%d nd=%d nx=%d ny=%d nz=%d",
		pt_Keys->DIM,
		pt_Grids != NULL ? pt_Grids->ND : -1,
		pt_Keys->NX,
		pt_Keys->NY,
		pt_Keys->NZ);


	if(pt_Keys->DIM == 2)
    {
		MaxX = pt_Keys->NX;
		MaxY = pt_Keys->NY;
		hosc_debug_log("InitialStruct4.valuemap2d.begin", "MaxX=%d MaxY=%d", MaxX, MaxY);
        ValueMap_2D = (double**)calloc(MaxX, sizeof(double*));
        if(ValueMap_2D == NULL)
		{
			hosc_debug_log("InitialStruct4.fail", "ValueMap_2D allocation failed");
            return -1;
		}
        for(int i = 0; i < MaxX; i++)
        {
            ValueMap_2D[i] = (double*)malloc(MaxY * sizeof(double));
            if(ValueMap_2D[i] == NULL)
			{
				hosc_debug_log("InitialStruct4.fail", "ValueMap_2D row allocation failed i=%d", i);
                return -1;
			}

			for (int j = 0; j < MaxY; ++j)
				ValueMap_2D[i][j] = -1.0;
        }
        for(int i = 0; i < pt_Grids->ND; i++)
        {
			if(pt_Grids->pt_Nodes[i].x < 0 || pt_Grids->pt_Nodes[i].y < 0 ||
			   pt_Grids->pt_Nodes[i].x >= MaxX || pt_Grids->pt_Nodes[i].y >= MaxY)
			{
				GsTLcerr << "pt_Grids->pt_Nodes[i].x" << pt_Grids->pt_Nodes[i].x << gstlIO::end;
				GsTLcerr << "pt_Keys->NX" << pt_Keys->NX << gstlIO::end;
				GsTLcerr << "pt_Grids->pt_Nodes[i].y" << pt_Grids->pt_Nodes[i].y << gstlIO::end;
				GsTLcerr << "pt_Keys->NY" << pt_Keys->NY << gstlIO::end;
				GsTLcerr << "valuemap outof range" << gstlIO::end;
				hosc_debug_log(
					"InitialStruct4.fail",
					"2D valuemap out of range i=%d node=(%d,%d) Max=(%d,%d)",
					i,
					pt_Grids->pt_Nodes[i].x,
					pt_Grids->pt_Nodes[i].y,
					MaxX,
					MaxY);
				return -1;
			}
            ValueMap_2D[pt_Grids->pt_Nodes[i].x][pt_Grids->pt_Nodes[i].y] = pt_Grids->pt_Nodes[i].val;

        }
		hosc_debug_log("InitialStruct4.valuemap2d.end", "done");
    }
    else
    {



        if(pt_Grids->pt_Nodes == NULL)
        {
            fprintf(stderr, "Error in function GetKeys: ");
            return -1;
        }
		MaxX = pt_Keys->NX;
		MaxY = pt_Keys->NY;
		MaxZ = pt_Keys->NZ;
		hosc_debug_log("InitialStruct4.valuemap3d.begin", "MaxX=%d MaxY=%d MaxZ=%d", MaxX, MaxY, MaxZ);
        ValueMap_3D = (double***)calloc(MaxX, sizeof(double**));
        if(ValueMap_3D == NULL)
		{
			hosc_debug_log("InitialStruct4.fail", "ValueMap_3D allocation failed");
            return -1;
		}
        for(int i = 0; i < MaxX; i++)
        {
            ValueMap_3D[i] = (double**)calloc(MaxY, sizeof(double*));
            if(ValueMap_3D[i] == NULL)
                return -1;
            for(int j = 0; j < MaxY; j++)
            {
                ValueMap_3D[i][j] = (double*)malloc(MaxZ * sizeof(double));
                if(ValueMap_3D[i][j] == NULL)
                    return -1;
                memset(ValueMap_3D[i][j], -1, MaxZ * sizeof(double));
            }
        }
		for(int i = 0; i < pt_Grids->ND; i++)
		{
			if(pt_Grids->pt_Nodes[i].x < 0 || pt_Grids->pt_Nodes[i].y < 0 || pt_Grids->pt_Nodes[i].z < 0 ||
			   pt_Grids->pt_Nodes[i].x >= MaxX || pt_Grids->pt_Nodes[i].y >= MaxY || pt_Grids->pt_Nodes[i].z >= MaxZ)
			{
				hosc_debug_log("InitialStruct4.fail", "3D valuemap out of range i=%d node=(%d,%d,%d) Max=(%d,%d,%d)", i, pt_Grids->pt_Nodes[i].x, pt_Grids->pt_Nodes[i].y, pt_Grids->pt_Nodes[i].z, MaxX, MaxY, MaxZ);
				return -1;
			}
			ValueMap_3D[pt_Grids->pt_Nodes[i].x][pt_Grids->pt_Nodes[i].y][pt_Grids->pt_Nodes[i].z] = pt_Grids->pt_Nodes[i].val;
        }
	}



	hosc_debug_log("InitialStruct4.end", "success");
	return 0;
}

