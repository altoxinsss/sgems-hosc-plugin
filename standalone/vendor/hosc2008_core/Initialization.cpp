#include "Initialization.h"


#define MaxStringLen 300
#define MaxNumLen 30


//for dll
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
    //读取KEYS文件时，假设KEYS是按照固定格式存储数据
    //下列均未对fgets进行出错检测
    //atoi未进行溢出检测
    char* pc_TempStr = (char*)calloc(MaxStringLen, sizeof(char));
    char* pc_TempParam = (char*)calloc(MaxNumLen, sizeof(char));
    if(pc_TempParam == NULL)
        return -1;
  
    //Get integer paragram
    int* p_TempKeys = (int* )pt_Keys;
    //涉及指针操作时需要考虑指针是否越界
    //10代表该结构体内有10个integer型成员
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
    
    //Get number of lags in each direction
    //无法采取指针偏移为结构体中的指针成员分配内存，因为calloc会返回新分配内存的地址
    pt_Keys->NLAG = (int* )calloc(pt_Keys->NDIR, sizeof(int));
    pt_Keys->LAG = (int* )calloc(pt_Keys->NDIR, sizeof(int));
    pt_Keys->LTOL = (int* )calloc(pt_Keys->NDIR, sizeof(int));
    pt_Keys->AZM = (int* )calloc(pt_Keys->NDIR, sizeof(int));
    //pt_Keys->ELEVATION = (int*) calloc(pt_Keys->NDIR, sizeof(int));
    //assert(pt_Keys->ELEVATION != NULL);
    pt_Keys->ATOL = (int* )calloc(pt_Keys->NDIR, sizeof(int));
    pt_Keys->BANDWH = (int* )calloc(pt_Keys->NDIR, sizeof(int));
    pt_Keys->DIP = (int* )calloc(pt_Keys->NDIR, sizeof(int));
    pt_Keys->DTOL = (int* )calloc(pt_Keys->NDIR, sizeof(int));
    pt_Keys->BANDWD = (int* )calloc(pt_Keys->NDIR, sizeof(int));
    bool TES = false;

	//FILE* newtest111 = fopen("tnmnmt1111.txt", "w");
	//fprintf(newtest111, "aaaaaaaaaaaaaaaaa");
	//fclose(newtest111);

    //涉及指针操作时需要考虑指针是否越界！！！
    //for(int k = 0; k< ((int*)pt_Keys+sizeof(Keys)-p_TempKeys)/sizeof(int*); k++) //由NDIR决定动态分配内存的指针有10个
    for(int k = 0; k< 3; k++)
    {
        //if(k < 4)
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
                    //if(i == pt_Keys->NDIR-1  && pc_TempParam[j] != '\0')
                        //fgets(pc_TempStr, MaxStringLen, FileKeys);
                    //pc_TempParam[j] = '\0';
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
        //for(int k = 0; k < 7; k++)
        for(int k = 0; k < 6; k++)
        {
            memset(pc_TempParam, 0, MaxNumLen*sizeof(char));
            for(int j = 0; j < MaxNumLen; j++)
            {
                pc_TempParam[j] = fgetc(FileKeys);
                //if(k == 6)
                    //pc_TempParam[j]=fgetc(FileKeys);
                if(pc_TempParam[j] == ' ' || pc_TempParam[j] == '\n')
                {
                    pc_TempParam[j] = '\0';
                        //fgets(pc_TempStr, MaxStringLen, FileKeys);
                    //pc_TempParam[j] = '\0';
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
    //Get ND and NVAR
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

    //读取点数据
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
    //free(pt_Keys);
    //pt_Keys = NULL;
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



    //printf("%p\n", FileKeys);

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
	//FILE* testGrid = fopen("testGrid2.txt", "w");
	//fprintf(testGrid, "%p\n", grid);
	
	if(pt_Grids == NULL)
	{
		return -1;
	}
	if(grid == NULL)
	{
		//fprintf(testGrid, "grid is null");
		return -1;
	}
	//fprintf(testGrid, "grid is not null");
	//fclose(testGrid);
	pt_Grids->ND = (int)grid[0][0];
	//GsTLcerr << pt_Grids->ND << gstlIO::end;
	pt_Grids->NVAR = (int)grid[0][1];
	TotalNodeNum = pt_Grids->ND;
	pt_Grids->pt_Nodes = (Nodes*)calloc(pt_Grids->ND, sizeof(Nodes));

	if(pt_Grids->pt_Nodes == NULL)
	{
		return -1;
	}
	double tempMean = 0;
	int TRANS = 0;
	for(int i = 1; i <= pt_Grids->ND; i++)
	{
		pt_Grids->pt_Nodes[i-1].x = (int)grid[i][0];
		pt_Grids->pt_Nodes[i-1].y = (int)grid[i][1];
		pt_Grids->pt_Nodes[i-1].z = (int)grid[i][2];
		pt_Grids->pt_Nodes[i-1].val = grid[i][3];
		tempMean += pt_Grids->pt_Nodes[i-1].val;
		pt_Grids->pt_Nodes[i-1].Flag_Used = false;
		if(i == 1 && pt_Grids->pt_Nodes[i-1].x == 1)
		{
			pt_Grids->pt_Nodes[i-1].x--;
			pt_Grids->pt_Nodes[i-1].y--;
			pt_Grids->pt_Nodes[i-1].z--;
			TRANS = 1;
		}
		if(i > 1 && TRANS == 1)
		{
			pt_Grids->pt_Nodes[i-1].x--;
			pt_Grids->pt_Nodes[i-1].y--;
			pt_Grids->pt_Nodes[i-1].z--;
		}
	}
	tempMean /= pt_Grids->ND;
	if(tempMean != 0)
	{
		for(int i = 0; i < pt_Grids->ND; i++)
		{
			pt_Grids->pt_Nodes[i].val -= tempMean;
			//会造成野指针
		}
	}
	#if 0
	else
	{
		for(int i = 0; i <= pt_Grids->ND; i++)
		{
			//会造成野指针
			delete[] grid[i];
		}
	}
	#endif
	for(int i = 0; i <= pt_Grids->ND; i++)
	{
		delete[] grid[i];
	}
	delete[] grid;
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
    f.flush();  // 立即写入！
}


int InitialStruct4(Keys* Key, double** grid)
{
	if(Key == NULL || grid == NULL)
		return -1;
	if(Key == NULL)
	{
		return -1;
	}
	pt_Keys = Key;
    pt_Grids = (Grid*)malloc(sizeof(Grid));

	if(GetGrid2(grid) == -1)
	{
		return -1;
	}

	//valuemap的大小需要设置成全局的
	if(pt_Keys->DIM == 2)
    {
		GetMaxElem(pt_Grids->pt_Nodes, pt_Grids->ND);
        ValueMap_2D = (double**)calloc(MaxX, sizeof(double*));
        if(ValueMap_2D == NULL)
            return -1;
        for(int i = 0; i < MaxX; i++)
        {
            ValueMap_2D[i] = (double*)malloc(MaxY * sizeof(double));
            if(ValueMap_2D[i] == NULL)
                return -1;
            //memset(ValueMap_2D[i], -1, pt_Keys->NY * sizeof(double));
			for (int j = 0; j < MaxY; ++j)
				ValueMap_2D[i][j] = -1.0;
        }
        for(int i = 0; i < pt_Grids->ND; i++)
        {
			if(pt_Grids->pt_Nodes[i].x >= MaxX || pt_Grids->pt_Nodes[i].y >= MaxY)
			{
				GsTLcerr << "pt_Grids->pt_Nodes[i].x" << pt_Grids->pt_Nodes[i].x << gstlIO::end;
				GsTLcerr << "pt_Keys->NX" << pt_Keys->NX << gstlIO::end;
				GsTLcerr << "pt_Grids->pt_Nodes[i].y" << pt_Grids->pt_Nodes[i].y << gstlIO::end;
				GsTLcerr << "pt_Keys->NY" << pt_Keys->NY << gstlIO::end;
				GsTLcerr << "valuemap outof range" << gstlIO::end;
				return -1;
			}
            ValueMap_2D[pt_Grids->pt_Nodes[i].x][pt_Grids->pt_Nodes[i].y] = pt_Grids->pt_Nodes[i].val;
			
        }
    }
    else
    {
        //MaxX = pt_Keys->NX;
        //MaxY = pt_Keys->NY;
        //MaxZ = pt_Keys->NZ;
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
                memset(ValueMap_3D[i][j], -1, MaxZ * sizeof(double));
            }
        }
        for(int i = 0; i < pt_Grids->ND; i++)
        {
            ValueMap_3D[pt_Grids->pt_Nodes[i].x][pt_Grids->pt_Nodes[i].y][pt_Grids->pt_Nodes[i].z] = pt_Grids->pt_Nodes[i].val;
        }
    }

	
	
	return 0;
}

