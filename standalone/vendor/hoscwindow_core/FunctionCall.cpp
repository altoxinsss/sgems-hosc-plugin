#include "FunctionCall.h"
#include "time.h"

static void ResetMomentWorkspaceForDLL()
{
    if(G_Moment != NULL)
    {
        free(G_Moment);
        G_Moment = NULL;
    }
    if(G_MomentIndex != NULL)
    {
        if(G_SeriesNum > 0)
        {
            Free2DVec_char(&G_MomentIndex, G_SeriesNum);
        }
        else
        {
            free(G_MomentIndex);
            G_MomentIndex = NULL;
        }
    }
    if(G_MomentHash != NULL)
    {
        free(G_MomentHash);
        G_MomentHash = NULL;
    }
    G_MomentHashSize = 0;
    G_MomentHashReady = false;
    G_MomentFull = false;
    G_SeriesNum = 0;
    G_CurrMomentPos = 0;
}




double** RecCalculation4DLL2(Keys* Key, double**& grid, int* getGridSize, int* crossection, int crossectionCount)
{

	////GsTLcerr << "RecCalculation4DLL2 start" << gstlIO::end;

	//FILE* testFunc = fopen("testFunc.txt", "w");
	//fprintf(testFunc, "start1\n");
	
	if(grid == NULL)
	{
		//fprintf(testFunc, "1");
		//GsTLcerr << "1" << gstlIO::end;
		return NULL;
	}
    if (InitialStruct4(Key, grid) == -1)
    {
		//fprintf(testFunc, "2");
		//GsTLcerr << "2" << gstlIO::end;
        return NULL;
    }
	////GsTLcerr << "3" << gstlIO::end;

	//fprintf(testFunc, "start2\n");

    TreeFounded = 0;
    int** LAGS;
    int LagNum;
    InitialLagMat(&LAGS, &LagNum);
    
	int gridSize = *getGridSize;
	////GsTLcerr << "gridSize"<< gridSize << gstlIO::end;

	int addCrossection = 0;
	if(crossectionCount > 0)
	{
		addCrossection = gridSize;
		for(int i = crossectionCount - 1; i >= 0; i--)
		{
			int tempAdd = gridSize * crossection[i];
			for(int j = i+1; j < crossectionCount; j++)
			{
				tempAdd *= Key->NLAG[j];
			}
			addCrossection += tempAdd;
		}
		addCrossection -= gridSize;
	}
	////GsTLcerr << addCrossection << gstlIO::end;

	//if(*getGridSize > gridSize)
		//*getGridSize = gridSize;


    double* CumulantResult = (double*)calloc(gridSize, sizeof(double));
    clock_t start, end;
    double cpu_time_used;
#if 0
    //初始化AVL根节点
    if(InitRoot2(pt_Keys->ORDER) == -1)
    {
        fprintf(stderr, "Error in function RecCumulant: Cant init Root node of AVLTree\n");
        return NULL;
    }
#endif


    start = clock();
    bool NotfoundTemplate = true;
	//fprintf(testFunc, "lets go\n");
	////GsTLcerr <<"pt_Keys->DIM:" <<pt_Keys->DIM << gstlIO::end;
    if(pt_Keys->DIM == 2)
    {
        //判断是否需要刷新AVL树
        //int RefreshSign1 = LAGS[gridSize-1][0];
        for(int i = gridSize-1; i >= 0; i--)
        {
			////mylog("start", i);
            Template* Templates;
            int TemplateNum = 0;
			////GsTLcerr << "TemplateConstruct_2D" << gstlIO::end;
            if(TemplateConstruct_2D(LAGS[i+addCrossection], pt_Keys->ORDER, &Templates, &TemplateNum, 0, pt_Grids->ND) != 0)
            {
                //fprintf(testFunc, "here %lf\n", CumulantResult[i]);
				////GsTLcerr << "continue " << i << gstlIO::end;
                continue;
            }
			////mylog("TemplateConstruct_2D", i);
            //SaveTemplateResult(Templates, pt_Keys->ORDER, TemplateNum);
            NotfoundTemplate = false;
			
			
				if(RootAVL2 != NULL)
				{
					for(int j = 1; j <= pt_Keys->ORDER; j++)
					{
						if(RootAVL2[j] != NULL)
						{
							FreeAVLTree(&RootAVL2[j]);
							RootAVL2[j] = NULL;
						}
					}
					free(RootAVL2);
					RootAVL2 = NULL;
				}
				//初始化AVL根节点
				if(InitRoot2(pt_Keys->ORDER) == -1)
				{
					fprintf(stderr, "Error in function RecCumulant: Cant init Root node of AVLTree\n");
					//GsTLcerr << "4" << gstlIO::end;
					return NULL;
				}

				////mylog("if(InitRoot2(pt_Keys->ORDER) == -1)", i);
			
			
            CumulantResult[i] = BetterRecCumulant3(Templates->index, pt_Keys->ORDER);
			////mylog("cumulant", i);
			////GsTLcerr << "CumulantResult[" <<i<<"]: "<< CumulantResult[i]<< gstlIO::end;
			


			RefreshTree(pt_Keys->ORDER);

			ResetMomentWorkspaceForDLL();
			////GsTLcerr << "FreeTemplate(&Templates);" << gstlIO::end;
            FreeTemplate(&Templates);
            Templates = NULL;

            
            if(CumulantResult[i] == myEmpty)
            {
                fprintf(stderr, "Error in main: Failed to calculate cumulant\n");
				//GsTLcerr << "5" << gstlIO::end;
                return NULL;
            }
			////GsTLcerr << "next while" << gstlIO::end;
        }
    }

    else if(pt_Keys->DIM == 3)
    {
        for(int i = gridSize-1; i >= 0; i--)
        {
            Template* Templates;
            int TemplateNum = 0;
            int SearchResult = TemplateConstruct_3D(LAGS[i+addCrossection], pt_Keys->ORDER, &Templates, &TemplateNum, 0, pt_Grids->ND);
            if(SearchResult != 0)
            {
                printf("TemplateConstruct_3D return: %d\n", SearchResult);
                continue;
            }
            NotfoundTemplate = false;
            
			
			
				if(RootAVL2 != NULL)
				{
					for(int j = 1; j <= pt_Keys->ORDER; j++)
					{
						if(RootAVL2[j] != NULL)
						{
							FreeAVLTree(&RootAVL2[j]);
							RootAVL2[j] = NULL;
						}
					}
					free(RootAVL2);
					RootAVL2 = NULL;
				}
				//初始化AVL根节点
				if(InitRoot2(pt_Keys->ORDER) == -1)
				{
					fprintf(stderr, "Error in function RecCumulant: Cant init Root node of AVLTree\n");
					return NULL;
				}
			
            CumulantResult[i] = BetterRecCumulant3(Templates->index, pt_Keys->ORDER);
        FoundCumulantin3D:
        #if 0
            for(int j = 0; j < pt_Keys->ORDER-3; j++)
            {
                free(G_Moments[j]);
                G_Moments[j] = NULL;
            }
            free(G_Moments);
            G_Moments = NULL;
        #endif
            //Free2DVec(&G_Moments, pt_Keys->ORDER+1);
            //G_Moments = NULL;
            //Free2DVec_int(&G_IndexSeries, G_SeriesNum);
            //G_IndexSeries = NULL;

			ResetMomentWorkspaceForDLL();
            FreeTemplate(&Templates);
            Templates = NULL;
            printf("%d: %f\n", i, CumulantResult[i]);
            
            if(CumulantResult[i] == myEmpty)
            {
                fprintf(stderr, "Error in main: Failed to calculate cumulant\n");
                return NULL;
            }
        }
    }

	////GsTLcerr << "if(NotfoundTemplate == true)" << gstlIO::end;

	if(NotfoundTemplate == true)
	{
		//GsTLcerr << "didnt find template"  << gstlIO::end;
	}

    end = clock();
	for(int i = 0; i <=  pt_Keys->ORDER; i++)
    {
        if(RootAVL2[i] != NULL)
        {
            free(RootAVL2[i]);
            RootAVL2[i] = NULL;
        }
    }
    free(RootAVL2);
    RootAVL2 = NULL;
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Program takes %f s\n", cpu_time_used);
    //printf("Tree Founded %d\n", TreeFounded);
	////GsTLcerr<<"Key->NDIR inner: "<<Key->NDIR<<gstlIO::end;
    double** targetGrid = SaveResult2Grid(CumulantResult, LAGS, LagNum, pt_Keys->NDIR, 0, gridSize-1);
	////GsTLcerr<<"Key->NDIR inner out: "<<Key->NDIR<<gstlIO::end;

    free(CumulantResult);
    Destruction();
    Free2DVec_int(&LAGS, LagNum);
    //if(NotfoundTemplate == false)
        //printf("Found some template\n");
	//fprintf(testFunc, "start6\n");
	//fclose(testFunc);

	//GsTLcerr << "RecCalculation4DLL2 end" << gstlIO::end;
    return targetGrid;
}

double** Moment4Dll(Keys* Key, double**& grid, int* getGridSize, int* crossection, int crossectionCount)
{
	//FILE* testFunc = fopen("testFunc.txt", "w");
	//fprintf(testFunc, "start1\n");
	
	if(grid == NULL)
	{
		//fprintf(testFunc, "1");
		return NULL;
	}
    if (InitialStruct4(Key, grid) == -1)
    {
		//fprintf(testFunc, "2");
        return NULL;
    }

	//fprintf(testFunc, "start2\n");

    TreeFounded = 0;
    int** LAGS;
    int LagNum;
    InitialLagMat(&LAGS, &LagNum);
    
	//int gridSize = *getGridSize;
	//gridSize = (gridSize > LagNum) ? LagNum : gridSize;
	//*getGridSize = gridSize;

	int gridSize = *getGridSize;

	int addCrossection = 0;
	if(crossectionCount > 0)
	{
		addCrossection = 1;
		for(int i = crossectionCount - 1; i >= 0; i--)
		{
			addCrossection = addCrossection * (crossection[i]-1) * gridSize;
		}
	}

    //printf("%d\n", LagNum);
	//fprintf(testFunc, "%d\n", gridSize);
    double* MomentResult = (double*)calloc(gridSize, sizeof(double));
    clock_t start, end;
    double cpu_time_used;




    start = clock();
    bool NotfoundTemplate = true;
	//fprintf(testFunc, "lets go\n");
    if(pt_Keys->DIM == 2)
    {
        for(int i = gridSize-1; i >= 0; i--)
        {
            Template* Templates;
            int TemplateNum = 0;
            if(TemplateConstruct_2D(LAGS[i+addCrossection], pt_Keys->ORDER, &Templates, &TemplateNum, 0, pt_Grids->ND) != 0)
            {
                //fprintf(testFunc, "here %lf\n", CumulantResult[i]);
                continue;
            }
            //SaveTemplateResult(Templates, pt_Keys->ORDER, TemplateNum);
            NotfoundTemplate = false;

			CumulantIndex* tempCumulant = (CumulantIndex*)calloc(1, sizeof(CumulantIndex));
			tempCumulant->index = Templates->index;
			tempCumulant->Order = pt_Keys->ORDER;
			int tempNumInFstFloor = InitialCoeff(tempCumulant);

             MomentResult[i] = G_Moment[0];

            //Free2DVec(&G_Moments, pt_Keys->ORDER+1);
            //G_Moments = NULL;
            //Free2DVec_int(&G_IndexSeries, G_SeriesNum);
            //G_IndexSeries = NULL;

			ResetMomentWorkspaceForDLL();
            FreeTemplate(&Templates);
            Templates = NULL;
        }
    }

    else if(pt_Keys->DIM == 3)
    {
        for(int i = gridSize-1; i >= 0; i--)
        {
            Template* Templates;
            int TemplateNum = 0;
            int SearchResult = TemplateConstruct_3D(LAGS[i+addCrossection], pt_Keys->ORDER, &Templates, &TemplateNum, 0, pt_Grids->ND);
            if(SearchResult != 0)
            {
                printf("TemplateConstruct_3D return: %d\n", SearchResult);
                continue;
            }
            NotfoundTemplate = false;

            MomentResult[i] = G_Moment[0];
        FoundCumulantin3D:
        
            //Free2DVec(&G_Moments, pt_Keys->ORDER+1);
            //G_Moments = NULL;
            //Free2DVec_int(&G_IndexSeries, G_SeriesNum);
            //G_IndexSeries = NULL;
			
			ResetMomentWorkspaceForDLL();
            FreeTemplate(&Templates);
            Templates = NULL;

        }
    }

	if(NotfoundTemplate == true)
	{
		//GsTLcerr << "zero template" << gstlIO::end;
	}


    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Program takes %f s\n", cpu_time_used);
    //printf("Tree Founded %d\n", TreeFounded);
    double** targetGrid = SaveResult2Grid(MomentResult, LAGS, LagNum, pt_Keys->NDIR, 0, gridSize-1);

    free(MomentResult);
    Destruction();
    Free2DVec_int(&LAGS, LagNum);
    //if(NotfoundTemplate == false)
        //printf("Found some template\n");
	//fprintf(testFunc, "start6\n");
	//fclose(testFunc);

    return targetGrid;
}

double** RecCalculation4DLL3(Keys* Key, double**& grid, int* getGridSize, int* crossection, int crossectionCount)
{

	//GsTLcerr << "start RecCalculation4DLL3" << gstlIO::end;
	if(grid == NULL)
	{
		return NULL;
	}
    if (InitialStruct4(Key, grid) == -1)
    {
        return NULL;
    }

    TreeFounded = 0;
    int** LAGS;
    int LagNum;
    InitialLagMat(&LAGS, &LagNum);
    
	//int gridSize = LagNum;
	//if(*getGridSize > gridSize)
		//*getGridSize = gridSize;


	int gridSize = *getGridSize;

	int addCrossection = 0;
	if(crossectionCount > 0)
	{
		addCrossection = 1;
		for(int i = crossectionCount - 1; i >= 0; i--)
		{
			addCrossection = addCrossection * (crossection[i]-1) * gridSize;
		}
	}

    double* CumulantResult = (double*)calloc(gridSize, sizeof(double));



    bool NotfoundTemplate = true;
    if(pt_Keys->DIM == 2)
    {
        for(int i = gridSize-1; i >= 0; i--)
        {
            Template* Templates;
            int TemplateNum = 0;
            if(TemplateConstruct_2D_old(LAGS[i+addCrossection], pt_Keys->ORDER, &Templates, &TemplateNum, 0, pt_Grids->ND) != 0)
            {
                continue;
            }
            NotfoundTemplate = false;
			
			double result = TrueCumulant(Templates->index, pt_Keys->ORDER);
			memset(g_all_partitions, 0, sizeof(Partition) * MAX_PARTITIONS);
			g_total_partitions = 0;

            CumulantResult[i] = result;

			//GsTLcerr << CumulantResult[i] << gstlIO::end;

            Free2DVec(&G_Moments, pt_Keys->ORDER+1);
            G_Moments = NULL;
            Free2DVec_int(&G_IndexSeries, G_SeriesNum);
            G_IndexSeries = NULL;
            FreeTemplate(&Templates);
            Templates = NULL;

            if(CumulantResult[i] == myEmpty)
            {
                fprintf(stderr, "Error in main: Failed to calculate cumulant\n");
				//fprintf(testFunc, "4");
                return NULL;
            }
        }
    }

    else if(pt_Keys->DIM == 3)
    {
        for(int i = gridSize-1; i >= 0; i--)
        {
            Template* Templates;
            int TemplateNum = 0;
            int SearchResult = TemplateConstruct_3D_old(LAGS[i+addCrossection], pt_Keys->ORDER, &Templates, &TemplateNum, 0, pt_Grids->ND);
            if(SearchResult != 0)
            {
                printf("TemplateConstruct_3D return: %d\n", SearchResult);
                continue;
            }
            NotfoundTemplate = false;
            
			double result = TrueCumulant(Templates->index, pt_Keys->ORDER);
			memset(g_all_partitions, 0, sizeof(Partition) * MAX_PARTITIONS);
			g_total_partitions = 0;
            CumulantResult[i] = result;

        FoundCumulantin3D:
            Free2DVec(&G_Moments, pt_Keys->ORDER+1);
            G_Moments = NULL;
            Free2DVec_int(&G_IndexSeries, G_SeriesNum);
            G_IndexSeries = NULL;
            FreeTemplate(&Templates);
            Templates = NULL;
            
            if(CumulantResult[i] == myEmpty)
            {
                fprintf(stderr, "Error in main: Failed to calculate cumulant\n");
                return NULL;
            }
        }
    }


    //printf("Tree Founded %d\n", TreeFounded);
	if(NotfoundTemplate == true)
	{
		//GsTLcerr << "zero template" << gstlIO::end;
	}
    double** targetGrid = SaveResult2Grid(CumulantResult, LAGS, LagNum, pt_Keys->NDIR, 0, gridSize-1);

    free(CumulantResult);
    Destruction();
    Free2DVec_int(&LAGS, LagNum);
    //if(NotfoundTemplate == false)
        //printf("Found some template\n");
	//fprintf(testFunc, "start6\n");
	//fclose(testFunc);

    return targetGrid;
}
