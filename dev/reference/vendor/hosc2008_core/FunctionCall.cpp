#include "FunctionCall.h"
#include "time.h"
#include <vector>

void ResetMomentWorkspaceForDLL()
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

void ClearAVLWorkspaceForDLL(int order)
{
    if(RootAVL2 == NULL)
    {
        return;
    }
    for(int i = 1; i <= order; i++)
    {
        if(RootAVL2[i] != NULL)
        {
            FreeAVLTree(&RootAVL2[i]);
            RootAVL2[i] = NULL;
        }
    }
    free(RootAVL2);
    RootAVL2 = NULL;
}

int InitAVLWorkspaceForDLL(int order)
{
    ClearAVLWorkspaceForDLL(order);
    if(InitRoot2(order) == -1)
    {
        fprintf(stderr, "Error in function RecCalculation4DLL2: failed to initialize AVL roots\n");
        return -1;
    }
    return 0;
}

static int EvaluateOneCumulantForDLL(Template* templates, int order, int dim, double* outValue)
{
    if(templates == NULL || outValue == NULL)
    {
        return -1;
    }

    if(PrecomputeMomentsForIndex(templates->index, order) == myEmpty)
    {
        fprintf(stderr, "Error in function RecCalculation4DLL2: PrecomputeMomentsForIndex failed\n");
        return -1;
    }

    if(InitAVLWorkspaceForDLL(order) != 0)
    {
        return -1;
    }

    double value = BetterRecCumulant3(templates->index, order);
    if(value == myEmpty)
    {
        return -1;
    }

    (void)dim;
    *outValue = value;
    return 0;
}

void ClearOneEvaluationWorkspaceForDLL(int order)
{
    ResetMomentWorkspaceForDLL();
    ClearAVLWorkspaceForDLL(order);
}

static bool LagMatchesCrossSectionForDLL(const Keys* key, const int* lag, const int* crossection, int crossectionCount)
{
    if(key == NULL || lag == NULL)
    {
        return false;
    }
    if(crossection == NULL || crossectionCount <= 0)
    {
        return true;
    }

    for(int i = 0; i < crossectionCount; i++)
    {
        if(crossection[i] < 0 || crossection[i] >= key->NLAG[i])
        {
            return false;
        }
        if(lag[i] != crossection[i] * key->LAG[i])
        {
            return false;
        }
    }
    return true;
}

static int CollectCrossSectionLagPositionsForDLL(
    const Keys* key,
    int** LAGS,
    int LagNum,
    const int* crossection,
    int crossectionCount,
    int expectedGridSize,
    std::vector<int>* selectedLagPositions)
{
    if(key == NULL || LAGS == NULL || selectedLagPositions == NULL)
    {
        return -1;
    }

    selectedLagPositions->clear();
    for(int i = 0; i < LagNum; i++)
    {
        if(LagMatchesCrossSectionForDLL(key, LAGS[i], crossection, crossectionCount) == true)
        {
            selectedLagPositions->push_back(i);
        }
    }

    if(expectedGridSize > 0 && static_cast<int>(selectedLagPositions->size()) != expectedGridSize)
    {
        fprintf(
            stderr,
            "Error in function RecCalculation4DLL2: cross-section selects %d lags, expected %d\n",
            static_cast<int>(selectedLagPositions->size()),
            expectedGridSize);
        return -1;
    }

    return 0;
}

static double** SaveSelectedResult2GridForDLL(
    const double* CumulantResult,
    int** LAGS,
    const std::vector<int>& selectedLagPositions,
    int DirectionNum)
{
    if(CumulantResult == NULL || LAGS == NULL)
    {
        return NULL;
    }

    const int rowsNum = static_cast<int>(selectedLagPositions.size());
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
            for(int j = 0; j < i; j++)
            {
                free(targetGrid[j]);
            }
            free(targetGrid);
            return NULL;
        }

        const int lagPosition = selectedLagPositions[static_cast<size_t>(i)];
        for(int j = 0; j < DirectionNum; j++)
        {
            targetGrid[i][j] = LAGS[lagPosition][j];
        }
        targetGrid[i][DirectionNum] = CumulantResult[i];
    }

    return targetGrid;
}




static double** RecCalculation4DLL2_Legacy(Keys* Key, double**& grid, int* getGridSize, int* crossection, int crossectionCount)
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
			//mylog("start", i);
            Template* Templates;
            int TemplateNum = 0;
			////GsTLcerr << "TemplateConstruct_2D" << gstlIO::end;
			//mylog("i+addCrossection=", i+addCrossection);
            if(TemplateConstruct_2D(LAGS[i+addCrossection], pt_Keys->ORDER, &Templates, &TemplateNum, 0, pt_Grids->ND) != 0)
            {
                //fprintf(testFunc, "here %lf\n", CumulantResult[i]);
				////GsTLcerr << "continue " << i << gstlIO::end;
                continue;
            }
			//mylog("TemplateConstruct_2D", i);
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

				//mylog("if(InitRoot2(pt_Keys->ORDER) == -1)", i);
			
			if(i >= gridSize)
			{
				//GsTLcerr <<i<< ">=" <<gridSize<< gstlIO::end;
			}
			////GsTLcerr << "BetterRecCumulant3" << gstlIO::end;
            CumulantResult[i] = BetterRecCumulant3(Templates->index, pt_Keys->ORDER);
			//mylog("cumulant", i);
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

double** RecCalculation4DLL2(Keys* Key, double**& grid, int* getGridSize, int* crossection, int crossectionCount)
{
    if(grid == NULL || getGridSize == NULL)
    {
        return NULL;
    }
    if(InitialStruct4(Key, grid) == -1)
    {
        return NULL;
    }

    TreeFounded = 0;
    int** LAGS = NULL;
    int LagNum = 0;
    if(InitialLagMat(&LAGS, &LagNum) != 0)
    {
        Destruction();
        return NULL;
    }

    const int gridSize = *getGridSize;
    std::vector<int> selectedLagPositions;
    if(CollectCrossSectionLagPositionsForDLL(
           pt_Keys,
           LAGS,
           LagNum,
           crossection,
           crossectionCount,
           gridSize,
           &selectedLagPositions) != 0)
    {
        Destruction();
        Free2DVec_int(&LAGS, LagNum);
        return NULL;
    }

    double* CumulantResult = (double*)calloc(selectedLagPositions.size(), sizeof(double));
    if(CumulantResult == NULL)
    {
        perror("Error in function RecCalculation4DLL2: calloc return NULL\n");
        Destruction();
        Free2DVec_int(&LAGS, LagNum);
        return NULL;
    }

    double** targetGrid = NULL;
    const clock_t start = clock();
    bool NotfoundTemplate = true;

    for(int i = static_cast<int>(selectedLagPositions.size()) - 1; i >= 0; i--)
    {
        Template* Templates = NULL;
        int TemplateNum = 0;
        const int lagPosition = selectedLagPositions[static_cast<size_t>(i)];
        int SearchResult = 0;

        if(pt_Keys->DIM == 2)
        {
            SearchResult = TemplateConstruct_2D(LAGS[lagPosition], pt_Keys->ORDER, &Templates, &TemplateNum, 0, pt_Grids->ND);
        }
        else if(pt_Keys->DIM == 3)
        {
            SearchResult = TemplateConstruct_3D(LAGS[lagPosition], pt_Keys->ORDER, &Templates, &TemplateNum, 0, pt_Grids->ND);
        }
        else
        {
            fprintf(stderr, "Error in function RecCalculation4DLL2: unsupported dimension %d\n", pt_Keys->DIM);
            free(CumulantResult);
            Destruction();
            Free2DVec_int(&LAGS, LagNum);
            return NULL;
        }

        if(SearchResult != 0)
        {
            continue;
        }
        NotfoundTemplate = false;

        if(EvaluateOneCumulantForDLL(Templates, pt_Keys->ORDER, pt_Keys->DIM, &CumulantResult[i]) != 0)
        {
            FreeTemplate(&Templates);
            ClearOneEvaluationWorkspaceForDLL(pt_Keys->ORDER);
            free(CumulantResult);
            Destruction();
            Free2DVec_int(&LAGS, LagNum);
            return NULL;
        }

        FreeTemplate(&Templates);
        ClearOneEvaluationWorkspaceForDLL(pt_Keys->ORDER);
    }

    const clock_t end = clock();
    const double cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Program takes %f s\n", cpu_time_used);

    if(NotfoundTemplate == false)
    {
        *getGridSize = static_cast<int>(selectedLagPositions.size());
        targetGrid = SaveSelectedResult2GridForDLL(CumulantResult, LAGS, selectedLagPositions, pt_Keys->NDIR);
    }

    free(CumulantResult);
    ClearOneEvaluationWorkspaceForDLL(pt_Keys->ORDER);
    Destruction();
    Free2DVec_int(&LAGS, LagNum);
    return targetGrid;
}

double** Moment4Dll(Keys* Key, double**& grid, int* getGridSize, int* crossection, int crossectionCount)
{
	//FILE* testFunc = fopen("testFunc.txt", "w");
	//fprintf(testFunc, "start1\n");
	mylog("start Moment4Dll", -1);
	
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
    int** LAGS = NULL;
    int LagNum = 0;
    if(InitialLagMat(&LAGS, &LagNum) != 0)
    {
        Destruction();
        return NULL;
    }
    
	const int gridSize = *getGridSize;
    std::vector<int> selectedLagPositions;
    if(CollectCrossSectionLagPositionsForDLL(
           pt_Keys,
           LAGS,
           LagNum,
           crossection,
           crossectionCount,
           gridSize,
           &selectedLagPositions) != 0)
    {
        Destruction();
        Free2DVec_int(&LAGS, LagNum);
        return NULL;
    }

    double* MomentResult = (double*)calloc(selectedLagPositions.size(), sizeof(double));
    if(MomentResult == NULL)
    {
        Destruction();
        Free2DVec_int(&LAGS, LagNum);
        return NULL;
    }
    clock_t start, end;
    double cpu_time_used;




    start = clock();
    bool NotfoundTemplate = true;
	//fprintf(testFunc, "lets go\n");
    for(int i = static_cast<int>(selectedLagPositions.size()) - 1; i >= 0; i--)
    {
        Template* Templates = NULL;
        int TemplateNum = 0;
        const int lagPosition = selectedLagPositions[static_cast<size_t>(i)];
        int SearchResult = 0;
        if(pt_Keys->DIM == 2)
        {
            SearchResult = TemplateConstruct_2D(LAGS[lagPosition], pt_Keys->ORDER, &Templates, &TemplateNum, 0, pt_Grids->ND);
        }
        else if(pt_Keys->DIM == 3)
        {
            SearchResult = TemplateConstruct_3D(LAGS[lagPosition], pt_Keys->ORDER, &Templates, &TemplateNum, 0, pt_Grids->ND);
        }
        else
        {
            fprintf(stderr, "Error in function Moment4Dll: unsupported dimension %d\n", pt_Keys->DIM);
            free(MomentResult);
            Destruction();
            Free2DVec_int(&LAGS, LagNum);
            return NULL;
        }

        if(SearchResult != 0)
        {
            if(pt_Keys->DIM == 3)
            {
                printf("TemplateConstruct_3D return: %d\n", SearchResult);
            }
            continue;
        }
        NotfoundTemplate = false;

        CumulantIndex* tempCumulant = (CumulantIndex*)calloc(1, sizeof(CumulantIndex));
        tempCumulant->index = Templates->index;
        tempCumulant->Order = pt_Keys->ORDER;
        int tempNumInFstFloor = InitialCoeff(tempCumulant);
        (void)tempNumInFstFloor;

        MomentResult[i] = G_Moment[0];
        mylog("MomentResult[i] = G_Moment[0];", MomentResult[i]);

        ResetMomentWorkspaceForDLL();
        free(tempCumulant);
        FreeTemplate(&Templates);
        Templates = NULL;
    }

	if(NotfoundTemplate == true)
	{
		//GsTLcerr << "zero template" << gstlIO::end;
	}


    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Program takes %f s\n", cpu_time_used);
    //printf("Tree Founded %d\n", TreeFounded);
    *getGridSize = static_cast<int>(selectedLagPositions.size());
    double** targetGrid = SaveSelectedResult2GridForDLL(MomentResult, LAGS, selectedLagPositions, pt_Keys->NDIR);

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
    int** LAGS = NULL;
    int LagNum = 0;
    if(InitialLagMat(&LAGS, &LagNum) != 0)
    {
        Destruction();
        return NULL;
    }
    
	const int gridSize = *getGridSize;
    std::vector<int> selectedLagPositions;
    if(CollectCrossSectionLagPositionsForDLL(
           pt_Keys,
           LAGS,
           LagNum,
           crossection,
           crossectionCount,
           gridSize,
           &selectedLagPositions) != 0)
    {
        Destruction();
        Free2DVec_int(&LAGS, LagNum);
        return NULL;
    }

    double* CumulantResult = (double*)calloc(selectedLagPositions.size(), sizeof(double));
    if(CumulantResult == NULL)
    {
        Destruction();
        Free2DVec_int(&LAGS, LagNum);
        return NULL;
    }


    bool NotfoundTemplate = true;
    for(int i = static_cast<int>(selectedLagPositions.size()) - 1; i >= 0; i--)
    {
        Template* Templates = NULL;
        int TemplateNum = 0;
        const int lagPosition = selectedLagPositions[static_cast<size_t>(i)];
        int SearchResult = 0;
        if(pt_Keys->DIM == 2)
        {
            SearchResult = TemplateConstruct_2D(LAGS[lagPosition], pt_Keys->ORDER, &Templates, &TemplateNum, 0, pt_Grids->ND);
        }
        else if(pt_Keys->DIM == 3)
        {
            SearchResult = TemplateConstruct_3D(LAGS[lagPosition], pt_Keys->ORDER, &Templates, &TemplateNum, 0, pt_Grids->ND);
        }
        else
        {
            fprintf(stderr, "Error in function RecCalculation4DLL3: unsupported dimension %d\n", pt_Keys->DIM);
            free(CumulantResult);
            Destruction();
            Free2DVec_int(&LAGS, LagNum);
            return NULL;
        }

        if(SearchResult != 0)
        {
            if(pt_Keys->DIM == 3)
            {
                printf("TemplateConstruct_3D return: %d\n", SearchResult);
            }
            continue;
        }
        NotfoundTemplate = false;
        
        double result = BellCumulantFromMultiplicityIndex(Templates->index, pt_Keys->ORDER);
        memset(g_all_partitions, 0, sizeof(Partition) * MAX_PARTITIONS);
        g_total_partitions = 0;
        CumulantResult[i] = result;

        if(G_Moments != NULL)
        {
            Free2DVec(&G_Moments, pt_Keys->ORDER+1);
            G_Moments = NULL;
        }
        if(G_IndexSeries != NULL && G_SeriesNum > 0)
        {
            Free2DVec_int(&G_IndexSeries, G_SeriesNum);
        }
        G_IndexSeries = NULL;
        G_SeriesNum = 0;
        FreeTemplate(&Templates);
        Templates = NULL;
        
        if(CumulantResult[i] == myEmpty)
        {
            fprintf(stderr, "Error in main: Failed to calculate cumulant\n");
            free(CumulantResult);
            Destruction();
            Free2DVec_int(&LAGS, LagNum);
            return NULL;
        }
    }


    //printf("Tree Founded %d\n", TreeFounded);
	if(NotfoundTemplate == true)
	{
		//GsTLcerr << "zero template" << gstlIO::end;
	}
    *getGridSize = static_cast<int>(selectedLagPositions.size());
    double** targetGrid = SaveSelectedResult2GridForDLL(CumulantResult, LAGS, selectedLagPositions, pt_Keys->NDIR);

    free(CumulantResult);
    Destruction();
    Free2DVec_int(&LAGS, LagNum);
    //if(NotfoundTemplate == false)
        //printf("Found some template\n");
	//fprintf(testFunc, "start6\n");
	//fclose(testFunc);

    return targetGrid;
}
