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

	int gridSize = *getGridSize;


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






    double* CumulantResult = (double*)calloc(gridSize, sizeof(double));
    clock_t start, end;
    double cpu_time_used;
#if 0

    if(InitRoot2(pt_Keys->ORDER) == -1)
    {
        fprintf(stderr, "Error in function RecCumulant: Cant init Root node of AVLTree\n");
        return NULL;
    }
#endif


    start = clock();
    bool NotfoundTemplate = true;


    if(pt_Keys->DIM == 2)
    {


        for(int i = gridSize-1; i >= 0; i--)
        {

            Template* Templates;
            int TemplateNum = 0;

            if(TemplateConstruct_2D(LAGS[i+addCrossection], pt_Keys->ORDER, &Templates, &TemplateNum, 0, pt_Grids->ND) != 0)
            {


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

				if(InitRoot2(pt_Keys->ORDER) == -1)
				{
					fprintf(stderr, "Error in function RecCumulant: Cant init Root node of AVLTree\n");

					return NULL;
				}




            CumulantResult[i] = BetterRecCumulant3(Templates->index, pt_Keys->ORDER);





			RefreshTree(pt_Keys->ORDER);

			ResetMomentWorkspaceForDLL();

            FreeTemplate(&Templates);
            Templates = NULL;


            if(CumulantResult[i] == myEmpty)
            {
                fprintf(stderr, "Error in main: Failed to calculate cumulant\n");

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



	if(NotfoundTemplate == true)
	{

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


    double** targetGrid = SaveResult2Grid(CumulantResult, LAGS, LagNum, pt_Keys->NDIR, 0, gridSize-1);


    free(CumulantResult);
    Destruction();
    Free2DVec_int(&LAGS, LagNum);






    return targetGrid;
}

double** Moment4Dll(Keys* Key, double**& grid, int* getGridSize, int* crossection, int crossectionCount)
{



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



    double* MomentResult = (double*)calloc(gridSize, sizeof(double));
    clock_t start, end;
    double cpu_time_used;




    start = clock();
    bool NotfoundTemplate = true;

    if(pt_Keys->DIM == 2)
    {
        for(int i = gridSize-1; i >= 0; i--)
        {
            Template* Templates;
            int TemplateNum = 0;
            if(TemplateConstruct_2D(LAGS[i+addCrossection], pt_Keys->ORDER, &Templates, &TemplateNum, 0, pt_Grids->ND) != 0)
            {

                continue;
            }

            NotfoundTemplate = false;

			CumulantIndex* tempCumulant = (CumulantIndex*)calloc(1, sizeof(CumulantIndex));
			tempCumulant->index = Templates->index;
			tempCumulant->Order = pt_Keys->ORDER;
			int tempNumInFstFloor = InitialCoeff(tempCumulant);

             MomentResult[i] = G_Moment[0];






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






			ResetMomentWorkspaceForDLL();
            FreeTemplate(&Templates);
            Templates = NULL;

        }
    }

	if(NotfoundTemplate == true)
	{

	}


    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Program takes %f s\n", cpu_time_used);

    double** targetGrid = SaveResult2Grid(MomentResult, LAGS, LagNum, pt_Keys->NDIR, 0, gridSize-1);

    free(MomentResult);
    Destruction();
    Free2DVec_int(&LAGS, LagNum);





    return targetGrid;
}

double** RecCalculation4DLL3(Keys* Key, double**& grid, int* getGridSize, int* crossection, int crossectionCount)
{


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



	if(NotfoundTemplate == true)
	{

	}
    double** targetGrid = SaveResult2Grid(CumulantResult, LAGS, LagNum, pt_Keys->NDIR, 0, gridSize-1);

    free(CumulantResult);
    Destruction();
    Free2DVec_int(&LAGS, LagNum);





    return targetGrid;
}
