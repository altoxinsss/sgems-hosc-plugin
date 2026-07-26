#include "FunctionCall.h"
#include "DebugLogger.h"
#include "time.h"
#include <vector>

namespace wmc = weighted_moment_core;

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

static void ClearAVLWorkspaceForDLL(int order)
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

static int InitAVLWorkspaceForDLL(int order)
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

static void ClearOneEvaluationWorkspaceForDLL(int order)
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

static double CalculateWeightedTemplateStatistic(
    Template* templates,
    int templateCount,
    int variableCount,
    const std::vector<int>& multiplicities,
    const wmc::WeightedMomentOptions& options,
    bool momentTarget)
{
    if(templates == NULL || templateCount <= 0 || variableCount <= 0)
    {
        return 0.0;
    }

    std::vector< std::vector<double> > sampleValues(static_cast<size_t>(templateCount));
    std::vector<wmc::WeightedSample> samples(static_cast<size_t>(templateCount));
    for(int i = 0; i < templateCount; ++i)
    {
        sampleValues[static_cast<size_t>(i)].resize(static_cast<size_t>(variableCount), 0.0);
        for(int j = 0; j < variableCount; ++j)
        {
            sampleValues[static_cast<size_t>(i)][static_cast<size_t>(j)] = templates[i].Value[j + 1];
        }
        samples[static_cast<size_t>(i)].values = &sampleValues[static_cast<size_t>(i)][0];
        samples[static_cast<size_t>(i)].value_count = variableCount;
        samples[static_cast<size_t>(i)].x = static_cast<double>(i);
        samples[static_cast<size_t>(i)].y = 0.0;
        samples[static_cast<size_t>(i)].z = 0.0;
    }

    if(momentTarget)
    {
        return wmc::weighted_moment(
            &samples[0], templateCount, &multiplicities[0], variableCount, options);
    }
    return wmc::weighted_cumulant_from_multiplicities(
        &samples[0], templateCount, &multiplicities[0], variableCount, options);
}

static int BuildWeightedMomentKeyForDLL(const int* index, int order, char* out, int outSize)
{
    if(index == NULL || out == NULL || outSize < order + 1)
    {
        return -1;
    }
    bool skipped = false;
    int zeroNum = 0;
    for(int i = order; i >= 1; --i)
    {
        if(index[i] == 0 && skipped == false)
        {
            zeroNum++;
            continue;
        }
        if(index[i] != 0 && skipped == false)
        {
            skipped = true;
        }
        out[i - 1] = static_cast<char>('0' + index[i]);
    }
    out[order - zeroNum] = '\0';
    return 0;
}

static double CalculateWeightedMomentForIndex(
    Template* templates,
    int templateCount,
    int variableCount,
    const int* index,
    const wmc::WeightedMomentOptions& options)
{
    if(templates == NULL || templateCount <= 0 || variableCount <= 0 || index == NULL)
    {
        return 0.0;
    }

    std::vector<int> multiplicities(static_cast<size_t>(variableCount), 0);
    for(int j = 0; j < variableCount; ++j)
    {
        multiplicities[static_cast<size_t>(j)] = index[j + 1];
    }

    std::vector< std::vector<double> > sampleValues(static_cast<size_t>(templateCount));
    std::vector<wmc::WeightedSample> samples(static_cast<size_t>(templateCount));
    for(int i = 0; i < templateCount; ++i)
    {
        sampleValues[static_cast<size_t>(i)].resize(static_cast<size_t>(variableCount), 0.0);
        for(int j = 0; j < variableCount; ++j)
        {
            sampleValues[static_cast<size_t>(i)][static_cast<size_t>(j)] = templates[i].Value[j + 1];
        }
        samples[static_cast<size_t>(i)].values = &sampleValues[static_cast<size_t>(i)][0];
        samples[static_cast<size_t>(i)].value_count = variableCount;
        samples[static_cast<size_t>(i)].x = static_cast<double>(i);
        samples[static_cast<size_t>(i)].y = 0.0;
        samples[static_cast<size_t>(i)].z = 0.0;
    }

    return wmc::weighted_moment(
        &samples[0], templateCount, &multiplicities[0], variableCount, options);
}

static int PrecomputeWeightedMomentsForAVL(
    Template* templates,
    int templateCount,
    const int* index,
    int order,
    const wmc::WeightedMomentOptions& options)
{
    if(templates == NULL || templateCount <= 0 || index == NULL || order <= 0)
    {
        return -1;
    }
    if(G_MomentFull == true)
    {
        return 0;
    }

    int* limits = (int*)calloc(order + 1, sizeof(int));
    if(limits == NULL)
    {
        return -1;
    }

    int momentCount = 1;
    for(int j = 1; j <= order; ++j)
    {
        limits[j] = index[j];
        momentCount *= (limits[j] + 1);
    }

    const int totalNum = momentCount;
    G_Moment = (double*)calloc(totalNum, sizeof(double));
    G_MomentIndex = (char**)calloc(totalNum, sizeof(char*));
    if(G_Moment == NULL || G_MomentIndex == NULL)
    {
        free(limits);
        return -1;
    }
    G_SeriesNum = totalNum;
    G_CurrMomentPos = 0;

    G_Moment[0] = CalculateWeightedMomentForIndex(templates, templateCount, order, index, options);
    G_MomentIndex[0] = (char*)calloc(order + 1, sizeof(char));
    if(G_MomentIndex[0] == NULL || BuildWeightedMomentKeyForDLL(index, order, G_MomentIndex[0], order + 1) != 0)
    {
        free(limits);
        return -1;
    }
    G_CurrMomentPos = 1;

    if(momentCount > 1)
    {
        int radix = 1;
        int* base = (int*)calloc(order + 1, sizeof(int));
        if(base == NULL)
        {
            free(limits);
            return -1;
        }
        for(int j = 1; j <= order; ++j)
        {
            base[j] = limits[j] + 1;
            radix *= base[j];
        }

        for(int code = 1; code < radix; ++code)
        {
            int* localExp = (int*)calloc(order + 1, sizeof(int));
            if(localExp == NULL)
            {
                free(base);
                free(limits);
                return -1;
            }
            localExp[0] = 0;
            int rem = code;
            for(int j = order; j >= 1; --j)
            {
                localExp[j] = rem % base[j];
                rem /= base[j];
            }

            G_Moment[code] = CalculateWeightedMomentForIndex(templates, templateCount, order, localExp, options);
            G_MomentIndex[code] = (char*)calloc(order + 1, sizeof(char));
            if(G_MomentIndex[code] == NULL ||
               BuildWeightedMomentKeyForDLL(localExp, order, G_MomentIndex[code], order + 1) != 0)
            {
                free(localExp);
                free(base);
                free(limits);
                return -1;
            }
            free(localExp);
        }
        G_CurrMomentPos = radix;
        free(base);
    }

    free(limits);
    G_MomentFull = true;
    G_MomentHashReady = false;
    if(G_MomentHash != NULL)
    {
        free(G_MomentHash);
        G_MomentHash = NULL;
    }
    G_MomentHashSize = 0;
    return 0;
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



			if(i >= gridSize)
			{

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

double** RecCalculation4DLL2(Keys* Key, double**& grid, int* getGridSize, int* crossection, int crossectionCount)
{
    hosc_debug_log("RecCalculation4DLL2.begin", "Key=%p grid=%p getGridSize=%p cross_count=%d", Key, grid, getGridSize, crossectionCount);
    if(grid == NULL || getGridSize == NULL)
    {
        hosc_debug_log("RecCalculation4DLL2.fail", "grid or getGridSize is NULL");
        return NULL;
    }
    if(InitialStruct4(Key, grid) == -1)
    {
        grid = NULL;
        hosc_debug_log("RecCalculation4DLL2.fail", "InitialStruct4 failed");
        return NULL;
    }
    grid = NULL;
    hosc_debug_log(
        "RecCalculation4DLL2.initialized",
        "gridSizeIn=%d dim=%d ndir=%d order=%d dims=%d,%d,%d",
        *getGridSize,
        pt_Keys != NULL ? pt_Keys->DIM : -1,
        pt_Keys != NULL ? pt_Keys->NDIR : -1,
        pt_Keys != NULL ? pt_Keys->ORDER : -1,
        pt_Keys != NULL ? pt_Keys->NX : -1,
        pt_Keys != NULL ? pt_Keys->NY : -1,
        pt_Keys != NULL ? pt_Keys->NZ : -1);

    TreeFounded = 0;
    int** LAGS = NULL;
    int LagNum = 0;
    if(InitialLagMat(&LAGS, &LagNum) != 0)
    {
        hosc_debug_log("RecCalculation4DLL2.fail", "InitialLagMat failed");
        Destruction();
        return NULL;
    }
    hosc_debug_log("RecCalculation4DLL2.lags", "LagNum=%d", LagNum);

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
        hosc_debug_log("RecCalculation4DLL2.fail", "CollectCrossSectionLagPositionsForDLL failed gridSize=%d", gridSize);
        Destruction();
        Free2DVec_int(&LAGS, LagNum);
        return NULL;
    }
    hosc_debug_log("RecCalculation4DLL2.cross", "selected=%d expected=%d", (int)selectedLagPositions.size(), gridSize);

    double* CumulantResult = (double*)calloc(selectedLagPositions.size(), sizeof(double));
    if(CumulantResult == NULL)
    {
        perror("Error in function RecCalculation4DLL2: calloc return NULL\n");
        hosc_debug_log("RecCalculation4DLL2.fail", "CumulantResult allocation failed selected=%d", (int)selectedLagPositions.size());
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
        if(i == static_cast<int>(selectedLagPositions.size()) - 1 || i == 0)
        {
            hosc_debug_log("RecCalculation4DLL2.loop", "i=%d lagPosition=%d", i, lagPosition);
        }

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
            hosc_debug_log("RecCalculation4DLL2.fail", "unsupported dimension %d", pt_Keys->DIM);
            free(CumulantResult);
            Destruction();
            Free2DVec_int(&LAGS, LagNum);
            return NULL;
        }

        if(SearchResult != 0)
        {
            if(i == static_cast<int>(selectedLagPositions.size()) - 1 || i == 0)
            {
                hosc_debug_log("RecCalculation4DLL2.loop.skip", "i=%d lagPosition=%d SearchResult=%d", i, lagPosition, SearchResult);
            }
            continue;
        }
        NotfoundTemplate = false;
        if(i == static_cast<int>(selectedLagPositions.size()) - 1 || i == 0)
        {
            hosc_debug_log("RecCalculation4DLL2.loop.template", "i=%d TemplateNum=%d Templates=%p", i, TemplateNum, Templates);
        }

        if(EvaluateOneCumulantForDLL(Templates, pt_Keys->ORDER, pt_Keys->DIM, &CumulantResult[i]) != 0)
        {
            hosc_debug_log("RecCalculation4DLL2.fail", "EvaluateOneCumulantForDLL failed i=%d lagPosition=%d", i, lagPosition);
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
        hosc_debug_log("RecCalculation4DLL2.save", "targetGrid=%p outSize=%d", targetGrid, *getGridSize);
    }
    else
    {
        hosc_debug_log("RecCalculation4DLL2.warning", "no template found");
    }

    free(CumulantResult);
    ClearOneEvaluationWorkspaceForDLL(pt_Keys->ORDER);
    Destruction();
    Free2DVec_int(&LAGS, LagNum);
    hosc_debug_log("RecCalculation4DLL2.end", "targetGrid=%p", targetGrid);
    return targetGrid;
}

double** Moment4Dll(Keys* Key, double**& grid, int* getGridSize, int* crossection, int crossectionCount)
{
    hosc_debug_log("Moment4Dll.begin", "Key=%p grid=%p getGridSize=%p cross_count=%d", Key, grid, getGridSize, crossectionCount);
    if(Key == NULL || grid == NULL || getGridSize == NULL)
    {
        hosc_debug_log("Moment4Dll.fail", "null input");
        return NULL;
    }
    if(InitialStruct4(Key, grid) == -1)
    {
        grid = NULL;
        hosc_debug_log("Moment4Dll.fail", "InitialStruct4 failed");
        return NULL;
    }
    grid = NULL;

    TreeFounded = 0;
    int** LAGS = NULL;
    int LagNum = 0;
    if(InitialLagMat(&LAGS, &LagNum) != 0)
    {
        hosc_debug_log("Moment4Dll.fail", "InitialLagMat failed");
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
        hosc_debug_log("Moment4Dll.fail", "CollectCrossSectionLagPositionsForDLL failed gridSize=%d", gridSize);
        Destruction();
        Free2DVec_int(&LAGS, LagNum);
        return NULL;
    }
    hosc_debug_log("Moment4Dll.cross", "selected=%d expected=%d", (int)selectedLagPositions.size(), gridSize);

    double* MomentResult = (double*)calloc(selectedLagPositions.size(), sizeof(double));
    if(MomentResult == NULL)
    {
        hosc_debug_log("Moment4Dll.fail", "MomentResult allocation failed selected=%d", (int)selectedLagPositions.size());
        Destruction();
        Free2DVec_int(&LAGS, LagNum);
        return NULL;
    }

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
            hosc_debug_log("Moment4Dll.fail", "unsupported dimension %d", pt_Keys->DIM);
            free(MomentResult);
            Destruction();
            Free2DVec_int(&LAGS, LagNum);
            return NULL;
        }

        if(SearchResult != 0)
        {
            if(i == static_cast<int>(selectedLagPositions.size()) - 1 || i == 0)
            {
                hosc_debug_log("Moment4Dll.loop.skip", "i=%d lagPosition=%d SearchResult=%d", i, lagPosition, SearchResult);
            }
            continue;
        }
        NotfoundTemplate = false;

        CumulantIndex* tempCumulant = (CumulantIndex*)calloc(1, sizeof(CumulantIndex));
        if(tempCumulant == NULL)
        {
            hosc_debug_log("Moment4Dll.fail", "tempCumulant allocation failed");
            FreeTemplate(&Templates);
            free(MomentResult);
            Destruction();
            Free2DVec_int(&LAGS, LagNum);
            return NULL;
        }
        tempCumulant->index = Templates->index;
        tempCumulant->Order = pt_Keys->ORDER;
        const int tempNumInFstFloor = InitialCoeff(tempCumulant);
        (void)tempNumInFstFloor;

        MomentResult[i] = G_Moment[0];
        ResetMomentWorkspaceForDLL();
        free(tempCumulant);
        FreeTemplate(&Templates);
        Templates = NULL;
    }

    const clock_t end = clock();
    const double cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Program takes %f s\n", cpu_time_used);

    double** targetGrid = NULL;
    if(NotfoundTemplate == false)
    {
        *getGridSize = static_cast<int>(selectedLagPositions.size());
        targetGrid = SaveSelectedResult2GridForDLL(MomentResult, LAGS, selectedLagPositions, pt_Keys->NDIR);
        hosc_debug_log("Moment4Dll.save", "targetGrid=%p outSize=%d", targetGrid, *getGridSize);
    }
    else
    {
        hosc_debug_log("Moment4Dll.warning", "no template found");
    }

    free(MomentResult);
    ResetMomentWorkspaceForDLL();
    Destruction();
    Free2DVec_int(&LAGS, LagNum);
    hosc_debug_log("Moment4Dll.end", "targetGrid=%p", targetGrid);
    return targetGrid;
}

double** RecCalculation4DLL3(Keys* Key, double**& grid, int* getGridSize, int* crossection, int crossectionCount)
{
    hosc_debug_log("RecCalculation4DLL3.begin", "Key=%p grid=%p getGridSize=%p cross_count=%d", Key, grid, getGridSize, crossectionCount);
    if(grid == NULL || getGridSize == NULL)
    {
        hosc_debug_log("RecCalculation4DLL3.fail", "grid or getGridSize is NULL");
        return NULL;
    }
    if(InitialStruct4(Key, grid) == -1)
    {
        grid = NULL;
        hosc_debug_log("RecCalculation4DLL3.fail", "InitialStruct4 failed");
        return NULL;
    }
    grid = NULL;
    hosc_debug_log(
        "RecCalculation4DLL3.initialized",
        "gridSizeIn=%d dim=%d ndir=%d order=%d dims=%d,%d,%d",
        *getGridSize,
        pt_Keys != NULL ? pt_Keys->DIM : -1,
        pt_Keys != NULL ? pt_Keys->NDIR : -1,
        pt_Keys != NULL ? pt_Keys->ORDER : -1,
        pt_Keys != NULL ? pt_Keys->NX : -1,
        pt_Keys != NULL ? pt_Keys->NY : -1,
        pt_Keys != NULL ? pt_Keys->NZ : -1);

    TreeFounded = 0;
    int** LAGS = NULL;
    int LagNum = 0;
    if(InitialLagMat(&LAGS, &LagNum) != 0)
    {
        hosc_debug_log("RecCalculation4DLL3.fail", "InitialLagMat failed");
        Destruction();
        return NULL;
    }
    hosc_debug_log("RecCalculation4DLL3.lags", "LagNum=%d", LagNum);

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
        hosc_debug_log("RecCalculation4DLL3.fail", "CollectCrossSectionLagPositionsForDLL failed gridSize=%d", gridSize);
        Destruction();
        Free2DVec_int(&LAGS, LagNum);
        return NULL;
    }
    hosc_debug_log("RecCalculation4DLL3.cross", "selected=%d expected=%d", (int)selectedLagPositions.size(), gridSize);

    double* CumulantResult = (double*)calloc(selectedLagPositions.size(), sizeof(double));
    if(CumulantResult == NULL)
    {
        hosc_debug_log("RecCalculation4DLL3.fail", "CumulantResult allocation failed selected=%d", (int)selectedLagPositions.size());
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
        if(i == static_cast<int>(selectedLagPositions.size()) - 1 || i == 0)
        {
            hosc_debug_log("RecCalculation4DLL3.loop", "i=%d lagPosition=%d", i, lagPosition);
        }

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
            hosc_debug_log("RecCalculation4DLL3.fail", "unsupported dimension %d", pt_Keys->DIM);
            free(CumulantResult);
            Destruction();
            Free2DVec_int(&LAGS, LagNum);
            return NULL;
        }

        if(SearchResult != 0)
        {
            if(i == static_cast<int>(selectedLagPositions.size()) - 1 || i == 0)
            {
                hosc_debug_log("RecCalculation4DLL3.loop.skip", "i=%d lagPosition=%d SearchResult=%d", i, lagPosition, SearchResult);
            }
            continue;
        }
        NotfoundTemplate = false;
        if(i == static_cast<int>(selectedLagPositions.size()) - 1 || i == 0)
        {
            hosc_debug_log("RecCalculation4DLL3.loop.template", "i=%d TemplateNum=%d Templates=%p", i, TemplateNum, Templates);
        }

        CumulantResult[i] = BellCumulantFromMultiplicityIndex(Templates->index, pt_Keys->ORDER);
        memset(g_all_partitions, 0, sizeof(Partition) * MAX_PARTITIONS);
        g_total_partitions = 0;
        ResetMomentWorkspaceForDLL();
        FreeTemplate(&Templates);
        Templates = NULL;

        if(CumulantResult[i] == myEmpty)
        {
            hosc_debug_log("RecCalculation4DLL3.fail", "Bell result is myEmpty i=%d lagPosition=%d", i, lagPosition);
            free(CumulantResult);
            Destruction();
            Free2DVec_int(&LAGS, LagNum);
            return NULL;
        }
    }

    double** targetGrid = NULL;
    if(NotfoundTemplate == false)
    {
        *getGridSize = static_cast<int>(selectedLagPositions.size());
        targetGrid = SaveSelectedResult2GridForDLL(CumulantResult, LAGS, selectedLagPositions, pt_Keys->NDIR);
        hosc_debug_log("RecCalculation4DLL3.save", "targetGrid=%p outSize=%d", targetGrid, *getGridSize);
    }
    else
    {
        hosc_debug_log("RecCalculation4DLL3.warning", "no template found");
    }

    free(CumulantResult);
    ResetMomentWorkspaceForDLL();
    Destruction();
    Free2DVec_int(&LAGS, LagNum);
    hosc_debug_log("RecCalculation4DLL3.end", "targetGrid=%p", targetGrid);
    return targetGrid;
}

double** WeightedMoment4Dll(
    Keys* Key,
    double**& grid,
    int* getGridSize,
    int* crossection,
    int crossectionCount,
    const wmc::WeightedMomentOptions& options,
    bool momentTarget)
{
    if(Key == NULL || grid == NULL || getGridSize == NULL)
    {
        return NULL;
    }
    if(InitialStruct4(Key, grid) == -1)
    {
        grid = NULL;
        return NULL;
    }
    grid = NULL;

    int** LAGS = NULL;
    int LagNum = 0;
    if(InitialLagMat(&LAGS, &LagNum) == -1)
    {
        Destruction();
        return NULL;
    }

    std::vector<int> selectedLagPositions;
    if(CollectCrossSectionLagPositionsForDLL(
           pt_Keys,
           LAGS,
           LagNum,
           crossection,
           crossectionCount,
           *getGridSize,
           &selectedLagPositions) != 0)
    {
        Free2DVec_int(&LAGS, LagNum);
        Destruction();
        return NULL;
    }

    std::vector<int> multiplicities(static_cast<size_t>(pt_Keys->ORDER), 0);
    multiplicities[0] = 1;
    for(int i = 1; i < pt_Keys->ORDER; ++i)
    {
        multiplicities[static_cast<size_t>(i)] = pt_Keys->NOD[i - 1];
    }

    double* resultValues = (double*)calloc(selectedLagPositions.size(), sizeof(double));
    if(resultValues == NULL)
    {
        Free2DVec_int(&LAGS, LagNum);
        Destruction();
        return NULL;
    }

    for(int row = static_cast<int>(selectedLagPositions.size()) - 1; row >= 0; --row)
    {
        Template* Templates = NULL;
        int TemplateNum = 0;
        const int lagPosition = selectedLagPositions[static_cast<size_t>(row)];
        int searchResult = 0;
        if(pt_Keys->DIM == 2)
        {
            searchResult = TemplateConstruct_2D(LAGS[lagPosition], pt_Keys->ORDER, &Templates, &TemplateNum, 0, pt_Grids->ND);
        }
        else if(pt_Keys->DIM == 3)
        {
            searchResult = TemplateConstruct_3D(LAGS[lagPosition], pt_Keys->ORDER, &Templates, &TemplateNum, 0, pt_Grids->ND);
        }
        else
        {
            free(resultValues);
            Free2DVec_int(&LAGS, LagNum);
            Destruction();
            return NULL;
        }

        if(searchResult != 0)
        {
            continue;
        }

        resultValues[row] = CalculateWeightedTemplateStatistic(
            Templates,
            TemplateNum,
            pt_Keys->ORDER,
            multiplicities,
            options,
            momentTarget);
        FreeTemplate(&Templates);
    }

    *getGridSize = static_cast<int>(selectedLagPositions.size());
    double** targetGrid = SaveSelectedResult2GridForDLL(resultValues, LAGS, selectedLagPositions, pt_Keys->NDIR);

    free(resultValues);
    Free2DVec_int(&LAGS, LagNum);
    Destruction();
    return targetGrid;
}

double** WeightedRecCalculation4DLL2(
    Keys* Key,
    double**& grid,
    int* getGridSize,
    int* crossection,
    int crossectionCount,
    const wmc::WeightedMomentOptions& options,
    bool momentTarget)
{
    hosc_debug_log("WeightedRecCalculation4DLL2.begin", "Key=%p grid=%p getGridSize=%p cross_count=%d momentTarget=%d", Key, grid, getGridSize, crossectionCount, momentTarget ? 1 : 0);
    if(Key == NULL || grid == NULL || getGridSize == NULL)
    {
        hosc_debug_log("WeightedRecCalculation4DLL2.fail", "null input");
        return NULL;
    }
    if(momentTarget)
    {
        hosc_debug_log("WeightedRecCalculation4DLL2.redirect", "moment target uses explicit weighted moment kernel");
        return WeightedMoment4Dll(Key, grid, getGridSize, crossection, crossectionCount, options, true);
    }
    if(InitialStruct4(Key, grid) == -1)
    {
        grid = NULL;
        hosc_debug_log("WeightedRecCalculation4DLL2.fail", "InitialStruct4 failed");
        return NULL;
    }
    grid = NULL;

    int** LAGS = NULL;
    int LagNum = 0;
    if(InitialLagMat(&LAGS, &LagNum) != 0)
    {
        hosc_debug_log("WeightedRecCalculation4DLL2.fail", "InitialLagMat failed");
        Destruction();
        return NULL;
    }

    std::vector<int> selectedLagPositions;
    if(CollectCrossSectionLagPositionsForDLL(
           pt_Keys,
           LAGS,
           LagNum,
           crossection,
           crossectionCount,
           *getGridSize,
           &selectedLagPositions) != 0)
    {
        hosc_debug_log("WeightedRecCalculation4DLL2.fail", "CollectCrossSectionLagPositionsForDLL failed gridSize=%d", *getGridSize);
        Free2DVec_int(&LAGS, LagNum);
        Destruction();
        return NULL;
    }
    hosc_debug_log("WeightedRecCalculation4DLL2.cross", "selected=%d expected=%d", (int)selectedLagPositions.size(), *getGridSize);

    double* resultValues = (double*)calloc(selectedLagPositions.size(), sizeof(double));
    if(resultValues == NULL)
    {
        hosc_debug_log("WeightedRecCalculation4DLL2.fail", "result allocation failed");
        Free2DVec_int(&LAGS, LagNum);
        Destruction();
        return NULL;
    }

    for(int row = static_cast<int>(selectedLagPositions.size()) - 1; row >= 0; --row)
    {
        Template* Templates = NULL;
        int TemplateNum = 0;
        const int lagPosition = selectedLagPositions[static_cast<size_t>(row)];
        int searchResult = 0;
        if(row == static_cast<int>(selectedLagPositions.size()) - 1 || row == 0)
        {
            hosc_debug_log("WeightedRecCalculation4DLL2.loop", "row=%d lagPosition=%d", row, lagPosition);
        }

        if(pt_Keys->DIM == 2)
        {
            searchResult = TemplateConstruct_2D(LAGS[lagPosition], pt_Keys->ORDER, &Templates, &TemplateNum, 0, pt_Grids->ND);
        }
        else if(pt_Keys->DIM == 3)
        {
            searchResult = TemplateConstruct_3D(LAGS[lagPosition], pt_Keys->ORDER, &Templates, &TemplateNum, 0, pt_Grids->ND);
        }
        else
        {
            hosc_debug_log("WeightedRecCalculation4DLL2.fail", "unsupported dim=%d", pt_Keys->DIM);
            free(resultValues);
            Free2DVec_int(&LAGS, LagNum);
            Destruction();
            return NULL;
        }

        if(searchResult != 0)
        {
            continue;
        }

        if(row == static_cast<int>(selectedLagPositions.size()) - 1 || row == 0)
        {
            hosc_debug_log("WeightedRecCalculation4DLL2.template", "row=%d TemplateNum=%d", row, TemplateNum);
        }

        ResetMomentWorkspaceForDLL();
        ClearAVLWorkspaceForDLL(pt_Keys->ORDER);
        if(PrecomputeWeightedMomentsForAVL(Templates, TemplateNum, Templates->index, pt_Keys->ORDER, options) != 0)
        {
            hosc_debug_log("WeightedRecCalculation4DLL2.fail", "PrecomputeWeightedMomentsForAVL failed row=%d", row);
            FreeTemplate(&Templates);
            free(resultValues);
            Free2DVec_int(&LAGS, LagNum);
            Destruction();
            return NULL;
        }

        if(InitAVLWorkspaceForDLL(pt_Keys->ORDER) != 0)
        {
            hosc_debug_log("WeightedRecCalculation4DLL2.fail", "InitAVLWorkspaceForDLL failed row=%d", row);
            FreeTemplate(&Templates);
            free(resultValues);
            Free2DVec_int(&LAGS, LagNum);
            Destruction();
            return NULL;
        }

        SetLowOrderCumulantMomentCacheMode(1);
        resultValues[row] = BetterRecCumulant3(Templates->index, pt_Keys->ORDER);
        SetLowOrderCumulantMomentCacheMode(0);
        if(row == 0)
        {
            std::vector<int> debugMultiplicities(static_cast<size_t>(pt_Keys->ORDER), 0);
            debugMultiplicities[0] = 1;
            for(int debugIndex = 1; debugIndex < pt_Keys->ORDER; ++debugIndex)
            {
                debugMultiplicities[static_cast<size_t>(debugIndex)] = pt_Keys->NOD[debugIndex - 1];
            }
            double debugBellValue = CalculateWeightedTemplateStatistic(
                Templates,
                TemplateNum,
                pt_Keys->ORDER,
                debugMultiplicities,
                options,
                false);
            hosc_debug_log(
                "WeightedRecCalculation4DLL2.compare",
                "row=0 avl=%.17g bell_ref=%.17g diff=%.17g templateNum=%d momentCount=%d",
                resultValues[row],
                debugBellValue,
                resultValues[row] - debugBellValue,
                TemplateNum,
                G_SeriesNum);
        }
        FreeTemplate(&Templates);
        ClearOneEvaluationWorkspaceForDLL(pt_Keys->ORDER);

        if(resultValues[row] == myEmpty)
        {
            hosc_debug_log("WeightedRecCalculation4DLL2.fail", "BetterRecCumulant3 returned myEmpty row=%d", row);
            free(resultValues);
            Free2DVec_int(&LAGS, LagNum);
            Destruction();
            return NULL;
        }
    }

    *getGridSize = static_cast<int>(selectedLagPositions.size());
    double** targetGrid = SaveSelectedResult2GridForDLL(resultValues, LAGS, selectedLagPositions, pt_Keys->NDIR);
    hosc_debug_log("WeightedRecCalculation4DLL2.save", "targetGrid=%p outSize=%d", targetGrid, *getGridSize);

    free(resultValues);
    Free2DVec_int(&LAGS, LagNum);
    Destruction();
    return targetGrid;
}
