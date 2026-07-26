

#ifndef __PLUGIN_h__
#define __PLUGIN_h__

#include <string>
#include <GsTLAppli/geostat/common.h>
#include <GsTLAppli/geostat/geostat_algo.h>

#include <GsTL/matrix_library/tnt_lib.h>
#include <GsTL/matrix_library/tnt/cmat.h>
#include <GsTL/utils/smartptr.h>
#include <GsTLAppli/grid/grid_model/geostat_grid.h>

#include <vector>

#include <GsTLAppli/geostat/filtersim_std/is_categorical.h>
#include <GsTLAppli/geostat/filtersim_std/TI_manipulation.h>


#include "FunctionCall.h"
#include "WeightedMomentCore.h"
#include <stdio.h>
#include <stdlib.h>


class RGrid;


int G_CalculateTarget;
int G_TargetDimension;

enum Algorithm {
    AVL,
    Bell
};

typedef struct
{
	int nx;
	int ny;
	int nz;
	double cx;
	double cy;
	double cz;
	double dx;
	double dy;
	double dz;
} GridShape;

class __declspec(dllexport) HOSC_WINDOW : public Geostat_algo
{
public:
    HOSC_WINDOW();
	~HOSC_WINDOW();

public:
	virtual bool initialize( const Parameters_handler* parameters, Error_messages_handler* errors );
	virtual int execute( GsTL_project* proj=0 );
	virtual std::string name() const { return "HOSC_WINDOW"; }
	static Named_interface* create_new_interface( std::string& ) { return new HOSC_WINDOW; }

public:
    bool get_Source_Data(const Parameters_handler* parameters, Error_messages_handler* error_mesgs );
    bool get_Target_Data(const Parameters_handler* parameters, Error_messages_handler* error_mesgs );

	bool get_Keys_Data(const Parameters_handler* parameters, Error_messages_handler* error_mesgs);
	bool get_Calculate_Target(const Parameters_handler* parameters, Error_messages_handler* error_mesgs);
	bool get_Calculate_Algorithm(const Parameters_handler* parameters, Error_messages_handler* error_mesgs);
	bool get_Weighted_Moment_Options(const Parameters_handler* parameters, Error_messages_handler* error_mesgs);
	bool get_My_Shape(const Parameters_handler* parameters, Error_messages_handler* error_mesgs);
	bool get_Sliding_Window(const Parameters_handler* parameters, Error_messages_handler* error_mesgs);
	bool get_Windows_Moment(const Parameters_handler* parameters, Error_messages_handler* error_mesgs);
	int calculateMomentWindow();
	void release_source_grid_cache();
	Geostat_grid* CreateSparseOutputPointSet(
		const char* gridName,
		const std::vector<Geostat_grid::location_type>& point_locations,
		const std::vector<double>& point_values);



	int set_Target_Data(double** targetGrid, int gridSize, int valuePos);

	RGrid* CreateMyGrid(GridShape* tempGridShape, const char* gridName);
	std::string build_output_grid_base_name(const GridShape* tempGridShape, bool pointset) const;
	int set_MyGrid_Data(double** targetGrid, int gridSize, int valuePos);



	void TestMyGrid();

public:

    Geostat_grid* source_grid_;
	RGrid* source_rgrid_;
    GsTLGridProperty* source_x;
	GsTLGridProperty* source_y;
	GsTLGridProperty* source_z;
	GsTLGridProperty* source_val;
	double** grid;
	int grid_alloc_rows_;
	int source_nx_;
	int source_ny_;
	int source_nz_;
	int source_min_x_;
	int source_min_y_;
	int source_min_z_;
	double source_origin_x_;
	double source_origin_y_;
	double source_origin_z_;
	double source_spacing_x_;
	double source_spacing_y_;
	double source_spacing_z_;
	std::vector<unsigned char> source_present_;


    RGrid* target_grid_;
	GsTLGridProperty* target_property_;
	string target_prop_name_;
	string target_dim_num;


	Geostat_grid* auto_grid_;
	GsTLGridProperty* auto_property_;
	string auto_prop_name_;
	string auto_grid_name_;
	bool auto_grid_registered_;
	GsTL_project* current_project_;
	GridShape myGridShape;
	std::vector<unsigned char> output_keep_mask_;



#if 0

	RGrid* getDefaultGrid();
	RGrid* default_grid_;
	GsTLGridProperty* default_property_name_;
#endif

	enum Algorithm al;
	int crossection[10];
	int crossectionNum;





	int momentWindowL;
	int momentWindowW;
	int momentWindowH;
	int momentOrder;
	int stepX;
	int stepY;
	int stepZ;
	GridShape myWindowShape;
	bool invertWindowOutput;
	bool enablePercentileClipping;
	double percentileClipMin;
	double percentileClipMax;
	bool useWeightedMoment;
	weighted_moment_core::WeightedMomentOptions weightedOptions;


	Keys* Key;
};



#endif
