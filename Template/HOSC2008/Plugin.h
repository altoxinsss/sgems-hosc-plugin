

#ifndef __PLUGIN_h__
#define __PLUGIN_h__

#include <string>
#include <GsTLAppli/geostat/common.h>
#include <GsTLAppli/geostat/geostat_algo.h>

#include <GsTL/matrix_library/tnt_lib.h>
#include <GsTL/matrix_library/tnt/cmat.h>
#include <GsTL/utils/smartptr.h>

#include <vector>

#include <GsTLAppli/geostat/filtersim_std/is_categorical.h>
#include <GsTLAppli/geostat/filtersim_std/TI_manipulation.h>


#include "FunctionCall.h"
#include "../../shared/WeightedMomentCore.h"
#include <stdio.h>
#include <stdlib.h>


class RGrid;
class Geostat_grid;



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

// This layout is frozen after direction reordering and before the calculator
// takes ownership of Key. Result writers must not read calculator-owned arrays
// because legacy and exact calculation paths release them before returning.
typedef struct
{
	bool valid;
	int displayed_direction_count;
	int row_column[3];
	int lag_step[3];
	int dimensions[3];
} TemplateOutputLayout;

class __declspec(dllexport) HOSC_TEMPLATE : public Geostat_algo
{
public:
    HOSC_TEMPLATE();
	~HOSC_TEMPLATE(){
		current_project_ = NULL;
	};

public:
	virtual bool initialize( const Parameters_handler* parameters, Error_messages_handler* errors );
	virtual int execute( GsTL_project* proj=0 );
	virtual std::string name() const { return "HOSC_TEMPLATE"; }
	static Named_interface* create_new_interface( std::string& ) { return new HOSC_TEMPLATE; }

public:
    bool get_Source_Data(const Parameters_handler* parameters, Error_messages_handler* error_mesgs );
    bool get_Target_Data(const Parameters_handler* parameters, Error_messages_handler* error_mesgs );

	bool get_Keys_Data(const Parameters_handler* parameters, Error_messages_handler* error_mesgs);
	bool get_Calculate_Target(const Parameters_handler* parameters, Error_messages_handler* error_mesgs);
	bool get_Calculate_Algorithm(const Parameters_handler* parameters, Error_messages_handler* error_mesgs);
	bool get_Weighted_Options(const Parameters_handler* parameters, Error_messages_handler* error_mesgs);
	bool get_Template_Sampling(const Parameters_handler* parameters, Error_messages_handler* error_mesgs);
	bool get_My_Shape(const Parameters_handler* parameters, Error_messages_handler* error_mesgs);
	bool get_Output_Slice_Options(const Parameters_handler* parameters, Error_messages_handler* error_mesgs);
	bool get_Output_Slice_Target(const Parameters_handler* parameters, Error_messages_handler* error_mesgs);
	bool get_Sliding_Window(const Parameters_handler* parameters, Error_messages_handler* error_mesgs);
	bool get_Windows_Moment(const Parameters_handler* parameters, Error_messages_handler* error_mesgs);
	int calculateMomentWindow();



	int set_Target_Data(double** targetGrid, int gridSize, int valuePos);

	RGrid* CreateMyGrid(GridShape* tempGridShape, const char* gridName);
	std::string build_output_grid_base_name(const GridShape* tempGridShape, const char* role) const;
	RGrid* CreateCartesianGridWithProperty(
		GridShape* tempGridShape,
		const char* gridName,
		const std::string& propertyName,
		GsTLGridProperty** out_property,
		std::string* created_name);
	bool create_3d_slice_output();
	bool create_3d_slice_output_from_selected_grid();
	int set_MyGrid_Data(double** targetGrid, int gridSize, int valuePos);
	bool write_Template_Output_Validity(double** targetGrid, int gridSize, int valuePos);
	bool write_Template_Exact_Support(double** targetGrid, int gridSize, int valuePos);



	void TestMyGrid();

public:

    Geostat_grid* source_geostat_grid_;
    RGrid* source_grid_;
    int source_node_count_;
    int source_nx_;
    int source_ny_;
    int source_nz_;
    GsTLGridProperty* source_x;
	GsTLGridProperty* source_y;
	GsTLGridProperty* source_z;
	GsTLGridProperty* source_val;
	double** grid;


    RGrid* target_grid_;
	GsTLGridProperty* target_property_;
	string target_prop_name_;
	string target_dim_num;


	RGrid* auto_grid_;
	GsTLGridProperty* auto_property_;
	string auto_prop_name_;
	string auto_grid_name_;
	RGrid* slice_grid_;
	GsTLGridProperty* slice_property_;
	string slice_grid_name_;
	GsTL_project* current_project_;
	GridShape myGridShape;
	TemplateOutputLayout output_layout_;



#if 0

	RGrid* getDefaultGrid();
	RGrid* default_grid_;
	GsTLGridProperty* default_property_name_;
#endif

	enum Algorithm al;
	bool weighted_enabled_;
	bool template_sampling_exact_;
	std::vector<ExactTemplateDirection> integer_direction_vectors_;
	weighted_moment_core::WeightedMomentOptions weighted_options_;
	int crossection[10];
	int crossectionNum;
	bool slice_output_enabled_;
	int slice_plane_mode_;
	int slice_position_;
	string slice_source_grid_name_;
	string slice_source_property_name_;




#if 0

	int momentWindowL;
	int momentWindowW;
	int momentWindowH;
	int momentOrder;
	int stepX;
	int stepY;
	int stepZ;
	GridShape myWindowShape;
#endif


	Keys* Key;
};



#endif
