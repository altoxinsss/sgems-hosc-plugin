/**********************************************************************
** Author: Nicolas Remy
** Copyright (C) 2002-2004 The Board of Trustees of the Leland Stanford Junior
**   University
** All rights reserved.
**
** This file is part of the "math" module of the Geostatistical Earth
** Modeling Software (GEMS)
**
** This file may be distributed and/or modified under the terms of the 
** license defined by the Stanford Center for Reservoir Forecasting and 
** appearing in the file LICENSE.XFREE included in the packaging of this file.
**
** This file may be distributed and/or modified under the terms of the
** GNU General Public License version 2 as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL included in the
** packaging of this file.
**
** This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
** WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
**
** See http://www.gnu.org/copyleft/gpl.html for GPL licensing information.
**
** Contact the Stanford Center for Reservoir Forecasting, Stanford University
** if any conditions of this licensing are not clear to you.
**
**********************************************************************/


#include <GsTLAppli/geostat/spatial_correlation_computer.h>
#include <GsTLAppli/math/correlation_measure.h>
#include <GsTLAppli/grid/grid_model/geostat_grid.h>
#include <GsTLAppli/grid/grid_model/gval_iterator.h>
#include <GsTLAppli/grid/grid_model/grid_property.h>

#include <GsTL/geometry/geometry_algorithms.h>
#include <GsTL/math/math_functions.h>


Variogram_parameters::Variogram_parameters() {
  lags_count = 0;
  lag_separation = 0;
  lag_tolerance = 0;
}


//==================================================

Spatial_correlation_computer::
Spatial_correlation_computer( Correlation_measure* measure )
  : prototype_( 0 ) {
  if( measure )
    prototype_ = measure->clone();
}


Spatial_correlation_computer::~Spatial_correlation_computer() {
  clean();
}

void Spatial_correlation_computer::clean() {
  for( unsigned int i = 0 ; i < measures_.size() ; i++ ) {
    delete measures_[i];
    measures_[i] = 0;
  } 
}

void Spatial_correlation_computer::resize( unsigned int i ) {
  clean();
  measures_.resize( i );
  for( unsigned int j = 0 ; j < measures_.size() ; j++ )
    measures_[j] = prototype_->clone();
}


void Spatial_correlation_computer::
set_correlation_measure( Correlation_measure* measure ) {
  delete prototype_;
  prototype_ = 0;
  if( measure )
    prototype_ = measure->clone();
}



void Spatial_correlation_computer::
compute( const Geostat_grid* grid, 
         const GsTLGridProperty* head_prop, const GsTLGridProperty* tail_prop, 
         const Variogram_parameters& parameters ) {


  // check that the properties belong to the grid
  if( !grid || !head_prop || !tail_prop ) return;
  if( !grid->property( head_prop->name() ) || 
      !grid->property( tail_prop->name() ) ) return;

  measures_.resize( parameters.lags_count );
  if( measures_.empty() ) return;


  typedef Geostat_grid::location_type::difference_type EuclideanVector;

  Geostat_grid::const_iterator head_it = grid->begin( head_prop ) ;
  Geostat_grid::const_iterator head_end = grid->end( head_prop ) ;

  for( ; head_it != head_end ; ++head_it ) {
    if( !head_it->is_informed() ) continue;

    Geostat_grid::const_iterator tail_it = grid->begin( tail_prop ) ;
    Geostat_grid::const_iterator tail_end = grid->end( tail_prop ) ;
    for( ; tail_it != tail_end ; ++tail_it ) {
      if( !tail_it->is_informed() ) continue;
  
      EuclideanVector vec( tail_it->location() - head_it->location() );
      if( !parameters.direction.is_colinear( vec ) ) continue;

      int lag = find_lag( parameters, euclidean_norm( vec ) );
      if( lag < 0 || lag >= measures_.size() ) continue;

      measures_[lag]->add_pair( head_it->property_value(), 
                                tail_it->property_value() );

    }
  }
}

  

std::vector<float> Spatial_correlation_computer::computed_values() {
  std::vector<float> result( measures_.size() );
  for( unsigned int i = 0 ; i < result.size() ; i++ )
    result[i] = measures_[i]->correlation();

  return result;
}



int Spatial_correlation_computer::
find_lag( const Variogram_parameters& parameters, float dist ) {
  int div = GsTL::ceil( dist / float( parameters.lag_separation ) );
  float lag_dist = div * parameters.lag_separation;

  if( lag_dist - parameters.lag_tolerance < dist && 
      dist < lag_dist + parameters.lag_tolerance ) return div;

  return -1;
}




//===========================================
void RGrid_correlation_computer::
compute( const Geostat_grid* grid, 
         const GsTLGridProperty* head_prop, const GsTLGridProperty* tail_prop, 
         const Variogram_parameters& parameters ) {


  // check that the properties belong to the grid
  if( !grid || !head_prop || !tail_prop ) return;
  if( !grid->property( head_prop->name() ) || 
      !grid->property( tail_prop->name() ) ) return;


}


