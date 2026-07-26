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

#ifndef __GSTLAPPLI_GEOSTAT_SPATIAL_CORRELATION_COMPUTER_H__
#define __GSTLAPPLI_GEOSTAT_SPATIAL_CORRELATION_COMPUTER_H__

#include <GsTLAppli/geostat/common.h>
#include <GsTLAppli/math/direction_3d.h>

#include <vector>


class Geostat_grid;
class GsTLGridProperty;
class Correlation_measure;



class Variogram_parameters {
public:
  Variogram_parameters();

  int lags_count;
  float lag_separation;
  float lag_tolerance;
  Direction_3d direction;
};



/** Computes a 1D variogram for un-structured sets of points
*/
class GEOSTAT_DECL Spatial_correlation_computer {
public:
  Spatial_correlation_computer( Correlation_measure* measure = 0 );
  ~Spatial_correlation_computer();

  void set_correlation_measure( Correlation_measure* measure );
  Correlation_measure* correlation_measure() { return prototype_; }


  virtual void compute( const Geostat_grid* grid, 
                        const GsTLGridProperty* head_prop,
                        const GsTLGridProperty* tail_prop, 
                        const Variogram_parameters& parameters ); 
  virtual std::vector<float> computed_values();


protected:
  int find_lag( const Variogram_parameters& parameters, float dist );
  void clean();
  void resize( unsigned int i );
  Spatial_correlation_computer( const Spatial_correlation_computer& ) {}


protected:
  std::vector<Correlation_measure*> measures_;
  Correlation_measure* prototype_;
};



/** Specialized for regular grids
*/
class GEOSTAT_DECL RGrid_correlation_computer : 
  public Spatial_correlation_computer {

public:
  virtual void compute( const Geostat_grid* grid, 
                        const GsTLGridProperty* head_prop,
                        const GsTLGridProperty* tail_prop, 
                        const Variogram_parameters& parameters ); 

};

#endif

