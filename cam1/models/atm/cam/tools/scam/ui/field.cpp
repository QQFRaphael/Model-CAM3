/*------------------------------------------------------------------------*
 * File: field.cpp 
 * $Author: jet $
 * $Id: field.cpp,v 1.1.6.1 2004/05/18 20:51:52 jet Exp $ *
 *------------------------------------------------------------------------*/

#ifndef lint
static char rcsid[] = "$Id: field.cpp,v 1.1.6.1 2004/05/18 20:51:52 jet Exp $";
#endif /* lint */

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "field.h"
#include "fortran.h"
#include "plot.h"
#include "manager.h"
#include "utils.h"
#include "ncfile.h"


// use the following defaults for min and max if the bldfld call does not specify them
const real_t DEFAULT_MIN = -10;
const real_t DEFAULT_MAX =  10;

using namespace ncfile;
//
// constructor
//
Field::Field ( const string& name, const string& longName, const string& plotUnits,
               const string& stdUnits, real_t plotMult, real_t minHint,         
               real_t  maxHint, bool isListed, bool isModifiable,    
               bool isAveraged, int numLevs, int fieldID ) 
    : _fieldID( fieldID ), _name( name ), _longName( longName ), _plotUnits( plotUnits ),
      _stdUnits( stdUnits ), _isModifiable( isModifiable ), _isAveraged( isAveraged ), 
      _isListed( isListed ), _minHint( minHint ), _maxHint( maxHint ), _plotMult( plotMult ),
      _numLevs( numLevs )
{

    if ( _minHint == _maxHint ) {
        _minHint = DEFAULT_MIN;
        _maxHint = DEFAULT_MAX;
    }
    
    _isSaved = false;
      // for single level fields, we'll use the data buffer
      //  for saving previous timepoints' data
      // (used for plotting and displaying points)
    if ( _numLevs == 1 ) 
        _data.resize( MAX_POINTS );
    else
        _data.resize( _numLevs );

    if ( _isModifiable )
        _restartData.resize( _data.size() );
    
      // 
      // note that _histVar is initialized in the History() constructor
      // 
    _histVar = 0;

}

Field::Field ( const NcVariable<real_t>& v )
    : _name( v.name() )
{
    _name = ( v.name() );
    _longName = ( NcCharAtt( v.attribute( "long_name" ))).value();
    _stdUnits = ( NcCharAtt( v.attribute( "units" ))).value();
    _plotUnits = ( NcCharAtt( v.attribute( "plot_units" ))).value();
    _numLevs = ( NcIntAtt( v.attribute( "num_levels" )));
    _plotMult =  ( NcFloatAtt( v.attribute( "plot_multiplier" )));

    _minHint = 0;
    _maxHint = 0;
    _isAveraged = false;
    _isListed = false;
    _isModifiable = false;
    _isSaved = false;
    _fieldID = -1;
    _histVar = 0;
}

Field::~Field()
{
    delete _histVar;
}

  
//
// copy the restart data, reset plots
//
void
Field::Reset()
{
    MANAGER.ResetField( *this );
      //
      // if this is a modifiable field, 
      // copy the restart data buffer to the data buffer...
    if ( IsModifiable() ) {
        _data = _restartData;
        MANAGER.WriteField( *this ); // ... and write it out
    }
      // ... else just zero out the diagnostic field's data buffer
    else 
        fill( _data.begin(), _data.end(), 0.0 );
    
    SetChanged();
    NotifyObservers();
    ClearChanged();
}
//
// copy the restart data, reset plots
//
void
Field::Zero()
{
      //
      // if this is a modifiable field, 
      // zero
    if ( IsModifiable() ) {
        fill( _data.begin(), _data.end(), 0.0 );
        MANAGER.WriteField( *this ); // ... and write it out
    }
      // ... else just zero out the diagnostic field's data buffer
    else 
        fill( _data.begin(), _data.end(), 0.0 );
    
    MANAGER.InitModel();
    SetChanged();
    NotifyObservers();
    ClearChanged();
}


//
// set element of field's data array - called when user wants
//  to modify a field through GUI
//
void
Field::SetDataPoint( real_t value, int index )
{

    _data[index] = value;
    MANAGER.WriteField( *this );     // change the value in the model
    SetChanged();
    NotifyObservers();
    ClearChanged();
}

//
// set all of field's data array 
//
void
Field::SetData( real_t* data )
{
    std::copy( data, data + _numLevs, _data.begin() );
    SetChanged();
    NotifyObservers();
    ClearChanged();
}

//
// copy  data to restart data for modifiable fields
//
void
Field::SaveRestartData()
{
    if ( IsModifiable() ) {
        MANAGER.ReadField( *this );
        _restartData = _data;
	//	printf("%f save restart****\n",_data[0]); 
    }
}

void 
Field::Update()
{
      // update is called in two situations: the model step was advanced or 
      // the model was reset. We only need to do something if the field
      // is currently being observed (by a plot or point list);
      // just need to read the new value from the model if we're being observed.
      // (the manager will handle reading the field if it is being saved)
    if ( NumObservers() > 0 ) {
        if ( !IsMultiLevel() ) {
            if ( _lastStepPlotted != MANAGER.CurrentStep() ) {
                  // shift values over by one in the data array
                  // because _data is a deque (which supports fast insertion at the front)
                  // this is actually pretty fast
                _data.push_front( 0 ); // insert 0 - overwritten by ReadField call below
                _data.pop_back();
                _lastStepPlotted = MANAGER.CurrentStep();
            }
        }
        MANAGER.ReadField( *this );
        NotifyObservers();
        ClearChanged();
    }
}

void
Field::WriteHistory( int record )
{
    MANAGER.ReadField( *this );
      // copy the data from the field to the history variable
    for ( int i=0; i<_numLevs; i++ ) 
        (*_histVar)[i] = _data[i];
      // _histVar has already been initialized in History() constructor
    (*_histVar).writeRecord( record );
    if ( _isAveraged ) 
        MANAGER.ResetField( *this );
}
