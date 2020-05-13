// $Id: dbgmodel.cpp,v 1.1.6.2 2004/08/19 15:05:46 jmccaa Exp $
// Implementation of DbgModel class

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>

#ifndef _RWSTD_NO_NAMESPACE
  using namespace std;
#endif

#include "dbgmodel.h"
#include "field.h"
#include "dataset.h"
#include "manager.h"
#include "msgdlg.h"
#include "fortran.h"
#include "ipc.h"                // for fieldinfo struct

DbgModel::DbgModel()
{
    cout << "DbgModel: entered DbgModel()" << endl;
    _dbgdata.resize( MAX_LEVELS );
}


DbgModel::~DbgModel()
{
    cout << "DbgModel: entered ~DbgModel()" << endl;
}


void
DbgModel::BuildFieldList()
{
    cout << "DbgModel: entered BuildFieldList() (300 fields)" << endl;
    char name[256];
    char longname[256];
    char units[256];
    char std_units[256];
    real_t mult;
    real_t min;
    real_t max;
    bool is_shown;
    bool is_modifiable;
    bool is_averaged;
    int size;

    for ( int i=0; i<100; i++ ) {
        if ( i%2 ) {
            sprintf( name, "FLD%d_multi", i );
            size = 18;
        }
        else {
            sprintf( name, "FLD%d_single", i );
            size = 1;
        }            
        if ( ! (i%5) ) {
            is_modifiable = true;
        }
        else {
            is_modifiable = false;
        }
        sprintf( longname, "FLD%d long name", i );
        sprintf( units, "FLD%d units", i );
        sprintf( std_units, "FLD%d std_units", i );
        mult = 1.0;
        min = -10;
        max = 10;
        is_shown = true;
        is_averaged = false;
        Field* f = new Field( name, longname, units,
                              std_units, mult, min,
                              max, is_shown, is_modifiable,
                              is_averaged, size, fields.size() );
        fields.insert( pair<string,Field*>( name, f ));
    }
}

bool
DbgModel::CheckDataset( const string& filename, int type )
{
    cout << "DbgModel: entered CheckDataset( " << filename << ", " << type <<endl;
    return Model::CheckDataset( filename, type );
}

int
DbgModel::Init(bool isRestart)
{
    cout << "DbgModel: entered Init()\n" 
         <<  "Initialization parameters:\n";
    for ( int i=0; i<NUM_SWITCHES; i++ ) 
        cout << "Switch #"<< i+1 << " " << switches[i];
    
    cout << "\nlat: " << latitude << ", lon: " << longitude 
         << "\nbasedate: " << baseDate << ", basesecs: " <<  baseSecs
         << "\nruntype: "<< runType << ",  steplen: " << stepLen << '\n';
    
    int i = 0;
    while( dataset[i] != 0 ) {
        cout << Dataset::TypeDescription( i ) << " : " 
             << dataset[i]->name() << '\n';
        i++;
    }
    cout << endl;
    
    currentStep = 0;
    initStatus = 0;
      //
      // get the model base pressure levels
      //
    currLevs = baseLevs = GetLevels();
    currILevs = baseILevs = GetILevels();

    return initStatus;
}

//
//  Retrieve field data from the model
//
void
DbgModel::ReadField( Field& f )
{

    cout << "DbgModel: entered ReadField( " << f.Name() << " )" << endl;
    if ( currentStep > 0 )
        for ( int i=0; i<MAX_LEVELS; i++ )
            _dbgdata[i] += (((rand()-(RAND_MAX/2.0))/RAND_MAX) / f.PlotMult()) / 10.0;   
    f.SetData( &_dbgdata[0] );
}

//
//  Retrieve the current model pressure levels
//
const vector<real_t>&
DbgModel::GetLevels()
{
    cout << "DbgModel: entered Levels()" << endl;
    numLevs = 18;
    currLevs.resize( numLevs );
    for( int i = 0; i<MANAGER.NumLevs(); i++ )
        currLevs[i] = i*1000/numLevs;
    return currLevs;
}

//
//  Retrieve the current model pressure levels
//
int
DbgModel::GetNumLevels()
{
    cout << "DbgModel: entered Levels()" << endl;
    numLevs = 18;
    return numLevs;
}

//
//  Retrieve the current model pressure levels
//
const vector<real_t>&
DbgModel::GetILevels()
{
    cout << "DbgModel: entered Levels()" << endl;
    numILevs = 19;
    currILevs.resize( numILevs );
    for( int i = 0; i<MANAGER.NumILevs(); i++ )
        currILevs[i] = i*1000/numILevs;
    return currILevs;
}

//
// Resets an averaged field
//
void
DbgModel::ResetField( const Field& f )
{
    cout << "DbgModel: entered ResetField( " << f.Name() << " )" << endl;
}

//
// Sets a field's value in the model
//
void
DbgModel::WriteField( const Field& f )
{
    cout << "DbgModel: entered WriteField( " << f.Name() << " )" << endl;
    for ( int i=0;i<f.NumLevs(); i++ )
        _dbgdata[i] = f[i];
}

void
DbgModel::Step()
{
    cout <<  "DbgModel: Step() : " << currentStep << endl;
    currentStep++;
}




