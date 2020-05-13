#include <sstream>

#include <cstdlib>
#include <cstring>

#include "model.h"
#include "dataset.h"
#include "manager.h"
#include "msgdlg.h"


//
// Model: constructor 
//
Model::Model()
{

      // model  initialization  parameters
    for ( int n=0; n<NUM_SWITCHES; n++ )
        switches[n] = 0;
    baseDate = 0;
    baseSecs = 0;
    runType = 0;
    latitude = 0;
    longitude = 0;
    initStatus = 0;           
    stepLen = 0;
    MANAGER.SetSeedVal(0);
    currentStep = 0;
    endStep = 0;
    restart = 0;
    for ( int n=0; n<NUM_INIT_FILES; n++ ) 
        dataset[n] = 0;
    
}

// 
//  checks out netCDF model initialization datasets, making sure that
//  all required variables are present, and that the latitude and 
//  longitudes match the initial conditions dataset
//  returns true if dataset checks out OK
//
bool
Model::CheckDataset( const string& filename, int type )
{
    string missingVars;
    try {
        Dataset data( filename, type );
        if ( !data.CheckVars( missingVars ) ){
            ShowMsg( __FILE__, __LINE__, 
                     "Dateset check failed for %s!\nRequired variables \"%s\" are missing",
                     filename.c_str(), missingVars.c_str() );
            return false;
        }
        if ( type == MODEL ) {
            initLatIdx = data.FindIdx( "lat", latitude ); 
            initLonIdx = data.FindIdx( "lon", longitude );
        }

          // don't need to check lats/lons in ozone dataset - they are 
          // interpolated
        if ( type == ANAL ) {
            if ( !data.CheckLat( latitude, initLatIdx ) ) {
                ShowMsg( __FILE__, __LINE__, 
                         "Latitudes in %s do not correspond to latitudes in initial conditions dataset %s", 
                         filename.c_str(), dataset[MODEL]->name().c_str()  );
                return false;
            }
            if ( !data.CheckLon( longitude, initLonIdx ) ) {
                ShowMsg( __FILE__, __LINE__, 
                         "Longitudes in %s do not correspond to longitudes in initial conditions dataset %s",
                         filename.c_str(), dataset[MODEL]->name().c_str() );
                return false;
            }
        }
        return true;
    } catch ( NcErr& e ) {
        ShowMsg( __FILE__, __LINE__, "ERROR: Unable to open dataset: %s\n", e.toString().c_str() );
        return false;
    }    
}

//
// return the state of the switches
//
bool          
Model::SwitchState( int index )
{
    return switches[index]; 
}

//
// return a pointer to an initialization dataset
//
const NcFile*
Model::GetDataset( int type ) 
{
      //
      // we need to have a dummy dataset available so 
      // that if a particular type of dataset hasn't been
      // set yet, we can still return a valid pointer. 
      // the dummy dataset has an invalid type, and it's 
      // name is ""
      //

    static Dataset dummy( -1 );

    if ( dataset[type] != 0 )
        return dataset[type];
    else
        return &dummy;
}

//
// return a pointer to an initialization dataset
//
const Ascii_Dataset*
Model::GetAsciiDataset( int type ) 
{
      //
      // we need to have a dummy dataset available so 
      // that if a particular type of dataset hasn't been
      // set yet, we can still return a valid pointer. 
      // the dummy dataset has an invalid type, and it's 
      // name is ""
      //

    static Ascii_Dataset dummy( -1 );

    if ( adataset[type] != 0 )
        return adataset[type];
    else
        return &dummy;
}

//
//  return a variable's data from one of the initialization datasets
//
NcCharAtt
Model::DatasetCharAtt( int datasetType, const string& attname ) throw (NcErr)
{
    NcCharAtt v( dataset[datasetType]->attribute( attname ) );
    return v;
}

//
//  return a variable's data from one of the initialization datasets
//
NcVariable<real_t>    
Model::DatasetRealVar( int datasetType, const string& varname ) throw (NcErr)
{
    NcVariable<real_t> v( dataset[datasetType]->variable( varname ) );
    v.read();
    return v;
}



//
//  return a variable's data from one of the initialization datasets
//
NcIntVar      
Model::DatasetIntVar( int datasetType, const string& varname ) throw (NcErr)
{
    NcIntVar v( dataset[datasetType]->variable( varname ) );
    v.read();
    return v;
}

// 
//  set a Model Dataset if it checks out OK
//
void
Model::SetDataset( int type, const string& filename ) throw (NcErr)
{
    string missingVars;

    Dataset* dset = new Dataset( filename, type );
//      if ( !dset->CheckVars( missingVars ) ) {
//          stringstream error;
//          error << "Required variable(s) \"" <<  missingVars
//                <<  "\" are missing from dataset";
//          throw NcErr( filename, missingVars, error.str(), dset->ncid() );
//      }
//        //
//        // check lats/lons against MODEL dataset or SIC dataset
//        //
      if ( type == MODEL ) {
          initLatIdx = dset->FindIdx( "lat", latitude ); 
          initLonIdx = dset->FindIdx( "lon", longitude );
      }
//      if ( type == ANAL   ) 
//          if ( !dset->CheckLat( latitude, initLatIdx ) ) 
//              throw NcErr( filename, "lat", "latitude does not correspond to initial dataset latitude", dset->ncid() );
//      if ( type == ANAL   ) 
//          if ( !dset->CheckLon( longitude, initLonIdx ) )
//              throw NcErr( filename, "lon", "longitude does not correspond to initial dsetset latitude", dset->ncid() );
    delete dataset[type];
    dataset[type] = dset;
}

// 
//  set a Model Dataset if it checks out OK
//
void
Model::SetAsciiDataset( int type, const string& filename )
{
    Ascii_Dataset* adset = new Ascii_Dataset( filename, type );
    adataset[type]=adset;
}

//
// set the latitude of the desired column
//
void
Model::SetLat( real_t lat)
{
      // must be in range of -90:90 degrees
    if ( lat >= -90 && lat <= 90 )
        latitude = lat;
    else {
        ShowMsg( __FILE__, __LINE__, "ERROR: Latitude %f is outside of valid range (-90:90)", lat );
        exit ( -1 );
    }
}

//
// set the longitude of the desired column
//
void
Model::SetLon( real_t lon )
{
      // must be in range of -180:180 degrees
    if ( lon >= -180 && lon <= 180 )
        longitude = lon;
    else {
        ShowMsg( __FILE__, __LINE__, "ERROR: Longitude %f is outside of valid range (-180:180)", lon );
        exit ( -1 );
    }
}

//
//  Set the state variable runType 
//
void
Model::SetRunType( int type )
{
    if ( type < 0 || type > MAX_RUNTYPE ) {
        cerr << " ERROR: "__FILE__"line, "<< __LINE__ 
             << ": Model::SetRunType(): unknown type " << type << endl;
        exit( -1 );
    }
    runType = type;
}

void
Model::SetSwitch( bool state, int index )
{
    switches[index] = state;
}


Model::~Model()
{
    for ( int i=0;i<NUM_INIT_FILES;i++ )
        delete dataset[i];
}         

