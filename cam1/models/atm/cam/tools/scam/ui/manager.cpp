/*------------------------------------------------------------------------*
 * File: manager.cpp 
 * $Author: jmccaa $
 * $Id: manager.cpp,v 1.1.6.3 2004/08/19 15:05:49 jmccaa Exp $ *
 *------------------------------------------------------------------------*/

#ifndef lint
static char rcsid[] = "$Id: manager.cpp,v 1.1.6.3 2004/08/19 15:05:49 jmccaa Exp $";
#endif /* lint */


#include <iostream>
#include <string>
//#include <strstream> /* this header will eventually be replaced with <sstream> */
#include <sstream> /* this header will eventually be replaced with <sstream> */

#include <climits>            
#include <cstdlib>
#include <unistd.h>            // getpid()
#include <ctype.h>
#include <math.h>

#include "manager.h"
#include "dataset.h"
#include "defaults.h"
#include "model.h"
#include "dbgmodel.h"
#include "field.h"
#include "fifomodel.h"
#include "rpcmodel.h"
#include "msgdlg.h"
#include "history.h"
#include "sicfile.h"
#include "runtype.h"
#include "utils.h"
#include "crm.h"

// initialize static member 
Manager* Manager::_instance = 0;

Manager::Manager()
{
      // model state initialization
    modelInited = false;
    histFile = 0;    
    defaultsFile = DEFAULTS_FILE;
}


Manager::~Manager()
{
      // remove the temporary history file in case it hasn't been saved
      // to avoid build up of temporary history files.
    stringstream tmp;
    tmp << userDataDir << "/.scamhist.tmp." << getpid();
    remove( tmp.str().c_str() );
    //    delete theModel;
    //    delete histFile;
}

void
Manager::SetDefaultsFile(const string& fileName)
{
  defaultsFile = fileName;
}

string 
Manager::GetDefaultsFile()
{
  return defaultsFile;
}

int Manager::BaseDate() 
{ 
    return theModel->BaseDate();
}

int Manager::BaseSecs() 
{ 
    return theModel->BaseSecs();
}

void
Manager::CreateModel( ModelType type )
{
  //    delete theModel;

    switch( type ) {
    case Model::FIFO:
        theModel = new FifoModel;
        break;
#ifdef SUN
    case Model::RPC:
        theModel = new RPCModel;
        break;
#endif        
    case Model::DBG:
        theModel = new DbgModel;
        break;
    default:
        cerr << "ERROR: "__FILE__":" << __LINE__
             << " Manager::SetModelType(): invalid type " 
             << int(type) << endl;
        exit( -1 );
    }
}


const NcFile&
Manager::GetDataset( int type ) 
{
    return *theModel->GetDataset( type ); 
}
const Ascii_Dataset&
Manager::GetAsciiDataset( int type ) 
{
    return *theModel->GetAsciiDataset( type ); 
}

NcVariable<real_t>
Manager::DatasetRealVar( int type, const string& varname ) throw (NcErr)
{
    return theModel->DatasetRealVar( type, varname ); 
}

NcCharAtt
Manager::DatasetCharAtt( int type, const string& attname ) throw (NcErr)
{
    return theModel->DatasetCharAtt( type, attname ); 
}

NcIntVar
Manager::DatasetIntVar( int type, const string& varname ) throw (NcErr)
{
    return theModel->DatasetIntVar( type, varname ); 
}

const FieldList&
Manager::GetFieldList() 
{ 
    return theModel->fields;
}


//
// return pointer to a field specified by name, or 0 if not found in the list
//
Field* 
Manager::FindField( const string& name ) { 
    Field* notFound = 0;
    FieldListItr it;
    if (( it = theModel->fields.find(name) )
        != theModel->fields.end() )
        return (*it).second;
    else
        return notFound;
}

string 
Manager::HistName()     
{ 
    return HistFile; 
}

//
// create and initialize the temp history file 
//
void
Manager::InitHistoryFile() throw (NcErr)
{

    stringstream tmp;
    tmp << userDataDir << ".scamhist.tmp." << getpid();
    
    tmpHistFile = tmp.str();

    delete histFile;
      
      // create a new history file using the value of tmpHistfile for the name
      // and the value of defaultsFile for the 'case' attribute.
    histFile = new History( tmpHistFile, defaultsFile, NcFile::CREATE, true );
}

//
// Initialize the model to starting state, initialize the history file
//
bool
Manager::InitModel( )
{
  if (IsModelInited()) {
    FieldListItr it = theModel->fields.begin(), end = theModel->fields.end();
    for (; it != end; ++it ) {
      Field* f = (*it).second;
      f->Reset();
    }
  }

    int initStatus = theModel->Init(IsModelInited());
    
      // Check for error conditions:
      // errors in Fortran initialization code will be returned in initStatus
      // (initStatus corresponds to filetypes where errors occurred)
    if ( initStatus != 0 ) {
        ShowMsg( __FILE__, __LINE__, 
                 "ERROR: Model initialization failed!\n Unable to load %s file: %s",
                 Dataset::TypeDescription(initStatus).c_str(), 
                 GetDataset(initStatus).name().c_str() );
        return ( modelInited = false );
    }
      // make sure that the number of levels that the model is defined on
      // is within allowable bounds
    if ( NumLevs() > MAX_LEVELS ) {
        ShowMsg( __FILE__, __LINE__, 
                 "ERROR: Number of model levels (%d) is greater \n than the maximum allowed(%d)",
                 NumLevs(), MAX_LEVELS) ;
        return ( modelInited = false );
    }
      // make sure that the pressure levels dataset has the same number of 
      // levels as the model was compiled for
    if ( NumLevs() > numDataLevs ) {
        ShowMsg( __FILE__, __LINE__, "ERROR:\n Number of levels in pressure levels dataset\n%s is %d\n Number of levels the model is compiled for is %d\n To use this dataset, you must recompile the model with PLEV set to %d", 
                        GetDataset(MODEL).name().c_str(), numDataLevs, NumLevs(), numDataLevs );
        modelInited = false;
        return false;
    }

    modelInited = true;

      //
      // build the field list (the field list may already have
      // been created, if this is an ensemble run, so make sure
      // it's empty first)
      //
    if ( NumFields() == 0 ) {
        theModel->BuildFieldList(); // Build the master list of fields
          // set saved fields from defaults file
        try { 
            SetDefaultSavedFields();
        } catch ( DfltsErr& e ) {
            ShowMsg( __FILE__, __LINE__,"Unable to set saved fields: %s\n", e.toString().c_str() );
        }
          // Add the fields to the list of observers 
        FieldListItr it = theModel->fields.begin(), end = theModel->fields.end();
        for (; it != end; ++it ) {
            AddObserver( (*it).second );
        }
	it = theModel->fields.begin();
	for (; it != end; ++it ) {
	  Field* f = (*it).second;
	  f->SaveRestartData();
	}
    }
      //---------------------------------------------------------------
      //
      // save field restart data
      //
      //---------------------------------------------------------------
      //---------------------------------------------------------------
      //
      // initialize the fields
      //
      //---------------------------------------------------------------
       else
	 if (IsModelInited()) {
	   FieldListItr it = theModel->fields.begin(), end = theModel->fields.end();
	   //	   for (; it != end; ++it ) {
	   //	     Field* f = (*it).second;
	   //	     f->SaveRestartData();
	   //	   }

	   //	   it = theModel->fields.begin();
	   for (; it != end; ++it ) {
	     Field* f = (*it).second;
	     f->Reset();
	     //	     f->print0();
	   }
	 }
      // initialize the history file
    try { 
      InitHistoryFile();
    } catch ( NcErr& e ) {
        ShowMsg( __FILE__, __LINE__, "ERROR: Unable to create the history file: %s\n", 
                 e.toString().c_str(), tmpHistFile.c_str() );
        exit( -1 );
    }

      // notify observers (i.e., GUI and fields) that the model state has changed
    SetChanged();
    NotifyObservers();
    ClearChanged();
    return true;
}

//
//  return number of seconds between the base date of the iop file, and the base
//  date of the model run
//
int Manager::IopStartOffset()
{
    return iopStartOffset; 
} 


bool Manager::IsShowingSettings()
{ 
    return showingSettings; 
}

bool          
Manager::IsValidDataset( const string& name, int type )
{
    try {
      // first see if the lat lon values exist in the initial dataset
      // this will also set initlat and initLon 
        theModel->CheckDataset(  GetDataset( MODEL ).name(), MODEL ); 
        return theModel->CheckDataset( name, type ); 
    } catch ( NcErr& e ) {
        return false;
    }
}

real_t Manager::Lat() 
{ 
    return theModel->Lat(); 
}

//
//  load an scam defaults file (including quickstart files)
//
void
Manager::LoadDefaults( const string& defaultsFile, bool isQuickstart ) throw ( IOErr )
{    
    this->defaultsFile =  defaultsFile;

      // set the default paths, variables, etc.    
    SetDefaults( defaultsFile, isQuickstart );
    
      // set the IOP max steps 
    if ( RunType() == IOP || RunType() == USER  )
        SetIopMaxSteps();

    if ( IsShowingSettings() )
        ShowSettings();

    if ( modelInited )
        SetDefaultSavedFields();
}

//
//  Load a Saved Initial Conditions (SIC) file
//
void
Manager::LoadInitialConditionsFile( const string& sicfilename ) throw ( NcErr )
{
      // Sicfile constructor will set all the 
      // necessary defaults
    Sicfile sic( sicfilename, NcFile::READ, true );
      // 
      // set the runtype to reflect that this is 
      //  a saved initial condition 
      //
    if (RunType()<SIC)
      theModel->SetRunType( RunType() + SIC );
    else
      theModel->SetRunType(SIC );
}

real_t Manager::Lon() 
{
    return theModel->Lon(); 
}

//
//  return iopMaxSteps if this is an IOP run, else INT_MAX
//
int
Manager::MaxStep()
{
    if ( RunType() == IOP || RunType() == USER )
        return iopMaxSteps;
    else
        return INT_MAX;         // defined in /usr/include/limits.h
}

//
//   Flush the netcdf output to disk and then rename it
//
void
Manager::SaveHistoryFile( const string& filename ) throw ( IOErr )
{    
      // make sure all data is written to disk
    histFile->flush();
      // rename the history file;
      // if it has already been saved, no need to rename it
      // ( in fact rename() will fail so don't try it )
    if ( filename == tmpHistFile )
        return;
    if (( rename( tmpHistFile.c_str(), filename.c_str() )) == -1 ) {
        stringstream errMsg;
        errMsg << "ERROR: Unable to rename history file (possibly trying to save across filesystems?) Output is in " << tmpHistFile;
        throw IOErr( errMsg.str().c_str(), filename );
    }
    else 
        HistFile = tmpHistFile = filename;

    return;
}

//
// Save the default model parameters
//
void
Manager::SaveDefaults( const string& filename, bool quickStart ) throw ( DfltsErr )
{
    Defaults d( filename, Defaults::WRITE );
      // quickstart files save a couple of extra parameters
    if ( quickStart ) {
        d.WriteDefault( "runtype", RunType() );
        d.WriteDefault( "basedate", BaseDate() );
        d.WriteDefault( "basesecs", BaseSecs() );
        d.WriteDefault( "iopstartoffset", IopStartOffset() );
    }
    d.WriteDefault( "globaldatadir", globalDataDir );
    d.WriteDefault( "iopdatadir", iopDataDir );
    d.WriteDefault( "boundarydatadir", boundaryDataDir );
    d.WriteDefault( "userdatadir", userDataDir );
    d.WriteDefault( "histfile", FileFromPath( HistFile ).c_str() );
    d.WriteDefault( "lat", Lat() );
    d.WriteDefault( "lon", Lon() );
    d.WriteDefault( "steplen", StepLen() );
    d.WriteDefault( "endstep", EndStep() );
    d.WriteDefault( "savefreq", saveFreq );
    d.WriteDefault( "timedisplayformat", TimeFormat() );
    d.WriteDefault( "showsettings", showingSettings );
    d.WriteDefault( "analysisfile", GetDataset( ANAL ).name() );
    d.WriteDefault( "modelfile", GetDataset( MODEL ).name() );
    d.WriteDefault( "userfile",  GetDataset( USER ).name() );
    d.WriteDefault( "iopfile", GetDataset( IOP ).name() );
    d.WriteDefault( "lsminifile", GetDataset( LSMINI ).name() );
    d.WriteDefault( "ozonfile", GetDataset( OZON ).name() );
    d.WriteDefault( "pressfile", GetDataset( PRES ).name() );
    d.WriteDefault( "absemsfile", GetDataset( ABSEMS ).name() );
    d.WriteDefault( "aeropticsfile", GetDataset( AEROPTICS ).name() );
    d.WriteDefault( "aermassfile", GetDataset( AERMASS ).name() );
    d.WriteDefault( "lsmsurffile", GetDataset( LSMSRF ).name() );
    d.WriteDefault( "lsmpftfile", GetAsciiDataset( LSMPFT ).name() );
    d.WriteDefault( "sstfile", GetDataset( SST ).name() );
    d.WriteDefault( "sicfile", GetDataset( SIC ).name() );
      // write out the list of saved fields
    string savefields;
    FieldListItr it = theModel->fields.begin(), end = theModel->fields.end();
    for (; it != end; ++it ) {
        Field* f = (*it).second;
        if ( f->IsSaved() ) 
             savefields += f->Name() + " ";
        
    }
    d.WriteDefault( "savefields", savefields );
    
      // write the switch descriptions
    for ( int i=0; i < NUM_SWITCHES; i++ ) {
        stringstream desc;
	//        desc << "switch_desc" << i+1 << ends;
        desc << "switch_desc" << i+1;
        d.WriteDefault( desc.str(), SwitchDesc( i ).c_str() );
        stringstream name;
	//        name << "switch" << i+1 << ends;
        name << "switch" << i+1;
        d.WriteDefault( name.str(), static_cast< int >( SwitchState( i )));
    }
}

//
//  Save an initial conditions file to disk
//
void
Manager::SaveInitialConditions( const string& filename ) throw ( NcErr )
{
    Sicfile  sic( filename, NcFile::CREATE, true );
}

void  
Manager::SetLat( real_t lat ) 
{ 
    theModel->SetLat( lat ); 
}

void  
Manager::SetLon( real_t lon )
{
    theModel->SetLon( lon); 
}

void  
Manager::SetRemoteHost( const string& host )
{
#ifdef SUN
    static_cast<RPCModel*>(theModel)->SetRemoteHost( host );
    if ( !static_cast<RPCModel*>(theModel)->Connect() ) {
        ShowMsg( __FILE__, __LINE__,"Couldn't connect to server");
        exit( -1 );
    }
#endif
}

void  Manager::SetRunType( int type ) 
{
    theModel->SetRunType( type ); 
    cout <<"Manager settting run type to "<<type<<"\n";
}

void  Manager::SetBaseDate( int bdate ) 
{
    theModel->SetBaseDate( bdate ); 
}

void  Manager::SetBaseSecs( int bsecs )
{
    theModel->SetBaseSecs( bsecs ); 
}

int  Manager::SaveFrequency()
{ 
    return saveFreq; 
}

void Manager::SetIopStartOffset( int secs ) 
{ 
    iopStartOffset = secs; 
}

void Manager::SetDataset( const string& name, int type ) throw ( NcErr )
{
    if ( ! name.empty() )
        theModel->SetDataset( type, name );
}
void Manager::SetAsciiDataset( const string& name, int type )
{
  if ( ! name.empty() ) {
    theModel->SetAsciiDataset( type, name );
  }
}

void Manager::SetDataNumLevs(int nlevs) 
{ 
    numDataLevs = nlevs; 
}
void Manager::SetDataNumILevs(int ilevs) 
{ 
    numDataILevs = ilevs; 
}

void Manager::SetSaveFreq( int newFreq ) 
{ 
    saveFreq = newFreq; 
}

void Manager::SetShowSettings( bool show ) 
{ 
    showingSettings = show;
}

void Manager::SetSwitch( bool state, int index ) 
{ 
    theModel->SetSwitch( state, index ); 
}

void Manager::SetTimeFormat( int format ) 
{
    timeFormat = static_cast< Time >( format ); 
}

//
// Set the default model parameters
//
void
Manager::SetDefaults( const string& filename, bool quickStart ) throw ( IOErr )
{
   
    Defaults d( filename, Defaults::READ );

    globalDataDir = d.GetStringDefault( "globaldatadir" );
    iopDataDir = d.GetStringDefault( "iopdatadir" );
    boundaryDataDir = d.GetStringDefault( "boundarydatadir" );
    userDataDir = d.GetStringDefault( "userdatadir" );
    HistFile = MakeAbsPath( d.GetStringDefault( "histfile" ), userDataDir );
    
    SetLat( d.GetRealDefault( "lat" ) );
    SetLon( d.GetRealDefault( "lon" ) );
    SetStepLen( d.GetIntDefault( "steplen" ) );
    SetEndStep( d.GetIntDefault( "endstep" ) );
    SetTimeFormat( d.GetIntDefault( "timedisplayformat" ) );
    SetSaveFreq( d.GetIntDefault( "savefreq" ) );
    SetShowSettings( d.GetIntDefault( "showsettings" ) );
      //
      // initial conditions file must be set first, or checking will fail
      //
    SetDataset( MakeAbsPath( d.GetStringDefault( "modelfile" ), globalDataDir ), MODEL );
    SetDataset(  MakeAbsPath( d.GetStringDefault( "analysisfile" ), globalDataDir ), ANAL );
    SetDataset( MakeAbsPath( d.GetStringDefault( "iopfile" ), iopDataDir), IOP );
    SetDataset( MakeAbsPath( d.GetStringDefault( "lsminifile" ), boundaryDataDir ), LSMINI );
    SetDataset( MakeAbsPath( d.GetStringDefault( "lsmsurffile" ), boundaryDataDir ), LSMSRF ); 
    SetDataset( MakeAbsPath( d.GetStringDefault( "ozonfile" ), boundaryDataDir ), OZON );
    try {
      SetDataset( MakeAbsPath( d.GetStringDefault( "pressfile" ), boundaryDataDir ), PRES );
    }catch ( DfltsErr& e ) {
      cerr << " Manager::SetDefaults(): No pressfile found in .scam_defaults." << endl;
      stringstream pname ;
      pname << "./data/boundary/press" << theModel->GetNumLevels() << ".nc" ;
      cerr << "Based on number of model levels, will try " << pname.str() << "." << endl;
      SetDataset( MakeAbsPath( pname.str(), boundaryDataDir ), PRES );
    }
    SetDataset( MakeAbsPath( d.GetStringDefault( "sstfile" ), boundaryDataDir ), SST ); 
    SetDataset( MakeAbsPath( d.GetStringDefault( "absemsfile" ), boundaryDataDir ), ABSEMS ); 
    SetDataset( MakeAbsPath( d.GetStringDefault( "aeropticsfile" ), boundaryDataDir ), AEROPTICS ); 
    SetDataset( MakeAbsPath( d.GetStringDefault( "aermassfile" ), boundaryDataDir ), AERMASS ); 
    SetAsciiDataset( MakeAbsPath( d.GetStringDefault( "lsmpftfile" ), boundaryDataDir ), LSMPFT ); 
    SetDataNumLevs( theModel->DatasetRealVar( PRES, "lev" ).size() );
    
      // switches in options dialog
    
    for ( int i=0; i < NUM_SWITCHES; i++ ) {
          // the "<< ends" is necessary because of a bug in
          // the stringstream class (otherwise gets junk at the end)
        stringstream desc;
	//        desc << "switch_desc" << i+1 << ends;
        desc << "switch_desc" << i+1;
        SetSwitchDesc( d.GetStringDefault( desc.str() ), i );
        stringstream name;
        name << "switch" << i+1;
        SetSwitch( bool( d.GetIntDefault( name.str() )), i );
    }
    
    if ( quickStart ) {
        SetBaseDate( d.GetIntDefault( "basedate" ));
        SetBaseSecs( d.GetIntDefault( "basesecs" ));
        SetIopStartOffset( d.GetIntDefault( "iopstartoffset" ));
        SetRunType( d.GetIntDefault( "runtype" ));
        if ( RunType() == USER ) 
            SetDataset( MakeAbsPath(d.GetStringDefault( "userfile" ),
                                    userDataDir ), USER );
        if ( RunType() == SIC ) 
            SetDataset( MakeAbsPath(d.GetStringDefault( "sicfile" ),
                                    globalDataDir ), SIC );
        if ( RunType() == IOP || RunType() == USER )
            SetIopMaxSteps();
          // issue a warning it the latitude and longitude specified
          // in the defaults file don't correspond to the values 
          // specified in the IOP dataset
        if ( RunType() == IOP || RunType() == USER )
	  //            if ( DatasetRealVar( RunType(), "lat" ) != Lat() ||
	  //                 DatasetRealVar( RunType(), "lon" ) != Lon() )
            if ( fabs(DatasetRealVar( RunType(), "lat" ) - Lat()) >.001 ||
                 fabs(DatasetRealVar( RunType(), "lon" ) - Lon()) > .002 )
                ShowMsg( __FILE__, __LINE__, "WARNING!: Latitude and longitude \
specified in \"%s\" (%f,%f),\ndo not match values specified in defaults file \"%s\" (%f,%f)",
                         GetDataset( RunType() ).name().c_str(), 
                         real_t(DatasetRealVar( RunType(), "lat" )),
                         real_t(DatasetRealVar( RunType(), "lon" )), 
                         filename.c_str(), Lat(), Lon() );
    }
}

//
//  Write out the data from the fields that are being saved to the history file
//
void
Manager::SaveFields()
{
      // compute the time in hours since the base date
      //  for writing to the history file.
    int  tsec = StepLen() * CurrentStep();
    double  thour = tsec / 3600.0;

    if ( RunType() == IOP ) {
        thour += iopStartOffset / 3600.0;
        tsec += iopStartOffset;
    }
    int index = CurrentStep() / saveFreq - 1;
    
    try {
        histFile->WriteTime( thour, tsec, index );
        histFile->WriteLevs( theModel->GetLevels(), index );
        histFile->WriteILevs( theModel->GetILevels(), index );
        for (FieldListItr it = theModel->fields.begin(); it != theModel->fields.end(); ++it ) {
            Field* f = (*it).second;
            if ( f->IsSaved() ) {
                f->WriteHistory( index );
            }
        }
    } catch ( NcErr& e ) {
        cerr <<  "ERROR: "__FILE__":" << __LINE__
             << " Manager::SaveFields(): can't write history: " << e << endl;
        exit( -1 );
    }
    if ( SwitchState( CRM_SW ) && CRM.IsCRMTextFile()) {
      histFile->flush();
      CRM.PostProcess( tmpHistFile);
    }
}

//
// Set the default save fields by looking at defaults file.
// This function must be called *after* the model has been 
// initialized or the fieldlist will be empty.
//
void
Manager::SetDefaultSavedFields() throw ( DfltsErr )
{
    if ( NumFields() == 0 ){
        cerr << "ERROR: "__FILE__":" << __LINE__
             << ": Manager::SetDefaultSavedFields() called before field list built!" << endl;
        exit( -1 );
    }

    Defaults d( defaultsFile, Defaults::READ );

    string deflt = d.GetStringDefault( "savefields" );
    size_t pos = 0;                // index 1 past word
    int prev_pos = 0;           // index of beginning of word
    vector< string > savedFields;
    
      // parse the default saved fields string, separating it into separate names
    while (( pos = deflt.find_first_of( ' ', pos )) != string::npos ) {
        savedFields.push_back( deflt.substr( prev_pos, pos-prev_pos ));
        prev_pos = ++pos;
    }
      // the last field may not have a trailing blank...
    //    savedFields.push_back( deflt.substr( prev_pos, pos-prev_pos ));
    
      // set all of the model fields to non-saved state
    FieldList mfields = theModel->fields;
    FieldListItr mfit = mfields.begin(), mfend = mfields.end();

      // iterate over the vector of saved field names, 
      // matching against the model fields
      // and setting the matched fields to saved state
    for ( size_t i=0; i<savedFields.size(); i++ ) {
        if (( mfit = mfields.find( savedFields[i] )) != mfend ) {
            Field* f = (*mfit).second;
            f->SetSave( true );
        }
    }
}


//
//  set the number of time slices left in iop dataset
//
void
Manager::SetIopMaxSteps() throw ( NcErr )
{
    NcIntVar tsec = theModel->DatasetIntVar( RunType(), "tsec"  );
    int ntime = tsec.size();
    iopMaxSteps = (( tsec[ntime-1] - tsec[0] ) - iopStartOffset ) / StepLen();
    
      // if there is only one timepoint in dataset, 
      // assume that it is a starting condition (e.g., USER file )
      // and there should be no limit on steps
    
    if ( ntime == 1 ) 
        iopMaxSteps = INT_MAX; // defined in /usr/include/limits.h

}

void
Manager::SetEndStep( int step )
{
    theModel->SetEndStep( step );
      // notify observers that end step is changed
    SetChanged();
    NotifyObservers();
    ClearChanged();
}

void 
Manager::SetStepLen( int len )
{
    if ( StepLen() != 0 ) {
        SetEndStep( (static_cast<real_t>( StepLen() ) * EndStep()) / len ) ;
    }
    theModel->SetStepLen( len );
    if ( RunType() == IOP || RunType() == USER )
        SetIopMaxSteps(); 
}

//
//  Set one of the model switch descriptions
//
void
Manager::SetSwitchDesc( const string& desc, int n )
{
    switchDesc[n] = desc;
}

//
// Show the current model initialization settings
//
void
Manager::ShowSettings()
{
    stringstream settings;
    
    settings << "Run Type: ";
    switch( RunType() ) {
    case MODEL :
        settings << "Global Model";
        break;
    case ANAL :
        settings << "Global Analysis\n" << "Analysis Dataset: " << GetDataset(ANAL).name() << '\n';
        break;
    case IOP :
        settings << "IOP\n" << "IOP Dataset: " << GetDataset(IOP).name() <<  '\n';
        break;
    case USER :
        settings << "USER\n" << "USER Dataset: " << GetDataset(USER).name() <<  '\n';
        break;
    }
    
    settings << "Global Initial Conditions Data: " << GetDataset(MODEL).name() <<  '\n'
             << "Pressure Level Dataset: " << GetDataset(PRES).name() <<  '\n'
             << "Latitude: " <<  Lat() <<  '\n'
             << "Longitude: " <<  Lon() <<  '\n'
             << "Switches:\n" ;
    for ( int i=0; i<NUM_SWITCHES; i++ ) 
        settings << i+1 << ") " << switchDesc[i] << ": " << SwitchState(i) << endl;
    ShowMsg( __FILE__, __LINE__, settings.str().c_str() );
}


//
//  Call the model step routine
//
void
Manager::StepModel()
{
    if ( CurrentStep() == 0 ) {
        FieldListItr it = theModel->fields.begin(), end = theModel->fields.end();
        for (; it != end; ++it ) {
            Field* f = (*it).second;
            f->SaveRestartData();
        }
    }
    
    theModel->Step();           
    
      // notify observers that model has been stepped
    SetChanged();
    NotifyObservers();
    ClearChanged();
    
      // make the call to SaveFields() last because it will reset
      // the number of outfield calls to zero for averaged fields.
      // (any subsequent calls to Model::ReadField() will return zero's
      //  until another outfld call is made)
    if ( ! ( CurrentStep() % saveFreq ))
        SaveFields();
}


string  Manager::SwitchDesc( int n ) 
{ 
    return switchDesc[n]; 
}

string  Manager::GetUserDataDir() 
{ 
  return userDataDir; 
}

bool Manager::SwitchState( int index ) 
{
    return theModel->SwitchState( index ); 
}

Time
Manager::TimeFormat()    
{
    return timeFormat;
}


