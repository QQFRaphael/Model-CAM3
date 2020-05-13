/*------------------------------------------------------------------------*
 * File: manager.h 
 * $Author: jmccaa $
 * $Id: manager.h,v 1.1.6.2 2004/08/19 15:05:50 jmccaa Exp $ *
 *------------------------------------------------------------------------*/
#ifndef manager_h
#define manager_h

#include <string>
#include <vector>

#include "realtype.h"
#include "defaults.h"
#include "max.h"
#include "model.h"
#include "timeconvert.h"
#include "history.h"
#include "field.h"
#include "ncfile.h"



// convenience symbol to access the static manager instance
#define MANAGER Manager::Instance()

class Manager : public Observable
{
public:
                 ~Manager();
    void          SetDefaultsFile( const string& fileName );
    string        GetDefaultsFile();
    int           BaseDate();
    const vector<real_t>& BaseLevels();         
    const vector<real_t>& BaseILevels();         
    int           BaseSecs();    
    void          CreateModel( ModelType type );
    int           CurrentStep();
    NcVariable<real_t>    DatasetRealVar( int datasetType, const string& varname ) throw (NcErr);
    NcCharAtt DatasetCharAtt( int datasetType, const string& attname  ) throw (NcErr);
    NcIntVar      DatasetIntVar( int datasetType, const string& varname ) throw (NcErr);
    int           EndStep();
    Field*        FindField( const string& name );
    const NcFile& GetDataset( int type );
    const Ascii_Dataset& GetAsciiDataset( int type );
    const FieldList& GetFieldList();
    string        HistName();
    static Manager& Instance();
    int           IopStartOffset();
    void          InitHistoryFile() throw (NcErr);
    bool          InitModel();
    bool          IsModelInited();
    bool          IsShowingSettings();
    bool          IsValidDataset( const string& name, int type );
    real_t         Lat();
    void          LoadDefaults( const string& filename, bool isQuickstart=false ) throw( IOErr );
    void          LoadInitialConditionsFile(const string& filename) throw (NcErr);
    real_t         Lon();
    int           MaxStep();
    int           NumFields();
    int           NumLevs();
    int           NumILevs();
    void          ReadField( Field& f );
    void          ResetField( const Field& f );
    void          RunModel();
    int           RunType();
    int           SaveFrequency();
    void          SaveDefaults( const string& filename, bool quickStart = false ) throw ( DfltsErr );
    void          SaveHistoryFile( const string& filename ) throw( IOErr );
    void          SaveInitialConditions( const string& filename ) throw( NcErr );
    void          SetBaseDate( int bdate );
    void          SetBaseSecs( int bsecs );
    void          SetDataset( const string& name, int type ) throw( NcErr );
    void          SetAsciiDataset( const string& name, int type );
    void          SetDataNumLevs(int nlevs);
    void          SetDataNumILevs(int plevs);
    void          SetDefaults( const string& filename, bool quickStart ) throw ( IOErr );
    void          SetEndStep( int step );
    void          SetIopMaxSteps() throw ( NcErr );
    void          SetIopStartOffset( int secs );
    void          SetLat( real_t lat );
    void          SetLon( real_t lon );
    void          SetRemoteHost( const string& host );
    void          SetRunType( int type );
    void          SetSaveFreq( int newFreq );
    void          SetSeedVal( int val ) { seedVal = val; }
    void          SetDefaultSavedFields() throw ( DfltsErr );
    void          SetShowSettings( bool show );
    void          SetStepLen( int len );
    void          SetSwitch( bool state, int index );
    void          SetSwitchDesc( const string& desc, int n );
    void          SetTimeFormat( int format );
    void          ShowSettings();
    int           StepLen();
    int           SeedVal() { return seedVal; }
    void          StepModel();
    void          StopModel();
    string        SwitchDesc( int n );
    bool          SwitchState( int index );
    Time          TimeFormat();
    void          WriteField( const Field& f );
    string        lsmpftfile;
    string        GetUserDataDir();
private:                        

      // private methods
                  Manager();
    void          SaveFields(); // write the saved fields to the history file

      // private fields
    static Manager*  _instance; // global static instance

    bool          modelInited;  // is model initialized
    int           saveFreq;     // freq of saves
    bool          showingSettings; // show settings at startup?
    string        tmpHistFile;  // temp history file
    string        HistFile;     // default history file name
    Time          timeFormat;   // format for displaying elapsed time

      // dataset info
    int           numDataLevs;  // levels in init dataset
    int           numDataILevs;  //interface levels in init dataset
    string        boundaryDataDir;
    string        globalDataDir;
    string        iopDataDir;
    string        userDataDir;

      // model state variables
    int           iopNumTime;   // time slices in iop dataset
    int           iopMaxSteps;  // steps available in iop dataset
    int           iopTimeInt;   // interval of iop dataset observations
    int           currentStep;  // current model step
    int           iopStartOffset; // offset from start of iop dataset
    string        switchDesc[NUM_SWITCHES];

    int           seedVal;   // time slices in iop dataset
    
    Model*        theModel;
    History*      histFile;
    string        defaultsFile; // startup defaults file
};

inline const vector<real_t>& Manager::BaseLevels() { return theModel->BaseLevels(); }
inline const vector<real_t>& Manager::BaseILevels() { return theModel->BaseILevels(); }
 
inline int  Manager::CurrentStep() { return theModel->CurrentStep(); }

inline int  Manager::EndStep() { return theModel->EndStep(); }

inline Manager& Manager::Instance() { 
    if ( _instance == NULL )  _instance = new Manager();
    return *_instance; }

inline bool Manager::IsModelInited() { return modelInited; }

inline int  Manager::NumFields() { return theModel->NumFields(); }

inline int  Manager::NumLevs() { return theModel->NumLevs(); }
inline int  Manager::NumILevs() { return theModel->NumILevs(); }

inline void Manager::ReadField( Field& f ) { theModel->ReadField( f ); }
inline void Manager::ResetField( const Field& f ) { theModel->ResetField( f ); }
inline int  Manager::RunType() { return theModel->RunType(); }

inline int  Manager::StepLen() { return theModel->StepLen(); }

inline void Manager::WriteField( const Field& f ) { theModel->WriteField( f ); }
#define DEFAULTS_FILE (".scam_defaults")


#endif /* manager_h */


    
