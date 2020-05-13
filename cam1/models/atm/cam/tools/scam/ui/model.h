#ifndef _MODEL_H
#define _MODEL_H

#include <map>
#include <string>
#include <vector>

#include "realtype.h"
#include "ncfile.h"
#include "ascii_dataset.h"
#include "max.h"
#include "runtype.h"
#include "field.h"

using ncfile::NcVariable;
using ncfile::NcIntVar;
using ncfile::NcFile;
using ncfile::NcCharAtt;

class Model
{
public:    
    enum ModelType{ FIFO, RPC, DBG };
      // constructor
    Model();
      // destructor
    virtual ~Model();
      // advance the model a single time step
    virtual void          Step() = 0;
      // initialize field list (ignore duplicates)
    virtual void          BuildFieldList() = 0;
      // check whether dataset is valid without setting it
    virtual bool          CheckDataset( const string& name, int type );
      // retrieve model dataset  variable
    virtual NcIntVar      DatasetIntVar( int datasetType, const string& varname  ) throw (NcErr);
      // retrieve model dataset  variable
    virtual NcVariable<real_t>    DatasetRealVar( int datasetType, const string& varname  ) throw (NcErr);
      // retrieve model dataset  attribute
    virtual NcCharAtt DatasetCharAtt( int datasetType, const string& attname  ) throw (NcErr);
      // retrieve model pressure levels
    virtual const vector<real_t>&  GetLevels() = 0;
      // retrieve model pressure levels
    virtual int           GetNumLevels() = 0;
      // retrieve model pressure levels
    virtual const vector<real_t>&  GetILevels() = 0;
      // initialize the actual model
    virtual int           Init(bool) = 0;
      // retrieve data for a field
    virtual void          ReadField( Field& ) = 0;
      // resets averaged fields in the model outfld buffer
    virtual void          ResetField( const Field& ) = 0;
      // set a model dataset
    virtual void          SetDataset( int type, const string& filename ) throw (NcErr);
    virtual void          SetAsciiDataset( int type, const string& filename );
      // set fields data in the model
    virtual void          WriteField( const Field& ) = 0;
    
    
    int           CurrentStep() { return currentStep; }
    int           EndStep() { return endStep; }
    const vector<real_t>& BaseLevels() { return baseLevs; }
    const vector<real_t>& BaseILevels() { return baseILevs; }
    int           BaseDate() { return baseDate; }
    int           BaseSecs() { return baseSecs; }
    const NcFile* GetDataset( int type );
    const Ascii_Dataset* GetAsciiDataset( int type );
    int           NumLevs() { return numLevs; }
    int           NumILevs() { return numILevs; }
    int           NumFields() { return fields.size(); }
    int           RunType() { return runType; }
    int           StepLen() { return stepLen; }
    real_t         Lat() { return latitude; }
    real_t         Lon() { return longitude; }
    bool          SwitchState( int index );
    void          SetBaseDate( int bdate ) { baseDate = bdate; }
    void          SetBaseSecs( int bsecs ) { baseSecs = bsecs; }
    void          SetEndStep( int step ) { endStep = step; }
    void          SetLat( real_t latitude );
    void          SetLon( real_t longitude );
    void          SetRunType( int type );
    void          SetStepLen( int seconds ) { stepLen = seconds; }
    void          SetSwitch( bool state, int index );

    FieldList  fields;          // the collection of fields corresponding 
                                // to those created in bldfld.F
    string lsmpftfile;

protected:
    
    int           numLevs;      // number of model pressure levels
    int           numILevs;      // number of model interface pressure levels
    vector<real_t> baseLevs;     // base model pressure levels (at initialization)
    vector<real_t> baseILevs;     // base model pressure levels (at initialization)
    vector<real_t> currLevs;     // current model pressure levels
    vector<real_t> currILevs;     // current model interface pressure levels
    int           endStep;      // ending step
    int           currentStep;  // current step
    int           initLatIdx, initLonIdx; // the lat and lon specified in the 
                                          // initial conditions dataset
      // parameters passed to the fortran model for initialization

    int           switches[NUM_SWITCHES];// model logical switches
    int           runType;      // model run type ( MODEL, IOP, ANALYSIS, etc)
    real_t         latitude;     // latitude of column
    real_t         longitude;    // longitude of column
    int           baseDate;     // date at start of model run (gmt)
    int           baseSecs;     // offset in seconds from base date
    int           initStatus;   // status returned from model upon initialization
    int           stepLen;      // timestep length in seconds
    int           restart;      // restart logical switch
    NcFile*       dataset[NUM_INIT_FILES]; // model initialization datasets
    Ascii_Dataset* adataset[NUM_INIT_FILES]; // model ascii initialization datasets
};

typedef Model::ModelType ModelType;


#endif // MODEL_H















