#ifndef _RPCMODEL_H
#define _RPCMODEL_H



#include "realtype.h"
#include "model.h"              // for Model class specification

class Field;
class Manager;
//
// RPCModel inherits from Model, thus keeping the same manager so that
//  they can be used almost interchangeably
//
class RPCModel : public Model
{
public:
    RPCModel(); // constructor
    virtual ~RPCModel();        // destructor
    void BuildFieldList();      // builds the list of fields    
    bool CheckDataset(const string& name, int type); // validates datasets
    bool Connect();             // connects to remote server
    void ReadField( Field& );    // retrieves field data from outfld buffer
    const vector<real_t>& GetLevels();  // retrieves pressure level data from model
    const vector<real_t>& GetILevels();  // retrieves pressure level data from model
    int  Init(bool);                // initializes the model
    void ResetField( const Field& );  // resets averaged fields in outfld buffer
    void WriteField( const Field& );    // sets field data in outfld buffer
    void SetRemoteHost( const string& host ); // set the server hostname
    void SetDataset( int type, const string& filename ) throw ( NcErr );
    void Step();                // time-steps model

    
    void ReadDatasetVar( int type, string& varname, int size, int* data );
    void ReadDatasetVar( int type, string& varname, int size, real_t* data );
    int  DatasetDimSize( int type, string& dimname );
private:
    string    _host;          // name of remote host
};

#endif // _RPCMODEL_H

