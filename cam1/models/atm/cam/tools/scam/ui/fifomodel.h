#ifndef _FIFOMODEL_H
#define _FIFOMODEL_H


#include "realtype.h"
#include "ipc.h"                // interprocess communication structs
#include "model.h"              // for Model class specification

class Field;

//
// FifoModel inherits from Model: only use the interface 
//  provided by Model in the public interface.
//
//
class FifoModel : public Model
{
public:
             FifoModel();       // constructor
    virtual ~FifoModel();       // destructor
    void BuildFieldList();      // construct the list of fields
    void ReadField( Field& );   // retrieves field data from outfld buffer
    const vector<real_t>& GetLevels();  // retrieves pressure level data from model
    int GetNumLevels();  // retrieves pressure level data from model
    const vector<real_t>& GetILevels();  // retrieves pressure level data from model
    int  Init(bool);                // initializes the model

    void ResetField( const Field& );  // resets averaged fields in outfld buffer
    void WriteField( const Field& );    // sets field data in outfld buffer
    void SetRemoteHost( const string& host ) { exit(-1); }  // should never be called
    void Step();                // time-steps model
    void WriteReq();
    void ReadReply();

private:
    srv_req          req;
    srv_reply        reply;
};


#endif // _FIFOMODEL_H






