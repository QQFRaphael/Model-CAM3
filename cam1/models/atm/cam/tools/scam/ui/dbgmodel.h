#ifndef _DBGMODEL_H
#define _DBGMODEL_H

#include "realtype.h"
#include "model.h"              // for Model class specification

class Field;

//
// DbgModel is used to debug the user interface
//
//
class DbgModel : public Model
{
public:
    DbgModel();                 // constructor
    virtual ~DbgModel();        // destructor
    void BuildFieldList();      // construct the list of fields
    bool CheckDataset( const string&, int );
    void ReadField( Field& );    // retrieves field data from outfld buffer
    const vector<real_t>&  GetLevels();  // retrieves pressure level data from model
    int GetNumLevels();  // retrieves pressure level data from model
    const vector<real_t>&  GetILevels();  // retrieves pressure level data from model
    int  Init(bool);                // initializes the model
    void ResetField( const Field& );  // resets averaged fields in outfld buffer
    void WriteField( const Field& );    // sets field data in outfld buffer
    void Step();                // time-steps model

private:
    vector<real_t> _dbgdata;
};

#endif // _DBGMODEL_H



