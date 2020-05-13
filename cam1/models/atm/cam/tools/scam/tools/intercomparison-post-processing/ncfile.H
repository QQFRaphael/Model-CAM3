#ifndef NCFILE_H
#define NCFILE_H

#include <netcdf.h>

class Field;
#define STR_LEN 256

class NCFile {
public:
    enum { READ, WRITE, CREATE };
          //: The three possible open modes for NCFiles


    NCFile( const char* filename, bool NCerrorHandling = false );
          //: constructor
          //
    virtual ~NCFile();
          //: destructor
          //
    size_t GetDimSize( const char* dimName ); 
          //: returns size of dimension
          //
    bool GetGlobalAttribute( const char* attName, double* value ); 
          //: returns double attribute value
          //
    bool GetAttribute( const char* varname, const char* attName, char* attValue ); 
          //: returns varname's text attribute "attName" value
          //
    bool GetAttribute( const char* varname, const char* attName, float* attValue ); 
          //: returns varname's float attribute "attName" value
          //
    void HandleError( int errorStatus );  
          //: outputs netcdf error message
          //
    bool Open( int mode = READ, bool clobber = false );
          //: opens a netcdf file
          //
    const char* Name() { return filename; }
          //: returns name of netcdf file
          //
    bool Read( const char* name, Field& f );
          //: read entire float variable
          //
    bool Read( const char* varName, int size, char* outText,
                        int index = 0 );
          //: read text variable
          //
    bool Read( const char* varName, int size, double* outData,
                        int index = 0 );
          //: read double variable
          //
    bool Read( const char* varName, int size, float* outData,
                        int index = 0 );
          //: read float variable
          //
    bool Read( const char* varName, int size, int* outData,
                        int index = 0 );
          //: read int variable
          //
    bool Read( const char* varName, float* outData );
          //: read entire float variable
          //
    bool Read( const char* varName, int* outData );
          //: read entire int variable
          //
    bool VarExists( const char* varName );  
          //: returns true if variable is in netcdf file
          //

protected:
    
    char* filename;		// name of opened netcdf file
    bool handleErrs;            // if true, show info about errors
    int  errStatus;             // error status returned from netcdf calls
    int  NCID;                  // id of netcdf file opened
    int  timeDimID;             // time dimension ID
    int  latDimID;              // latitude dimension ID
    int  lonDimID;              // longitude dimension ID
    int  levDimID;              // level dimension ID
    int  stringDimID;           // string dimension ID

    bool          Read( const char* varName, int size, void* outData,
                        int index, nc_type dataType );
          //: read variable of type dataType
          //
    bool          Write( const char* varName, int size, const void* inData,
                         int index, nc_type dataType );
          //: write variable of type dataType
};

#endif /* NCFILE_H */
 
