#ifndef IOEXCEPT_H
#define IOEXCEPT_H

// I/O exception classes

#include <sstream>
#include <string>
#include <cstdio>
#include <netcdf.h>
 
#ifndef _RWSTD_NO_NAMESPACE
  using namespace std;
#endif

class IOErr {
public:
    IOErr( const string& error, const string& filename, 
           const string& srcfile = "", int srcline = 0 ) 
        : _error( error ), _filename( filename ), _srcfile( srcfile ), _srcline( srcline )  {}
    virtual const string& toString(); // returns description of the exception
    string   srcInfo();         // returns information about where the exception occurred 
                                //  or empty string if no info available
    virtual ~IOErr() {}

protected:
    string     _error;          // description of the error 
    string     _filename;       // name of file being processed when exception was generated
    string     _srcfile;        // name of source file in which exception was generated
    int        _srcline;        // line in source file where exception was generated
    string     _errmsg;         // used to generate error message in toString() method
};

class DfltsErr: public IOErr {
public:
    enum ErrType { IO_ERR, MISSING_VAL };
    DfltsErr( const string& error, const string& filename, ErrType type,
              const string& srcfile = "", int srcline = 0 )  
        : IOErr( error, filename, srcfile, srcline ), _type( type ) {}

    const string& toString();   

protected:    
    ErrType _type;
};

class NcErr: public IOErr {
public:
    NcErr( const string& filename, const string& varname, const string& error, 
           int ncid,  const string& srcfile = "", int srcline = 0 ) 
        : IOErr( error, filename ), _varname( varname ), _ncid( ncid ) {}

    NcErr( const string& filename, const string& error, 
           int ncid, const string& srcfile = "", int srcline = 0  ) 
        : IOErr( error, filename, srcfile, srcline ), _varname( "GLOBAL" ), _ncid( ncid ) {} 

    const string& toString();

protected:
    string _varname;
    int   _ncid;
};

inline string
IOErr::srcInfo() {
#ifdef NDEBUG
    return string("");
#endif
    if ( _srcfile.empty() ) 
        return string("");  // we don't have any information; return empty string
    stringstream info;
    info << "( exception occurred at" << _srcfile << ":" << _srcline;
    return info.str(); 
}

inline const string& 
IOErr::toString() 
{
    _errmsg = "I/O Error: in file \"" + _filename + "\": " + _error + srcInfo();
    return _errmsg;
}

inline ostream& operator<<( ostream& o, IOErr& e ) 
{
    o << e.toString();
    return o;
}

inline const string& 
DfltsErr::toString() 
{
    switch ( _type ) {
    case DfltsErr::MISSING_VAL:
        _errmsg =  "SCAM Defaults Error: in defaults file \"" + _filename + "\": value for \"" + _error
            + "\" not found." + srcInfo();
            break;
    case DfltsErr::IO_ERR:
        _errmsg =  "I/O Error: while processing defaults file \"" + _filename + "\": " + _error + srcInfo();
        break;
    }
    return _errmsg;
}

inline ostream& operator<<( ostream& o, DfltsErr& e )
{ 
    o << e.toString();
    return o;
}

inline const string&
NcErr::toString()
{
    if ( _varname == "GLOBAL" ) {
        _errmsg =  "NetCDF Error: in file \"" + _filename + "\": " + _error + srcInfo();
    } else {
        _errmsg = "NetCDF Error: in file \"" + _filename + "\", variable \"" + _varname 
          + "\": " + _error + srcInfo();
    }
    return _errmsg;
}

inline ostream& operator<<( ostream& o, NcErr& e ) 
     {
    o << e.toString();
    return o;
}

#endif // IOEXCEPT_H



