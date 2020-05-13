/*------------------------------------------------------------------------*
 * File: defaults.h 
 * $Author: jmccaa $
 * $Id: defaults.h,v 1.1.6.2 2004/08/19 15:05:46 jmccaa Exp $ *
 *------------------------------------------------------------------------*/

#ifndef DEFAULTS_H
#define DEFAULTS_H

#include <string>
#include <vector>
#include <fstream>
#include "realtype.h"
#include "max.h"
#include "ioerr.h"
#define DEFAULTS_FILE (".scam_defaults")

class Defaults
{
public:
    enum OpenMode { READ, WRITE };

    Defaults( const string& fileName, OpenMode mode = READ ) throw ( DfltsErr );
    ~Defaults(); 
    
    string  GetStringDefault(  const string& defaultname ) throw ( DfltsErr );
    int     GetIntDefault(  const string& defaultname ) throw ( DfltsErr );
    real_t   GetRealDefault(  const string& defaultname ) throw ( DfltsErr );

    void  WriteDefault( const string& defaultname, int    defaultvalue ) throw ( DfltsErr );
    void  WriteDefault( const string& defaultname, real_t  defaultvalue ) throw ( DfltsErr );
    void  WriteDefault( const string& defaultname, const string&  defaultvalue ) throw ( DfltsErr );

private:

    string  _filename;
    OpenMode _mode;
    vector< string > defaultsText;
    ofstream *output;
};



#endif /* DEFAULTS_H */
 

