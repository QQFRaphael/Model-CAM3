/*------------------------------------------------------------------------*
 * File: defaults.cpp 
 * $Author: jet $
 * $Id: defaults.cpp,v 1.1.6.1 2004/05/18 20:51:52 jet Exp $ *
 *------------------------------------------------------------------------*/

#ifndef lint
static char rcsid[] = "$Id: defaults.cpp,v 1.1.6.1 2004/05/18 20:51:52 jet Exp $";
#endif /* lint */

#include <ctype.h>
#include <fstream>
#include <sstream>
#include <string>
#include "msgdlg.h"
#include "defaults.h"

#ifndef _RWSTD_NO_NAMESPACE
  using namespace std;
#endif

Defaults::Defaults( const string& filename, OpenMode mode ) throw ( DfltsErr ) :
    _filename( filename ), _mode( mode )
{
    if ( _mode == READ ) {
        ifstream in( _filename.c_str() );
        if ( in == 0 ) 
            throw DfltsErr( "Unable to open file for reading", _filename.c_str(), DfltsErr::IO_ERR );
          //
          // read in the defaults file, storing the lines in a vector
          //
        string textline;
        while ( getline( in, textline ) ) {
              // get rid of any preceding whitespace
            while( isspace( textline[0] ) ) 
                textline.erase( 0, 1 ); 
            defaultsText.push_back( textline );
            if ( in.eof() )
                break;
        }
    }
    else {
        if (( output = new ofstream( _filename.c_str() )) == 0 ) 
            throw DfltsErr( "Unable to open file for writing", _filename.c_str(), DfltsErr::IO_ERR );
        *output << "# SCAM Defaults File #\n" ;
          // check for error conditions after the first write
          // such as caused by write permission problems
        if ( ! output->good() ) {
            delete output;
            throw DfltsErr( "Unable to write to file (check permissions)", _filename.c_str(), DfltsErr::IO_ERR );
        }
    }
    
}

Defaults::~Defaults()
{
    if ( _mode == WRITE ) {
        output->close();
        delete output;
    }
}

int
Defaults::GetIntDefault(  const string& defaultName ) throw ( DfltsErr )
{
    int val;
    val = atoi( GetStringDefault( defaultName ).c_str() );
    return val;
}

real_t
Defaults::GetRealDefault(  const string& defaultName ) throw ( DfltsErr )
{
    real_t val;
    val = atof( GetStringDefault( defaultName).c_str() );
    return val;
}

string
Defaults::GetStringDefault(  const string& defaultName ) throw ( DfltsErr )
{
      //
      // look for the target default
      //
    vector< string >::iterator it( defaultsText.begin() );
    
    string value;
    const string delimiter( "\"\'" );
    string::size_type valstart = 0;
    string::size_type valend = 0;

    for ( ; it != defaultsText.end(); ++it ) {
          // get substring of same size as target and compare
        string sub( *it, 0, defaultName.size() );
        
          // see if we have a match
        if ( sub == defaultName ) {
              // now get the substring between the delimiters
            if (( valstart = (*it).find_first_of( delimiter )) 
                == string::npos )  
                continue; // didn't find beginning delimiter
            if (( valend = (*it).find_last_of( delimiter ))
                == string::npos )  
                continue; // didn't find ending delimiter
              // got the string, remove the delimiters
            value = string( *it, valstart+1, (valend-valstart) -1 );
            break;
        }
    }
      // throw exception if we didn't find the default
    if ( value.empty() ) 
        throw DfltsErr( defaultName.c_str(), _filename.c_str(), DfltsErr::MISSING_VAL );

    return value;
}

void
Defaults::WriteDefault(  const string& defaultName, const string& defaultValue ) throw ( DfltsErr )
{
    *output << defaultName << "=\"" << defaultValue << "\"\n";
}

void
Defaults::WriteDefault(  const string& defaultName, real_t defaultValue ) throw ( DfltsErr )
{
    *output << defaultName << "=\"" << defaultValue << "\"\n";
}

void
Defaults::WriteDefault(  const string& defaultName, int defaultValue ) throw ( DfltsErr )
{
    *output << defaultName << "=\"" << defaultValue << "\"\n";
}










