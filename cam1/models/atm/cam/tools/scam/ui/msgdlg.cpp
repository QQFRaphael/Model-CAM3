/*------------------------------------------------------------------------*
 * File: msgdlg.cpp 
 * $Author: jet $
 * $Id: msgdlg.cpp,v 1.1.6.1 2004/05/18 20:52:00 jet Exp $ *
 *------------------------------------------------------------------------*/

#ifndef lint
static char rcsid[] = "$Id: msgdlg.cpp,v 1.1.6.1 2004/05/18 20:52:00 jet Exp $";
#endif /* lint */

#include <string>
#include <sstream>
#include <cstdarg>
#include <cstdio>
#include <iostream>

#ifndef _RWSTD_NO_NAMESPACE
  using namespace std;
#endif

#include <qmsgbox.h>
#include <qfont.h>

#include "msgdlg.h"

  //
  //  Display a modal dialog with the message and
  //   an OK button
  //
  // if the symbol DBUG is defined, also display
  // the filename and line number where the message originated.
  // This is accomplished by using the preprocessor symbols 
  // __FILE__ and __LINE__ which are replaced by cpp with a character
  // string and integer respectively representing the current file and line number
  //
static bool useGUI = TRUE;
void
ShowMsg( string msg )
{
    cerr << "\nSCAM: " << msg << '\n' << endl;;
    if ( useGUI )
      QMessageBox::warning( NULL, "SCAM", msg.c_str());
}

void
ShowMsg( const char* filename, int line_number, const char* format ...  )
{
    va_list args;
    char fmtmsg[2000];

    va_start( args, format );
    vsprintf( fmtmsg, format, args );
    va_end( args );

    stringstream message;
    message << fmtmsg;
#ifndef NDEBUG
    message << "\n (" << filename << ":" << line_number << ")";
#endif
    //    message << ends;
    cerr << "\nSCAM:  " << message.str() << '\n' << endl;;
    
    if ( useGUI )
      QMessageBox::warning( NULL, "SCAM", message.str().c_str());
    //        QMessageBox::warning( NULL, "SCAM", message.str() );
}
 
void
ShowMsg( const char* filename, int line_number, string msg )
{
    stringstream message;
    message << msg;
#ifndef NDEBUG
    message << "\n (" << filename << ":" << line_number << ")" ;
#endif
    //    message << ends;
    cerr << message << '\n' << endl;
    if ( useGUI )
        QMessageBox::warning( NULL, "SCAM", message.str().c_str() );
}

void
ShowMsg( const char* format ...  )
{
    va_list args;
    char message[2000];

    va_start( args, format );
    vsprintf( message, format, args );
    va_end( args );

    cerr << "\nSCAM:  "  << message << '\n' << endl;

    if ( useGUI )
        QMessageBox::warning( NULL, "SCAM", message );
}
 
void
SetGraphicalMode( bool _useGUI )   
{
    useGUI = _useGUI;
}




