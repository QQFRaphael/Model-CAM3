/*------------------------------------------------------------------------*
 * File: msgdlg.h 
 * $Author: jet $
 * $Id: msgdlg.h,v 1.1.6.1 2004/05/18 20:52:00 jet Exp $ *
 *------------------------------------------------------------------------*/
#ifndef MSGDLG_H
#define MSGDLG_H

#include <string>
#ifndef _RWSTD_NO_NAMESPACE
  using namespace std;
#endif


  //
  //  Displays a modal dialog with the message and
  //   an OK button
  //
  // if the symbol DBUG is defined, also display
  // the filename and line number where the message originated.
  // This is accomplished by using the preprocessor symbols 
  // __FILE__ and __LINE__ which are replaced by cpp with a character
  // string and integer respectively representing the current file and line number
  //

void ShowMsg( const char* filename, int line_number, const char* format ... );

void ShowMsg( const char* filename, int line_number, const string& msg );

// if true, display messages in dialogs, otherwise only to stdout
void SetGraphicalMode ( bool useGUI );  

#endif /* MSGDLG_H */



