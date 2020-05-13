/*------------------------------------------------------------------------*
 * File: sicfile.h 
 * $Author: jet $
 * $Id: sicfile.h,v 1.1.6.1 2004/05/18 20:52:03 jet Exp $
 *------------------------------------------------------------------------*/

#ifndef SICFILE_H
#define SICFILE_H

#include "realtype.h"
#include "ncfile.h"

class Sicfile: public NcFile
{
public:
    Sicfile( const string& name, OpenMode mode, bool clobber );
    ~Sicfile();
};

#endif // SICFILE_H

