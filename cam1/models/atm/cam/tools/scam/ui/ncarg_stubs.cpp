/*------------------------------------------------------------------------*
 * File: ncarg_stubs.cpp 
 * $Author: jet $
 * $Id: ncarg_stubs.cpp,v 1.1.6.1 2004/05/18 20:52:00 jet Exp $ *
 *------------------------------------------------------------------------*/

#ifndef lint
static char rcsid[] = "$Id: ncarg_stubs.cpp,v 1.1.6.1 2004/05/18 20:52:00 jet Exp $";
#endif /* lint */

#ifdef NO_NCARG

#include "ncarg.h" 

extern "C" {
int
CreatePlot( const char* filename, const char* varname,
            const char* longname, const char* units,
            real_t field_multiplier,
            int startIdx, int endIdx, int timeFormat,
            real_t steps_per_hour, int basedate, bool average,
            int outputDeviceID )
{ return 0; }

int
OpenOutputDevice(int,int)
{ return 0; }

void
CloseOutputDevice(int)
{}

void
ClearOutputDevice(int)
{}

void
ActivateOutputDevice(int)
{}

void
DeactivateOutputDevice(int)
{}

void
OpenGKS()
{}

void
CloseGKS()
{}

} // extern "C"

#endif // defined NO_NCARG









