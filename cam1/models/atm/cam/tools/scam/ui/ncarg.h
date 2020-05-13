/*------------------------------------------------------------------------*
 * File: ncarg.h 
 * $Author: jet $
 * $Id: ncarg.h,v 1.1.6.1 2004/05/18 20:52:00 jet Exp $ *
 *------------------------------------------------------------------------*/
#ifndef NCARG_H
#define NCARG_H

#include "realtype.h"
#include "max.h"

extern "C" {
int
CreatePlot( const char* filename, const char* varname,
            const char* longname, const char* units,
            real_t field_multiplier,
            int startIdx, int endIdx, int timeFormat,
            real_t steps_per_hour, int basedate, bool average,
            int outputDeviceID )
;

int
OpenOutputDevice(int,int)
;

void
CloseOutputDevice(int)
;

void
ClearOutputDevice(int)
;

void
ActivateOutputDevice(int)
;

void
DeactivateOutputDevice(int)
;

void
OpenGKS()
;

void
CloseGKS()
;

} // extern "C"
 
#define POSTSCRIPT 20
#define CGM 1
#define XWIND 8
#define MAX_GKS_XWIN 14
#define OUTPUT_FILE_ID MAX_FIELDS + 1 // can't conflict with field id

#endif // NCARG_H










