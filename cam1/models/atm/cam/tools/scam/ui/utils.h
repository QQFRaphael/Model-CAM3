/*------------------------------------------------------------------------*
 * File: utils.h 
 * $Author: jet $
 * $Id: utils.h,v 1.1.6.1 2004/05/18 20:52:03 jet Exp $ *
 *------------------------------------------------------------------------*/
#ifndef UTILS_H
#define UTILS_H

#include <string>
#include "realtype.h"

#ifndef _RWSTD_NO_NAMESPACE
  using namespace std;
#endif

string
MakeAbsPath( const string& inputfile, const string& inputdir );

string
FileFromPath( const string& pathname );

// returns closest number of same
// magnitude that begins with 1, 2, 5, or 10. 
real_t 
Nicenum( real_t x, bool round );

// return rounded number of same magnitude
real_t 
RoundNum( real_t x ); 

void
FindMinMax( const real_t *data, int size, real_t* min,real_t* max);

real_t
PressToHeight( real_t millibars );
    
real_t
HeightToPress( real_t height );    

#endif /* UTILS_H */




