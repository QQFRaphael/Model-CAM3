/*------------------------------------------------------------------------*
 * File: utils.cpp 
 * $Author: jet $
 * $Id: utils.cpp,v 1.1.6.1 2004/05/18 20:52:03 jet Exp $ *
 *------------------------------------------------------------------------*/

#ifndef lint
static char rcsid[] = "$Id: utils.cpp,v 1.1.6.1 2004/05/18 20:52:03 jet Exp $";
#endif /* lint */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "max.h"

const real_t LOGBASE = 2.7182818284590452354;

string
MakeAbsPath( const string& inputfile, const string& inputdir )
{
    string abspath;
      // if the inputfile begins with a '.' or '/' assume it's a path
    if ( inputfile[0]=='/' || inputfile[0]=='.' )
         return  inputfile;

    abspath = inputdir;
    abspath += inputfile;

    return abspath;
}

string
FileFromPath( const string& pathname )
{
    char  filename[MAX_PATH_LEN];
    char path[MAX_PATH_LEN];
    
    char* ptr;
    strcpy( path, pathname.c_str() );
    
    if ( strtok( path, "/" ) == NULL )
        return pathname;
    while ( ( ptr = strtok( NULL, "/" ) ) != NULL ) {
        strcpy( filename, ptr );
     }

    return filename;
}

real_t 
Nicenum( real_t x, bool round )
{
    real_t exp;
    real_t f, nicex;
    bool neg = false;
    
    if ( x < 0 ) {
        x = -x;
        neg = true;
    }
    
    if ( x == 0 )
        exp = 0;
    else
        exp = floor( log10( x ) );
    
    f = RoundNum( x / pow(10., exp ) );	/* integer between 1 and 9 */
    
    if ( round ) {
        if ( f < 2 )
            f = 1;
        else if ( f < 4 )
            f = 2;
        else if ( f < 8 )
            f = 5;
        else
            f = 10;
    }
    
    nicex =  f * pow( 10., exp );

    if ( neg )
        nicex = -nicex;
    
    return nicex;
}

real_t 
RoundNum( real_t x )
{
    real_t c, f;
    
    c = ceil( x );
    f = floor( x );
    
    if ( fabs(x-c) < fabs(x-f) )
        return c;
    else
        return f;
}

real_t
PressToHeight( real_t millibars )
{
    if( millibars <= 1 )
        millibars = 1;
    return fabs( -7 * log( millibars/1000 ));
}
    
real_t
HeightToPress( real_t height )
{
    return( 1000 * pow( LOGBASE, height/(-7)));
}





