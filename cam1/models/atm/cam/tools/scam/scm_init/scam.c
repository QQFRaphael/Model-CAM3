/*
 * $Id: scam.c,v 1.1.6.2 2004/10/18 17:40:40 jmccaa Exp $
 *  main() routine for fifo based client-server version
 */

/*#include <rpc/rpc.h>*/
#include <sys/time.h>

#include "fortran.h"
#include "scam_fifo.h"
#ifndef HP
/*#include "scam_rpc.h"*/
#endif

#ifndef lint
static char rcsid[] = "$Id: scam.c,v 1.1.6.2 2004/10/18 17:40:40 jmccaa Exp $";
#endif 

#ifndef DEBUG_INIT
int
main( int argc, char* argv[] )
{
    server_init( argc, argv );
    return server_run();
}
    
#endif









