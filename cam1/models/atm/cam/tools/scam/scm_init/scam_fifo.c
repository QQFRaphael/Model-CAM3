#ifndef lint
static char rcsid[] = "$Id: scam_fifo.c,v 1.1.6.2 2004/08/19 15:05:42 jmccaa Exp $";
#endif
#define _POSIX_C_SOURCE   1      /* needed on HP */
#include <errno.h>
#include <signal.h>             /* signals */
#include <stdio.h>              /* printf, et al. */
#include <stdlib.h>             /* exit */
#include <string.h>             /* strerror */
#include <sys/types.h>          /* various type defs */
#include <sys/stat.h>           /* mkfifo */
#include <sys/wait.h>           /* wait() */
#include <sys/signal.h>  
#include <unistd.h>

#include <fcntl.h>              /* file io */

#include "fortran.h"            /* for FORTRAN mangled function names */
#include "c_outfld.h"             /* function protos for get_field, set_field, etc */
#include "ipc.h"                /* structs used for interprocess communication */
#include "max.h"

#if (defined USE_4BYTE_REAL )
typedef float   real;
#else
typedef double   real;
#endif

/* globals */
static char      req_fifo[MAX_PATH_LEN];
static int       req_fd;   /*  file descriptor */
static char      reply_fifo[MAX_PATH_LEN];
static int       reply_fd;  /*   file descriptor */
static int       gui_pid; 
static int       sigterm_status;

/*--------------------------------------------------------
 * cleanup: cleans up fifos, exits
 *------------------------------------------------------*/
static void
cleanup()
{
    /*
     *  don't bother to check success on these calls because
     *  the fifos may have already been removed by the gui.
     */
    close( req_fd );
    close( reply_fd );
    remove( req_fifo );
    remove( reply_fifo );
}


/*--------------------------------------------------------
 * handle_interrupt: handles interrupts to exit gracefully
 *             input: number of the signal caught - unused
 *--------------------------------------------------------*/
static void
handle_interrupt( int sig )
{
    int child_status;
    
    if ( sig == SIGTERM ) {
      sigterm_status=1;
      return;
    }
    if ( sig == SIGCHLD ) {
      waitpid( gui_pid, &child_status, WNOHANG );
      if ( !WIFSIGNALED( child_status ) || !WIFEXITED( child_status )) 
	return;             /* gui only suspended, i.e., with ^z */
    }
    cleanup();
    kill( gui_pid, SIGTERM );
    exit( 0 );
}


/*-----------------------------------------------------------------------------
 * init_model_svc: calls init_model with pararameter passed from client,
 *                     sets a timer.
 * input: pointer to request buffer, reply buffer(unused)
 * returns: status 
 *----------------------------------------------------------------------------*/
int
init_model_svc( void* req_buf, void* reply_buf )
{
    int status;
    init_data* init;
    real latitude, longitude;
    
    init = req_buf;
    
    /*
     * call the model initialization function (FORTRAN)
     */
                                /* init_model expects double ??? pointers for lat,lon */
    latitude = init->lat;
    longitude = init->lon;
    init_model( init->sw,   &latitude, &longitude,
                &init->bd,  &init->bs,  &init->rt,
                &init->sl,  &init->rst, &init->ie, &status , init->ana, 
                init->ini,init->iop, init->lsmini,init->lsmsurf ,
                init->ozo,init->prs, init->sic, 
                init->sst, init->abs,init->optics,init->mass, 
		init->lsmpft, init->usr,
                (int)strlen( init->ana ),
                (int)strlen( init->ini ),
                (int)strlen( init->iop ),   
                (int)strlen( init->lsmini ), 
                (int)strlen( init->lsmsurf ), 
                (int)strlen( init->ozo ),
                (int)strlen( init->prs ),
                (int)strlen( init->sic ),
                (int)strlen( init->sst ),
                (int)strlen( init->abs ),
                (int)strlen( init->optics ),
                (int)strlen( init->mass ),
                (int)strlen( init->lsmpft ),
                (int)strlen( init->usr ) );
    return status;
}

/*-----------------------------------------------------------------------------
 * step_svc: calls model time-stepping routine
 * input: pointer to request buffer(unused), reply buffer (unused)
 *----------------------------------------------------------------------------*/
int
step_svc( void* req_buf, void* reply_buf )
{
    stepon();
    
      /* currently stepon doesn't return any status */
    return 0;
}

/*-----------------------------------------------------------------------------
 * get_field_svc: retrieves field data for a layer field
 * input: pointer to request buffer, reply buffer
 * returns: status 
 *----------------------------------------------------------------------------*/
int
get_field_svc( void* req_buf, void* reply_buf )
{
    int* fid;                   /* field id */
    fid = (int*) req_buf;
    get_field( *fid, (real_t*)reply_buf );
    return 0;
}

/*-----------------------------------------------------------------------------
 * reset_field_svc: resets averaged field data
 * input: pointer to request buffer, reply buffer(unused)
 * returns: status 
 *----------------------------------------------------------------------------*/
int
reset_field_svc( void* req_buf, void* reply_buf )
{
    int* fid;                    /* field id */

    fid = (int*) req_buf;
    
    /*
     * call reset_field  to reset the field's averaging mechanism
     */
    reset_field( *fid ); 

    return 0;
}


/*-----------------------------------------------------------------------------
 * set_field_svc: sets field data
 * input: pointer to request buffer, reply buffer(unused)
 * returns: status 
 *----------------------------------------------------------------------------*/
int
set_field_svc( void* req_buf, void* reply_buf )
{
    field_data* field;
    field = (field_data*) req_buf;
    /*
     * call set_field  to retrieve field's values
     */
    set_field( field->id, field->data );
    return 0;
}

/*-----------------------------------------------------------------------------
 * init_field_svc: initializes a field in the outfld buffer
 * input: pointer to request buffer, reply buffer
 * returns: status 
 *----------------------------------------------------------------------------*/
int
get_field_info_svc( void* req_buf, void* reply_buf )
{
    int* id;
    field_info* fi;

    id = (int*) req_buf;
    fi = (field_info*) reply_buf;

    return get_field_info( *id, fi ); 
}

/*-----------------------------------------------------------------------------
 * get_levels_svc: retrieves model pressure levels
 * input: pointer to request buffer, reply buffer
 * returns: status 
 *----------------------------------------------------------------------------*/
int
get_levels_svc( void* req_buf, void* reply_buf )
{
    real model_levels[MAX_LEVELS]; /* levels filled in by  get_levels(FORTRAN) */
    real model_ilevels[MAX_LEVELS]; /* levels filled in by  get_levels(FORTRAN) */
    int nlev;                   /* number of levels */
    int ilev;                   /* number of interface levels */
    int i;
    level_data* lev;
    int* number_only;

    number_only = (int*) req_buf;
    lev = (level_data*) reply_buf;
    /*
     * call get_levels (FORTRAN)  to retrieve the pressure level values
     */

    get_levels( &nlev, model_levels,&ilev, model_ilevels, number_only ); 
    
      /*
       * fill in the reply struct 
       */
    lev->nlev = nlev;
    lev->ilev = ilev;
    for ( i=0; i<nlev; i++ )
        lev->levels[i] = model_levels[i];
    for ( i=0; i<ilev; i++ )
        lev->ilevels[i] = model_ilevels[i];
    
    return 0;                   /* always succeeds ! */
} 

/*------------------------------------------------------------------------
 * server_init: starts the user-interface process running, sets up signal
 *              handling, creates the FIFOs
 *------------------------------------------------------------------------*/
void
server_init( int argc, char* argv[] )
{
    struct sigaction in;  
    
      /* register the cleanup function to be called at exit */
    if ( atexit( cleanup ) != 0 )
        fprintf(stderr, "Warning: can't register exit function \"cleanup()\n" );
    
    
    
      /*
       *    Set up a handler for catching signal interrupts
       */
    in.sa_handler = handle_interrupt;
    in.sa_flags = 0;
    if ( sigemptyset( &in.sa_mask ) == -1 ) {
        fprintf( stderr, "Error: cannot clear interrupt handler: %s\n",
                 strerror( errno ));
        exit( -1 );
    }
    if ( sigaction( SIGCHLD, &in, NULL ) == -1 ) {
        fprintf( stderr, "Warning: cannot install interrupt handler: %s\n",
                 strerror( errno ));
    }
    if ( sigaction( SIGPIPE, &in, NULL ) == -1 ) {
        fprintf( stderr, "Warning: cannot install interrupt handler: %s\n",
                 strerror( errno ));
    }
    if ( sigaction( SIGINT, &in, NULL ) == -1 ) {
        fprintf( stderr, "Warning: cannot install timeout handler: %s\n",
                 strerror( errno ));
    }
    if ( sigaction( SIGTERM, &in, NULL ) == -1 ) {
        fprintf( stderr, "Warning: cannot install timeout handler: %s\n",
                 strerror( errno ));
    }
    sigterm_status=0;
#if ( ! defined HP ) && ( ! defined IRIX64 )
    if ( sigaction( SIGBUS, &in, NULL ) == -1 ) {
        fprintf( stderr, "Warning: cannot install interrupt handler: %s\n",
                 strerror( errno ));
    }
#endif
    if ( sigaction( SIGFPE, &in, NULL ) == -1 ) {
        fprintf( stderr, "Warning: cannot install interrupt handler: %s\n",
                 strerror( errno ));
    }
    if ( sigaction( SIGSEGV, &in, NULL ) == -1 ) {
        fprintf( stderr, "Warning: cannot install interrupt handler: %s\n",
                 strerror( errno ));
    }
    if ( sigaction( SIGHUP, &in, NULL ) == -1 ) {
        fprintf( stderr, "Warning: cannot install interrupt handler: %s\n",
                 strerror( errno ));
    }

      /*
       * create request and reply fifos
       */ 
    sprintf( req_fifo, "/tmp/scam_req%d", (int)getpid() );
    if ( mkfifo( req_fifo, S_IRWXU ) < 0 ) {
        fprintf( stderr, "Error: cannot create fifo: %s\n",
                 strerror( errno ));
        exit( -1 );
    }
    sprintf( reply_fifo, "/tmp/scam_reply%d", (int)getpid() );
    if ( mkfifo( reply_fifo, S_IRWXU ) < 0 ) {
        fprintf( stderr, "Error: cannot create fifo: %s\n",
                 strerror( errno ));
        exit( -1 );
    }
    
      /*
       * fork and exec the user interface process
       */
    if ( ( gui_pid = fork() ) < 0 ) {
        fprintf( stderr, "Error: fork error: %s\n",
                 strerror( errno ));
        cleanup();
        exit( -1 );
    }
    if ( gui_pid == 0 ) { /* child process */
        if ( execv( "./scamgui", argv ) < 0 ) {
            fprintf( stderr, "Error: cannot start GUI (scamgui): %s\n",
                     strerror( errno ));
            cleanup();
        }
        exit( -1 );             /* shouldn't get here */
    }

      /*
       * open the fifo for reading/writing 
       */
    
    sprintf( req_fifo, "/tmp/scam_req%d", (int)getpid() );
    sprintf( reply_fifo, "/tmp/scam_reply%d", (int)getpid() );
    
    if (( req_fd = open( req_fifo, O_RDONLY ) ) < 0 ) {
        fprintf( stderr, "Server error: cannot open request fifo: %s : %s\n",
                 req_fifo, strerror( errno ));
        exit( -1 );
    }
    
    if (( reply_fd = open( reply_fifo, O_WRONLY ) ) < 0 ) {
        fprintf( stderr, "Server error: cannot open reply fifo: %s : %s\n",
                 reply_fifo,strerror( errno ));
        exit( -1 );
    }

}

/*---------------------------------------------------------
 * server_run: blocks on read from request pipe,
 *      calls appropriate function to handle requests
 *---------------------------------------------------------*/
int
server_run()
{
    srv_req   req;
    srv_reply reply;
    
    while (sigterm_status==0) {
      again:
        if ( read( req_fd, &req, sizeof( req )) != sizeof( req ) ) {
                                /* check for interrupted system call */
            if (errno == EINTR && sigterm_status==0) goto again;
            return 0;
        }
        
        switch ( req.id ) {
        case INIT_MODEL_REQ_ID:
            reply.status = init_model_svc( req.buf, reply.buf );
            break;
        case FIELD_INFO_REQ_ID:
            reply.status = get_field_info_svc( req.buf, reply.buf );
            break;
        case GET_FIELD_REQ_ID:
            reply.status = get_field_svc( req.buf, reply.buf );
            break;
        case GET_LEVELS_REQ_ID:
            reply.status = get_levels_svc( req.buf, reply.buf );
            break;
        case SET_FIELD_REQ_ID:
            reply.status = set_field_svc( req.buf, reply.buf );
            break;
        case RESET_FIELD_REQ_ID:
            reply.status = reset_field_svc( req.buf, reply.buf );
            break;
        case STEP_REQ_ID:
            reply.status = step_svc( req.buf, reply.buf );
            break;
        default:
            fprintf( stderr, "Server error: invalid request id: %d\n", req.id );
            break;
        }
        
        if ( write( reply_fd, &reply, sizeof( reply ) ) < sizeof( reply ) ) {
	    fprintf( stderr, "Server error: WriteRequest() - write: %s\n",strerror( errno ));
	    return -1;
	}
    }
    return 0;
}
