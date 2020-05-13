// $Id: fifomodel.cpp,v 1.1.6.2 2004/08/19 15:05:47 jmccaa Exp $
// Implementation of FifoModel class


#if ( defined IRIX ) || ( defined sgi )
  // need to define this to fix signal.h sa_handler definition 
#define _LANGUAGE_C_PLUS_PLUS
#endif

#ifndef lint
static char rcsid[]="$Header: /fs/cgd/csm/models/CVS.REPOS/tools/atm_tools/ccm/scam/ui/Attic/fifomodel.cpp,v 1.1.6.2 2004/08/19 15:05:47 jmccaa Exp $";
#endif

#include <cerrno>
#include <csignal>              // signals
#include <cstring>              // strncpy
#include <sys/types.h>          // system types
#include <unistd.h>             // getpid()
#include <fcntl.h>              // file ops

#include "max.h"                // definition of NUM_SWITCHES
#include "ipc.h"                // structs used in interprocess communication
#include "dataset.h"
#include "field.h"              // declaration of Field class
#include "fifomodel.h"          // declaration of FifoModel class
#include "msgdlg.h"
#include "manager.h"

#define TIMEOUT (100)           // freq of alarms in seconds

static  char reply_fifo[MAX_PATH_LEN];
static  char req_fifo[MAX_PATH_LEN];
static  int  reply_fd;   // file descriptor
static  int  req_fd;    // file descriptor
static  int  server_pid; // server process id

//
// The following routines have to be declared as C, because of some
//  screwups in the system header file /usr/include/sys/signal.h on
//  some platforms.  
//
extern "C" {
//------------------------------------------------------------------
// Cleanup:
//------------------------------------------------------------------
static void
cleanup()
{  
    close( req_fd );
    close( reply_fd );
    remove( req_fifo );
    remove( reply_fifo );
    kill( getppid(), SIGTERM );
}

//--------------------------------------------------------
// handle_intrpt: handles interrupts to exit gracefully
//               input: number of the signal caught - unused
//--------------------------------------------------------
static void
handle_intrpt( int sig )
{
    int  ppid;                  /* parent process id */
    
    if ( sig == SIGALRM ) {
          //
          // check if parent process is still alive
          // if it has terminated and we've been inherited by init
          // (process id = 1), then we should exit
          //
        ppid = getppid();
        if ( ppid != 1 ){
              //
              // parent still alive; reset the timeout alarm
              //
            alarm(TIMEOUT);
            return;
        }
    }
      //
      // All signals besides SIGALRM will cause exit
      // 
      // print out the signal number unless it's just SIGTERM

    if ( sig != SIGTERM )
        cerr << "scamgui: caught signal " << sig <<", exiting.\n";
    cleanup();
    exit( 0 );
}

void
install_sig_handler()
{    
    //
    // Set up a handler for catching various signals
    // that will allow us to exit gracefully, without
    // leaving open pipes and the server running
    //
    struct sigaction in;        // signal for interrupt 

    in.sa_handler = handle_intrpt;
    in.sa_flags = 0;
    if ( sigemptyset( &in.sa_mask ) == -1 ) {
        cerr << "WARNING "__FILE__":"<<__LINE__
             << " : cannot clear interrupt handler: " << strerror( errno ) << endl;
    }
    if ( sigaction( SIGINT, &in, NULL ) == -1 ) {
        cerr << "WARNING "__FILE__":"<<__LINE__
             << " : cannot install interrupt handler: " << strerror( errno )<< endl;
    }
    if ( sigaction( SIGALRM, &in, NULL ) == -1 ) {
        cerr << "WARNING "__FILE__":"<<__LINE__
             << " : cannot install interrupt handler: " << strerror( errno )<< endl;
    }
    if ( sigaction( SIGTERM, &in, NULL ) == -1 ) {
        cerr << "WARNING "__FILE__":"<<__LINE__
             << " : cannot install interrupt handler: " << strerror( errno )<< endl;
    }
    if ( sigaction( SIGHUP, &in, NULL ) == -1 ) {
        cerr << "WARNING "__FILE__":"<<__LINE__
             << " : cannot install interrupt handler: " << strerror( errno )<< endl;
    }
    if ( sigaction( SIGPIPE, &in, NULL ) == -1 ) {
        cerr << "WARNING "__FILE__":"<<__LINE__
             << " : cannot install interrupt handler: " << strerror( errno )<< endl;
    }
    if ( sigaction( SIGBUS, &in, NULL ) == -1 ) {
        cerr << "WARNING "__FILE__":"<<__LINE__
             << " : cannot install interrupt handler: " << strerror( errno )<< endl;
    }
    if ( sigaction( SIGFPE, &in, NULL ) == -1 ) {
        cerr << "WARNING "__FILE__":"<<__LINE__
             << " : cannot install interrupt handler: " << strerror( errno )<< endl;
    }
    if ( sigaction( SIGSEGV, &in, NULL ) == -1 ) {
        cerr << "WARNING "__FILE__":"<<__LINE__
             << " : cannot install interrupt handler: " << strerror( errno )<< endl;
    }
}

} // extern C

//------------------------------------------------------
//  FifoModel: constructor
//------------------------------------------------------

FifoModel::FifoModel()
        : Model()
{
    install_sig_handler();
    
      // register the cleanup function to be called at exit
    if ( atexit( cleanup ) != 0 )
        cerr << "WARNING "__FILE__":"<<__LINE__
             << " FifoModel::FifoModel() - can't register exit function \"cleanup()\n";

    server_pid = getppid();


    alarm(TIMEOUT);             // set the timeout

      // 
      // open the  request fifo for writing, and the reply fifo for reading
      //
      //  
    sprintf( req_fifo, "/tmp/scam_req%d", (int)getppid() );
    sprintf( reply_fifo, "/tmp/scam_reply%d", (int)getppid() );


    if (( req_fd = open( req_fifo, O_WRONLY ) ) < 0 ) {
        ShowMsg( __FILE__, __LINE__, "ERROR: Cannot connect to scam (not running?)\nScamgui is executed automatically by scam.");
        cleanup();
        exit( -1 );
        ;
    }
    if (( reply_fd = open( reply_fifo, O_RDONLY ) ) <  0 ) {
        ShowMsg( __FILE__, __LINE__, "ERROR: Cannot connect to scam (not running?)");
        cleanup();
        exit( -1 );
    }

}
    
//-----------------------------------------------------
//  ~FifoModel: destructor, cleans up client connection
//-----------------------------------------------------

FifoModel::~FifoModel()
{
    cleanup();
}


//---------------------------------------------------------------------------
// BuildFieldList: builds the list of fields by retrieving information stored 
//                  in the model's outfld buffer 
//---------------------------------------------------------------------------
void
FifoModel::BuildFieldList()
{
    field_info* fi;             // summary of field information from outfld buffer 
    int* fid = (int*)req.buf;   // id of field to retrieve
    *fid = 0;
    req.id = FIELD_INFO_REQ_ID;

    fi = (field_info*)req.buf;

    while (1) {
        reply.status =  -1;
        WriteReq();
        ReadReply();
        fi = (field_info*)reply.buf;
        if ( reply.status == 0 ) {
              // check for duplicates
            if ( fields.find( fi->name ) != fields.end() )
                break;
            Field* f = new Field( fi->name, fi->longname, fi->units,
                                  fi->std_units, fi->mult, fi->min,
                                  fi->max, bool(fi->is_shown), bool(fi->is_modifiable),
                                  bool(fi->is_averaged), fi->size, fields.size() );
              // insert the field into the map (see map.h)
            fields.insert( fieldpair( f->Name(), f ));
        }
        else
            break;
        ++*fid;                 // get next field
    }
}

//-----------------------------------------------------
//  Init: calls the model initialization function
//         inputs: model initialization pararameters
//-----------------------------------------------------
int
FifoModel::Init(bool isRestart)
{
    init_data* ini;        // init data struct
    req.id = INIT_MODEL_REQ_ID;
    ini = (init_data*)req.buf; // map init data struct to request buffer

      //
      // Fill in the init data struct with the data used for initialization
      //

    for ( int i=0; i<NUM_SWITCHES; i++ )
        ini->sw[i] = switches[i];
    ini->lat = latitude;
    ini->lon = longitude;
    ini->bd = baseDate;
    ini->bs = baseSecs;
    ini->rt = runType;
    ini->sl = stepLen;
    ini->rst = isRestart;
    ini->ie = MANAGER.SeedVal();
    strcpy( ini->ana, dataset[ANAL]->name().c_str() );
    strcpy( ini->ini, dataset[MODEL]->name().c_str() );
    strcpy( ini->iop, dataset[IOP]->name().c_str() );
    strcpy( ini->lsmini, dataset[LSMINI]->name().c_str() );
    strcpy( ini->ozo, dataset[OZON]->name().c_str() );
    strcpy( ini->prs, dataset[PRES]->name().c_str() );
    strcpy( ini->abs, dataset[ABSEMS]->name().c_str() );
    strcpy( ini->optics, dataset[AEROPTICS]->name().c_str() );
    strcpy( ini->mass, dataset[AERMASS]->name().c_str() );
    strcpy( ini->sst, dataset[SST]->name().c_str() );
    strcpy( ini->lsmpft, adataset[LSMPFT]->name().c_str() );
    strcpy( ini->lsmsurf, dataset[LSMSRF]->name().c_str() );
    if ( dataset[SIC] ) strcpy( ini->sic, dataset[SIC]->name().c_str() );
    if ( dataset[USER] ) strcpy( ini->usr, dataset[USER]->name().c_str() );
    
    WriteReq();
    ReadReply();

    if ( reply.status != 0 ) // error on initialization
        return  reply.status;
    
    currLevs = baseLevs = GetLevels(); // get the model pressure levels    
    currILevs = baseILevs = GetILevels(); // get the model interface pressure levels    

    currentStep = 0;
    return 0;
}

//-----------------------------------------------------
//  ReadField: Retrieves a field's value from the model
//           input: Field object to get data for
//-----------------------------------------------------
void
FifoModel::ReadField( Field& f )
{
    int* field_id;

    req.id = GET_FIELD_REQ_ID;
    field_id = (int*)req.buf;
    reply.status =  -1;
    
    *field_id = f.FieldID();
    WriteReq();
    ReadReply();
    
    if ( reply.status != 0 ) {
        cerr << "ERROR: "__FILE__"line, " << __LINE__
             << " : ReadField() - bad server reply status: " << reply.status <<endl;
        cleanup();
        exit( -1 );
    }
    
    f.SetData( (real_t*)reply.buf );
}


//-----------------------------------------------------
//  GetNumLevels: Retrieves the models pressure levels data
//-----------------------------------------------------

int
FifoModel::GetNumLevels()
{
  level_data* ld;
  int* number_only;

  number_only = (int*)req.buf;
  *number_only = 1;
  req.id = GET_LEVELS_REQ_ID;
  reply.status =  -1;
  
  WriteReq();
  ReadReply();
  
  if ( reply.status != 0 ) {
    cerr << "ERROR: "__FILE__":" << __LINE__
	 << " : FifoModel::GetNumLevels() - bad server reply status" << endl;;
    cleanup();
    exit( -1 );
  }
  
  
  ld = (level_data*)reply.buf;
  numLevs = ld->nlev;
  
  return numLevs;
}

//-----------------------------------------------------
//  GetLevels: Retrieves the models pressure levels data
//-----------------------------------------------------

const vector<real_t>&
FifoModel::GetLevels()
{
    level_data* ld;
    int* number_only;

    number_only = (int*)req.buf;
    *number_only = 0;
    req.id = GET_LEVELS_REQ_ID;
    reply.status =  -1;

    WriteReq();
    ReadReply();

    
    if ( reply.status != 0 ) {
        cerr << "ERROR: "__FILE__":" << __LINE__
             << " : FifoModel::GetLevels() - bad server reply status" << endl;;
        cleanup();
        exit( -1 );
    }
    
    
    ld = (level_data*)reply.buf;
    numLevs = ld->nlev;
    currLevs.resize( numLevs );
    std::copy( ld->levels, ld->levels+numLevs, currLevs.begin() );

    return currLevs;
}
//-----------------------------------------------------
//  GetILevels: Retrieves the model interface pressure levels data
//-----------------------------------------------------

const vector<real_t>&
FifoModel::GetILevels()
{
    level_data* ld;
    int* number_only;

    number_only = (int*)req.buf;
    *number_only = 0;
    req.id = GET_LEVELS_REQ_ID;
    reply.status =  -1;

    WriteReq();
    ReadReply();

    
    if ( reply.status != 0 ) {
        cerr << "ERROR: "__FILE__":" << __LINE__
             << " : FifoModel::GetILevels() - bad server reply status" << endl;;
        cleanup();
        exit( -1 );
    }
    
    
    ld = (level_data*)reply.buf;
    numILevs = ld->ilev;
    currILevs.resize( numILevs );
    std::copy( ld->ilevels, ld->ilevels+numILevs, currILevs.begin() );

    return currILevs;
}
//---------------------------------------------------------------------------
// Resets an averaged field
//---------------------------------------------------------------------------
void
FifoModel::ResetField( const Field& f )
{
    int* field_id;

    req.id = RESET_FIELD_REQ_ID;
    field_id = (int*)req.buf;
    *field_id = f.FieldID();
    reply.status =  -1;

    WriteReq();
    ReadReply();
    
    if ( reply.status != 0 ) {
        cerr << "ERROR: "__FILE__":" << __LINE__
             << " : FifoModel::ResetField() - bad server reply status for field " << f.Name() << endl;;
        cleanup();
        exit( -1 );
    }
}


//---------------------------------------------------------------------------
// WriteField: Sets a field's value in the model
//          input: Field object containing data to be set in the model
//---------------------------------------------------------------------------
void
FifoModel::WriteField( const Field& f )
{
    field_data* fd;

    req.id = SET_FIELD_REQ_ID;
    fd = (field_data*)req.buf;
    reply.status =  -1;

    fd->id = f.FieldID();
    copy( f.Data().begin(), f.Data().begin()+f.NumLevs(), fd->data );
    WriteReq();
    ReadReply();
        
    if ( reply.status != 0 ) {
        cerr << "ERROR: "__FILE__":" << __LINE__ 
             << " FifoModel::WriteField() : bad server reply status: " << reply.status << endl;
        cleanup();
        exit( -1 );
    }
}

//-----------------------------------------------------
//  Step: calls the model time stepping function
//-----------------------------------------------------
void
FifoModel::Step()
{
    req.id = STEP_REQ_ID;
    reply.status =  -1;
        
    WriteReq();
    ReadReply();
    if ( reply.status != 0 ){
        cerr<< "ERROR: "__FILE__"line, " << __LINE__ 
            << " : FifoModel::Step() - bad server reply status: " << reply.status<<endl;
        cleanup();
        exit( -1 );
    }
    currentStep++;    
}

void
FifoModel::WriteReq()
{
    if ( write( req_fd, &req, sizeof( req ) ) < (int)sizeof( req ) ) {
        cerr << "ERROR: "__FILE__":"<<__LINE__
             <<" : FifoModel::WriteReq(): " << strerror( errno ) << endl;
        cleanup();
        exit( -1 );
    }
}

void
FifoModel::ReadReply()
{
  int temp;
  again:
  temp = read( reply_fd, &reply, sizeof( reply ) );
    if ( temp != sizeof( reply ) ) {
      //    if ( read( reply_fd, &reply, sizeof( reply ) ) != sizeof( reply ) ) {
        if ( errno == EINTR )   
              // interrupted system call, ignore and restart read()
            goto again;
        if ( errno == 0 )  { 
              // read didn't fail but got fewer bytes than expected:
              // this indicates that the server has gone down
            cerr << "ERROR: "__FILE__":" << __LINE__ 
                 << " FifoModel::ReadReply() - Model aborted" << endl;
	    cerr << "Expected " << sizeof( reply ) << " bytes." << endl;
	    cerr << "Received " << temp << " bytes." << endl;
            cleanup();
            exit( -1 );
        }
    }
}


