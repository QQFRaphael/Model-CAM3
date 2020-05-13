/*------------------------------------------------------------------------*
 * File: main.cpp 
 * $Author: jmccaa $
 * $Id: main.cpp,v 1.1.6.2 2004/08/19 15:05:49 jmccaa Exp $ *
 *------------------------------------------------------------------------*/

#ifndef lint
static char rcsid[] = "$Id: main.cpp,v 1.1.6.2 2004/08/19 15:05:49 jmccaa Exp $";
#endif /* lint */

//#include <strstream>
#include <sstream>
#include <cstdlib>
#include <unistd.h>
#include <ctime>
#include <qapp.h>

#include "MainWndImpl.h"
#include "manager.h"
#include "model.h"
#include "msgdlg.h"
#include "observer.h"
#include "crm.h"

void
PrintUsage();

int
main( int argc, char **argv )
{
    string remoteHost("localhost"); 
    char validOptions[] = "i:c:n:o:r:t:h:d";
    bool useGUI = true; 
    ModelType  modelType = Model::FIFO; // run fifo model by default
    int c;                      // option flag
    string startupfile = "";    // file containing startup parameters
    string savefile = "";       // file in which to save model output
    string crm_inputfile = ""; // ascii input file for crm
    string crm_outputfile = ""; // ascii output file for crm
    int  steps = 0;             // number of model steps to run
    int  reps = 1;              // number of repetions
    bool crm = false;           // flag for CRM model

      //-----------------------------------------------------------
      // parse command line options
      //-----------------------------------------------------------
    if ( argc > 1 ) {
        while ( (c = getopt( argc, argv, validOptions ) ) != EOF ) {
            switch (c) {
	    case 'i': // -i <crm_inputfile>; ascii input file for crm
	        if ( ! optarg ) { PrintUsage(); exit(1); }
		crm_inputfile = optarg;
		crm = true;
		break;
	    case 'c': // -c <crm_outputfile>; ascii output file for crm
	        if ( ! optarg ) { PrintUsage(); exit(1); }
		crm_outputfile = optarg;
		break;
            case 'n':               // -ng = no gui
                if ( strcmp( optarg, "g" ) ) {
                    cerr << "ERROR: "__FILE__":"<<__LINE__
                         << "main() - Unrecognized option -" 
                         << optopt <<" " << optarg << endl;
                    PrintUsage();
                }
                else
                    useGUI = false;
                break;
            case 'o':               // name for saving history output
                if ( ! optarg )
                    PrintUsage();
                savefile = optarg;
                break;
            case 'r':               // reps
                if ( ! optarg )
                    PrintUsage();
                reps = atoi( optarg );
                if ( reps < 1 ) {
                    cerr << "ERROR: "__FILE__", line" << __LINE__ 
                         << ": main() - Invalid number of reps: " << reps << endl;
                    PrintUsage();
                }
                break;
            case 'd':               // debug model
                 modelType = Model::DBG;
                break;
            case 't':               // number of timesteps to run
                steps = atoi( optarg );
                if ( steps < 0 ) {
                    cerr << "ERROR: "__FILE__",line" << __LINE__
                         << " main(): invalid number of steps: " << steps << endl;
                    PrintUsage();
                }
                break;
            case 'h':           // specify remote hostname
                if ( ! optarg )
                    PrintUsage();
                remoteHost =  optarg;
                modelType = Model::RPC;
                break;
            case '?':
                PrintUsage();
                break;
            default:
                cerr << "ERROR: "__FILE__":" << __LINE__
                     << " main(): Unrecognized option \"" << optopt << "\"\n" ;
                PrintUsage();                
            }
        }
    }
      // check if there is an argument that is not an option
      // - this would be a quickstart file
    if ( optind < argc )
        startupfile = argv[optind++];
      // set graphical mode to false until GUI is initialized (see msgdlg.cpp)
    SetGraphicalMode( false );

    // CRM creates and loads its own defaults file
    if ( crm ) {
      if ( startupfile.empty() ) startupfile = DEFAULTS_FILE;
      startupfile = CRM.Setup(startupfile,crm_inputfile);
      if ( !crm_outputfile.empty() ) CRM.SaveCRMTextFile( crm_outputfile.data() );
    }
    // the first reference to MANAGER actually creates the single global static
    // instance (i.e., it is a singleton pattern)
    // If specified, the startup file becomes the default file
    if ( startupfile.empty() ) {
      MANAGER.SetDefaultsFile( DEFAULTS_FILE );
    } else {
      MANAGER.SetDefaultsFile( startupfile );
    }

    MANAGER.CreateModel( modelType );
      
    if ( modelType == Model::RPC )
        MANAGER.SetRemoteHost( remoteHost );
    
      //-----------------------------------------------------------
      // START - GUI
      //-----------------------------------------------------------
    if ( useGUI ) { 
          // create and start up the GUI.
      //        sleep (10);
        QApplication theApp( argc, argv );
        SetGraphicalMode( true );
        MainWndImpl* gui = new MainWndImpl(0,"MainWndImpl",FALSE,0,startupfile.c_str());
        (void) &gui; // suppress unused variable warning
        return theApp.exec();   // give control to Qt's event loop
    }
      //-----------------------------------------------------------
      // START - NO GUI
      //-----------------------------------------------------------
    if ( startupfile.empty() ) {
        cerr << "Start-up filename is mandatory with -ng option!\n";
        PrintUsage();
    }
    
    MANAGER.SetShowSettings( true );

    try {
        MANAGER.LoadDefaults( startupfile, true );
    } catch ( IOErr& e ) {
        cerr << "ERROR: "__FILE__", line " << __LINE__
             << " main(): While loading startupfile " << startupfile << ": " << e << endl;
        exit ( -1 );
    }     
    if ( steps >  MANAGER.MaxStep() ) {
        cerr << "ERROR: "__FILE__":" << __LINE__
             << " main(): Max number of steps available in dataset: " << MANAGER.MaxStep() << endl;;
        exit ( -1 );
    }
    if ( steps > 0 && steps <= MANAGER.MaxStep() )
        MANAGER.SetEndStep( steps );
      //
      // do <reps> repetitions of the run; useful if 
      // doing something with stochastic processes
      //
    MANAGER.InitModel();
    for ( int r=0; r<reps; r++ ) {
        cout << "Model step: 0" ;
        cout.flush();
        while ( MANAGER.CurrentStep() < MANAGER.EndStep() ) {
              // advance the model. use the manager to do this so
              // that the fields will get saved properly
            MANAGER.StepModel();
              // print out the current step every 100 steps
            if (( MANAGER.CurrentStep() % 100 ) == 0 ) 
                cout << MANAGER.CurrentStep();
            else
                cout << ".";
            cout.flush();
        }
        
        if ( savefile.empty() ) { // no savefile name given
	  //            stringstream s;
	  //            s << "./userdata/scam_" << time(NULL) << ".nc";
	  //            savefile = s.str();
	  savefile = MANAGER.HistName();
        }
        try { 
            if ( reps > 1 ) {
                char filename[256];
                sprintf( filename, "%s.%d", savefile.c_str(), r );
                MANAGER.SaveHistoryFile( filename );
            }
            else {
                MANAGER.SaveHistoryFile( savefile );
                cout << "\nSaved model output in " <<  savefile << endl;
            }
        }
        catch ( NcErr& e ) {
            cerr << "ERROR: "__FILE__":" << __LINE__ 
                 << " main(): Can't save history: \n" << e << endl;
            exit( -1 );
        }
    MANAGER.InitModel();
    }
    cout << "\nscam: done." << endl;
    return 0;
}


void
PrintUsage()
{
    cerr << "usage: scam [-ng [-o <outputfile>] [-t <timesteps>] [-r <repetions>]] [-h <server>] [<start-up file>]" 
         << endl;;

    exit( -1 );
}






