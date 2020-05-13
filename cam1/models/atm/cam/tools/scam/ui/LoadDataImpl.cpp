#ifndef lint
static char rcsid[] = "$Id: LoadDataImpl.cpp,v 1.1.6.2 2004/08/19 15:05:43 jmccaa Exp $";
#endif /* lint */


#include <stdlib.h>
#include <stdio.h>
#include <qfiledlg.h>
#include <qtooltip.h>
#include <iostream>

#ifndef _RWSTD_NO_NAMESPACE
  using namespace std;
#endif

#include "defaults.h"
#include "msgdlg.h"
#include "MainWndImpl.h"
#include "ncfile.h"
#include "manager.h"
#include "runtype.h"
#include "LoadDataImpl.h"
#include "crm.h"
#include "string.h"

#define QS 100  // the quickstart button ID
#define CRM_BUTTON 101 // the CRM button ID

/* 
 *  Constructs a LoadDataImpl which is a child of 'parent', with the 
 *  name 'name' and widget flags set to 'f' 
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
LoadDataImpl::LoadDataImpl( QWidget* parent,  const char* name, bool modal, WFlags fl )
    : SelectDataForm( parent, name, modal, fl )
{
    theGUI = ( MainWndImpl* )parent;
    //    connect( CancelPB, SIGNAL( clicked() ), this, SLOT( hide() ) );
}

/*  
 *  Destroys the object and frees any allocated resources
 */
LoadDataImpl::~LoadDataImpl()
{
    // no need to delete child widgets, Qt does it all for us
}

/* 
 * protected slot
 */
void LoadDataImpl::SetRunType(int runtype)
{
  //    qWarning( "LoadDataImpl::SetRunType(int) not yet implemented!" ); 
    QString filename;
    QFileInfo  fi( MANAGER.GetDataset( MANAGER.RunType() ).name().c_str());
    string filestring;
    string defaultsfile;

    switch( runtype ) {
    case MODEL:
        hide();
        MANAGER.SetRunType( MODEL );
        theGUI->ShowGlobalDlg();
        break;
    case ANAL:
        hide();
        MANAGER.SetRunType( ANAL );
        theGUI->ShowGlobalDlg();
        break;
    case IOP:
        hide();
        MANAGER.SetRunType( IOP );
        theGUI->ShowIopDlg();
        break;
    case QS:
        filename = QFileDialog::getOpenFileName( ".", "*.scm" );
        if ( filename.isNull() )
            break;
        hide();

        try {
            MANAGER.LoadDefaults( filename.data(), TRUE );
        } catch ( IOErr& e ) {
            ShowMsg( __FILE__, __LINE__, "ERROR: in QuickStart file: %s\n",  e.toString().c_str() );
            return;
        }            
        
        MANAGER.InitModel();
        break;
    case CRM_BUTTON:
        filename = QFileDialog::getOpenFileName( ".", "*.in" );
        if ( filename.isNull() )
	  break;
        hide();
	filestring = filename.data();
	defaultsfile = MANAGER.GetDefaultsFile();
	
	cout << "Reading from " << filestring << endl;
	defaultsfile = CRM.Setup(defaultsfile,filestring);

        try {
	  MANAGER.LoadDefaults( defaultsfile, TRUE );
        } catch ( IOErr& e ) {
	  ShowMsg( __FILE__, __LINE__, "ERROR: in QuickStart file: %s\n",  e.toString().c_str() );
	  return;
        }            

        MANAGER.InitModel();
        break;
    case USER:
        filename = QFileDialog::getOpenFileName( fi.dirPath(), "*.nc", this );
        if ( filename.isNull() )
            break;
        hide();
        try {
            MANAGER.SetDataset( filename.data(), USER );
            MANAGER.SetRunType( USER );
            theGUI->ShowIopDlg();
        } catch( NcErr& e ) {
            ShowMsg( __FILE__, __LINE__, "ERROR: - in User dataset: %s\n",  e.toString().c_str() );
            break;
        }
        break;
    case SIC:
        hide();
        filename = QFileDialog::getOpenFileName( fi.dirPath(), "*.nc" );
        if ( filename.isNull() )
            break;
        try {
            MANAGER.LoadInitialConditionsFile( filename.data() );
            MANAGER.InitModel();
        } catch ( NcErr& e ) {
            ShowMsg( __FILE__, __LINE__, "ERROR: Unable to load Saved Initial Conditions file: %s\n", e.toString().c_str() ); 
        }
          //  the runtype is set in LoadInitialConditionsFile() above
        break;
    default:
        cerr << "ERROR: "__FILE__":" << __LINE__ 
             << " LoadDlg::SetRunType() - unknown type " << runtype << endl;
        exit( -1 );
    }

}

