
#ifndef lint
static char rcsid[] = "$Id: MainWndImpl.cpp,v 1.1.6.2 2004/08/19 15:05:43 jmccaa Exp $";
#endif /* lint */


#include <sstream>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>

#ifndef _RWSTD_NO_NAMESPACE
  using namespace std;
#endif

#include <qapp.h>
#include <qfiledlg.h>
#include <qfont.h>
#include <qmsgbox.h>
#include <qpixmap.h>
#include <qpoint.h>
#include <qpushbutton.h>
#include <qcombobox.h>
#include <qtooltip.h>

#include "defaults.h"
#include "SelectGlobalDataDlgImpl.h"
#include "manager.h"
#include "IOPSelectDateDlgImpl.h"
#include "ncfile.h"
#include "OptionsDlgImpl.h"
#include "plot.h"
#include "PlotDlgImpl.h"
#include "PostPlottingDlgImpl.h"
#include "LoadDataImpl.h"
#include "msgdlg.h"
#include "timeconvert.h"
#include "fieldlistbox.h"
#include "MainWndImpl.h"

/* 
 *  Constructs a MainWndImpl which is a child of 'parent', with the 
 *  name 'name' and widget flags set to 'f' 
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
MainWndImpl::MainWndImpl( QWidget* parent,  const char* name, bool modal, WFlags fl, const char* startupfile )
    : MainWnd( parent, name, modal, fl )
{
      // watch manager
    MANAGER.AddObserver( this );
    
    
      // create the various dialogs associated with the main window buttons    
    IopDlg = new IOPSelectDateDlgImpl(this );

    IopDlg->hide();
    LoadDlg = new LoadDataImpl( this );    
    LoadDlg->hide();
    GlobalDlg = new SelectGlobalDataDlgImpl( this );
    GlobalDlg->hide();
    OptionDlg = new OptionsDlgImpl( this );
    OptionDlg->Init();
    OptionDlg->hide();
    
      // stepTimer repeatedly generates a 'timeout' signal
      // when all window events have been processed, the signal
      // is connected to the Step() routine.
    stepTimer = new QTimer( this );	

    connect( stepTimer, 
             SIGNAL( timeout() ),        SLOT( StepBtnClicked() ));

      // set the step button to autorepeat when held down
    StepPB->setAutoRepeat( TRUE );
    
      // set the initial state for the main window buttons/widgets
    isRunning = FALSE;
    currentState = STARTING;
    
      // if quickstart provided, load it and initialize the model
    if ( startupfile[0] != 0 ) {
        try {
            MANAGER.LoadDefaults( startupfile, TRUE );
            MANAGER.InitModel( );
        } catch ( IOErr& e ) {
            ShowMsg( __FILE__, __LINE__,"ERROR: while trying to load quickstart file \"%s\":\n %s\n", startupfile, e.toString().c_str() );
        }
    } else {
      // load defaults file
      
      try {
        MANAGER.LoadDefaults( DEFAULTS_FILE );
      } catch ( IOErr& e ) {
        ShowMsg( __FILE__, __LINE__,"ERROR: while trying to load startup defaults file \"%s\":\n %s\n", DEFAULTS_FILE, e.toString().c_str() );
        exit( -1 );
      }
    }

      // set the time widgets
      // note that timeFormat is set when LoadDefaults() is called
    theTimeFormatPulldown->setCurrentItem( MANAGER.TimeFormat() );
    
#ifdef NO_NCARG
    QToolTip::add( PostPlotPB, "Post-plotting disabled" );
    EnablePostPlotting( FALSE );
#else
    if ( ! getenv( "NCARG_ROOT" ) ) {
        ShowMsg( __FILE__, __LINE__, "Warning: NCARG_ROOT not defined!\nYou must set the environment variable NCARG_ROOT\n to the scam root directory. Post plotting will be disabled." );
        QToolTip::add( PostPlotPB, "Post-plotting disabled" );
        EnablePostPlotting( FALSE );
    }else{
        QToolTip::add( PostPlotPB, "Post-plotting" );
        EnablePostPlotting( TRUE );
    }
#endif

      // set the main window to call exit() when closed
    //    parent->setMainWidget( this );
      // don't allow resizing
    setFixedSize( this->size() ); 
    Update();
    show();
}

/*  
 *  Destroys the object and frees any allocated resources
 */
MainWndImpl::~MainWndImpl()
{
    // no need to delete child widgets, Qt does it all for us
}

//
// enable or disable post plotting button
//
void
MainWndImpl::EnablePostPlotting( bool enable )
{
    PostPlotPB->setEnabled( postPlottingIsEnabled = enable );
}

//
//  add fields to the list in the main window 
// 
void
MainWndImpl::FillFieldListBox()
{
    char fieldName[100];

    if ( FieldListBoxCW->count() > 0 )
        return;                 // already filled

    FieldListBoxCW->setAutoUpdate( FALSE );
      // setting autoupdate to false speeds up fill    
    FieldListBoxCW->clear();   // remove old entries
    FieldListConstItr it = MANAGER.GetFieldList().begin();
    for ( ; it != MANAGER.GetFieldList().end(); ++it ) {
        Field* f = (*it).second;
        if ( f->IsListed() ) {
            sprintf( fieldName, "%s - %s" , f->Name().c_str(),  f->LongName().c_str() );
              // FieldListBoxText does the coloring of modifiable fields
            FieldListBoxCW->insertItem( new FieldListBoxText( fieldName, f->IsModifiable() ) );
        }
    }
    FieldListBoxCW->setAutoUpdate( TRUE );
    FieldListBoxCW->repaint();
}

void
MainWndImpl::Reposition( QWidget* w)
{
      //
      // set the positions of the dialogs relative to the desktop:
      // if they are partially off the desktop, move them so they
      // will be entirely within the desktop
      //
    QWidget* desktop = QApplication::desktop();
    int dr, db, wl, wr, wt, wb, ww, wh;
    
    dr = desktop->width();
    db = desktop->height();
    wl = w->x();
    ww = w->width();
    wr = wl + ww; 
    wt = w->y();
    wh = w->height();
    wb = wt + wh;
    
      // check left - right position
    if ( wl < 20 )
        wl = 20;
    else if ( wr > dr )
        wl = (dr - ww) - 20;
      // check top -  botton position
    if ( wt < 20 )
        wt = 20;
    else if ( wb > db ) 
        wt = (db - wh) - 20;
    
    w->move( wl, wt );
}
/* 
 * protected slot
 */
//
// FieldSelected()
//  slot that receives signal when a fieldListBox item is double clicked
//        - index is position in list
void MainWndImpl::FieldSelected(int index)
{
    if ( index < 0 || !FieldListBoxCW->isEnabled() )
        return;

    Field* f = MANAGER.FindField( FieldListBoxCW->FieldName( index ) );
    PlotDlgImpl* p = new PlotDlgImpl( f );
    f->Update();
    p->show();
}
//
// set the current time in the mainwnd
//
void
MainWndImpl::SetCurrentTimeDisplay( int steps )
{
    TimeConverter tc( steps );
    CurrentLE->setText( tc.TimeToString().c_str() );
}


//
// set the end time in the mainwnd
//
void
MainWndImpl::SetEndTimeDisplay( int endStep )
{
    TimeConverter tc( endStep );
    EndLE->setText( tc.TimeToString().c_str() );
}

//
/* 
 * protected slot
 */
void MainWndImpl::LoadBtnClicked()
{
    Reposition( LoadDlg );
    LoadDlg->show();
    //    qWarning( "MainWndImpl::LoadBtnClicked() not yet implemented!" ); 
}
/* 
 * protected slot
 */
void MainWndImpl::OptionBtnClicked()
{
    // the Init() call is made when the dialog is shown (see OptionDlg::OptionDlg())
    OptionDlg->Init();
    Reposition( OptionDlg );
    OptionDlg->show();
    //    qWarning( "MainWndImpl::OptionBtnClicked() not yet implemented!" ); 
}
/* 
 * protected slot
 */
//
// post button callback - creates a new instance of the 
//   post dialog and initializes it. (it automatically becomes visible
//   when it is created)
//
void MainWndImpl::PostBtnClicked()
{
    PostDlg = new PostPlottingDlgImpl(); 
    if ( PostDlg->Init() ) {
        Reposition( PostDlg );
        PostDlg->show();
    }
    //    qWarning( "MainWndImpl::PostBtnClicked() not yet implemented!" ); 
}
/* 
 * protected slot
 */
//
// quit button callback: quits the application
//
void MainWndImpl::QuitBtnClicked()
{
    qApp->quit();
    delete &MANAGER;
    //    delete stepTimer;
    //    qWarning( "MainWndImpl::QuitBtnClicked() not yet implemented!" ); 
}
/* 
 * protected slot
 */
//
// Restart button callback
//
void MainWndImpl::RestartBtnClicked()
{
    MANAGER.InitModel( );
    //    qWarning( "MainWndImpl::RestartBtnClicked() not yet implemented!" ); 
}
/* 
 * protected slot
 */
//
// Run button callback: starts or stops the model
//
void MainWndImpl::RunBtnClicked()
{
    if ( isRunning == TRUE ) {
        stepTimer->stop();
        isRunning = FALSE;
    }
    else {
        stepTimer->start( 0, FALSE );
        isRunning = TRUE;
    }
    Update();
}
/* 
 * protected slot
 */
//
// Save button callback
//
void MainWndImpl::SaveBtnClicked()
{
    QString filename = QFileDialog::getSaveFileName( MANAGER.HistName().c_str(), "*.nc", NULL );
    if ( filename.isNull() ) 
        return;
    try {
        MANAGER.SaveHistoryFile( filename.data() );
    }
    catch ( NcErr& e ) {
        ShowMsg( __FILE__, __LINE__, "ERROR: Unable to save history file: %s\n", 
                 e.toString().c_str() );
    }
    //    qWarning( "MainWndImpl::SaveBtnClicked() not yet implemented!" ); 
}
/* 
 * protected slot
 */
//
// make the global data selection dialog visible
//
void MainWndImpl::ShowGlobalDlg()
{
    if ( GlobalDlg->Init() ) {
        Reposition( GlobalDlg );
        GlobalDlg->show();
    }
  //    qWarning( "MainWndImpl::ShowGlobalDlg() not yet implemented!" ); 
}
/* 
 * protected slot
 */
//
// make the iop date selection dialog visible
//
void MainWndImpl::ShowIopDlg()
{
    IopDlg->Init();
    Reposition( IopDlg );
    IopDlg->show();
    //    qWarning( "MainWndImpl::ShowIopDlg() not yet implemented!" ); 
}
/* 
 * protected slot
 */
//
// step button callback; also connected to stepTimer signal:
//  steps the model, also calls the save button callback 
//  after the last step
//
void MainWndImpl::StepBtnClicked()
{
    if ( MANAGER.CurrentStep() < MANAGER.EndStep() )
        MANAGER.StepModel();
    else {
        isRunning = FALSE;
        stepTimer->stop();
        SaveBtnClicked();
    }
    Update();
    //    qWarning( "MainWndImpl::StepBtnClicked() not yet implemented!" ); 
}
/* 
 * protected slot
 */
//
//  slot connected to timeformat pulldown menu:
//    changes the time format
//
void MainWndImpl::TimeFormatChanged(int format)
{
    MANAGER.SetTimeFormat( format );
    SetEndTimeDisplay( MANAGER.EndStep() );
    SetCurrentTimeDisplay( MANAGER.CurrentStep() );
  //    qWarning( "MainWndImpl::TimeFormatChanged(int) not yet implemented!" ); 
}
void MainWndImpl::EndLEClicked(const QString & timestr)
{
  EndLE->setFocus();
}

void MainWndImpl::SetEndTimeStep()
{
  QString timestr=EndLE->text();
  SetEndTimeStep(timestr);
}
void MainWndImpl::SetEndTimeStep(const QString & timestr)
{
    MANAGER.SetEndStep( timestr.toInt() );
    SetEndTimeDisplay( MANAGER.EndStep() );
}



//
//  implements a state machine to control the state of the gui
//
void
MainWndImpl::Update()
{
      //
      // determine whether to transition to a different
      // state based on the current state and the state of 
      // the model
      //
    switch ( currentState ) {
    case STARTING:
        if ( MANAGER.IsModelInited() == TRUE )
            currentState = INITIALIZED;
        break;
    case INITIALIZED:
        if ( isRunning )
            currentState = RUNNING;
        else if ( MANAGER.CurrentStep() > 0 )
            currentState = STOPPED;
        break;
    case RUNNING:
        if ( ! isRunning ) {
            if ( MANAGER.CurrentStep() != MANAGER.EndStep()  )
                currentState = STOPPED;
            else
                currentState = FINISHED;
        }
        break;
    case STOPPED:
        if ( isRunning )
            currentState = RUNNING;
        else if ( MANAGER.CurrentStep() == 0 )
            currentState = INITIALIZED;
        break;
    case FINISHED:
        if ( MANAGER.CurrentStep() == 0 )
            currentState = INITIALIZED;
        else if ( MANAGER.CurrentStep() != MANAGER.EndStep()  )
                currentState = STOPPED;
        break;
    }

      //
      //  set the new state of the widgets 
      //
    switch ( currentState ) {
    case STARTING:
        LoadDataPB->setEnabled( TRUE );
        LoadDataPB->setFocus();
        SaveDataPB->setEnabled( FALSE );
        OptionsPB->setEnabled( FALSE );
        PostPlotPB->setEnabled( postPlottingIsEnabled );
        QuitPB->setEnabled( TRUE );
        RunPB->setEnabled( FALSE );
        RunPB->setText( "&Run" );
        StepPB->setEnabled( FALSE );
        RestartPB->setEnabled( FALSE );
        FieldListBoxCW->setEnabled( FALSE );
        SetCurrentTimeDisplay( MANAGER.CurrentStep() );
        SetEndTimeDisplay( MANAGER.EndStep() );
        break;
    case INITIALIZED:
        LoadDataPB->setEnabled( TRUE );
        SaveDataPB->setEnabled( FALSE );
        OptionsPB->setEnabled( TRUE );
        PostPlotPB->setEnabled( postPlottingIsEnabled );
        QuitPB->setEnabled( TRUE );
        RunPB->setEnabled( TRUE );
        RunPB->setText( "&Run" );
        StepPB->setEnabled( TRUE );
        RestartPB->setEnabled( FALSE );
        FillFieldListBox();
        FieldListBoxCW->setEnabled( TRUE );
        FieldListBoxCW->setFocus(); 
        SetCurrentTimeDisplay( MANAGER.CurrentStep() );
        SetEndTimeDisplay( MANAGER.EndStep() );
        break;
    case RUNNING:
        LoadDataPB->setEnabled( FALSE );
        SaveDataPB->setEnabled( FALSE );
        OptionsPB->setEnabled( FALSE );
        PostPlotPB->setEnabled( FALSE );
        QuitPB->setEnabled( TRUE );
        RunPB->setEnabled( TRUE );
        RunPB->setText( "Stop" );
        StepPB->setEnabled( FALSE );
        RestartPB->setEnabled( FALSE );
        FieldListBoxCW->setEnabled( TRUE );
        FieldListBoxCW->setFocus(); 
        SetCurrentTimeDisplay( MANAGER.CurrentStep() );
          // don't call SetEndTimeDisplay() here: it causes flickering
          // because it's a line edit (slow updates).
        break;
    case STOPPED:
        LoadDataPB->setEnabled( FALSE );
        SaveDataPB->setEnabled( TRUE );
        OptionsPB->setEnabled( FALSE );
        PostPlotPB->setEnabled( postPlottingIsEnabled );
        QuitPB->setEnabled( TRUE );
        RunPB->setEnabled( TRUE );
        RunPB->setText( "&Run" );
        StepPB->setEnabled( TRUE );
        RestartPB->setEnabled( TRUE );
        FieldListBoxCW->setEnabled( TRUE );
        FieldListBoxCW->setFocus(); 
        SetCurrentTimeDisplay( MANAGER.CurrentStep() );
        SetEndTimeDisplay( MANAGER.EndStep() ); 
        break;
    case FINISHED:
        LoadDataPB->setEnabled( FALSE );
        SaveDataPB->setEnabled( TRUE );
        OptionsPB->setEnabled( FALSE );
        PostPlotPB->setEnabled( postPlottingIsEnabled );
        QuitPB->setEnabled( TRUE );
        RunPB->setEnabled( FALSE );
        RunPB->setText( "&Run" );
        StepPB->setEnabled( FALSE );
        RestartPB->setEnabled( TRUE );
        FieldListBoxCW->setEnabled( TRUE );
        FieldListBoxCW->setFocus(); 
        SetCurrentTimeDisplay( MANAGER.CurrentStep() );
        SetEndTimeDisplay( MANAGER.EndStep() );
        break;        
    default:
        cerr <<  "ERROR: "__FILE__ << ":"<<__LINE__ 
             <<" : MainWndImpl::Update() - illegal state: "<< int( currentState ) << endl;
        exit( -1 );
    }

}
