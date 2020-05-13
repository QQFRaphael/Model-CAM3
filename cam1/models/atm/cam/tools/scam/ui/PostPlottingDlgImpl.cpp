
#ifndef lint
static char rcsid[] = "$Id: PostPlottingDlgImpl.cpp,v 1.1.6.1 2004/05/18 20:51:49 jet Exp $";
#endif /* lint */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <qfiledlg.h>
#include <qtooltip.h>
#include <qpushbutton.h>
#include <qcombobox.h>
#include <qcheckbox.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <iostream>

#ifndef _RWSTD_NO_NAMESPACE
  using namespace std;
#endif

#include "fieldlistbox.h"
#include "MainWndImpl.h"
#include "max.h"
#include "manager.h"
#include "model.h"
#include "msgdlg.h"
#include "history.h"
#include "ncarg.h"
#include "timeconvert.h"
#include "PostPlottingDlgImpl.h"


#define POSTSCRIPT_BTN_ID 0
#define CGM_BTN_ID 1

#define WKSID 0
#define FID   1

using namespace ncfile;

int PostPlottingDlgImpl::numPlots = 0;
int PostPlottingDlgImpl::nextGksId = 0;

/* 
 *  Constructs a PostPlottingDlgImpl which is a child of 'parent', with the 
 *  name 'name' and widget flags set to 'f' 
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
PostPlottingDlgImpl::PostPlottingDlgImpl( QWidget* parent,  const char* name, bool modal, WFlags fl )
    : PostPlotDlg( parent, name, modal, fl )
{
    for ( int i=0; i<MAX_GKS_XWIN; i++ ) 
        plotId[i][FID] = -1;
    FieldListBoxCW->setMultiSelection( TRUE );
}

/*  
 *  Destroys the object and frees any allocated resources
 */
PostPlottingDlgImpl::~PostPlottingDlgImpl()
{
    ClosePlots();
    delete plotFile;
    // no need to delete child widgets, Qt does it all for us
}

void
PostPlottingDlgImpl::closeEvent( QCloseEvent* e )
{
      // overrides the default QWidget implementation
      // to delete this widget, rather than just hiding it
      // so that we don't keep on creating new dialogs while
      // just hiding the old ones.
    delete this;
}
/* 
 * protected slot
 */
//
// this exists so that the dismiss button can be connected to a slot
//
void PostPlottingDlgImpl::CloseDlg()
{
    close();              // kill the widget
}

void
PostPlottingDlgImpl::ClosePlots()
{
    for ( int i=0; i < (int)FieldListBoxCW->count(); i++ ) {
        ClosePlot( i );
    }
    CloseGKS();    
}

void
PostPlottingDlgImpl::ClosePlot( int fieldId )
{
    for ( int i=0; i<MAX_GKS_XWIN; i++ )
        if ( plotId[i][FID] == fieldId ) {
            CloseOutputDevice( plotId[i][WKSID] );
            plotId[i][FID] = -1;
            numPlots--;
            return;
        }
}

bool
PostPlottingDlgImpl::Init()
{
    NcDoubleVar hours;
    NcIntAtt stepLen;
    ShowPlotsPB->setEnabled( TRUE );
    
      // let user choose file to display
    inputFile = QFileDialog::getOpenFileName( MANAGER.HistName().c_str(), "*.nc", this );
    if ( inputFile.isNull() ) {
        close();
        return FALSE;
    }
    try {
        plotFile = new History( inputFile.data() );
    } catch ( NcErr& e ) {
        ShowMsg( __FILE__, __LINE__, "ERROR: Couldn't open history file: %s\n", e.toString().c_str() );
        close();
        return FALSE;
    }
    timeFormat = MANAGER.TimeFormat();
    TimeFormatCB->setCurrentItem( timeFormat );

      // fill in the values for start and end time
    try {

        if ( plotFile->hasVariable("bdate" ))
            bdate = plotFile->variable("bdate");
        else
            bdate = plotFile->variable("basedate");
        hours = plotFile->variable("time");

        bdate.read();
        hours.read();

        stepLen = NcIntAtt(plotFile->attribute( "time_step_length"));
    } catch ( NcErr& e ) {
        ShowMsg( __FILE__, __LINE__, "ERROR: Can't read history file: %s\n", e.toString().c_str() );
        close();
        return FALSE;
    }

      // fill in the filename label
    QFileInfo fi( inputFile );
    FileNameTL->setText( fi.fileName() ); 

      // get the list of fields that are present in the dataset
      // and fill the listbox 
    FieldListBoxCW->setAutoUpdate( FALSE );
    FieldListBoxCW->clear();
    
    FieldList fl = plotFile->PlottableFields();
    for ( FieldListItr it = fl.begin(); it != fl.end(); ++it )
        FieldListBoxCW->insertItem( (*it).first.c_str() );
    
    FieldListBoxCW->setAutoUpdate( TRUE );
    FieldListBoxCW->setCurrentItem( 0 );
    FieldListBoxCW->repaint();
    FieldListBoxCW->setFocus();

    steps_per_hour = int(3600.0/stepLen[0]);

    double minHours = hours[0];
    double maxHours = hours[hours.size()-1];

    if( maxHours == 0 ) {       // no data!
        close();
        return FALSE;
    }    
    minSteps = (int)(minHours * steps_per_hour);
    maxSteps = (int)(maxHours * steps_per_hour);


    startStep = minSteps;
    endStep = maxSteps;

    TimeConverter tc( bdate, endStep );
    EndingStepLE->setText( tc.TimeToString( timeFormat ).c_str() );
    tc.SetStep( startStep );
    StartingStepLE->setText( tc.TimeToString( timeFormat ).c_str() );

    TimeAvgCB->setChecked( FALSE );
    OpenGKS();
    return TRUE;
}
/* 
 * protected slot
 */
void PostPlottingDlgImpl::SavePlots(int btnID)
{
    char    outputfile[MAX_PATH_LEN];
    char    filter[10];
    int     outputType;         // CGM = 1 or PostScript = 20

      //
      // NCAR Graphics names postscript output files "gmetaXX.ps"
      // where XX is the postscript workstation ID.
      // It names CGM files "gmeta"
      //
    switch (btnID) {
    case CGM_BTN_ID:
        outputType = CGM;
        strcpy( outputfile, "gmeta" );
        strcpy( filter, "*.cgm" );
        break;
    case POSTSCRIPT_BTN_ID:
        outputType = POSTSCRIPT;
        sprintf( outputfile, "gmeta%d.ps", OUTPUT_FILE_ID );
        strcpy( filter, "*.ps" );
         break;
    default:
        cerr << "ERROR:"__FILE__":"<< __LINE__ 
             << "  PostPlottingDlgImpl::SavePlots(): invalid button id : " << btnID << endl;
        exit( -1 );
    }

    SetStartStep();
    SetEndStep();

    OpenOutputDevice( OUTPUT_FILE_ID, outputType );
    ActivateOutputDevice( OUTPUT_FILE_ID );

    char fieldname[256];

    for ( int i=0; i < (int)FieldListBoxCW->count(); i++ ) {
        if ( FieldListBoxCW->isSelected( i ) ) {
            strcpy( fieldname, FieldListBoxCW->text(i) );
            
            try {
                Field f( NcVariable<real_t>(plotFile->variable(fieldname)) ); 
                CreatePlot( inputFile, f.Name().c_str(), f.LongName().c_str(),
                            f.PlotUnits().c_str(), f.PlotMult(), 
                            startStep, endStep,
                            TimeFormatCB->currentItem(),
                            steps_per_hour, bdate,
                            TimeAvgCB->isChecked(), OUTPUT_FILE_ID );
                
            } catch ( NcErr& e ) {
                ShowMsg( __FILE__, __LINE__,"Can't save plot: %s\n", e.toString().c_str() );
                DeactivateOutputDevice( OUTPUT_FILE_ID );
                CloseOutputDevice( OUTPUT_FILE_ID );
                return;
            }
            
        }
    }

    DeactivateOutputDevice( OUTPUT_FILE_ID );
    CloseOutputDevice( OUTPUT_FILE_ID );
    
    QString filename =   QFileDialog::getSaveFileName( ".", filter, this );

    if ( filename.isNull() ) 
        remove( outputfile );
    else 
        rename( outputfile, filename.data() );

}
/* 
 * protected slot
 */
void PostPlottingDlgImpl::SetEndStep()
{
    TimeConverter tc( bdate, 0 );
    tc.SetTime( EndingStepLE->text().ascii(), timeFormat );
    endStep = tc.Steps();

    if ( endStep > maxSteps || endStep <= minSteps ) {
        endStep = maxSteps;
    }

    tc.SetStep( endStep );
    EndingStepLE->setText( tc.TimeToString( timeFormat).c_str() );
}
/* 
 * protected slot
 */
void PostPlottingDlgImpl::SetStartStep()
{
    TimeConverter tc( bdate, 0 );
    tc.SetTime( StartingStepLE->text().ascii(), timeFormat );

      // make sure we have a legitimate start step
    startStep = tc.Steps();
    if ( startStep < minSteps || startStep >= maxSteps) {
        startStep = minSteps;
    }
    tc.SetStep( startStep );

    StartingStepLE->setText( tc.TimeToString( timeFormat).c_str() );
}
/* 
 * protected slot
 */
void PostPlottingDlgImpl::ShowPlots()
{
    int i;
    int wksid;                  // GKS X window ID
    int status = 0;             // status returned from ncarg call
    
      //
      // check the validity of the input
      //
    SetStartStep();
    SetEndStep();
    
    
      // Plot the fields
    
      // first close any plots that were open, but have now
      // been deselected
    for ( i=0; i<(int)FieldListBoxCW->count(); i++ ) {
        if ( ! FieldListBoxCW->isSelected( i ) ) 
            ClosePlot( i );
    }
      // next, open all of the selected plots
    char fieldname[256];

    for ( i=0; i<(int)FieldListBoxCW->count(); i++ ) {
        if ( FieldListBoxCW->isSelected( i ) ) {
            strcpy( fieldname, FieldListBoxCW->text(i) );
            if ( (wksid = OpenPlot( i ) ) < 0 ) {
                ShowMsg( __FILE__, __LINE__, "Can't display field %s.\nMaximum number of open plot windows (%d) exceeded.", fieldname, MAX_GKS_XWIN );
                continue;
            }
            ActivateOutputDevice( wksid );
            ClearOutputDevice( wksid );
            try {
                  // initialize a field from the plotfile variable
                Field f( NcVariable<real_t>(plotFile->variable(fieldname) ));
                status =  CreatePlot( inputFile, f.Name().c_str(), f.LongName().c_str(),
                                      f.PlotUnits().c_str(), f.PlotMult(), 
                                      startStep, endStep,
                                      TimeFormatCB->currentItem(),
                                      steps_per_hour, bdate,
                                      TimeAvgCB->isChecked(), wksid );
            } catch ( NcErr& e ) {
                ShowMsg( __FILE__, __LINE__,"Can't show %s plot: %s\n", fieldname, e.toString().c_str() );
                DeactivateOutputDevice( wksid );
                return;
            }
            DeactivateOutputDevice( wksid );
            if ( status < 0 ) {
                ClosePlot( i );
                FieldListBoxCW->setSelected( i, FALSE );
            }
        }
    }
}
int
PostPlottingDlgImpl::OpenPlot( int fieldId )
{
    for ( int i=0; i<MAX_GKS_XWIN; i++ )
        if ( plotId[i][FID] == fieldId ) // plot is already open
            return plotId[i][WKSID];
      // get the next available GKS workstation id
    for ( int i=0; i<MAX_GKS_XWIN; i++ )
        if ( plotId[i][FID] == -1 ) {
            plotId[i][FID] = fieldId;
            plotId[i][WKSID] = nextGksId++;
            OpenOutputDevice ( plotId[i][WKSID], XWIND );
            numPlots++;
            return plotId[i][WKSID];
        }
    return -1;               
}

/* 
 * protected slot
 */
//
// This version is called by double clicking in the list of fields
//
void PostPlottingDlgImpl::ShowPlots(int)
{
    ShowPlots();
}
/* 
 * protected slot
 */
void PostPlottingDlgImpl::SetTimeFormat(int format)
{
    timeFormat = (Time) format;  
    TimeConverter tc( bdate, 0 );
    tc.SetStep( endStep );
    EndingStepLE->setText( tc.TimeToString( timeFormat ).c_str() );
    tc.SetStep( startStep );
    StartingStepLE->setText( tc.TimeToString( timeFormat ).c_str() );
}

