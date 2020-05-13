/****************************************************************************
** Form implementation generated from reading ui file 'PostPlottingDlg.ui'
**
** Created: Wed May 16 15:25:56 2001
**      by:  The User Interface Compiler (uic)
**
** WARNING! All changes made in this file will be lost!
****************************************************************************/
#include "PostPlottingDlg.h"

#include <qbuttongroup.h>
#include <qcheckbox.h>
#include <qcombobox.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qpushbutton.h>
#include "fieldlistbox.h"
#include <qlayout.h>
#include <qvariant.h>
#include <qtooltip.h>
#include <qwhatsthis.h>

/* 
 *  Constructs a PostPlotDlg which is a child of 'parent', with the 
 *  name 'name' and widget flags set to 'f' 
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
PostPlotDlg::PostPlotDlg( QWidget* parent,  const char* name, bool modal, WFlags fl )
    : QDialog( parent, name, modal, fl )
{
    if ( !name )
	setName( "PostPlotDlg" );
    resize( 406, 490 ); 
    setCaption( tr( "Post Plotting" ) );

    TextLabel15 = new QLabel( this, "TextLabel15" );
    TextLabel15->setGeometry( QRect( 140, 260, 20, 20 ) ); 
    TextLabel15->setText( tr( "to" ) );

    TimeAvgCB = new QCheckBox( this, "TimeAvgCB" );
    TimeAvgCB->setGeometry( QRect( 20, 300, 190, 21 ) ); 
    TimeAvgCB->setText( tr( "Plot Time-Averaged Data" ) );

    ClosePB = new QPushButton( this, "ClosePB" );
    ClosePB->setGeometry( QRect( 210, 410, 190, 60 ) ); 
    ClosePB->setText( tr( "Close" ) );

    StartingStepLE = new QLineEdit( this, "StartingStepLE" );
    StartingStepLE->setGeometry( QRect( 20, 260, 110, 25 ) ); 

    EndingStepLE = new QLineEdit( this, "EndingStepLE" );
    EndingStepLE->setGeometry( QRect( 160, 260, 100, 25 ) ); 
    EndingStepLE->setSizePolicy( QSizePolicy( (QSizePolicy::SizeType)1, (QSizePolicy::SizeType)0, EndingStepLE->sizePolicy().hasHeightForWidth() ) );

    TimeFormatCB = new QComboBox( FALSE, this, "TimeFormatCB" );
    TimeFormatCB->insertItem( tr( "Steps" ) );
    TimeFormatCB->insertItem( tr( "Hours" ) );
    TimeFormatCB->insertItem( tr( "Days" ) );
    TimeFormatCB->insertItem( tr( "Date (yymmdd)" ) );
    TimeFormatCB->setGeometry( QRect( 280, 260, 110, 29 ) ); 
    QToolTip::add(  TimeFormatCB, tr( "Change the units of the time range" ) );

    ShowPlotsPB = new QPushButton( this, "ShowPlotsPB" );
    ShowPlotsPB->setGeometry( QRect( 10, 410, 200, 60 ) ); 
    ShowPlotsPB->setText( tr( "Show Plots" ) );
    ShowPlotsPB->setDefault( TRUE );
    QToolTip::add(  ShowPlotsPB, tr( "(Re)display plots of selected fields" ) );

    FieldListBoxCW = new FieldListBox( this, "FieldListBoxCW" );
    FieldListBoxCW->setGeometry( QRect( 20, 50, 370, 200 ) ); 
    QToolTip::add(  FieldListBoxCW, tr( "Select fields to display/save" ) );

    FileNameTL = new QLabel( this, "FileNameTL" );
    FileNameTL->setGeometry( QRect( 20, 10, 370, 30 ) ); 
    FileNameTL->setFrameShape( QLabel::Box );
    FileNameTL->setFrameShadow( QLabel::Plain );
    FileNameTL->setLineWidth( 2 );
    FileNameTL->setText( tr( "" ) );
    QToolTip::add(  FileNameTL, tr( "Dataset being plotted" ) );

    SavePlotBG = new QButtonGroup( this, "SavePlotBG" );
    SavePlotBG->setGeometry( QRect( 0, 320, 401, 81 ) ); 
    SavePlotBG->setFrameShape( QButtonGroup::NoFrame );
    SavePlotBG->setFrameShadow( QButtonGroup::Plain );
    SavePlotBG->setTitle( tr( "" ) );

    SavePSPB = new QPushButton( SavePlotBG, "SavePSPB" );
    SavePSPB->setGeometry( QRect( 10, 10, 390, 35 ) ); 
    SavePSPB->setText( tr( "Save Plots to PostScript File" ) );
    QToolTip::add(  SavePSPB, tr( "Save plots of selected fields to disk in Postscript format" ) );

    SaveCGMPB = new QPushButton( SavePlotBG, "SaveCGMPB" );
    SaveCGMPB->setGeometry( QRect( 10, 40, 390, 35 ) ); 
    SaveCGMPB->setText( tr( "Save Plots to CGM file" ) );
    QToolTip::add(  SaveCGMPB, tr( "Save plots of selected fields to disk in CGM format" ) );

    // signals and slots connections
    connect( StartingStepLE, SIGNAL( returnPressed() ), this, SLOT( SetStartStep() ) );
    connect( EndingStepLE, SIGNAL( returnPressed() ), this, SLOT( SetEndStep() ) );
    connect( ClosePB, SIGNAL( clicked() ), this, SLOT( CloseDlg() ) );
    connect( ShowPlotsPB, SIGNAL( clicked() ), this, SLOT( ShowPlots() ) );
    connect( SavePlotBG, SIGNAL( clicked(int) ), this, SLOT( SavePlots(int) ) );
    connect( FieldListBoxCW, SIGNAL( selected(int) ), this, SLOT( ShowPlots(int) ) );
    connect( TimeFormatCB, SIGNAL( activated(int) ), this, SLOT( SetTimeFormat(int) ) );
}

/*  
 *  Destroys the object and frees any allocated resources
 */
PostPlotDlg::~PostPlotDlg()
{
    // no need to delete child widgets, Qt does it all for us
}

void PostPlotDlg::CloseDlg()
{
    qWarning( "PostPlotDlg::CloseDlg(): Not implemented yet!" );
}

void PostPlotDlg::SavePlots(int)
{
    qWarning( "PostPlotDlg::SavePlots(int): Not implemented yet!" );
}

void PostPlotDlg::SetEndStep()
{
    qWarning( "PostPlotDlg::SetEndStep(): Not implemented yet!" );
}

void PostPlotDlg::SetStartStep()
{
    qWarning( "PostPlotDlg::SetStartStep(): Not implemented yet!" );
}

void PostPlotDlg::ShowPlots()
{
    qWarning( "PostPlotDlg::ShowPlots(): Not implemented yet!" );
}

void PostPlotDlg::ShowPlots(int)
{
    qWarning( "PostPlotDlg::ShowPlots(int): Not implemented yet!" );
}

void PostPlotDlg::SetTimeFormat(int)
{
    qWarning( "PostPlotDlg::SetTimeFormat(int): Not implemented yet!" );
}

