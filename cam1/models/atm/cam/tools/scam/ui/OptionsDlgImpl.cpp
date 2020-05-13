#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include <qkeycode.h>
#include <qfiledlg.h>
#include <qlineedit.h>
#include <qtabdialog.h>
#include <qbuttongroup.h>
#include <qlabel.h>
#include <qscrollbar.h>
#include <qlayout.h>
#include <stdio.h>
#include <qstring.h>
#include <qpushbutton.h>

#include "manager.h"
#include "MainWndImpl.h"
#include "OptionsDlgImpl.h"
#include "msgdlg.h"
#include "ncfile.h"
#include "model.h"
#include "fieldlistbox.h"
#include "dataset.h"
#include "defaults.h"
#include "crm.h"

const int MAX_STEPLEN = 4800;
const int MIN_STEPLEN = 0;
/* 
 *  Constructs a OptionsDlgImpl which is a child of 'parent', with the 
 *  name 'name' and widget flags set to 'f' 
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
OptionsDlgImpl::OptionsDlgImpl( QWidget* parent,  const char* name, bool modal, WFlags fl )
    : OptionsDlg( parent, name, modal, fl )
{
    theMainWnd = (MainWndImpl*) parent;

      //------------------------------------------------------
      //
      // set up page one of the tab dialog - fields to save
      //
      //------------------------------------------------------

    MANAGER.SetSwitchDesc( SpecifiedSurfPropCB->text().ascii(), SRFPROP_SW );
    MANAGER.SetSwitchDesc( ForcingCB->text().ascii(), FRC3D_SW );
    MANAGER.SetSwitchDesc( RelaxationCB->text().ascii(), RELAX_SW );
    MANAGER.SetSwitchDesc( FixCamIOPCB->text().ascii(), FIX_DIV3D_SW );
    MANAGER.SetSwitchDesc( DiurnalAvgCB->text().ascii(), DIURNAL_AVG_SW );
    MANAGER.SetSwitchDesc( CRMCB->text().ascii(), CRM_SW );
    

      //------------------------------------------------------
      //
      // set up page two of the tab dialog - switches
      //
      //------------------------------------------------------

    QString s;
    
    QGroupBox* g = new QGroupBox( "Set Model Logical Switches", tab_2 );

    for ( int i=0; i<NUM_USER_SWITCHES; i++ ) {
        theSwitches[i] = new QCheckBox( g );
        theSwitches[i]->setText( s.setNum( i+1 )  );
        theSwitches[i]->setMinimumSize( theSwitches[i]->sizeHint() );
        theSwitchDesc[i] = new QLineEdit( g );
        theSwitchDesc[i]->setMinimumSize( theSwitchDesc[i]->sizeHint() );
    }
    
      // place everything into a grid layout
    QGridLayout* g2 = new QGridLayout( g, NUM_USER_SWITCHES/2 + 1, 4, 20, 10 );
    
      // split up the switches into 2 columns
    for ( int i=0; i<NUM_USER_SWITCHES/2; i++ ) {
        g2->addWidget( theSwitches[i], i+1, 0 );
        g2->addWidget( theSwitchDesc[i], i+1, 1 );
    }
    for ( int i=NUM_USER_SWITCHES/2; i<NUM_USER_SWITCHES; i++ ) {
        g2->addWidget( theSwitches[i], i-NUM_USER_SWITCHES/2+1, 2 );
        g2->addWidget( theSwitchDesc[i], i-NUM_USER_SWITCHES/2+1, 3 );
    }

      // let all the stretch be taken up by the line edits in columns 1,3
    g2->setColStretch( 1, 1 );
    g2->setColStretch( 3, 1 );
    g2->activate();

    QVBoxLayout* v = new QVBoxLayout( tab_2, 15 );
    v->addWidget( g );
    v->addSpacing( 30 );
    v->activate();

}

/*  
 *  Destroys the object and frees any allocated resources
 */
OptionsDlgImpl::~OptionsDlgImpl()
{
    // no need to delete child widgets, Qt does it all for us
}

/* 
 * public slot
 */
void OptionsDlgImpl::Init()
{
    LineEdit1->setText(MANAGER.GetDataset(ANAL).name().c_str());
    LineEdit2->setText(MANAGER.GetDataset(MODEL).name().c_str());
    LineEdit3->setText(MANAGER.GetDataset(LSMINI).name().c_str());
    LineEdit4->setText(MANAGER.GetDataset(LSMSRF).name().c_str());
    LineEdit5->setText(MANAGER.GetDataset(OZON).name().c_str());
    LineEdit6->setText(MANAGER.GetDataset(PRES).name().c_str());
    LineEdit7->setText(MANAGER.GetDataset(SST).name().c_str());
    LineEdit8->setText(MANAGER.GetDataset(ABSEMS).name().c_str());
    LineEdit9->setText(MANAGER.GetDataset(AEROPTICS).name().c_str());
    LineEdit10->setText(MANAGER.GetDataset(AERMASS).name().c_str());
    LineEdit11->setText(MANAGER.GetDataset(IOP).name().c_str());
    datasetLbl[ANAL]=LineEdit1;
    datasetLbl[MODEL]=LineEdit2;
    datasetLbl[LSMINI]=LineEdit3;
    datasetLbl[LSMSRF]=LineEdit4;
    datasetLbl[OZON]=LineEdit5;
    datasetLbl[PRES]=LineEdit6;
    datasetLbl[SST]=LineEdit7;
    datasetLbl[ABSEMS]=LineEdit8;
    datasetLbl[AEROPTICS]=LineEdit9;
    datasetLbl[AERMASS]=LineEdit10;
    datasetLbl[IOP]=LineEdit11;

    OutputFrequencySB->setValue( MANAGER.SaveFrequency() );
    ModelStepLengthSB->setValue( MANAGER.StepLen() );

    FieldListBoxCW->setMultiSelection( TRUE );
    for (int i=0; i<10; i++ )
        FieldListBoxCW->insertItem(""); // placeholders to get size right

    //    if ( FieldListBoxCW->count() > 0 )
    //        return;                 // already filled

    FieldListBoxCW->setAutoUpdate( FALSE );
      // setting autoupdate to false speeds up fill    
    FieldListBoxCW->clear();   // remove old entries
    FieldListConstItr it = MANAGER.GetFieldList().begin();
    char fieldName[100];
    for ( ; it != MANAGER.GetFieldList().end(); ++it ) {
        Field* f = (*it).second;
        if ( f->IsListed() ) {
            sprintf( fieldName, "%s - %s" , f->Name().c_str(),  f->LongName().c_str() );
              // FieldListBoxText does the coloring of modifiable fields
            FieldListBoxCW->insertItem( new FieldListBoxText( fieldName, f->IsModifiable() ) );
        }
    }

    for ( int i=0; i<(int)FieldListBoxCW->count(); i++ ) {
      Field* f = MANAGER.FindField( FieldListBoxCW->FieldName( i ) );
      if ( f->IsSaved() != FieldListBoxCW->isSelected( i ) ) {
	FieldListBoxCW->setSelected( i,TRUE ) ;
      }
    }
    FieldListBoxCW->setAutoUpdate( TRUE );
    FieldListBoxCW->setFocus();
    FieldListBoxCW->repaint();

      // set the state of the switches 
    for ( int i=0; i<NUM_USER_SWITCHES; i++) {
        theSwitches[i]->setChecked( MANAGER.SwitchState( i ) );
        theSwitchDesc[i]->setText( MANAGER.SwitchDesc( i ).c_str() );
    }

    CRMCB->setChecked( MANAGER.SwitchState( CRM_SW ) );
    SpecifiedSurfPropCB->setChecked(MANAGER.SwitchState(SRFPROP_SW));
    ForcingCB->setChecked(MANAGER.SwitchState(FRC3D_SW ));
    RelaxationCB->setChecked(MANAGER.SwitchState(RELAX_SW ));
    FixCamIOPCB->setChecked(MANAGER.SwitchState(FIX_DIV3D_SW ));
    DiurnalAvgCB->setChecked(MANAGER.SwitchState(DIURNAL_AVG_SW ));
    
}

void OptionsDlgImpl::ChooseInitialDataset( int type )
{
  ChooseDataset(type);
}
void OptionsDlgImpl::SetSeedValue()
{
  if (SeedValLE->text() != 'None')
    MANAGER.SetSeedVal(SeedValLE->text().toInt());
}
void OptionsDlgImpl::SetSeedValue(int)
{

}

void OptionsDlgImpl::ChooseBoundaryDataset( int type)
{
  type+=LSMINI;
  ChooseDataset(type);

}

void OptionsDlgImpl::ChooseIOPDataset( int type)
{
  type+=IOP;
  ChooseDataset(type);
}

/* 
 * protected slot
 */
void OptionsDlgImpl::ApplyChanges()
{

    bool reinitialize = FALSE;  // whether the model needs to be reinitialized based on 
                                // changes made to the options


    // set the state of the switches 
    for ( int i=0; i<NUM_USER_SWITCHES; i++) {
          // check whether the state has changed; if yes, need to reinitialize 
        if ( MANAGER.SwitchState( i ) != theSwitches[i]->isChecked() ) {
            MANAGER.SetSwitch( theSwitches[i]->isChecked(), i );
            reinitialize = TRUE;
        }
        MANAGER.SetSwitchDesc( theSwitchDesc[i]->text().ascii(), i );
    }
    if ( PerturbInitialConditionsCB->isChecked() != MANAGER.SwitchState( PERT_INIT_SW ) ) {
        MANAGER.SetSwitch( PerturbInitialConditionsCB->isChecked(), PERT_INIT_SW );
        reinitialize = TRUE;
    }
    if ( PerturbForcingCB->isChecked() != MANAGER.SwitchState( PERT_FRC_SW ) ) {
        MANAGER.SetSwitch( PerturbForcingCB->isChecked(), PERT_FRC_SW );
        reinitialize = TRUE;
    }
    if ( DiurnalAvgCB->isChecked() != MANAGER.SwitchState( DIURNAL_AVG_SW ) ) {
        MANAGER.SetSwitch( DiurnalAvgCB->isChecked(), DIURNAL_AVG_SW );
        reinitialize = TRUE;
    }
    if ( SpecifiedSurfPropCB->isChecked() != MANAGER.SwitchState( SRFPROP_SW ) ) {
        MANAGER.SetSwitch( SpecifiedSurfPropCB->isChecked(), SRFPROP_SW );
        reinitialize = TRUE;
    }
    if ( RelaxationCB->isChecked() != MANAGER.SwitchState( RELAX_SW ) ) {
        MANAGER.SetSwitch( RelaxationCB->isChecked(), RELAX_SW );
        reinitialize = TRUE;
    }
    if ( FixCamIOPCB->isChecked() != MANAGER.SwitchState( FIX_DIV3D_SW ) ) {
        MANAGER.SetSwitch( FixCamIOPCB->isChecked(), FIX_DIV3D_SW );
        reinitialize = TRUE;
    }
    if ( ForcingCB->isChecked() != MANAGER.SwitchState( FRC3D_SW ) ) {
        MANAGER.SetSwitch( ForcingCB->isChecked(), FRC3D_SW );
        reinitialize = TRUE;
    }
    if ( CRMCB->isChecked() != MANAGER.SwitchState( CRM_SW ) ) {
      MANAGER.SetSwitch( CRMCB->isChecked(), CRM_SW );
      reinitialize = TRUE;
      if ( CRMCB->isChecked() ) {
	MANAGER.SetEndStep( 1 );
	bool savedFieldsHaveChanged = FALSE;        
	for ( int i=0; i<(int)FieldListBoxCW->count(); i++ ) {
	  if ( FieldListBoxCW->FieldName( i ) == "FUL" || 
	       FieldListBoxCW->FieldName( i ) == "FDL" || 
	       FieldListBoxCW->FieldName( i ) == "QRS" || 
	       FieldListBoxCW->FieldName( i ) == "QRL" || 
	       FieldListBoxCW->FieldName( i ) == "FUS" || 
	       FieldListBoxCW->FieldName( i ) == "FDS" ) {
	    FieldListBoxCW->setSelected( i,TRUE ) ;
	    Field* f = MANAGER.FindField( FieldListBoxCW->FieldName( i ) );
	    if ( f->IsSaved() != FieldListBoxCW->isSelected( i ) ) {
	      savedFieldsHaveChanged = TRUE;
	      f->SetSave( FieldListBoxCW->isSelected( i ) );
	    }
	  }
	}
      }
    }
     
    MANAGER.SetShowSettings( ShowSettingCB->isChecked() );
    
      // set the model step length
    if ( MANAGER.StepLen() != ModelStepLengthSB->value() ) {
        MANAGER.SetStepLen( ModelStepLengthSB->value() );
        theMainWnd->SetEndTimeDisplay( MANAGER.EndStep() );
        reinitialize = TRUE;
    }

      // set the save frequency
    MANAGER.SetSaveFreq( OutputFrequencySB->value() );

      // set the save fields

    bool savedFieldsHaveChanged = FALSE;        
    for ( int i=0; i<(int)FieldListBoxCW->count(); i++ ) {
        Field* f = MANAGER.FindField( FieldListBoxCW->FieldName( i ) );
        if ( f->IsSaved() != FieldListBoxCW->isSelected( i ) ) {
            savedFieldsHaveChanged = TRUE;
            f->SetSave( FieldListBoxCW->isSelected( i ) );
        }
    }

    try {
          // set the datasets
      const char* entry = LineEdit1->text();
      if ( MANAGER.GetDataset( ANAL ).name() != entry ){
	MANAGER.SetDataset( entry, ANAL );
	reinitialize = TRUE;
      }
      entry = LineEdit2->text();
      if ( MANAGER.GetDataset( MODEL ).name() != entry ){
	MANAGER.SetDataset( entry, MODEL );
	reinitialize = TRUE;
      }
      entry = LineEdit3->text();
      if ( MANAGER.GetDataset( LSMINI ).name() != entry ){
	MANAGER.SetDataset( entry, LSMINI );
	reinitialize = TRUE;
      }
      entry = LineEdit4->text();
      if ( MANAGER.GetDataset( LSMSRF ).name() != entry ){
	MANAGER.SetDataset( entry, LSMSRF );
	reinitialize = TRUE;
      }
      entry = LineEdit5->text();
      if ( MANAGER.GetDataset( OZON ).name() != entry ){
	MANAGER.SetDataset( entry, OZON );
	reinitialize = TRUE;
      }
      entry = LineEdit6->text();
      if ( MANAGER.GetDataset( PRES ).name() != entry ){
	MANAGER.SetDataset( entry, PRES );
	reinitialize = TRUE;
      }
      entry = LineEdit7->text();
      if ( MANAGER.GetDataset( SST ).name() != entry ){
	MANAGER.SetDataset( entry, SST );
	reinitialize = TRUE;
      }
      entry = LineEdit8->text();
      if ( MANAGER.GetDataset( ABSEMS ).name() != entry ){
	MANAGER.SetDataset( entry,ABSEMS );
	reinitialize = TRUE;
      }
      entry = LineEdit9->text();
      if ( MANAGER.GetDataset(AEROPTICS ).name() != entry ){
	MANAGER.SetDataset( entry,AEROPTICS );
	reinitialize = TRUE;
      }
      entry = LineEdit10->text();
      if ( MANAGER.GetDataset( AERMASS ).name() != entry ){
	MANAGER.SetDataset( entry,AERMASS );
	reinitialize = TRUE;
      }
      entry = LineEdit11->text();
      if ( MANAGER.GetDataset( IOP ).name() != entry ){
	MANAGER.SetDataset( entry, IOP );
	reinitialize = TRUE;
      }
        
      // save CRM text file
      if ( CRMTextFileCB->isChecked() ) {
	CRM.SaveCRMTextFile( crmTxtFile );
	// Certain fields MUST be saved to produce a CRM text file
	Field* f;
	f =  MANAGER.FindField( "QRS" );
	f->SetSave( TRUE );
	f =  MANAGER.FindField( "QRL" );
	f->SetSave( TRUE );
	f =  MANAGER.FindField( "FUL" );
	f->SetSave( TRUE );
	f =  MANAGER.FindField( "FDL" );
	f->SetSave( TRUE );
	f =  MANAGER.FindField( "FUS" );
	f->SetSave( TRUE );
	f =  MANAGER.FindField( "FDS" );
	f->SetSave( TRUE );
	savedFieldsHaveChanged = TRUE;
      }
      
          // if the saved fields have changed, reinitialize the history file
        if( savedFieldsHaveChanged )
            MANAGER.InitHistoryFile();  
        
          // save default settings file
        if ( SetDefaultCB->isChecked() ) 
            MANAGER.SaveDefaults( DEFAULTS_FILE );
        
          // save quick-start file
        if ( CreateQSCB->isChecked() ) 
            MANAGER.SaveDefaults( quickStartFile, TRUE );
        
          // save SIC file
        if ( CreateInitialConditionsCB->isChecked() ) 
            MANAGER.SaveInitialConditions( sicFile );
        
        if ( reinitialize )
            MANAGER.InitModel( );

    } catch ( IOErr& e ) {
        ShowMsg( __FILE__, __LINE__, "ERROR: Unable to set option: %s\n", e.toString().c_str() );
    }
	  //    qWarning( "OptionsDlgImpl::ApplyChanges() not yet implemented!" ); 
}
/* 
 * protected slot
 */
//
//  close the dialog without making any changes
//
void OptionsDlgImpl::CancelChanges()
{
     hide();
}
/* 
 * protected slot
 */
void OptionsDlgImpl::ChooseDataset(int type)
{
  //    qWarning( "OptionsDlgImpl::ChooseDataset(int) not yet implemented!" ); 
    bool isValidSelection = false;
      // display file open dialog set to directory of current dataset
    QFileInfo fi( MANAGER.GetDataset( type ).name().c_str() );

    while ( ! isValidSelection ) { 
        QString qstr = QFileDialog::getOpenFileName( fi.dirPath() , "*.nc", this );
        if ( qstr.isNull() ) 
            return;
        try { 
	  if ( isValidSelection = MANAGER.IsValidDataset( qstr.data(), type ) ) 
	    datasetLbl[type]->setText( qstr );
	  else 
	    ShowMsg( __FILE__, __LINE__,"%s is not a valid %s dataset!\n", 
                        qstr.data(), Dataset::TypeDescription(type).c_str());
        } catch ( NcErr& e ) {
	  ShowMsg( __FILE__, __LINE__,"%s is not a valid dataset!\n", qstr.data(), e.toString().c_str());
        }
    }
}
/* 
 * protected slot
 */
void OptionsDlgImpl::QuickStartFileName( bool on)
{
    if ( on ) {
        QString qstr = QFileDialog::getSaveFileName( ".", "*.scm", this );
        if ( ! qstr.isNull() ) { 
            strcpy( quickStartFile, qstr );
        }
        else
            CreateQSCB->setChecked( FALSE );
    }

}
/* 
 * protected slot
 */
void OptionsDlgImpl::ShowCurrentSettings()
{
    MANAGER.ShowSettings();
}
/* 
 * protected slot
 */
void OptionsDlgImpl::SicFileName(bool on)
{
    if ( on ) {
        QFileInfo fi( MANAGER.GetDataset( SIC ).name().c_str() );
        QString qstr = QFileDialog::getSaveFileName( fi.dirPath(), "*.nc", this );
        if ( ! qstr.isNull() ) {  // isNull is True if no file selected
            strcpy( sicFile, qstr );
        }
        else
            CreateInitialConditionsCB->setChecked( FALSE );
    }
}
/* 
 * protected slot
 */
void OptionsDlgImpl::CRMTxtFileName(bool on)
{
  if ( on ) {
    QString qstr = QFileDialog::getSaveFileName( MANAGER.GetUserDataDir().data(), "*.txt", this );
    if ( ! qstr.isNull() ) {  // isNull is True if no file selected
      crmTxtFile = qstr;
    }
    else
      CRMTextFileCB->setChecked( FALSE );
  }
}

