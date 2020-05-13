#ifndef lint
static char rcsid[] = "$Id: IOPSelectDateDlgImpl.cpp,v 1.1.6.1 2004/05/18 20:51:45 jet Exp $";
#endif /* lint */


#include <stdio.h>
#include <stdlib.h>
#include <limits.h>          
#include <qfiledlg.h>
#include <qlabel.h>
#include <qstring.h>
#include <qlistbox.h>
#include <qtooltip.h>

#include "dataset.h"
#include "manager.h"
#include "msgdlg.h"
#include "ncfile.h"
#include "timeconvert.h"
#include "IOPSelectDateDlgImpl.h"
using namespace ncfile;
/* 
 *  Constructs a IOPSelectDateDlgImpl which is a child of 'parent', with the 
 *  name 'name' and widget flags set to 'f' 
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
IOPSelectDateDlgImpl::IOPSelectDateDlgImpl( QWidget* parent,  const char* name, bool modal, WFlags fl )
    : IOPDateSelection( parent, name, modal, fl )
{
}

/*  
 *  Destroys the object and frees any allocated resources
 */
IOPSelectDateDlgImpl::~IOPSelectDateDlgImpl()
{
    // no need to delete child widgets, Qt does it all for us
}

/* 
 * protected slot
 */
void IOPSelectDateDlgImpl::CommitSelection()
{
    int index = SelectDateLB->currentItem();
    CommitSelection( index );
}
void IOPSelectDateDlgImpl::CommitSelection( int index )
{
    int    outDate;             // date converted from seconds
    int    outSecs;             // seconds after conversion

    try {
        NcIntVar tsec = MANAGER.DatasetIntVar( MANAGER.RunType(), "tsec" );
        NcIntVar bdate = MANAGER.DatasetIntVar( MANAGER.RunType(), "bdate" );
        NcVariable<real_t> lat = MANAGER.DatasetRealVar( MANAGER.RunType(), "lat" );
        NcVariable<real_t> lon = MANAGER.DatasetRealVar( MANAGER.RunType(), "lon" );
/*  
 * If this is a cam generated IOP file then adjust longitude to be
 * -180 to 180
 */
        if ( MANAGER.GetDataset( MANAGER.RunType()).hasAttribute( "CAM_GENERATED_FORCING" )){
//        NcCharAtt camiop = MANAGER.DatasetCharAtt( MANAGER.RunType(), "CAM_GENERATED_FORCING" );
//	if (string(camiop.value()).compare("create SCAM IOP dataset")== 0) {
	  for (size_t i=0;i<lon.size(); i++ )
            if ( lon[i] > 180.0 )
	      lon[i] -= 360.0;
	}
        TimeConverter::SecondsToDate( tsec[index], bdate, outDate, outSecs ); 
        MANAGER.SetBaseDate( outDate );
        MANAGER.SetBaseSecs( outSecs );
        MANAGER.SetIopStartOffset(  tsec[index] - tsec[0] );
        MANAGER.SetLat( lat );
        MANAGER.SetLon( lon );
        MANAGER.SetIopMaxSteps(); 
        MANAGER.SetEndStep( MANAGER.MaxStep() ); // ensures valid endstep
    } catch ( NcErr& e ) {
        ShowMsg( __FILE__, __LINE__, "ERROR: Unable to set dataset:\n%s\n", e.toString().c_str() );
        return;
    }
    MANAGER.InitModel();
    hide();
}

void IOPSelectDateDlgImpl::FillDateList()
{
    int    outDate;             // date converted from seconds
    int    outSecs;             // seconds after conversion
    char   dateString[100];     // string form of date

    NcIntVar tsec, bdate;

    try {
        tsec = MANAGER.DatasetIntVar( MANAGER.RunType(), "tsec" );
        bdate = MANAGER.DatasetIntVar( MANAGER.RunType(), "bdate" );
    } catch ( NcErr& e ) {
        ShowMsg( __FILE__, __LINE__, "ERROR: Can't get variable:\n%s", e.toString().c_str() );
        return;
    }

    SelectDateLB->clear();       // remove all of the old items
    SelectDateLB->setAutoUpdate( FALSE ); // to speed up inserts into list

    for ( unsigned i = 0; i < tsec.size(); i++ ) {
        TimeConverter::SecondsToDate( tsec[i], bdate, outDate, outSecs ); 
        sprintf( dateString, "%4d - %d %.2d:%.2d",
                 i, outDate , outSecs/3600, (int)(outSecs%3600) / 60 );
        SelectDateLB->insertItem( dateString );
          // set the current item based on the basedate currently set in the model
        if ( MANAGER.BaseDate() > outDate || 
             (  MANAGER.BaseDate() == outDate && MANAGER.BaseSecs() >= outSecs ) )
            SelectDateLB->setCurrentItem( i );
    }
      // make sure the current item has been set
      // if  MANAGER.BaseDate() is less than any date in the dataset
      // it won't have gotten set just above.
    if ( SelectDateLB->currentItem() < 0 )
        SelectDateLB->setCurrentItem( 0 );
    SelectDateLB->centerCurrentItem();
    SelectDateLB->setAutoUpdate( TRUE );
    SelectDateLB->repaint();
}

void IOPSelectDateDlgImpl::Init()
{
    // either IOP or USER 
    
    while ( ! MANAGER.IsValidDataset( MANAGER.GetDataset( MANAGER.RunType() ), MANAGER.RunType() )) 
        SelectFile();
    
      // use a QFileInfo object to easily obtain basename
    QFileInfo fi( MANAGER.GetDataset( MANAGER.RunType() ).name().c_str() );
    FilenameLbl->setText( fi.fileName() );
    FillDateList();
}
/* 
 * protected slot
 */
void IOPSelectDateDlgImpl::SelectFile()
{
  //    qWarning( "IOPSelectDateDlgImpl::SelectFile() not yet implemented!" ); 
    bool isValidSelection = false;

    while ( !isValidSelection ) { 
        QFileInfo fi( MANAGER.GetDataset( MANAGER.RunType() ).name().c_str() );
        QString qstr = QFileDialog::getOpenFileName( fi.dirPath(), "*.nc", this );
        if ( qstr.isNull() )   // isNull is True if no file selected
            return;
        if ( isValidSelection = MANAGER.IsValidDataset( qstr.data(), MANAGER.RunType() ) ) {
            MANAGER.SetDataset( qstr.data(), MANAGER.RunType() );
            fi.setFile( qstr );  
            FilenameLbl->setText( fi.fileName() );
            Init();
            SelectDateLB->setCurrentItem( 0 );
            SelectDateLB->centerCurrentItem();
        }
        else 
            ShowMsg( __FILE__, __LINE__,"ERROR: %s is not a valid %s dataset", qstr.data(), 
                    Dataset::TypeDescription(MANAGER.RunType()).c_str());
    }
}

