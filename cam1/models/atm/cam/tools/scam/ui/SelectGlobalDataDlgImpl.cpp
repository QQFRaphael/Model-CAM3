#ifndef lint
static char rcsid[] = "$Id: SelectGlobalDataDlgImpl.cpp,v 1.1.6.1 2004/05/18 20:51:50 jet Exp $";
#endif /* lint */

#include <math.h>
#include <qfiledlg.h>
#include <qfileinf.h>
#include <qmsgbox.h>
#include <qpixmap.h>
#include <qstring.h>
#include <stdlib.h>
#include <stdio.h>
#include <qcombobox.h>
#include <qlineedit.h>

#include "dataset.h"
#include "globalmap.h"
#include "msgdlg.h"
#include "manager.h"
#include "runtype.h"
#include "SelectGlobalDataDlgImpl.h"

/* 
 *  Constructs a SelectGlobalDataDlgImpl which is a child of 'parent', with the 
 *  name 'name' and widget flags set to 'f' 
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
SelectGlobalDataDlgImpl::SelectGlobalDataDlgImpl( QWidget* parent,  const char* name, bool modal, WFlags fl )
    : SelectGlobalDatadlg( parent, name, modal, fl )
{
    char dateString[20];
    int i;
    
    
    theWorldMap = new Map( this );
    char* months[] = { "January", "February", "March", "April", "May", "June",
                      "July", "August", "September", "October", "November",
                      "December" };
    
 
     // fill in the date pulldowns
    MonthCB->clear();
    for ( i=0; i < 12; i++ ) {
        sprintf( dateString,"%s", months[i] );
        MonthCB->insertItem( dateString );
    }
    DayCB->clear();
    for ( i=1; i <= 31; i++ ) {
        sprintf( dateString,"%d", i );
        DayCB->insertItem( dateString );
    }
    TimeCB->clear();
    for ( i=0; i < 24; i++ ) {
        sprintf( dateString,"%d:00", i );
        TimeCB->insertItem( dateString );
    }

}

/*  
 *  Destroys the object and frees any allocated resources
 */
SelectGlobalDataDlgImpl::~SelectGlobalDataDlgImpl()
{
    // no need to delete child widgets, Qt does it all for us
}
//
// Find the latitude and longitude in the input dataset that are closest
//  to the input latitude and longitude;
//
bool
SelectGlobalDataDlgImpl::FindNearestLatLon( real_t inLat, real_t inLon,
                              real_t* outLat, real_t* outLon )
{
    real_t newVal;
    int initType;
    
    MANAGER.RunType() == MODEL ? initType=MODEL : initType=ANAL;

    if ( numDatasetLats < 1 || numDatasetLons < 1 )
        return FALSE;
        // error checking
    if ( inLat < -90  || inLat > 90 ) {
        ShowMsg( __FILE__, __LINE__,  "ERROR: Latitude is out of range in %s %s: %f",  
                 Dataset::TypeDescription( initType ).c_str(), 
                 MANAGER.GetDataset( initType ).name().c_str(), inLat );
        return FALSE;
    }
    else if ( inLon < -180 || inLon > 180 ) {
        ShowMsg( __FILE__, __LINE__,  "ERROR: Longitude is out of range in %s %s: %f",  
                 Dataset::TypeDescription( initType ).c_str(), 
                 MANAGER.GetDataset( initType ).name().c_str(), inLon );
        return FALSE;
    }

    // search for the nearest latitude
    *outLat = lats[0];
    int i = 1;
    
    while ( i < numDatasetLats ) {
        newVal = lats[i];
        if ( fabs( inLat - newVal ) < fabs( inLat - *outLat )) {
            *outLat = newVal;
        }
        i++;
    }
     
    // search for the nearest longitude
    // note that range of longitude in dataset is 0 - 360
    inLon = ( inLon > 0 ) ? inLon : inLon + 360;
    *outLon = 0;
    i = 1;
    
    while ( i < numDatasetLons ) {
        newVal = lons[i];
        if ( fabs( inLon - newVal ) < fabs( inLon - *outLon )) {
            *outLon = newVal;
        }
        i++;
    }
    *outLon = ( *outLon <= 180 ) ? *outLon : *outLon - 360;
    return TRUE;
}
bool
SelectGlobalDataDlgImpl::Init()
{
    int status = false;
    int initType;

    
    MANAGER.RunType() == MODEL ? initType=MODEL : initType=ANAL;
    try { 
        NcIntVar dates = MANAGER.DatasetIntVar( initType, "date" );
        NcIntVar secs = MANAGER.DatasetIntVar( initType, "datesec" );
        
        MANAGER.SetBaseDate( dates[0] );
        MANAGER.SetBaseSecs( secs[0] );
        
        year =  dates[0] / 10000;
        month = ( dates[0] % 10000 ) / 100;
        day =  dates[0] % 100;
        hour = secs[0] /3600;
        
          // set the date pulldowns
        MonthCB->setCurrentItem( month - 1 );
        MonthCB->repaint();
        DayCB->setCurrentItem( day -1 );
        DayCB->repaint();
        TimeCB->setCurrentItem( hour );
        TimeCB->repaint();
        
        SetLatsLons();
        status = true;
    } catch ( NcErr& e ) {
        ShowMsg( __FILE__, __LINE__, "ERROR: While reading %s: %s", 
                 Dataset::TypeDescription( initType ).c_str(),
                 e.toString().c_str() );
        status = false;
    }

    return status;
}

/* 
 * public slot
 */
void SelectGlobalDataDlgImpl::InitModel()
{
    hide();
    MANAGER.InitModel();
}
void
SelectGlobalDataDlgImpl::SetLatsLons()
{
    lats = MANAGER.DatasetRealVar( MODEL, "lat" );
    lons = MANAGER.DatasetRealVar( MODEL, "lon" );

    numDatasetLats = lats.size();
    numDatasetLons = lons.size();

    theWorldMap->SetNumLatsLons( numDatasetLats, numDatasetLons );
    
    // initialize the lat and lon line edits
    char latString[10], lonString[10];

    sprintf( latString, "%4.1f", MANAGER.Lat() );
    sprintf( lonString, "%4.1f",  MANAGER.Lon() );
    LatitudeLE->setText( latString );
    LongitudeLE->setText( lonString );
}
/* 
 * public slot
 */
void SelectGlobalDataDlgImpl::SetRunType(int)
{
    qWarning( "SelectGlobalDataDlgImpl::SetRunType(int) not yet implemented!" ); 
}

/* 
 * protected slot
 */
void SelectGlobalDataDlgImpl::SetDay(int day_)
{
     int bdate;
    day = day_ + 1;
    bdate = year * 10000 + month * 100 + day;
    MANAGER.SetBaseDate( bdate );    
}
/* 
 * protected slot
 */
void SelectGlobalDataDlgImpl::SetHour(int hour_)
{
    hour = hour_ ;
    MANAGER.SetBaseSecs( hour * 3600 );
}
/* 
 * protected slot
 */
void SelectGlobalDataDlgImpl::SetMonth(int month_)
{
    int bdate;
    month = month_ + 1;
    UpdateDayPulldown();        // call before calculating bdate
    bdate = year * 10000 + month * 100 + day;
    MANAGER.SetBaseDate( bdate );

}
/* 
 * protected slot
 */
void SelectGlobalDataDlgImpl::SetNearestLatLon()
{
    char latString[10], lonString[10];

    real_t nearestLat, nearestLon;

    real_t requestedLat = atof( LatitudeLE->text() );
    real_t requestedLon = atof( LongitudeLE->text() );

 
    if ( FindNearestLatLon( requestedLat, requestedLon,
                            &nearestLat, &nearestLon ) ) {
        sprintf( latString, "%4.1f", nearestLat );
        sprintf( lonString, "%4.1f", nearestLon );
        LatitudeLE->setText( latString );
        LongitudeLE->setText( lonString );
        MANAGER.SetLat( nearestLat );
        MANAGER.SetLon( nearestLon );
        theWorldMap->SetNumLatsLons( numDatasetLats, numDatasetLons );
        theWorldMap->DrawHighlightClmn( (bool)FALSE );  // erase previous column
        theWorldMap->SetHighlightClmn( nearestLat, nearestLon );
        theWorldMap->DrawHighlightClmn();  
    }
    else {
        // reset to old value if selection is out of range
        sprintf( latString, "%4.1f", MANAGER.Lat() );
        sprintf( lonString, "%4.1f", MANAGER.Lon() );
        LatitudeLE->setText( latString );
        LongitudeLE->setText( lonString );
    }


}
/* 
 * protected slot
 */
void SelectGlobalDataDlgImpl::UpdateDayPulldown()
{
    static int days_per_month[] = {31,28,31,30,31,30,31,31,30,31,30,31};
    char dateString[20];
    
    DayCB->clear();
    for ( int i=1; i <= days_per_month[month-1]; i++ ) {
        sprintf( dateString,"%d", i );
        DayCB->insertItem( dateString );
    }
    if ( day > days_per_month[month-1] )
        day = days_per_month[month-1];
    DayCB->setCurrentItem( day - 1 );
    DayCB->repaint();

}

