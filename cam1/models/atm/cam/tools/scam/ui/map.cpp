/*------------------------------------------------------------------------*
 * File: map.cpp 
 * $Author: jet $
 * $Id: map.cpp,v 1.1.6.1 2004/05/18 20:51:57 jet Exp $ *
 *------------------------------------------------------------------------*/

#ifndef lint
static char rcsid[] = "$Id: map.cpp,v 1.1.6.1 2004/05/18 20:51:57 jet Exp $";
#endif /* lint */

// File: map.cpp

#include <stdio.h>
#include <qimage.h>
#include <qlineedit.h>
#include "manager.h"
#include "globalmap.h"
#include "SelectGlobalDataDlgImpl.h"

#include "ncfile.h"
#include "map.xbm"


Map::Map(  SelectGlobalDataDlgImpl* globalDlg, const char* name )
        : QFrame( globalDlg, name )
{
    
    this->theGlobalDlg = globalDlg;
    theMapPM = new QPixmap();
    theMapPM->loadFromData( world_xbm_data, world_xbm_len, "XBM", QPixmap::Color);
    setFrameStyle( QFrame::Panel | QFrame::Raised );
    frameBorder = 2;
    setLineWidth( frameBorder );
    QPainter painter;
    painter.begin( this );	
    painter.drawPixmap( 0, 0, *theMapPM ); 
    painter.end();
    setGeometry( 10, 60, theMapPM->width() + 4, theMapPM->height() + 4 );
}


//
// inherited from QFrame
//
void
Map::DrawContents( QPainter* painter )
{

    SetIncrements();
    SetHighlightClmn( MANAGER.Lat(), MANAGER.Lon() );
    bitBlt( this, frameBorder, frameBorder, theMapPM ); 
    DrawLatLon( painter );
    DrawHighlightClmn( painter );
}

void
Map::DrawHighlightClmn( QPainter* painter )
{
    painter->setBrush( QColor( red ) );
    painter->setPen( QColor( darkGray ) );
    painter->drawRect( xHCol, yHCol, horizIncr+1, vertIncr+1 );
}


void
Map::DrawHighlightClmn( bool highlight )
{
    QPainter painter;
    painter.begin( this );
    if ( highlight ) {
        painter.setBrush( QColor( red ) );
        painter.setPen( QColor( darkGray ) );
        painter.drawRect( xHCol, yHCol, horizIncr+1, vertIncr+1 );
    }
    else // fill in the old column location with pixmap contents
        bitBlt( this, xHCol+1, yHCol+1, theMapPM, xHCol-1, yHCol-1, horizIncr-1,vertIncr - 1 );
    painter.end();
    
}

void
Map::DrawLatLon( QPainter* painter )
{
    int mapWidth = theMapPM->width();
    int mapHeight = theMapPM->height();

    painter->setPen( QColor( darkGray) );
        // draw the latitude lines
    int i;
    for ( i=frameBorder; i <= mapHeight+frameBorder; i+=vertIncr ) {
        painter->drawLine( frameBorder, i, mapWidth, i );
    }

        // draw the longitude lines
    for ( i=frameBorder; i <= mapWidth+frameBorder; i+=horizIncr ) {
        painter->drawLine(  i, frameBorder, i, mapHeight );
    }

}

void
Map::mousePressEvent( QMouseEvent* theEvent )
{ 
    // determine which lat and lon were pressed
    QPoint theMousePos = theEvent->pos();
    int xpos = theMousePos.x() - frameBorder;
    int ypos = theMousePos.y() - frameBorder;
    // ignore clicks in the frame
    if ( xpos < 0 || xpos > theMapPM->width() ||
         ypos < 0 || ypos > theMapPM->height() )
        return;
    // convert x and y positions to lat and lon
    real_t lon = ( real_t( xpos ) / theMapPM->width() ) * 360;
    if ( lon > 180 )
        lon -= 360;
    real_t lat = ( real_t ( ypos ) / theMapPM->height() ) * -180 + 90;
    // erase old column and draw new column
    DrawHighlightClmn( (bool)FALSE );
    SetHighlightClmn( lat, lon );
    DrawHighlightClmn( (bool)TRUE ); 
    // update the global dialog and the mi
    real_t nearestLat, nearestLon;
    theGlobalDlg->FindNearestLatLon( lat, lon,
                                     &nearestLat, &nearestLon );
    char latString[10], lonString[10];
    sprintf( latString, "%4.1f", nearestLat );
    sprintf( lonString, "%4.1f", nearestLon );
    theGlobalDlg->LatitudeLE->setText( latString );
    theGlobalDlg->LongitudeLE->setText( lonString );
    MANAGER.SetLat( nearestLat );
    MANAGER.SetLon( nearestLon );
}


//
//  needed to add this because map lat and lon lines weren't
//   getting redrawn properly when window was moved off the
//   edge of the screen ( only partially redrawn );
//
void
Map::paintEvent( QPaintEvent* )
{
    QPainter p;
    p.begin( this );

    DrawContents( &p );
    p.end();
}

void
Map::SetHighlightClmn( real_t lat, real_t lon )
{
        // highlight  column corresponding to the lat, lon
        // lat and lon should be in the range of -90:90, -180:180
        // respectively
    
    xHCol = int ( ( ( ( lon >= 0 ) ? lon : lon + 360 )
                   / 360 ) * numLons ) * horizIncr + frameBorder;
    yHCol = int ( ( ( lat -  90 ) /  -180 ) * numLats ) * vertIncr + frameBorder;
}
 
void
Map::SetIncrements()
{
    if ( numLats < 1 || numLons < 1 ) {
        vertIncr = theMapPM->height();
        horizIncr = theMapPM->width();
    }
    else {
        vertIncr = theMapPM->height() / numLats;
        horizIncr = theMapPM->width() / numLons;
    }
}
