/*------------------------------------------------------------------------*
 * File: map.h 
 * $Author: jet $
 * $Id: globalmap.h,v 1.1.6.1 2004/05/18 20:51:54 jet Exp $ *
 *------------------------------------------------------------------------*/

#ifndef Map_included
#define Map_included

#include <qframe.h>
#include <qpixmap.h>
#include <qpainter.h>
#include <realtype.h>

class SelectGlobalDataDlgImpl;
class Manager;

class Map : public QFrame
{
    friend class SelectGlobalDataDlgImpl;

public:
               Map( SelectGlobalDataDlgImpl* globalDlg, const char* name = NULL );
    void       DrawContents( QPainter* );
    void       mousePressEvent( QMouseEvent* );
    void       DrawLatLon( QPainter* painter );
    void       DrawHighlightClmn( bool highlight = TRUE );
    void       SetHighlightClmn( real_t lat, real_t lon );
    void       SetNumLatsLons( int nlat, int nlon ) { numLats = nlat,
                                                          numLons = nlon; }
private:
    
    void       paintEvent( QPaintEvent * );
    SelectGlobalDataDlgImpl* theGlobalDlg;
    QPixmap*   theMapPM;
    int        xHCol;
    int        yHCol;
    int        numLats, numLons;
    int        vertIncr;
    int        horizIncr;
    int        frameBorder;   // thickness of the frame border
    void       DrawHighlightClmn( QPainter* painter );
    void       SetIncrements();
    
};
    
#endif // Map_included

