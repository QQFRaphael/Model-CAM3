/*------------------------------------------------------------------------*
 * File: plot.h 
 * $Author: mvr $
 * $Id: plot.h,v 1.1.6.1.28.1 2006/02/01 21:29:26 mvr Exp $ *
 *------------------------------------------------------------------------*/

#ifndef Plot_included
#define Plot_included

#include <qwidget.h>
#include <qpntarry.h>
#include <qwmatrix.h>
#include <qpaintdevice.h>
#include <stdlib.h>

#include "realtype.h"
#include "observer.h"
#include "manager.h"


class Field;
class PlotDlgImpl;
class ChangeAxisScaleDlgImpl;

const int   MAX_POINTS = 504;   // max points to draw for 1-D plots
const int   MAX_TICKS = 50;     // max ticks to draw on plot axis line

enum PlotScaleType { HEIGHT, PRESSURE };
enum Axis { X, Y };

class Plot : public QWidget, public Observer
{
    Q_OBJECT

    friend class ChangeAxisScaleDlgImpl;
    friend class PlistDlgImpl;
    friend class PlotDlgImpl;

public:
                Plot( Field* field, PlotDlgImpl* parent );
                ~Plot();
    int         NumPoints() { return numPoints; }
    void        DrawPlot( QPaintDevice* pd );
    void        Update();

public slots:
    void        AutoScale( Axis a );

private:
    real_t       xMin, xMax, yMin, yMax;  // mins, maxs of the plot axes
                                         // in world coordinates
    real_t       xIncrement, yIncrement;
    int         numPoints;      // number of points to display
    QPointArray pointArray;     // point array
    int         pointSize;	// diameter of the points
    int         currentPoint;	// index of selected point
    QColor      colors[2];      // color array
    QWMatrix    transMatrix, inverseTransMatrix;
    PlotScaleType  verticalScale;

      // static members to ensure that the ranges of all the plots 
      // will be the same, i.e., if you change one plot, the change
      // will be reflected in all the plots.
    static PlotScaleType globalVerticalScale;
    static real_t  minPressure, maxPressure, minHeight, maxHeight, minTime, maxTime;

    
    void        DrawAxes( QPaintDevice* pd );
    void        DrawLabels( QPaintDevice* pd );
    void        DrawPoints( QPaintDevice* pd );
    void        MovePoint( const QPoint& theMousePos );
    int         PlotWidth();
    int         PlotLeft();
    int         PlotRight();
    int         PlotHeight();
    int         PlotTop();
    int         PlotBottom();
    void        SelectPoint( QPoint cursorPos );
    void        SetTransformation(); // world to screen transform
      // overridden Qt virtual functions
    void        paintEvent( QPaintEvent * );
    void        resizeEvent( class QResizeEvent * );
    void        mousePressEvent( QMouseEvent *);
    void        mouseMoveEvent( QMouseEvent *);
    void        mouseDoubleClickEvent( QMouseEvent* theMouseEvent );   

    ChangeAxisScaleDlgImpl*   theScaleDlg;
    Field*      theField;       // the field object that created this
    int         xTickPos[MAX_TICKS];
    int         yTickPos[MAX_TICKS];
    char        xTickLabels[MAX_TICKS][10];  
    char        yTickLabels[MAX_TICKS][10];  
    PlotDlgImpl*    parent;
};



#endif // Plot_included

