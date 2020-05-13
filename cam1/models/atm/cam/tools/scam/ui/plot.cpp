#ifndef lint
static char rcsid[] = "$Id: plot.cpp,v 1.1.6.2 2004/08/19 15:05:50 jmccaa Exp $";
#endif /* lint */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include <qlayout.h>
#include <qfiledlg.h>
#include <qpixmap.h>
#include <qpainter.h>
#include <qpaintdevice.h>
#include <qtooltip.h>

#include "plot.h"
#include "PlotDlgImpl.h"
#include "manager.h"
#include "PlistDlgImpl.h"
#include "ChangeAxisScaleDlgImpl.h"
#include "field.h"
#include "model.h"
#include "utils.h"

const int   DEFAULT_POINT_SIZE = 6;
const int   TICKSIZE = 4;
const int   TICK_SPACING = 30;   // preferred number of pixels between ticks

real_t  Plot::minPressure = 0;
real_t  Plot::maxPressure = 1000;
real_t  Plot::minHeight = 0;
real_t  Plot::maxHeight = 36;
PlotScaleType Plot::globalVerticalScale = PRESSURE;
real_t  Plot::minTime = 0;
real_t  Plot::maxTime = 24;

Plot::Plot( Field*  theField, PlotDlgImpl* parent ) 
    : QWidget( (QWidget*)parent, "" )
{
    verticalScale = globalVerticalScale;
    
    this->theField = theField;
      // make the plot an observer of the field to monitor changes
    theField->AddObserver( this );
    this->parent = parent;

    pointSize = DEFAULT_POINT_SIZE;
    numPoints = 0;
    colors[0] = QColor( black );
    colors[1] = QColor( red );


    // ----------- 1 - dimensional fields -------------------
    if ( !theField->IsMultiLevel() ) {
        yMin = theField->MinHint();
        yMax = theField->MaxHint();
        xMin = minTime;
        xMax = maxTime;
    }
    // ----------- 2 - dimensional fields -------------------
    else {
        xMin = theField->MinHint();
        xMax = theField->MaxHint();
        if ( verticalScale == PRESSURE ) {
            yMin = minPressure;  
            yMax = maxPressure;
        } 
        else {
            yMin = PressToHeight( maxPressure );
            yMax = PressToHeight( minPressure );
        }
    }
    
    if ( !theField->IsMultiLevel() ) {
        pointArray.resize( MAX_POINTS );
        numPoints = 1;
    }
    else {
        numPoints = theField->NumLevs();
        pointArray.resize( numPoints );
    }
}

Plot::~Plot()
{
    theField->DeleteObserver(this);
}


//
// Automatically scales the X or Y axis to appropriate range
//     for the current data values
//
void
Plot::AutoScale( Axis axis )
{
    real_t min, max;
   
      // determine what the max and min should be from the data

    if ( !theField->IsMultiLevel() || axis == X ) {
        min = *min_element( theField->Data().begin(), theField->Data().end() ) * theField->PlotMult();
        max = *max_element( theField->Data().begin(), theField->Data().end() ) * theField->PlotMult();

          // if min == max (all values same), we need to set an arbitrary range
        if ( min == max ) {
            if ( min == 0 ) {
                min = -1.0;
                max = 1.0;
            }
            else {
                min -= min*.2*theField->PlotMult();
                max += max*.2*theField->PlotMult();
            }
        }
    }
    else {
      if ( theField->NumLevs() == MANAGER.NumLevs() ) {
        min = *std::min_element( MANAGER.BaseLevels().begin(), MANAGER.BaseLevels().end() ); 
        max = *std::max_element( MANAGER.BaseLevels().begin(), MANAGER.BaseLevels().end() ); 
      } else {
        min = *std::min_element( MANAGER.BaseILevels().begin(), MANAGER.BaseILevels().end() ); 
        max = *std::max_element( MANAGER.BaseILevels().begin(), MANAGER.BaseILevels().end() ); 
      }
    }
    
    if ( axis == X ) {
        xMin = min;
        xMax = max;
    }
    
    if ( axis == Y ) {
        if ( !theField->IsMultiLevel() || verticalScale == PRESSURE ) {
            yMin = min;
            yMax = max;
        }
        else {
            yMin = PressToHeight( max );
            yMax = PressToHeight( min );
        }
    }

    DrawPlot( this );
}

//
//
//
void
Plot::SetTransformation()
{


      //
      //  First figure out what the range of the axes will be
      //

      //  X Axis 
    xIncrement  = Nicenum( (xMax-xMin)/(PlotWidth()/TICK_SPACING), TRUE );
    xMin = floor( xMin / xIncrement ) * xIncrement;
    xMax = ceil( xMax / xIncrement ) * xIncrement;

      // Y Axis
    yIncrement  = Nicenum( (yMax - yMin) /(PlotHeight()/TICK_SPACING), TRUE );
    yMin = floor( yMin / yIncrement ) * yIncrement;
    yMax = ceil( yMax / yIncrement ) * yIncrement;

      //
      // now set the transformation matrix
      //
    transMatrix.reset();
    real_t xScaleFactor, yScaleFactor;
    xScaleFactor = PlotWidth()/(xMax - xMin);
    yScaleFactor = PlotHeight()/(yMax - yMin);
    
    if ( theField->IsMultiLevel() ) {
        if ( verticalScale == PRESSURE ) {
            transMatrix.translate( -xMin*xScaleFactor, -yMin*yScaleFactor );
            transMatrix.scale( xScaleFactor*theField->PlotMult(), yScaleFactor );
        }
        else {                 
            transMatrix.translate( -xMin*xScaleFactor, yMax*yScaleFactor );
            transMatrix.scale( xScaleFactor*theField->PlotMult(), -yScaleFactor );
        }
    }
    else {                      // single level fields
        transMatrix.translate( -xMin*xScaleFactor, yMax*yScaleFactor );    
        transMatrix.scale( xScaleFactor, -yScaleFactor*theField->PlotMult() );
    }
    inverseTransMatrix = transMatrix.invert();
}


//
// Handles mouse press events for the plot widget.
//
void
Plot::mousePressEvent( QMouseEvent *e )
{
    QPoint theMousePos = e->pos();
    SelectPoint( theMousePos );
}
    
// Handles point move events for the plot widget.
// changes the data value at a point;
// 
void 
Plot::mouseMoveEvent( QMouseEvent *theMouseEvent)
{
    MovePoint( theMouseEvent->pos() );
}


void
Plot::MovePoint( const QPoint& cursorPos )
{
    if ( MANAGER.CurrentStep() == 0 && currentPoint >= 0 && theField->IsModifiable() ) {	
        
          // update the field's data array to reflect change
          // made in the point plot;
          // the new value is in screen coordinates,
          // so we must convert back to world coordinates.
          // subtract the plot left and plot top distances
          // to make the click relative to the plot drawing area
        double worldX, worldY;
        
        inverseTransMatrix.map( double(cursorPos.x()-PlotLeft()), double(cursorPos.y()-PlotTop()), &worldX, &worldY );
        
        if ( !theField->IsMultiLevel() ) 
            theField->SetDataPoint( worldY, 0 );
        else 
            theField->SetDataPoint( worldX, currentPoint );
        DrawPoints( this );
    }
}

//
// Handles paint events for the plot widget.
//
void 
Plot::paintEvent( QPaintEvent * )
{
    DrawPlot( this );           // draw entire plot
}

void
Plot::DrawPlot( QPaintDevice* pd ) 
{    
    if (pd == this)
        erase();                

    SetTransformation();
    DrawAxes( pd );
    DrawLabels( pd );
    DrawPoints( pd );
}


//
// draw the plot axes
//
void 
Plot::DrawAxes( QPaintDevice* pd )
{
    real_t x,y;                  // position of the tick marks
    real_t xscale, yscale;       // world to screen scale factors

      // ---------- calculate placement of ticks ------
      // note that the x and y increments were calculated in 
      // SetTransformation()

      //
      //  x ticks
      //
    assert( xMax != xMin );

    xscale = PlotWidth() / ( xMax - xMin );        

    int numXticks =  int( RoundNum( ( xMax - xMin )/xIncrement ))*2 + 1;
 
    x = xMin;
    for ( int i = 0; i <= numXticks/2; i++ ) {
        xTickPos[i] = (int)( ( x - xMin) * xscale) + PlotLeft();
        if ( fabs(xIncrement) > fabs( x ) * 10000 ) // fix roundoff error around zero
            x=0;
        sprintf( xTickLabels[i], "%3g", x );
        x += xIncrement;
    }

      // add the minor ticks
    x = xMin + xIncrement/2;
    for ( int i=numXticks/2+1; i < numXticks; i++ ) {
        xTickPos[i] = (int)( ( x - xMin) * xscale) + PlotLeft();
        strcpy( xTickLabels[i], "" ); // no label
        x += xIncrement;
    }

      //
      //  y ticks
      //
    assert( yMax != yMin );
    yscale = PlotHeight() / ( yMax - yMin );
    
    int numYticks = int( RoundNum( ( yMax - yMin )/yIncrement ))*2 + 1;
    y = yMin;
    for ( int i=0; i <= numYticks/2; i++ ) {
          // ------ surface fields -------
        if ( !theField->IsMultiLevel() || verticalScale == HEIGHT ) 
            yTickPos[i] = PlotBottom() - int((y - yMin) * yscale);
          // ------ layer fields -------
        else 
            yTickPos[i] = int ((y - yMin) * yscale) + PlotTop();
        
        if ( fabs(yIncrement) > fabs( y ) * 10000 )
            y=0;
        sprintf( yTickLabels[i], "%3g", y );
        y += yIncrement;
    }        

      // add the minor ticks
    y = yMin + yIncrement/2;

    for( int i=numYticks/2+1; i<numYticks ; i++ ) {
        if ( !theField->IsMultiLevel() || verticalScale == HEIGHT ) 
            yTickPos[i] = PlotBottom() - int((y - yMin) * yscale);
          // ------ layer fields -------
        else 
            yTickPos[i] = int ((y - yMin) * yscale) + PlotTop();
        strcpy( yTickLabels[i], "" ); // no label
        y += yIncrement;
    }

    
    int textWidth = 0;          // set these to 0 so that text is
    int textHeight = 0;         // positioned accurately 
    int i;
    
    QPainter painter;
    QPen     pen( black, 2 );   // black, width 2
    
    painter.begin( pd );
    painter.setPen( pen );
    
      // draw the axis lines
      // the adjustment by 2 is to allow for the thickness of the lines
    painter.moveTo( PlotLeft() - 2, PlotTop() ); 
    painter.lineTo( PlotLeft() - 2, PlotBottom() + 2 );
    painter.lineTo( PlotRight(), PlotBottom() + 2 );
    
      // draw the tick marks and tick labels
    pen.setWidth( 1 );
    painter.setPen( pen );
    painter.setFont( QFont ( "helvetica", 10 ) );
    int format = ( AlignTop | AlignHCenter | DontClip );
      // draw the X-Axis ticks
    for ( i = 0; i < numXticks; i++ ) {
        painter.moveTo( xTickPos[i], PlotBottom() );
        painter.lineTo( xTickPos[i], PlotBottom() + TICKSIZE );
        painter.drawText( xTickPos[i], PlotBottom() + TICKSIZE,
                          textWidth, textHeight, format,
                          xTickLabels[i], strlen( xTickLabels[i] ) );
    }
      // draw the Y-Axis ticks
    format = ( AlignRight | AlignVCenter | DontClip );
    for ( i = 0; i < numYticks; i++ ) {
        painter.moveTo( PlotLeft(), yTickPos[i] );
        painter.lineTo( PlotLeft() - TICKSIZE , yTickPos[i]);
        painter.drawText( PlotLeft() - TICKSIZE - 3, yTickPos[i],
                          textWidth, textHeight,
                          format, yTickLabels[i], strlen( yTickLabels[i] ) );
    }
    painter.end();
}

void
Plot::DrawLabels( QPaintDevice* pd )
{
    int textWidth = 1;          // set these to 1 so that text is
    int textHeight = 1;         // positioned accurately 

    QPainter painter;
    painter.begin( pd );


      // draw the title
    painter.setFont( QFont ( "helvetica", 12, QFont::Bold ) );
    int format = ( AlignTop | AlignHCenter | DontClip );
    painter.drawText( PlotWidth() / 2 + PlotLeft(), PlotTop()/3, textWidth , textHeight,
                      format,  theField->LongName().c_str(), theField->LongName().size() );
    
    char label[256];
      // draw the X-axis label
    
    if ( !theField->IsMultiLevel() )
        strcpy( label, "hours" );
    else
        strcpy( label, theField->PlotUnits().c_str() );

    painter.setFont( QFont ( "helvetica", 10 ) );
    format = ( AlignTop | AlignHCenter | DontClip );
    painter.drawText( ( PlotRight() + PlotLeft() ) / 2 , 
                      PlotBottom() + 15 , textWidth, textHeight,
                      format, label, strlen( label ) ); 
    
      // draw the Y-axis label
    
    if ( !theField->IsMultiLevel() )
        strcpy( label, theField->PlotUnits().c_str() );
    else {
        if ( verticalScale == PRESSURE )
            strcpy( label, "millibars" );
        else
            strcpy( label, "km (approx)" );
    }

    painter.setFont( QFont ( "helvetica", 10 ) );
    format = ( AlignTop | AlignHCenter | DontClip );
    painter.drawText( PlotLeft(), PlotTop() - 15, textWidth, textHeight,
                      format, label, strlen( label ) );
    
    painter.end();
}

void
Plot::DrawPoints( QPaintDevice* pd ) 
{
      //
      // initialize the array of points
      //

    
    for ( int i=0; i<numPoints; i++ ) {
        double x, y;
        if ( !theField->IsMultiLevel() ) {
            double steps_per_hour = 3600.0/MANAGER.StepLen();
            transMatrix.map( i/steps_per_hour, (*theField)[i], &x, &y );
        }
        else {
	  if ( theField->NumLevs() == MANAGER.NumLevs() ) {
            if ( verticalScale == PRESSURE )
	      transMatrix.map( (*theField)[i], MANAGER.BaseLevels()[i], &x, &y );
            else
	      transMatrix.map( (*theField)[i], PressToHeight(MANAGER.BaseLevels()[i]), &x, &y );
	  } else {
            if ( verticalScale == PRESSURE )
	      transMatrix.map( (*theField)[i], MANAGER.BaseILevels()[i], &x, &y );
            else
	      transMatrix.map( (*theField)[i], PressToHeight(MANAGER.BaseILevels()[i]), &x, &y );
	  }
        }
        pointArray[i] = QPoint(x,y);
    }
    
      //  use an offscreen pixmap to avoid flickering
    QPixmap offscreenPM( PlotWidth(), PlotHeight()  );
    offscreenPM.setOptimization( QPixmap::BestOptim );
    offscreenPM.fill( colorGroup().background()  ); // fill background color
    
    QPainter painter;
    painter.begin( &offscreenPM );	
    painter.setPen( colors[0] );
    painter.drawPolyline( pointArray, 0, numPoints );	

      // if points are very small, draw 
      // boundary in  same color as center
    if( pointSize < 5 )
      painter.setPen( colors[1] ); 

    painter.setBrush( colors[1] );
    for ( int i = 0; i < numPoints; i++ ) { // draw the data points
 	painter.drawEllipse( pointArray[i].x()-pointSize/2, pointArray[i].y()-pointSize/2, 
                             pointSize, pointSize );
    }
    painter.end();
      // copy the pixmap back onto the paint device
    bitBlt( pd,                 // destination
            PlotLeft(), PlotTop(), // x,y position in destination
            &offscreenPM,       // source
            0, 0,               // x,y position in  source
            -1, -1 );           // copy whole source
            
}

//
//  handles resizing of the plot when user resizes window
//
void
Plot::resizeEvent( QResizeEvent* )
{
    pointSize = PlotHeight()/80 + PlotWidth()/80;
    pointSize = std::min( pointSize, 8 );
    pointSize = std::max( pointSize, 3 );

    DrawPlot( this );
}

//
// Find the index of a point described by x,y in
//  the pointArray, set currentPoint
//
void
Plot::SelectPoint( QPoint cursor )
{
      // make cursor position relative to plot drawing area
    cursor -= QPoint( PlotLeft(), PlotTop() );
      // max distance at which mouse press will select point
    QPoint slop( 1.5 * pointSize, 1.2 * pointSize );

    QRect grabRegion( cursor - slop, cursor + slop );

    currentPoint = -1;    
    
      // iterate over all points, looking for point close to mouse press
    for ( int i=0; i < numPoints; i++ ) {
	if ( grabRegion.contains( pointArray[i] ) )
	    currentPoint = i;
    }
}


//
// Handles mouse double click events for the Axes
//
void 
Plot::mouseDoubleClickEvent( QMouseEvent* theMouseEvent )
{
    QPoint cursorPosition = theMouseEvent->pos();

      // rects define a region near each of the axes where double click 
      // will activate the scale dialog box

    QRect yAxisRegion( PlotLeft() - 30, PlotTop(), 60, PlotHeight() );
    QRect xAxisRegion( PlotLeft(), PlotBottom() - 30 , PlotWidth(), 60 );

    if ( yAxisRegion.contains( cursorPosition ) )
        parent->ShowScaleDlg( Y );
    else if ( xAxisRegion.contains( cursorPosition ) )
        parent->ShowScaleDlg( X );
}


int
Plot::PlotWidth()
{
    return PlotRight() - PlotLeft();
}

int
Plot::PlotLeft()
{
    return int( width() * 0.15 );
}

int
Plot::PlotRight()
{
    return int( width() * 0.9 );
}

int
Plot::PlotHeight()
{
    return PlotBottom() - PlotTop();
}

int
Plot::PlotTop()
{
    return int( height() * 0.15 );
}

int
Plot::PlotBottom()
{
    return int( height() * 0.85 );
}


//
// implements Observer update actions
//
void
Plot::Update()
{
      // in the following mess of logic, the purpose is to avoid having to
      // redraw the entire plot, unless required because the global ranges have
      // changed

      // range has been changed ?
    if ( theField->IsMultiLevel() ) {
          // vertical scale has changed from pressure to height
        if ( verticalScale != globalVerticalScale ) {
            verticalScale = globalVerticalScale;
              // note the max and min will reverse when switching from pressure to height
              // and vice versa
            if ( verticalScale == PRESSURE ) {
                yMax = HeightToPress( yMin );
                yMin = HeightToPress( yMax );
            }
            if ( verticalScale == HEIGHT ) {
                yMax = PressToHeight( yMin );
                yMin = PressToHeight( yMax );
            }
            DrawPlot( this );
        }
        else {  // max or min has changed ?
            if ( verticalScale == PRESSURE ) {
                if ( yMax != maxPressure || yMin != minPressure ) {
                    yMax = maxPressure;
                    yMin = minPressure;
                    DrawPlot( this );
                }
            }
            if ( verticalScale == HEIGHT ) {
                if ( yMax != maxHeight || yMin != minHeight ) {
                    yMax = maxHeight;
                    yMin = minHeight;
                    DrawPlot( this );
                }
            }
        }
        
    }
      // the number of points to plot will change for single level fields
    if ( ! theField->IsMultiLevel() ) {
        double steps_per_hour = 3600.0/MANAGER.StepLen();
        numPoints = std::min( int((xMax-xMin)*steps_per_hour), MAX_POINTS );
        numPoints = std::min( numPoints, MANAGER.CurrentStep()+1 );
        
        if ( xMax != maxTime || xMin != minTime ) {
            xMin = minTime;
            xMax = maxTime;
            DrawPlot( this );
        }
    }
    
    
    DrawPoints( this );          // only redraw the points (quicker than redrawing the whole plot)
}    





