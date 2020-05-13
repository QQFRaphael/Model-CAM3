
#ifndef lint
static char rcsid[] = "$Id: ChangeAxisScaleDlgImpl.cpp,v 1.1.6.1 2004/05/18 20:51:44 jet Exp $";
#endif /* lint */

#include <stdio.h>
#include <stdlib.h>
#include <qlineedit.h>
#include <qcheckbox.h>
#include <qpushbutton.h>

#include "field.h"
#include "plot.h"
#include "manager.h"
#include "PlotDlgImpl.h"
#include "utils.h"


#include "ChangeAxisScaleDlgImpl.h"

/* 
 *  Constructs a ChangeAxisScaleDlgImpl which is a child of 'parent', with the 
 *  name 'name' and widget flags set to 'f' 
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
ChangeAxisScaleDlgImpl::ChangeAxisScaleDlgImpl( Plot* thePlot, QWidget* parent,  const char* name, bool modal, WFlags fl )
    : ChangeAxisScale( parent, name, modal, fl )
{
  this->thePlot = thePlot;
  //  connect( DismissPB, SIGNAL( clicked() ),this, SLOT( reject() ));
}

/*  
 *  Destroys the object and frees any allocated resources
 */
ChangeAxisScaleDlgImpl::~ChangeAxisScaleDlgImpl()
{
    // no need to delete child widgets, Qt does it all for us
}

/* 
 * protected slot
 */
void ChangeAxisScaleDlgImpl::ApplyChanges()
{
    if ( thePlot->theField->IsMultiLevel() ) {
        if ( theAxis == X ) {
            thePlot->xMin = atof( MinLE->text() );
            thePlot->xMax = atof( MaxLE->text() );
        }
        if ( theAxis == Y ) {
            thePlot->yMin = atof( MinLE->text() );
            thePlot->yMax = atof( MaxLE->text() );
        }
    }

    if ( ! thePlot->theField->IsMultiLevel() ) {
        if ( theAxis == X ) {
            thePlot->xMin = atof( MinLE->text() );
            thePlot->xMax = atof( MaxLE->text() );
        }
        if ( theAxis == Y ) {
            thePlot->yMin = atof( MinLE->text() );
            thePlot->yMax = atof( MaxLE->text() );
        } 
    }    
    
    if ( MultiplierLE->isEnabled() )
        thePlot->theField->SetPlotMult( atof( MultiplierLE->text() ) );
    if ( UnitsLE->isEnabled() )
        thePlot->theField->SetPlotUnits( UnitsLE->text().ascii() );

    thePlot->DrawPlot( thePlot );

    PropagateChanges();

    hide();
}
/* 
 * protected slot
 */
void ChangeAxisScaleDlgImpl::AutoScale()
{
    ApplyChanges();
    thePlot->AutoScale( theAxis );
    PropagateChanges();
    hide();
}
//
// propagates changes made in one plot to all the plots
//
void ChangeAxisScaleDlgImpl::PropagateChanges() 
{
    
    if ( thePlot->theField->IsMultiLevel() && theAxis == Y ) {
        if ( Plot::globalVerticalScale == PRESSURE ) {
            Plot::minPressure = thePlot->yMin;
            Plot::maxPressure = thePlot->yMax;
        } 
        if ( Plot::globalVerticalScale == HEIGHT ) {
            Plot::minHeight = thePlot->yMin;
            Plot::maxHeight = thePlot->yMax;
        }             
    }

    if ( ! thePlot->theField->IsMultiLevel() &&  theAxis == X ) {
            Plot::minTime = thePlot->xMin;
            Plot::maxTime = thePlot->xMax;
    }
    
      // make sure changes are propagated to all plots
    FieldListItr it = const_cast<FieldList&>(MANAGER.GetFieldList()).begin();
    for ( ; it != MANAGER.GetFieldList().end(); ++it )
        (*it).second->Update();

}

//
//  Sets which axis to edit and fills in the line edits
//       with the appropriate values.  Always call this
//       function before "show()"ing the dialog.
//
void ChangeAxisScaleDlgImpl::SetAxis( Axis axis )
{
    PressureCB->hide();
    HeightCB->hide();
    MinLE->setEnabled( TRUE );
    MultiplierLE->setEnabled( TRUE );
    UnitsLE->setEnabled( TRUE );
    AutoscalePB->setEnabled( TRUE );
    theAxis = axis;

    char str[50];
    real_t min, max;
    
    if ( theAxis == X ) {
        min = thePlot->xMin;
        max = thePlot->xMax;
    }

        // ------ Y axis ------
    else {
        min = thePlot->yMin;
        max = thePlot->yMax;
        if (  thePlot->theField->IsMultiLevel() ) {
            PressureCB->show();
            HeightCB->show();
            PressureCB->setChecked( Plot::globalVerticalScale == PRESSURE );
            HeightCB->setChecked( Plot::globalVerticalScale != PRESSURE );
        }
    }
    sprintf( str, "%5g" , min );
    MinLE->setText( str );
    sprintf( str, "%5g" , max );
    MaxLE->setText( str );
     
    if ( theAxis == X && !thePlot->theField->IsMultiLevel() ) {
        AutoscalePB->setEnabled( FALSE );
        MinLE->setEnabled( FALSE );
        MultiplierLE->setText( "N/A" );
        MultiplierLE->setEnabled( FALSE );
        UnitsLE->setText( "N/A" );
        UnitsLE->setEnabled( FALSE );
    }
    else if ( theAxis == Y && thePlot->theField->IsMultiLevel() ) {
        MultiplierLE->setText( "N/A" );
        MultiplierLE->setEnabled( FALSE );
        UnitsLE->setEnabled( FALSE );
        if ( Plot::globalVerticalScale == HEIGHT )
            UnitsLE->setText( "kilometers" );
        else
            UnitsLE->setText( "millibars" );
    }
    else {
        sprintf( str, "%5g" , thePlot->theField->PlotMult() );
        MultiplierLE->setText( str );
        UnitsLE->setText( thePlot->theField->PlotUnits().c_str() );
    }
}
/* 
 * protected slot
 */
void ChangeAxisScaleDlgImpl::SetScaleType(int type)
{
        // note that the ID's of the buttons in the dialog
        // are set to the defaults 0 and 1
        // when they are created in the constructor
    
    char str[100];

      // set the static member of Plot so that all plots will
      //  have the same scale

    Plot::globalVerticalScale = (PlotScaleType)type;

    if ( type == PRESSURE ) {
        UnitsLE->setText( "millibars" );
        thePlot->yMin = HeightToPress( atof( MaxLE->text() ) );
        thePlot->yMax = HeightToPress( atof( MinLE->text() ) );
    }
    else {
        UnitsLE->setText( "kilometers" );
        thePlot->yMin = PressToHeight( atof( MaxLE->text() ) );
        thePlot->yMax = PressToHeight( atof( MinLE->text() ) );
    }
        
    thePlot->DrawPlot( thePlot );

    sprintf( str, "%5g" , thePlot->yMin );
    MinLE->setText( str );
    sprintf( str, "%5g" , thePlot->yMax );
    MaxLE->setText( str );

    PropagateChanges();
}

