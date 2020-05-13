
#ifndef lint
static char rcsid[] = "$Id: PlistDlgImpl.cpp,v 1.1.6.1 2004/05/18 20:51:48 jet Exp $";
#endif /* lint */


#include <stdio.h>
#include <stdlib.h>

#include "field.h"
#include "PlotDlgImpl.h"
#include "manager.h"
#include "utils.h"
#include "PlistDlgImpl.h"
#include "qlistbox.h"
#include "qlabel.h"
#include "qlineedit.h"

/* 
 *  Constructs a PlistDlgImpl which is a child of 'parent', with the 
 *  name 'name' and widget flags set to 'f' 
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
PlistDlgImpl::PlistDlgImpl( Plot* thePlot, QWidget* parent,  const char* name, bool modal, WFlags fl )
    : PlistLayout( parent, name, modal, fl )
{
    theField = thePlot->theField;    

    theField->AddObserver( this );

    char dummyString[30];


      // create space in the list

    int numEntries;
    if ( !theField->IsMultiLevel() )
        numEntries = int(3600.0*24.0/MANAGER.StepLen());
    else
        numEntries = theField->NumLevs();
    for ( int i = 0; i < numEntries; i++ )
        theValueList->insertItem( dummyString );
    if ( !theField->IsMultiLevel() ) {
        theUnitsLabel->setText( "hours" );
        theLevelLbl->hide();
    }
    else {
        if ( Plot::globalVerticalScale == PRESSURE )
            theUnitsLabel->setText( "pressure" );
        else
            theUnitsLabel->setText( "km (approx)" );
    }
    
    currentItem = 0;
    setFixedSize( size() );
    setCaption( theField->Name().c_str() );    
}

/*  
 *  Destroys the object and frees any allocated resources
 */
PlistDlgImpl::~PlistDlgImpl()
{
    // no need to delete child widgets, Qt does it all for us
    theField->DeleteObserver( this );
}

/* 
 * public slot
 */
void PlistDlgImpl::HideSelf()
{
    close( TRUE );              // destroy the widget
}

/* 
 * protected slot
 */
void PlistDlgImpl::EditPointValue(int index)
{
    if ( ! theField->IsModifiable() )
        return;
    
    if ( theField->IsMultiLevel() || index == 0 ) {
        currentItem = index;
        char number[20];
        sprintf( number, "%f", (*theField)[index] * theField->PlotMult());
        theValueEditorLE->setText( number );
        theValueEditorLE->selectAll();
        theValueEditorLE->setFocus();
    } 

}
/* 
 * protected slot
 */
void PlistDlgImpl::SetPointValue()
{
    real_t newValue = atof( theValueEditorLE->text() ) / theField->PlotMult();
    theField->SetDataPoint( newValue, currentItem );
    theValueList->repaint();
    theValueList->setFocus();
}
/* 
 * protected slot
 */
void PlistDlgImpl::UpdateValueEditor(int index)
{
    char number[20];
    currentItem = index;
    sprintf( number, "%f", (*theField)[index] * theField->PlotMult());
    theValueEditorLE->setText( number );

}
void
PlistDlgImpl::Update()
{
    char number[50];

    theValuesLabel->setText( theField->PlotUnits().c_str() );

    theValueList->setAutoUpdate( FALSE );

      // erase old entries if we're starting at step zero
      // for a surface field
    if ( MANAGER.CurrentStep() == 0 && !theField->IsMultiLevel() )
        for ( uint i=0; i<theValueList->count(); i++ )
            theValueList->changeItem( "", i );

      // fill in the new values
    real_t hours = 0;
    real_t hours_per_step = MANAGER.StepLen()/3600.0;
    
    if ( !theField->IsMultiLevel() ) {
        int numValsToShow = std::min( int(theValueList->count()), MANAGER.CurrentStep()+1 );
        for ( int i=0; i < numValsToShow; i++ ){
            sprintf( number, "         %5.1f  - %10g", hours, (*theField)[i] * theField->PlotMult() );        
            hours += hours_per_step;
            theValueList->changeItem( number, i );
        }
    }
    else {
        for ( uint i=0; i < theValueList->count(); i++ ){
            if (  Plot::globalVerticalScale == HEIGHT ) 
                sprintf( number, "%4d  - %7.2f  - %10g", i+1, PressToHeight(MANAGER.BaseLevels()[i]), (*theField)[i] );
            else
                sprintf( number, "%4d  - %7.2f  - %10g", i+1, MANAGER.BaseLevels()[i], (*theField)[i] * theField->PlotMult() );
            theValueList->changeItem( number, i );
        }

    }
    theValueList->setAutoUpdate( TRUE );
        theValueList->repaint();
}

