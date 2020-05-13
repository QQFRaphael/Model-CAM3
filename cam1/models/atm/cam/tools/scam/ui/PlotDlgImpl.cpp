#ifndef lint
static char rcsid[] = "$Id: PlotDlgImpl.cpp,v 1.1.6.1 2004/05/18 20:51:49 jet Exp $";
#endif /* lint */

#include <qlayout.h>
#include <qfiledlg.h>
#include <qpainter.h>
#include <qtooltip.h>
#include <qgroupbox.h>
#include <qpixmap.h>

#include "manager.h"
#include "PlistDlgImpl.h"
#include "ChangeAxisScaleDlgImpl.h"
#include "field.h"
#include "model.h"
#include "utils.h"

#include "PlotDlgImpl.h"

/* 
 *  Constructs a PlotDlgImpl which is a child of 'parent', with the 
 *  name 'name' and widget flags set to 'f' 
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
PlotDlgImpl::PlotDlgImpl( Field*  theField, QWidget* parent,  const char* name, bool modal, WFlags fl )
    : PlotDlg( parent, name, modal, fl )
{
    thePlot = new Plot( theField, this );
    this->theField = theField;
    thePlot->setMinimumHeight( 50 );
      // Add the plot widget
      // Set vertical stretch factor to 10 to let the plot stretch
      // vertically. It will stretch horizontally because there are no
      // widgets beside it to take up horizontal stretch.
    PlotDlgLayout->setDirection( QBoxLayout::Up );
    PlotDlgLayout->insertWidget( 1,thePlot, 10 );
    setCaption( theField->Name().c_str() );

    theScaleDlg = new  ChangeAxisScaleDlgImpl( thePlot );
    theScaleDlg->hide();

}

/*  
 *  Destroys the object and frees any allocated resources
 */
PlotDlgImpl::~PlotDlgImpl()
{
    // no need to delete child widgets, Qt does it all for us
}

/* 
 * protected slot
 */
void PlotDlgImpl::AutoScale()
{
    if ( !theField->IsMultiLevel() )
        thePlot->AutoScale( Y );       
    else
        thePlot->AutoScale( X );
}
/* 
 * protected slot
 */
void PlotDlgImpl::HideSelf()
{
    close( TRUE );              // kill the widget
}
/* 
 * protected slot
 */
void PlotDlgImpl::Reset()
{
    if ( MANAGER.CurrentStep() != 0 )
        return;
    theField->Reset();
}
/* 
 * protected slot
 */
void PlotDlgImpl::Zero()
{
    if ( MANAGER.CurrentStep() != 0 )
        return;
    theField->Zero();
}
/* 
 * protected slot
 */
void PlotDlgImpl::SavePlot()
{
     QString qstr = QFileDialog::getSaveFileName( ".", "*.bmp", this );
     if ( qstr.isNull() )
         return;

     QPixmap plotPM( thePlot->size() );
     plotPM.fill( colorGroup().background()  );      // fill in background color

     thePlot->DrawPlot( (QPaintDevice*)&plotPM ); // draw plot onto the pixmap
     
     plotPM.save( qstr, "BMP" );
}
/* 
 * protected slot
 */
void PlotDlgImpl::SetScaleType()
{
    qWarning( "PlotDlgImpl::SetScaleType() not yet implemented!" ); 
}
/* 
 * protected slot
 */
void PlotDlgImpl::ShowPlist()
{
    PlistDlgImpl* p = new PlistDlgImpl( thePlot );
    theField->Update();
    p->show();
}
//
//  Popup plot scaling dialog
//
void PlotDlgImpl::ShowScaleDlg( Axis axis )
{
    theScaleDlg->SetAxis( axis );
    theScaleDlg->show();
}

