#ifndef CHANGEAXISSCALEDLGIMPL_H
#define CHANGEAXISSCALEDLGIMPL_H
#include "realtype.h"
#include "ChangeAxisScaleDlg.h"
#include "plot.h"

class QButtonGroup;
class QLineEdit;

class ChangeAxisScaleDlgImpl : public ChangeAxisScale
{ 
    Q_OBJECT

public:
    ChangeAxisScaleDlgImpl( Plot* thePlot= 0, QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
    ~ChangeAxisScaleDlgImpl();
    void      SetAxis( Axis whichAxis ); // set the axis to modify

protected slots:
    void ApplyChanges();
    void AutoScale();
    void SetScaleType(int);

private:
    void      PropagateChanges();
    Plot*       thePlot; 
    Axis        theAxis;     // which axis to modify, X or Y


};

#endif // CHANGEAXISSCALEDLGIMPL_H
