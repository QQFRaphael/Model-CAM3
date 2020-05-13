#ifndef PLOTDLGIMPL_H
#define PLOTDLGIMPL_H
#include "realtype.h"
#include "plot.h"
#include "PlotDlg.h"

class Field;
class Manager;

class PlotDlgImpl : public PlotDlg
{ 
    Q_OBJECT
    friend class ChangeAxisScaleDlgImpl;
    friend class MainWndImpl;

public:
    PlotDlgImpl(  Field*  theField = 0,QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
    ~PlotDlgImpl();

public slots:
    void        ShowScaleDlg( Axis axis );

protected slots:
    void AutoScale();
    void HideSelf();
    void Reset();
    void Zero();
    void SavePlot();
    void SetScaleType();
    void ShowPlist();

private:
    Field*      theField;
    Plot*       thePlot;
    ChangeAxisScaleDlgImpl*   theScaleDlg;
};

#endif // PLOTDLGIMPL_H
