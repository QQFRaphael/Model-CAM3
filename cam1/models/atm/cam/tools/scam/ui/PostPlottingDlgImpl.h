#ifndef POSTPLOTTINGDLGIMPL_H
#define POSTPLOTTINGDLGIMPL_H
#include "PostPlottingDlg.h"
#include <qframe.h>
#include <qbttngrp.h>
#include "realtype.h"
#include "max.h"
#include "ncarg.h"
#include "timeconvert.h"
#include "history.h"

class Field;
class MainWndImpl;

class PostPlottingDlgImpl : public PostPlotDlg
{ 
    Q_OBJECT

public:
    PostPlottingDlgImpl( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
    ~PostPlottingDlgImpl();
    bool     Init();

protected slots:
    void CloseDlg();
    void SavePlots(int);
    void SetEndStep();
    void SetStartStep();
    void ShowPlots();
    void ShowPlots(int);
    void SetTimeFormat(int);

private:
    void    closeEvent( QCloseEvent * );
    void    ClosePlot( int fieldId );
    void    ClosePlots();
    int     OpenPlot( int fieldId );
    QString inputFile;
    History* plotFile;
    int     startStep;
    int     endStep;
    double  steps_per_hour;
    int     minSteps;
    int     maxSteps;

    Time     timeFormat;
    NcIntVar  bdate;
    MainWndImpl*  parent;
    static int numPlots;
    static int nextGksId;
    int     plotId[MAX_GKS_XWIN][2];
};

#endif // POSTPLOTTINGDLGIMPL_H
