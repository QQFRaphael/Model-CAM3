#ifndef MAINWNDIMPL_H
#define MAINWNDIMPL_H
#include "MainWnd.h"
#include <qlist.h>
#include <qtimer.h>

#include "realtype.h"
#include "observer.h"
#include "list.h"

//#include "IOPSelectDataDlgImpl.h"
class Defaults;
class IOPSelectDateDlgImpl; 
class SelectGlobalDataDlgImpl;
class LoadDataImpl;
class OptionsDlgImpl;
class PostPlottingDlgImpl;
class Model;

class MainWndImpl : public MainWnd,public Observer
{ 
    Q_OBJECT

    friend class Manager;
    friend class LoadDataImpl;
    friend class QWidget;

public:
    MainWndImpl( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 , const char* startupfile = 0);
    ~MainWndImpl();
    void EnablePostPlotting( bool enable );
    void FillFieldListBox();
    bool PostPlottingIsEnabled() { return postPlottingIsEnabled; }
    void SetCurrentTimeDisplay( int step );
    void SetEndTimeDisplay( int endStep );
    void Update();

public slots:
    void EndLEClicked(const QString&);
    void FieldSelected(int);
    void LoadBtnClicked();
    void OptionBtnClicked();
    void PostBtnClicked();
    void QuitBtnClicked();
    void RestartBtnClicked();
    void RunBtnClicked();
    void SaveBtnClicked();
    void ShowGlobalDlg();
    void ShowIopDlg();
    void StepBtnClicked();
    void TimeFormatChanged(int);
    void SetEndTimeStep(const QString&);
    void SetEndTimeStep();

 private:
    enum State         { STARTING, INITIALIZED, RUNNING, STOPPED, FINISHED };
    State              currentState;
    void               Reposition( QWidget * );
    bool               postPlottingIsEnabled;
    bool               isRunning;
    IOPSelectDateDlgImpl* IopDlg;
    LoadDataImpl* LoadDlg ;
    SelectGlobalDataDlgImpl* GlobalDlg;
    OptionsDlgImpl* OptionDlg;
    PostPlottingDlgImpl* PostDlg;
    QTimer* stepTimer;
};

#endif // MAINWNDIMPL_H
