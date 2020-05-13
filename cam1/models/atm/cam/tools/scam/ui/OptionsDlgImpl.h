#ifndef OPTIONSDLGIMPL_H
#define OPTIONSDLGIMPL_H
#include "OptionsDlg.h"
#include <qbttngrp.h>
#include <qlined.h>
#include <qcheckbox.h>
#include <qspinbox.h>
#include <qtabdlg.h>
#include <qlabel.h>

#include "realtype.h"
#include "MainWndImpl.h"
#include "runtype.h"
#include "max.h"
class Manager;
class FieldListBox;

class OptionsDlgImpl : public OptionsDlg
{ 
    Q_OBJECT

public:
    OptionsDlgImpl( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
    ~OptionsDlgImpl();

public slots:
    void Init();

protected slots:
    void ApplyChanges();
    void CancelChanges();
    void ChooseDataset(int);
    void ChooseInitialDataset(int);
    void ChooseBoundaryDataset(int);
    void ChooseIOPDataset(int);
    void QuickStartFileName( bool);
    void ShowCurrentSettings();
    void SicFileName(bool);
    void CRMTxtFileName(bool);
    void SetSeedValue(int);
    void SetSeedValue();


private:

    MainWndImpl*   theMainWnd;
    QLineEdit*     theSwitchDesc[NUM_USER_SWITCHES];
    QSpinBox*      theSaveFreqSB;
    QSpinBox*      theStepLenSB;
    QCheckBox*     theSwitches[NUM_USER_SWITCHES];
    QCheckBox*     theSrfprpCB;
    QCheckBox*     theCreateQSCB;
    QCheckBox*     theCreateSICCB;
    QCheckBox*     theRelaxCB;
    QCheckBox*     theFrc3dCB;
    QCheckBox*     thePertInitCB;
    QCheckBox*     thePertFrcCB;
    QCheckBox*     theShowStartCB;
    QCheckBox*     theSaveDefaultsCB;
    FieldListBox*  theSavedFieldsLB;
    QPushButton*   chooseBtn[NUM_INIT_FILES];
    QLineEdit*        datasetLbl[NUM_INIT_FILES];
    
    char           sicFile[MAX_PATH_LEN];
    char           quickStartFile[MAX_PATH_LEN];
    QString crmTxtFile;

};

#endif // OPTIONSDLGIMPL_H
