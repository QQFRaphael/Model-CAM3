/****************************************************************************
** Form interface generated from reading ui file 'PostPlottingDlg.ui'
**
** Created: Wed May 16 15:25:54 2001
**      by:  The User Interface Compiler (uic)
**
** WARNING! All changes made in this file will be lost!
****************************************************************************/
#ifndef POSTPLOTDLG_H
#define POSTPLOTDLG_H

#include <qvariant.h>
#include <qdialog.h>
class QVBoxLayout; 
class QHBoxLayout; 
class QGridLayout; 
class FieldListBox;
class QButtonGroup;
class QCheckBox;
class QComboBox;
class QLabel;
class QLineEdit;
class QPushButton;

class PostPlotDlg : public QDialog
{ 
    Q_OBJECT

public:
    PostPlotDlg( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
    ~PostPlotDlg();

    QLabel* TextLabel15;
    QCheckBox* TimeAvgCB;
    QPushButton* ClosePB;
    QLineEdit* StartingStepLE;
    QLineEdit* EndingStepLE;
    QComboBox* TimeFormatCB;
    QPushButton* ShowPlotsPB;
    FieldListBox* FieldListBoxCW;
    QLabel* FileNameTL;
    QButtonGroup* SavePlotBG;
    QPushButton* SavePSPB;
    QPushButton* SaveCGMPB;

protected slots:
    virtual void CloseDlg();
    virtual void SavePlots(int);
    virtual void SetEndStep();
    virtual void SetStartStep();
    virtual void ShowPlots();
    virtual void ShowPlots(int);
    virtual void SetTimeFormat(int);

};

#endif // POSTPLOTDLG_H
