#ifndef PLISTDLGIMPL_H
#define PLISTDLGIMPL_H

#include <qlined.h>

#include "realtype.h"
#include "observer.h"
#include "PlotDlgImpl.h"
#include "field.h"
#include "manager.h"
#include "PlistDlg.h"

class PlistDlgImpl : public PlistLayout, public Observer
{ 
    Q_OBJECT

public:
    PlistDlgImpl( Plot* thePlot, QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
    ~PlistDlgImpl();

    void     Update();

public slots:
    void HideSelf();

protected slots:
    void EditPointValue(int);
    void SetPointValue();
    void UpdateValueEditor(int);
private:
    int            currentItem;
    int            numItems;          // number of items to display in list
    Field*         theField;
    QLineEdit*     theEditBox;
    const real_t*   values;            // values to display in list

};

#endif // PLISTDLGIMPL_H
