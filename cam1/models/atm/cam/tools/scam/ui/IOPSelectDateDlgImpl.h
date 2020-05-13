#ifndef IOPSELECTDATEDLGIMPL_H
#define IOPSELECTDATEDLGIMPL_H
#include "realtype.h"
#include "max.h"
#include "IOPSelectDateDlg.h"

class IOPSelectDateDlgImpl : public IOPDateSelection
{ 
    Q_OBJECT

public:
    IOPSelectDateDlgImpl( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
    ~IOPSelectDateDlgImpl();
    void    Init();

protected slots:
    void CommitSelection();
    void CommitSelection(int);
    void SelectFile();

private:
    void    FillDateList();     // insert dates from file into list

    char              currFilename[MAX_PATH_LEN];  // file currently being displayed
    real_t             lat, lon;           // dataset latitude and longitude

};

#endif // IOPSELECTDATEDLGIMPL_H
