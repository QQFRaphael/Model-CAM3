#ifndef LOADDATAIMPL_H
#define LOADDATAIMPL_H
#include "LoadData.h"

class MainWndImpl;
class Manager;

class LoadDataImpl : public SelectDataForm
{ 
    Q_OBJECT

public:
    LoadDataImpl( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
    ~LoadDataImpl();

private:
    MainWndImpl*  theGUI;

protected:
    void SetRunType(int type);


};

#endif // LOADDATAIMPL_H
