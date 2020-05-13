#ifndef SELECTGLOBALDATADLGIMPL_H
#define SELECTGLOBALDATADLGIMPL_H
#include "SelectGlobalDataDlg.h"

#include "realtype.h"
#include "ncfile.h"
#include <qcombobox.h>
#include <qlayout.h>

class Map;

using ncfile::NcVariable;

class SelectGlobalDataDlgImpl : public SelectGlobalDatadlg
{ 
    Q_OBJECT

public:
    SelectGlobalDataDlgImpl( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
    ~SelectGlobalDataDlgImpl();

public slots:
    void InitModel();
    void SetRunType(int);
    bool       Init();
    bool       FindNearestLatLon( real_t inLat,   real_t inLon,
                                  real_t* outLat, real_t* outLon );

protected slots:
    void SetDay(int);
    void SetHour(int);
    void SetMonth(int);
    void SetNearestLatLon();
    void UpdateDayPulldown();

private:
    Map*       theWorldMap;
    int        numDatasetLats, numDatasetLons;
    int        year, month, day, hour;
    NcVariable<real_t> lats, lons;     // initial conditions file latitudes
                                       // and longitudes variables
    void       SetLatsLons();
    QHBoxLayout* theGlobalMapFrameLayout;
};

#endif // SELECTGLOBALDATADLGIMPL_H
