/*------------------------------------------------------------------------*
 * File: history.h 
 * $Author: jet $
 * $Id: history.h,v 1.1.6.1 2004/05/18 20:51:54 jet Exp $
 *------------------------------------------------------------------------*/

#ifndef HISTORY_H
#define HISTORY_H

#include <string>
#include <vector>
#include "realtype.h"
#include "max.h"
#include "field.h"
#include "ncfile.h"

using ncfile::NcFile;
using ncfile::NcVariable;
using ncfile::NcDoubleVar;
using ncfile::NcCharVar;
using ncfile::NcIntVar;


class History: public NcFile
{

 public:
    History( const string& name, const string& case_name = "", OpenMode mode = READ, bool clobber = false );
    ~History();
    void RemoveHistVar();
    void WriteTime( double hours, int seconds, int index );
    void WriteLevs( const vector<real_t>& levels, int index );
    void WriteILevs( const vector<real_t>& levels, int index );
    FieldList PlottableFields();

private:
    NcDoubleVar timeVar; 
    NcVariable<real_t> presVar, ipresVar, latVar, lonVar, levVar, ilevVar;
    NcCharVar swDescVar;
    NcIntVar tsecVar, switchVar,  bdateVar;
};

#endif // HISTORY_H












