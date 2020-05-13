/*------------------------------------------------------------------------*
 * File: field.h 
 * $Author: jet $
 * $Id: field.h,v 1.1.6.1 2004/05/18 20:51:53 jet Exp $ *
 *------------------------------------------------------------------------*/
#ifndef _Field_h
#define _Field_h

#include <map>
#include <deque>
#include "realtype.h"
#include "observer.h"
#include "ncfile.h"


using ncfile::NcVariable;


class Field: public Observer, public Observable
{
    friend class History;       // give History class access to the _histVar member

public:
    Field( const string& name, const string& longName, const string& plotUnits,
           const string& stdUnits, real_t plotMult, real_t  minHint,         
           real_t maxHint, bool isListed, bool isModifiable,    
           bool isAveraged, int numLevs, int fieldID );
    Field( const NcVariable<real_t>& v );
    virtual       ~Field();
    int           FieldID()      const { return _fieldID; }
    const deque<real_t>& Data()   const { return _data; }
    const NcVariable<real_t>& HistVar()    const { return *_histVar; }
    bool          IsAveraged()   const { return _isAveraged; }
    bool          IsListed()     const { return _isListed; }
    bool          IsModifiable() const { return _isModifiable; }
    bool          IsMultiLevel() const { return ( _numLevs > 1 ); }
    bool          IsSaved()      const { return _isSaved; }
    const string& LongName()     const { return _longName; }
    real_t         MaxHint()      const { return _maxHint; }
    real_t         MinHint()      const { return _minHint; }
    const string& Name()         const { return _name; }
    int           NumLevs()      const { return _numLevs; }
    real_t         PlotMult()     const { return _plotMult; }
    const string& PlotUnits()    const { return _plotUnits; }
    void print0() const {std::printf("%2f \n",_data[0]);}
    void          Reset();      // set data and plots to starting state
    void          Zero();      // set data and plots to starting state
    void          SaveRestartData(); // copy contents of data to restart data
    void          SetPlotUnits( const string& plotUnits ) { _plotUnits= plotUnits; }
    void          SetPlotMult( real_t plotMult ) { _plotMult = plotMult; }
    void          SetMaxHint( real_t maxHint ) { _maxHint = maxHint; }
    void          SetMinHint( real_t minHint ) { _minHint = minHint; }
    void          SetSave( bool isSaved ) { _isSaved = isSaved; }
    void          SetDataPoint( real_t value, int index ); // set a single data point
    void          SetData( real_t* data );
    void          ShowPlot();      // Create field's plot widget
    const string& Units() const { return _stdUnits; }
    void          Update();    // Update field's plot widget
    void          WriteHistory( int record ); // write field data to history file
    real_t         operator[]( int i ) const { return _data[i]; }
    real_t&        operator[]( int i ) { return _data[i]; }
protected:
    int           _fieldID;     // unique field number corresponding to id in model
    string        _name;        // short name for field
    string        _longName;    // descriptive field name
    string        _plotUnits;   // the units data are displayed in
    string        _stdUnits;    // the units data are saved in (SI) in history
    deque<real_t>  _data;        // data buffer (supports efficient insertion at front for 1-D fields)
    deque<real_t>  _restartData; // initial data for restarting         
    bool          _isModifiable; // Modifiable or diagnostic field?    
    bool          _isAveraged;  // TRUE if saving average over timesteps
    bool          _isListed;    // Field shows up in field list?
    bool          _isSaved;     // Field being saved?     
    real_t         _minHint;     // minimum value hint 
    real_t         _maxHint;     // maximum value hint
    real_t         _plotMult;    // multiplication factor for conversion of units in plots
    int           _lastStepPlotted;  // last step that has been plotted
    int           _numLevs;     // number of pressure levels 
    NcVariable<real_t>*     _histVar;     // for output to history file 
};

typedef map< string, Field* > FieldList;

typedef FieldList::iterator FieldListItr;
typedef FieldList::const_iterator FieldListConstItr;
typedef pair<string, Field*> fieldpair;


#endif // Field_h










