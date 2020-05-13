/*------------------------------------------------------------------------*
 * File: ncfile.h 
 * $Author: mvr $
 * $Id: ncfile.h,v 1.1.6.1.28.1 2006/02/01 21:29:25 mvr Exp $ *
 *------------------------------------------------------------------------*/
#ifndef NCFILE_H
#define NCFILE_H

#include <netcdf.h>
#include <string>
#include <vector>

#include <algorithm>
#include <numeric>

#include <cstdio>
#include <cstdlib>
#include <limits.h>
#include <iostream>

#include "realtype.h"
#include "ioerr.h"              // for i/o exception classes
namespace ncfile {

using namespace std;

typedef signed char byte;


//
//  NcDimension class represents a dimension of a netCDF file. It has
//    name, id, and size member variables. All netCDF files have a number 
//    of dimensions associated with them, which are used to define 
//    the shapes of the variables in netCDF files. Therefore, you can add
//    a dimension to a netCDF file, and then use that dimension to define
//    a variable. Conversely, you can get the dimensions that are defined
//    in a netCDF file, and the dimensions that a particular variable has. 
//

class NcDimension {
    friend class NcFile;
public:
      // Constructor: creates a NcDimension object with the specified name,
      //       and size. The third argument is specified to make 
      //       this an unlimited (i.e., time/record) dimension.
      //       The id member variable gets initialized when the 
      //       dimension is added to a file.
    NcDimension( const string& name, size_t size,  bool isUnlimDim = false ) 
        : _name( name ), _size( size ), _id( -1 ), _isUnlimDim( isUnlimDim ) {
        if (_size == NC_UNLIMITED) _isUnlimDim = true; /* just to be safe */ } 

      // Copy constructor: creates a copy of the NcDimension object passed as argument
      //   also used when an NcDimension object is initialized with another
    NcDimension( const NcDimension& d )
        : _name( d._name ), _size( d._size ), _id( d._id ), _isUnlimDim( d._isUnlimDim ) {}

      // Default constructor: creates a NcDimension object with size = 0,
      //  (normally this would be assigned to later)
    NcDimension()
        : _size( 0 ), _id( -1 ), _isUnlimDim( false ) {}

      // id: returns the id that was assigned to this dimension in the netcdf file
    int               id() const { return _id; }

      // isUnlimDim: returns true if this is an unlimited (i.e., time/record) dimension
    bool              isUnlimDim() const { return _isUnlimDim; }

      // size: returns the length of this dimension
    size_t            size() const { return _size; }

      // assignment operator: copies the input NcDimension object
    NcDimension&      operator=( const NcDimension& d );
    
protected:
    string            _name;    
    size_t            _size;    
    mutable int       _id;      
    mutable bool      _isUnlimDim;
};


//
// NcAttBase is the base class for NcAttribute. This class contains the 
//   metadata for an attribute. 
//
class NcAttBase {
    friend class NcFile;
    friend class NcVarBase;

public:
      // destructor:
    virtual             ~NcAttBase() {}

      // name: returns name of the attribute
    const string&       name() const { return _name; }

      // type: returns the type of a variable (e.g., NC_CHAR, NC_INT, etc )
    virtual  nc_type    type() const = 0;

      // clone: returns a pointer to a copy of the derived attribute (uses "new", so
      //        user is responsible for deletion)
    virtual NcAttBase*  clone() const = 0;

      // assignment operator
    NcAttBase&          operator=( const NcAttBase& a );

protected:
      // Note that all of the constructors are protected, i.e., this class
      // is only used internally; objects of NcAttBase should not be created 
      // by the user
    NcAttBase() {}
    NcAttBase( const NcAttBase& a ) 
        : _name( a._name ) {}
    NcAttBase( const string& name )
        : _name( name ) {}
    
    string              _name;
};
 
//
// NcAttribute class represents attributes in a netCDF file. The values of the 
//   attribute can be any one of the enumerated nc_types (i.e, char, int ,double, etc)
//   so this is a template class. To retrieve an attribute of type int, you use
//   an NcAttribute<int> object, etc. Since the amount of data associated with
//   an attribute is typically quite small (as opposed to a variable), the attribute
//   is initialized with it's value when it is created. To retrieve the attribute's 
//   value, one can use any of several methods: value() returns a pointer to the
//   value; there are overloaded subscript operators; and there are overloaded 
//   conversion operators, so that you can, for instance, pass an NcAttribute<int>
//   object to a function expecting an int, and the attribute's value will be passed. 
//   
template <class T>
class NcVariable;
 
template <class T>
class NcAttribute : public NcAttBase {
    friend class NcFile;
    template <class U> friend class NcVariable;
public:
                       NcAttribute( const string& name, size_t size, const T* value );
                         // convenience for defining text attributes
                       NcAttribute( const string& name, const char* value ); 
                       NcAttribute( const string& name, const T& value );
                       NcAttribute() {}
                       NcAttribute( const NcAttribute& a );
                       NcAttribute( const NcAttBase& a );
                       ~NcAttribute() {};
    const T*           value() const { return &_value[0]; }
    virtual nc_type    type() const; // defined by template specializations
    virtual size_t     size() const { return _value.size(); }
    virtual NcAttribute* clone() const { return new NcAttribute( *this ); }
    NcAttribute&       operator=( const NcAttribute& a );
    NcAttribute&       operator=( const NcAttBase& a );
    T                  operator[] ( int i ) const { return _value[i]; }  
    T&                 operator[] ( int i ) { return _value[i]; }
                       operator vector< T > () { return _value; }
                       operator const vector< T > () const { return _value; }
                       operator T () const { return _value[0]; }
    
protected:
    void               readValue();
    void               writeValue() const;
    vector< T >        _value; 
    void               copyValue( const NcAttBase& a );

};

 class NcFile;
class NcVarBase {
    friend class NcFile;

public:
    virtual            ~NcVarBase();
    void               addAttribute( const NcAttBase& a );
    const NcAttBase&   attribute( const string& name ) const throw (NcErr);
    const NcAttBase&   attribute( int attnum ) const throw (NcErr);
    const vector< NcAttBase* >& atts() const { return _atts; }
    const vector< NcDimension >& dims() const { return _dims; }
    const string&      name() const { return _name; } 
    bool               hasAttribute( const string& name ) const;
    void               rename( const string& name ) throw (NcErr);
    void               setAttributes( const NcAttBase *const atts[], int natts );
    void               setAttributes( const vector< NcAttBase* >& atts );
    int                id() const { return _id; }
    int                numAtts() const { return _atts.size(); }
    int                numDims() const { return _dims.size(); }
    const int*         dimids() const;
    size_t             size() const; 
    virtual nc_type    type() const { return nc_type( -1 ); }
    virtual NcVarBase* clone() const { return (NcVarBase*) 0; };
                       NcVarBase( const NcVarBase& v );
                       NcVarBase()
                           : _id( -1 ), _parent( 0 ) {}
                       NcVarBase( const string& name, const NcDimension dims[], int ndims );
                       NcVarBase( const string& name, const vector< NcDimension >& dims );
    
protected:    
    void               setDimensions( const NcDimension dims[], int ndims );
    void               setDimensions( const vector< NcDimension >& dims );
    NcVarBase&         operator=( const NcVarBase& v );

    string             _name;
      //  _id and _parent are mutable so that we can 
      //  call NcFile::addVariable() on const Variables
    mutable int        _id;
    mutable NcFile*    _parent;  

    vector< NcAttBase* > _atts;
    vector< NcDimension > _dims;
};

template <class T>
class NcVariable : public NcVarBase {
    friend class NcFile;
    friend class NcVarBase;
public:
                       NcVariable( const string& name, const vector< NcDimension >& dims );
                       NcVariable( const string& name, const NcDimension dims[], int ndims );
                       NcVariable( const NcVariable& v );
                       NcVariable( const NcVarBase& v );
                       NcVariable(); 
                       ~NcVariable() {}
    typename vector<T>::iterator begin()             { return _data.begin(); }
    typename vector<T>::const_iterator begin() const { return _data.begin(); }
    NcVariable*        clone() const        { return new NcVariable( *this ); }
    const T*           data() const         { return &_data[0]; }
    typename vector<T>::iterator end()               { return _data.end(); }
    typename vector<T>::const_iterator end() const   { return _data.end(); }
    T                  min() const;
    T                  max() const;
    T                  mean() const;
    T                  pop_back();
    void               push_back( const T& val ) { _data.push_back( val ); }
    void               push_front( const T& val ); 
    void               read() throw ( NcErr );
    void               readRecord( int index ) throw ( NcErr );
    T                  range() const;
    nc_type            type() const;  // defined by template specializations
    void               write() const throw ( NcErr );
    void               writeRecord( int index ) const throw ( NcErr );
    NcVariable&        operator=( const NcVariable& v );
    NcVariable&        operator=( const NcVarBase& v );
    NcVariable&        operator=( const T& t );
    NcVariable&        operator=( const vector< T >& t );
    const T&           operator[]( int i ) const { return _data[i]; }
    T&                 operator[]( int i ) { return _data[i]; }
    void               operator*=( const T& val );
    void               operator*=( const NcVariable& v );
    void               operator/=( const T& val );
    void               operator/=( const NcVariable& v );
    void               operator-=( const T& val );
    void               operator-=( const NcVariable& v );
    void               operator+=( const T& val );
    void               operator+=( const NcVariable& v );
                       operator T () const { return *_data.begin(); }
                       operator const T* () const { return _data[0]; }

protected:

    int                read( const size_t* start, const size_t* count );
    int                write( const size_t* start, const size_t* count ) const ;
    vector< T >        _data;
};


class NcFile {
public:
    enum OpenMode      { READ, WRITE, CREATE };
                        NcFile( const string& name, OpenMode mode=READ, bool clobber=false )
                            throw ( NcErr );
                       NcFile() {}
                       NcFile( const NcFile& f );
    virtual            ~NcFile();
    void               addAttribute( const NcAttBase& att, int varid = NC_GLOBAL ) throw ( NcErr );
    void               addAttribute( const vector< NcAttBase* >& atts, int varid = NC_GLOBAL ) throw ( NcErr );
    void               addAttribute( const NcAttBase* atts[], int natts, int varid = NC_GLOBAL ) throw ( NcErr );
    void               addDimension( const NcDimension& dim );
    void               addVariable( const NcVarBase& var ) throw ( NcErr ) ;
    const NcAttBase&   attribute( const string& name, int id=NC_GLOBAL ) throw ( NcErr );
    NcDimension        dimension( const string& name ) throw ( NcErr );
    void               endDefineMode();
    void               enterDefineMode();
    void               flush() throw ( NcErr );            
    bool               hasAttribute( const string& attName ) const;
    bool               hasDimension( const string& dimName ) const;
    bool               hasUnlimDimension() const;
    bool               hasVariable( const string& varName ) const;
    const string&      name() const { return _name; }
    int                ncid() const { return _ncid; }
    int                numAtts() const;
    int                numDims() const;    
    int                numVars() const;
    NcDimension        unlimDimension() const throw ( NcErr );
    NcVarBase          variable( const string& name ) throw ( NcErr );
    NcVarBase          variable( int id ) throw ( NcErr );
                       operator const string& () const { return _name; }
    const NcFile&      operator=( const NcFile& f );

protected:

    string      _name;          // name of opened netcdf file
    int         _status;        // error status returned from netcdf calls
    int         _ncid;          // id of netcdf file opened
    OpenMode    _mode;          // file read/write/create mode 
};

typedef NcFile::OpenMode OpenMode;

typedef vector<NcAttBase*>::iterator        att_itr;
typedef vector<NcAttBase*>::const_iterator  const_att_itr;
typedef vector<NcVarBase*>::iterator        var_itr;
typedef vector<NcVarBase*>::const_iterator  const_var_itr;
typedef vector<NcDimension>::iterator       dim_itr;
typedef vector<NcDimension>::const_iterator const_dim_itr;

typedef NcAttribute<char>         NcCharAtt;
typedef NcAttribute<byte>         NcByteAtt;
typedef NcAttribute<short>        NcShortAtt;
typedef NcAttribute<int>          NcIntAtt;
typedef NcAttribute<float>        NcFloatAtt;
typedef NcAttribute<double>       NcDoubleAtt;

typedef NcVariable<char>          NcCharVar;
typedef NcVariable<byte>          NcByteVar;
typedef NcVariable<short>         NcShortVar;
typedef NcVariable<int>           NcIntVar;
typedef NcVariable<float>         NcRealVar;
typedef NcVariable<double>        NcDoubleVar;

typedef vector<char>::iterator    NcCharVarItr;
typedef vector<byte>::iterator    NcByteVarItr;
typedef vector<short>::iterator   NcShortVarItr;
typedef vector<int>::iterator     NcIntVarItr;
typedef vector<float>::iterator   NcRealVarItr;
typedef vector<double>::iterator  NcDoubleVarItr;


//**********************************************************************
//
//  START OF ATTRIBUTE CLASS MEMBER FUNCTION DEFINITIONS
//
//**********************************************************************

template <class T>
NcAttribute<T>::NcAttribute( const string& name, size_t size, const T* value ):
    NcAttBase( name )
{ 
    for ( size_t i=0; i<size; i++ )
        _value.push_back( value[i] );
}


template <class T>
NcAttribute<T>::NcAttribute( const string& name, const char* text ):
    NcAttBase( name )
{ 
    for ( size_t i=0; i<strlen(text)+1; i++ )
        _value.push_back( text[i] );
}

template <class T>
NcAttribute<T>::NcAttribute( const string& name, const T& value ) :
    NcAttBase( name )
{ 
    _value.push_back( value );
}

template <class T>
NcAttribute<T>::NcAttribute( const NcAttBase& a ) :
    NcAttBase( a )
{ 
    copyValue( a );
}

template <class T>
NcAttribute<T>::NcAttribute( const NcAttribute& a ) :
    NcAttBase( a )
{ 
    _value = a._value;
}

template <class T>
NcAttribute<T>& NcAttribute<T>::operator=( const NcAttBase& a )
{
    NcAttBase::operator=( a );
    copyValue( a );
    return *this;
}

template <class T>
NcAttribute<T>& NcAttribute<T>::operator=( const NcAttribute& a )
{
    NcAttBase::operator=( a );
    _value = a._value;
    return *this;
}

//
// copies data from one NcAttribute to another performing type conversion 
//
template <class T>
inline void NcAttribute<T>::copyValue( const NcAttBase& source )
{
    _value.clear();
    
    switch( source.type() ) {
    case NC_CHAR:
    {
        const NcCharAtt* a = static_cast<const NcCharAtt*>( &source );
        for ( size_t i=0; i<a->size(); i++ ) 
            _value.push_back( (*a)[i] ); 
        break;
    }
    case NC_BYTE:
    {
        const NcByteAtt* a = static_cast<const NcByteAtt*>( &source );
        for ( size_t i=0; i<a->size(); i++ ) 
            _value.push_back( (*a)[i] ); 
        break;
    }
    case NC_SHORT:
    {
        const NcShortAtt* a = static_cast<const NcShortAtt*>( &source );
        for ( size_t i=0; i<a->size(); i++ ) 
            _value.push_back( (*a)[i] ); 
        break;
    }
    case NC_INT: 
    {
        const NcIntAtt* a = static_cast<const NcIntAtt*>( &source );
        for ( size_t i=0; i<a->size(); i++ ) 
            _value.push_back( (*a)[i] ); 
        break;
    }
    case NC_FLOAT:
    {
        const NcFloatAtt* a = static_cast<const NcFloatAtt*>( &source );
        for ( size_t i=0; i<a->size(); i++ ) 
            _value.push_back( (*a)[i] ); 
        break;
    }
    case NC_DOUBLE:
    {
        const NcDoubleAtt* a = static_cast<const NcDoubleAtt*>( &source );
        for ( size_t i=0; i<a->size(); i++ ) 
            _value.push_back( (*a)[i] ); 
        break;
    }
    default:    
        cerr << "ERROR: "__FILE__", line " << __LINE__
             << ":  NcAttribute::copyValue: unknown type -\n" << int(source.type()) << endl;
        exit( -1 );
    }
}

//
// Template specializations to NcAttribute::type()
//
template<>  inline nc_type NcAttribute<char>::type() const {
     return NC_CHAR;
}
template<>  inline nc_type NcAttribute<byte>::type() const {
     return NC_BYTE;
}
template<>  inline nc_type NcAttribute<short>::type() const {
     return NC_SHORT;
}
template<>  inline nc_type NcAttribute<int>::type() const {
     return NC_INT;
}
template<>  inline nc_type NcAttribute<float>::type() const {
     return NC_FLOAT;
}
template<>  inline nc_type NcAttribute<double>::type() const {
     return NC_DOUBLE;
}



//**********************************************************************
//
//  START OF VARIABLE CLASS MEMBER FUNCTION DEFINITIONS
//
//**********************************************************************

inline size_t  
NcVarBase::size() const {
      // all variables have a size of at least 1
    size_t sz = 1;

    for ( size_t i=0; i<_dims.size(); i++ ) {
          // ignore dimensions of size 0 (e.g., record dimension can be zero initially )
        if ( _dims[i].size() != 0 ) 
            sz *= _dims[i].size();
    }
    return sz;
}

template <class T>
inline NcVariable<T>&  NcVariable<T>::operator=( const NcVariable& v )
{
    NcVarBase::operator=( v );
    _data = v._data;
    return *this;
}

template <class T>
inline NcVariable<T>&  NcVariable<T>::operator=( const NcVarBase& v )
{
    NcVarBase::operator=( v );
    _data.resize( size() );
    return *this;
}


template <class T>
inline NcVariable<T>&  NcVariable<T>::operator=( const T& t )
{
    _data.clear();
    _data.push_back( t );
    return *this;
}

template <class T>
inline NcVariable<T>&  NcVariable<T>::operator=( const vector< T >& t )
{
    _data = t;
    return *this;
}


template <class T>
inline NcVariable<T>::NcVariable( const string& name, const NcDimension dims[], int ndims )
    : NcVarBase( name, dims, ndims ) 
{
    _data.resize( size() );
}

template <class T>
inline NcVariable<T>::NcVariable( const string& name, const vector< NcDimension >& dims )
    : NcVarBase( name, dims ) 
{
    _data.resize( size() );
}

template <class T>
inline NcVariable<T>::NcVariable( const NcVariable& v ) 
        : NcVarBase( v )
{ 
    _data = v._data;
}

template <class T>
inline NcVariable<T>::NcVariable( const NcVarBase& v ) 
        : NcVarBase( v )
{ 
    _data.resize( size() );
}

template <class T>
inline NcVariable<T>::NcVariable() 
{
}

// read an entire variable

template <class T>
inline void NcVariable<T>::read()  throw ( NcErr )
{
    size_t count[NC_MAX_VAR_DIMS];
    size_t start[NC_MAX_VAR_DIMS];
    
    for ( size_t i = 0; i<_dims.size(); i++ ) {
        start[i] = 0;
        count[i] = _dims[i].size();
    }
    
    int status = read( start, count );
    if ( status != NC_NOERR )
        throw NcErr(  _name.c_str(), nc_strerror( status ), _parent->ncid() );
}


// read a single record of a variable

template <class T>
inline void NcVariable<T>::readRecord( int record )  throw ( NcErr )
{
    size_t count[NC_MAX_VAR_DIMS];
    size_t start[NC_MAX_VAR_DIMS];
    
    for ( size_t i = 0; i<_dims.size(); i++ ) {
        if ( _dims[i].isUnlimDim() ) {
            start[i] = record;
            count[i] = 1;
        }
        else {
            start[i] = 0;
            count[i] = _dims[i].size();
        } 
    }
    
    int status = read( start, count );
    if ( status != NC_NOERR )
        throw NcErr(  _name.c_str(), nc_strerror( status ), _parent->ncid() );
}


// write an entire variable

template <class T>
inline void NcVariable<T>::write() const throw ( NcErr )
{
    size_t count[NC_MAX_VAR_DIMS];
    size_t start[NC_MAX_VAR_DIMS];
    
    for ( size_t i = 0; i<_dims.size(); i++ ) {
        start[i] = 0;
        count[i] = _dims[i].size();
    }
    
    int status = write( start, count );
    if ( status != NC_NOERR )
        throw NcErr(  _name.c_str(), nc_strerror( status ), _parent->ncid() );
}


// write a single record of a variable

template <class T>
inline void NcVariable<T>::writeRecord( int record ) const throw ( NcErr )
{
    size_t count[NC_MAX_VAR_DIMS];
    size_t start[NC_MAX_VAR_DIMS];
    
    for ( size_t i = 0; i<_dims.size(); i++ ) {
        if ( _dims[i].isUnlimDim() ) {
            start[i] = record;
            count[i] = 1;
        }
        else {
            start[i] = 0;
            count[i] = _dims[i].size();
        } 
    }
    
    int status = write( start, count );
    if ( status != NC_NOERR )
        throw NcErr(  _name.c_str(), nc_strerror( status ), _parent->ncid() );
}



template <class T>
inline void NcVariable<T>::operator*=( const T& val ) 
{
    for ( size_t i=0; i<val.size(); i++ )
        _data[i] *= val;
}

template <class T>
inline void NcVariable<T>::operator*=( const NcVariable& v ) 
{
    for ( size_t i=0; i<v.size(); i++ )
        _data[i] *= v[i];
}

template <class T>
inline void NcVariable<T>::operator/=( const T& val ) 
{
    for ( size_t i=0; i<val.size(); i++ )
        _data[i] /= val;
}

template <class T>
inline void NcVariable<T>::operator/=( const NcVariable& v ) 
{
    for ( size_t i=0; i<v.size(); i++ )
        _data[i] /= v[i];
}

template <class T>
inline void NcVariable<T>::operator+=( const T& val ) 
{
    for ( size_t i=0; i<val.size(); i++ )
        _data[i] += val;
}

template <class T>
inline void NcVariable<T>::operator+=( const NcVariable& v ) 
{
    for ( size_t i=0; i<v.size(); i++ )
        _data[i] += v[i];
}

template <class T>
inline void NcVariable<T>::operator-=( const T& val ) 
{
    for ( size_t i=0; i<val.size(); i++ )
        _data[i] -= val;
}

template <class T>
inline void NcVariable<T>::operator-=( const NcVariable& v ) 
{
    for ( size_t i=0; i<v.size(); i++ )
        _data[i] -= v[i];
}

template <class T>
inline T NcVariable<T>::mean() const
{
    return accumulate( _data.begin(), _data.end(), T(0) )/static_cast< double >( _data.size() );
}


template <class T>
inline T NcVariable<T>::min() const
{
    return *min_element( _data.begin(), _data.end() );
}

template <class T>
inline T NcVariable<T>::max() const
{
    return *max_element( _data.begin(), _data.end() );
}

template <class T>
inline T NcVariable<T>::range() const
{
    return max() - min();
}

template <class T>
inline void NcVariable<T>::push_front( const T& val )
{
    _data.insert( _data.begin(), 1, val );
}

template <class T>
inline T NcVariable<T>::pop_back()
{
    T back = *_data.rbegin();
    _data.pop_back();
    return back;
}

//
// template specializations for NcVariable::type()
//
template<>  inline nc_type NcVariable<char>::type() const {
     return NC_CHAR;
}
 
template<>  inline nc_type NcVariable<byte>::type() const {
     return NC_BYTE;
}
 
template<>  inline nc_type NcVariable<short>::type() const {
     return NC_SHORT;
}
 
template<>  inline nc_type NcVariable<int>::type() const {
     return NC_INT;
}
 
template<>  inline nc_type NcVariable<float>::type() const {
     return NC_FLOAT;
}
 
template<>  inline nc_type NcVariable<double>::type() const {
     return NC_DOUBLE;
}
 
//
// template specializations for NcVariable::read()
//

template<>  inline int NcVariable<char>::read( const size_t* start, const size_t* count ) {
    return  nc_get_vara_text( _parent->ncid(), _id, start, count, &_data[0] );
}

template<>  inline int NcVariable<byte>::read( const size_t* start, const size_t* count ) {
    return  nc_get_vara_schar( _parent->ncid(), _id, start, count, &_data[0] );
}

template<>  inline int NcVariable<short>::read( const size_t* start, const size_t* count ) {
    return  nc_get_vara_short( _parent->ncid(), _id, start, count, &_data[0] );
}

template<>  inline int NcVariable<int>::read( const size_t* start, const size_t* count ) {
    return  nc_get_vara_int( _parent->ncid(), _id, start, count, &_data[0] );
}

template<>  inline int NcVariable<float>::read( const size_t* start, const size_t* count ) {
    return  nc_get_vara_float( _parent->ncid(), _id, start, count, &_data[0] );
}

template<>  inline int NcVariable<double>::read( const size_t* start, const size_t* count ) {
    return  nc_get_vara_double( _parent->ncid(), _id, start, count, &_data[0] );
}

//
// template specializations for NcVariable::write()
//

template<>  inline int NcVariable<char>::write( const size_t* start, const size_t* count ) const {
    return  nc_put_vara_text( _parent->ncid(), _id, start, count, &_data[0] );
}

template<>  inline int NcVariable<byte>::write( const size_t* start, const size_t* count ) const {
    return  nc_put_vara_schar( _parent->ncid(), _id, start, count, &_data[0] );
}

template<>  inline int NcVariable<short>::write( const size_t* start, const size_t* count ) const {
    return  nc_put_vara_short( _parent->ncid(), _id, start, count, &_data[0] );
}

template<>  inline int NcVariable<int>::write( const size_t* start, const size_t* count ) const {
    return  nc_put_vara_int( _parent->ncid(), _id, start, count, &_data[0] );
}

template<>  inline int NcVariable<float>::write( const size_t* start, const size_t* count ) const {
    return  nc_put_vara_float( _parent->ncid(), _id, start, count, &_data[0] );
}

template<>  inline int NcVariable<double>::write( const size_t* start, const size_t* count ) const {
    return  nc_put_vara_double( _parent->ncid(), _id, start, count, &_data[0] );
}

} // end namespace ncfile

#endif /* NCFILE_H */
 

