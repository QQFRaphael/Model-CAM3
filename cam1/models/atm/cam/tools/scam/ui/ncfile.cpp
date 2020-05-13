/*------------------------------------------------------------------------*
 * File: ncfile.cpp 
 * $Author: jet $
 * $Id: ncfile.cpp,v 1.1.6.1 2004/05/18 20:52:01 jet Exp $ *
 *------------------------------------------------------------------------*/

#ifndef lint
static char rcsid[] = "$Id: ncfile.cpp,v 1.1.6.1 2004/05/18 20:52:01 jet Exp $";
#endif /* lint */

#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "ncfile.h"


namespace ncfile {
    
//
//: Constructor: opens a netcdf file
//  Inputs: 
//     file - name of the netcdf file;

NcFile::NcFile( const string& name, OpenMode mode, bool clobber ) throw ( NcErr )
    : _name( name ), _mode( mode )
{
    int status;
    switch ( _mode ) {
    case READ:
        status = nc_open( _name.c_str(), NC_NOWRITE, &_ncid );
        break;
    case WRITE:
        status = nc_open( _name.c_str(), NC_WRITE, &_ncid );
        break;
    case CREATE:
        status = nc_create( _name.c_str(), 
                            ( clobber ? NC_CLOBBER : NC_NOCLOBBER ), &_ncid );
        break;
    default:
        cerr << "ERROR: "__FILE__":" << __LINE__
             << ": NcFile::NcFile() - invalid open mode: " 
             << int(_mode) << endl;
        exit( -1 );
        
    }
    if ( status != NC_NOERR ) 
        throw NcErr( _name.c_str(), nc_strerror( status ), 
                     -1, __FILE__, __LINE__ );

}


//
//  destructor
//    closes the netcdf file 
//
NcFile::~NcFile()
{
    nc_close( _ncid );
}

//
//  Add a vector of attributes to the file 
//
void 
NcFile::addAttribute( const vector< NcAttBase* >& atts, int varid ) throw ( NcErr )
{
    for ( const_att_itr it = atts.begin(); it != atts.end(); ++it ) 
        addAttribute( *(*it), varid );
}

//
//  Add a single attribute to the file 
//
void 
NcFile::addAttribute( const NcAttBase& att, int varid ) throw ( NcErr )
{
    int status;    
    switch( att.type() ) {
    case NC_CHAR:
    {
        const NcCharAtt* a = dynamic_cast< const NcCharAtt* >(&att);
        status = nc_put_att_text( _ncid, varid, a->name().c_str(), 
                                  a->size(), a->value() );
    }
    break;
    case NC_BYTE:
    {
        const NcByteAtt* a = dynamic_cast< const NcByteAtt* >(&att);
        status = nc_put_att_schar( _ncid, varid, a->name().c_str(), 
                                   a->type(), a->size(), a->value() );
    }
    break;
    case NC_SHORT:
    {
        const NcShortAtt* a = dynamic_cast< const NcShortAtt* >(&att);
        status = nc_put_att_short( _ncid, varid, a->name().c_str(), 
                                   a->type(), a->size(), a->value() );
    }
    break;
    case NC_INT:
    {
        const NcIntAtt* a = dynamic_cast< const NcIntAtt* >(&att);
        status = nc_put_att_int( _ncid, varid, a->name().c_str(), 
                                 a->type(), a->size(), a->value() );
    }
    break;
    case NC_FLOAT:
    {
        const NcFloatAtt* a = dynamic_cast< const NcFloatAtt* >(&att);
        status = nc_put_att_float( _ncid, varid, a->name().c_str(), 
                                   a->type(), a->size(), a->value() );
    } 
    break;
    case NC_DOUBLE:
    {
        const NcDoubleAtt* a = dynamic_cast< const NcDoubleAtt* >(&att);
        status = nc_put_att_double( _ncid, varid, a->name().c_str(), 
                                    a->type(), a->size(), a->value() );
    }
    break;
    default:
        cerr << "ERROR: "__FILE__":"<< __LINE__ 
             << ":  NcFile::addAttribute() - invalid type: " 
             << int(att.type()) << endl;
        exit( -1 );
    }
    if ( status != NC_NOERR )
        throw NcErr( name().c_str(), (&att)->name().c_str(), nc_strerror( status ), 
                     _ncid, __FILE__, __LINE__ );
    
}

void
NcFile::addVariable( const NcVarBase& var ) throw ( NcErr )
{
    var._parent = this;
    int status = 0;
    if (( status = nc_def_var( _ncid, var.name().c_str(), var.type(), var.numDims(),
                                var.dimids(), &var._id )) != NC_NOERR ) {
        cerr << "ERROR: "__FILE__":"<< __LINE__ 
             << " NcFile::AddVariable() : " 
             << name() << ", " << var.name() << ", "
             << nc_strerror( status ) << endl;
        exit( -1 );
    }
    
    addAttribute( var.atts() , var.id() );
}

// 
// numDims:
//  return the number of dimensions used by the file
//
int NcFile::numDims() const
{
    int ndims, nvars, natts, unlim;
    nc_inq( _ncid, &ndims, &nvars, &natts, &unlim );
    return ndims;
}


// 
// numVars:
//  return the number of variables defined in the file
//
int NcFile::numVars() const
{
    int ndims, nvars, natts, unlim;
    nc_inq( _ncid, &ndims, &nvars, &natts, &unlim );
    return nvars;
}


// 
// numVars:
//  return the number of global attributes associated with the file
//
int NcFile::numAtts() const
{
    int ndims, nvars, natts, unlim;
    nc_inq( _ncid, &ndims, &nvars, &natts, &unlim );
    return natts;
}

//
// variable:
// retrieves the metadata associated with a variable specified by name.
//
NcVarBase
NcFile::variable( const string& varname ) throw ( NcErr )
{
    int id, status;
      // get the var id

    if (( status = nc_inq_varid( _ncid, varname.c_str(), &id )) != NC_NOERR ) 
        throw NcErr( _name.c_str(), varname.c_str(), nc_strerror( status ), 
                     _ncid, __FILE__, __LINE__ );
    
    return variable( id );
}  

//
// variable:
// retrieves the metadata associated with a variable specified by id.
//
NcVarBase
NcFile::variable( int id ) throw ( NcErr )
{
    nc_type type;
    int status;
    int natts, ndims;
    int dimids[NC_MAX_DIMS];
    char varname[NC_MAX_NAME];
    vector< NcDimension > dims;
    
    
      // get the variable metadata
    
    if (( status = nc_inq_var( _ncid, id, varname, &type, &ndims, dimids, &natts )) 
        != NC_NOERR ) 
        throw NcErr( _name.c_str(), varname, nc_strerror( status ), 
                     _ncid, __FILE__, __LINE__ );
    
      // get the dimensions and attributes
    char dimname[NC_MAX_NAME];
    for ( int i=0; i<ndims; i++ ) {
        nc_inq_dimname( _ncid, dimids[i], dimname );
        dims.push_back( dimension( dimname ));
    }

    NcVarBase v( varname, dims );
    
    
      // set the id's
    v._id = id;                

    v._parent = this;
    
      // add the attributes
    for ( int i=0; i<natts; i++ ) {
        char attname[NC_MAX_NAME];
        nc_inq_attname( _ncid, v._id, i, attname );
        v.addAttribute(  attribute( attname, v._id ));
    }

    return v;
}

//
// addDimension:
// add a dimension to the file
//
void
NcFile::addDimension( const NcDimension& dim )
{
    int status;
      // define the dimension
    if (( status = nc_def_dim( _ncid, dim._name.c_str(), dim._size, &dim._id )) 
        != NC_NOERR ) {
        cerr << "ERROR: "__FILE__":" << __LINE__
             << " NcFile::AddDimension(): "<<  dim._name << " " << nc_strerror( status );
        exit( -1 );
    }
}

//
// unlimDimension:
// returns the unlimited dimension of the file if one is defined, otherwise
//  throws an NcErr exception (a call to hasUnlimDimension() can determine
//     if one is defined)
//
NcDimension
NcFile::unlimDimension() const throw( NcErr )
{
    int status;
    int dimid;
    char dimname[NC_MAX_NAME];
    size_t dimsize;

    if (( status = nc_inq_unlimdim( _ncid, &dimid ))
        != NC_NOERR )  {
        throw ( NcErr( _name.c_str(), "no unlimited dimension", 
                       _ncid, __FILE__, __LINE__ ));  
    }
    
    if (( status = nc_inq_dim( _ncid, dimid, dimname, &dimsize ))
        != NC_NOERR ) {
        cerr << "ERROR: "__FILE__":" << __LINE__
             << ": NcFile::unlimDimension(): " << nc_strerror( status );
        exit( -1 );
    }

    NcDimension d(dimname, dimsize, true);
    d._id = dimid;
    return d;
}

//
//  hasDimension:
//       returns true if the named dimension is defined in the file
//
bool
NcFile::hasDimension( const string& dimName ) const
{
    int dimid;                  // placeholder
    return ( nc_inq_dimid( _ncid, dimName.c_str(), &dimid )  == NC_NOERR );
}

//
//  hasUnlimDimension:
//       returns true if the file has an unlimited (record) dimension
//
bool
NcFile::hasUnlimDimension() const
{
    int dimid;
    return ( nc_inq_unlimdim( _ncid, &dimid )  == NC_NOERR );
}


//
//  hasVariable:
//       returns true if the named variable is defined in the file
//
bool
NcFile::hasVariable( const string& varName ) const
{
    int varID;
    return ( nc_inq_varid( _ncid, varName.c_str(), &varID )  == NC_NOERR );
}

bool
NcFile::hasAttribute( const string& name ) const
{
    int attnum;                 // placeholder
    return ( nc_inq_attid( _ncid, NC_GLOBAL, name.c_str(), &attnum ) 
             == NC_NOERR );
}

void
NcFile::endDefineMode()
{
   nc_enddef( _ncid );
}

void
NcFile::enterDefineMode()
{
    nc_redef( _ncid );
}

//
// synchronize the file ( i.e. flush buffers to disk )
//
void
NcFile::flush() throw ( NcErr )
{
    int status;
    if (( status = nc_sync( _ncid )) != NC_NOERR ) 
        throw NcErr( _name.c_str(), nc_strerror( status ), 
                     _ncid, __FILE__, __LINE__ );
}

//
// attribute:
// retrieves the metadata associated with a attribute specified by name.
//  returns a const reference to 
// a static local attribute which should not be used directly.
//
const NcAttBase&
NcFile::attribute( const string& attname, int varId ) throw ( NcErr )
{
    nc_type type;
    size_t size;

    static NcCharAtt ca;
    static NcByteAtt ba;
    static NcShortAtt sa;
    static NcIntAtt ia;
    static NcFloatAtt fa;
    static NcDoubleAtt da;

    NcAttBase* att = 0;
    int status;

    if((status = nc_inq_att(_ncid, varId, attname.c_str(), &type, &size)) != NC_NOERR )
        throw NcErr( _name.c_str(), attname.c_str(), nc_strerror( status ), 
                     _ncid, __FILE__, __LINE__ );

    switch ( type ) {
    case NC_CHAR: 
    {
        char* val = new char[size];
        status = nc_get_att_text( _ncid, varId, attname.c_str(), val );
        ca = NcCharAtt( attname, size, val );
        att = &ca;
        delete[] val;
    }
    break;
    case NC_BYTE:
    {
        signed char* val = new signed char[size];
        status = nc_get_att_schar( _ncid, varId, attname.c_str(), val );
        ba = NcByteAtt( attname, size, val );
        att = &ba;
        delete[] val;
    }
    break;
    case NC_SHORT:
    {
        short* val = new short[size];
        status = nc_get_att_short( _ncid, varId, attname.c_str(), val );
        sa = NcShortAtt( attname, size, val );
        att = &sa;
        delete[] val;
    }
    break;
    case NC_INT:
    {
        int* val = new int[size];
        status = nc_get_att_int( _ncid, varId, attname.c_str(), val );
        ia = NcIntAtt( attname, size, val );
        att = &ia;
        delete[] val;
    }
    break;
    case NC_FLOAT:
    {
        float* val = new float[size];
        status = nc_get_att_float( _ncid, varId, attname.c_str(), val );
        fa = NcFloatAtt( attname, size, val );
        att = &fa;
        delete[] val;
    }
    break;
    case NC_DOUBLE:
    {
        double* val = new double[size];
        status = nc_get_att_double( _ncid, varId, attname.c_str(), val );
        da = NcDoubleAtt( attname, size, val );
        att = &da;
        delete[] val;
    }
    break;
    }

    if ( status != NC_NOERR ) 
        throw NcErr( _name.c_str(), attname.c_str(), nc_strerror( status ), 
                     _ncid, __FILE__, __LINE__ );
    
    return *att;
}

//
// dimension:
//
NcDimension
NcFile::dimension( const string& dimname ) throw ( NcErr )
{
    char name[NC_MAX_NAME];
    int id;
    size_t size;
    int unlimDimId;
    int status;

    if (( status = nc_inq_dimid( _ncid, dimname.c_str(), &id )) != NC_NOERR )
        throw NcErr( _name.c_str(), dimname.c_str(), nc_strerror( status ), 
                     _ncid, __FILE__, __LINE__ );
 
    if (( status = nc_inq_dim( _ncid, id, name, &size )) != NC_NOERR ) 
        throw NcErr( _name.c_str(), dimname.c_str(), nc_strerror( status ), 
                     _ncid, __FILE__, __LINE__ );

    if (( status = nc_inq_unlimdim( _ncid, &unlimDimId )) != NC_NOERR ) 
        throw NcErr( _name.c_str(), dimname.c_str(), nc_strerror( status ), 
                     _ncid, __FILE__, __LINE__ );
    

    NcDimension d( name, size, ( id == unlimDimId ));
    d._id = id;
    return d;
}




NcAttBase&  
NcAttBase::operator=( const NcAttBase& a )
{   
    _name = a._name; 
    return *this;
}


//======================================================================
//======================================================================
//======================================================================
 

NcVarBase::NcVarBase( const NcVarBase& v ) 
    : _name( v._name ), _id( v._id ), _parent( v._parent )
{

      // copy dimensions and attributes
    setDimensions( v._dims );
    setAttributes( v._atts );
}

NcVarBase::NcVarBase( const string& name, const NcDimension dims[], int ndims ) 
    : _name( name ), _id( -1 )
{
    setDimensions( dims, ndims );
}

NcVarBase::NcVarBase( const string& name, const vector< NcDimension >& dims ) 
    : _name( name ), _id( -1 )
{
    setDimensions( dims );
}

NcVarBase&
NcVarBase::operator=( const NcVarBase& v )
{
    _name = v._name;
    _id = v._id;
    _parent = v._parent;
    setDimensions( v._dims );
    setAttributes( v._atts );
    
    return *this;
}

 
NcVarBase::~NcVarBase()
{
    
    for ( att_itr it = _atts.begin(); it != _atts.end(); ++it )
            delete *it;
}

//
//  Get a variable attribute.
//   
const NcAttBase& NcVarBase::attribute( const string& attname ) const
    throw (NcErr)
{
    for ( const_att_itr it = _atts.begin(); it != _atts.end(); ++it )
        if ( attname == (*it)->name() )
            return **it;
      // didn't find it
    char errmsg[256];
    sprintf( errmsg, "attribute \"%s\" not found", attname.c_str() );
    throw NcErr( attname.c_str(), errmsg, 
                 _parent->ncid(), __FILE__, __LINE__ );
}  


bool NcVarBase::hasAttribute( const string& attname ) const
{
    for ( const_att_itr it = _atts.begin(); it != _atts.end(); ++it )
        if ( attname == (*it)->name() )
            return true;
    return false;
}

void NcVarBase::rename( const string& name ) throw (NcErr)
{
    int status = nc_rename_var( _parent->ncid(), _id, name.c_str() );
    if ( status != NC_NOERR ) 
        throw NcErr( _parent->name().c_str(),_name.c_str(), nc_strerror( status ),
                     _parent->ncid(), __FILE__, __LINE__ );
    _name = name;
}

void              
NcVarBase::setDimensions( const NcDimension dims[], int ndims )
{
      // clear the current dimensions
    _dims.clear();

      // copy the new dimensions
    for ( int i=0; i<ndims; i++ )
        _dims.push_back( dims[i] );
}

void              
NcVarBase::setDimensions( const vector< NcDimension >& dims )
{
      // clear the current dimensions
    _dims.clear();

      // copy the new dimensions
    for ( const_dim_itr it = dims.begin(); it != dims.end(); ++it ) 
        _dims.push_back( *it );
}
  
void              
NcVarBase::setAttributes( const vector< NcAttBase* >& atts )
{
      // clear the current dimensions
    for ( att_itr it = _atts.begin(); it != _atts.end(); ++it )
        delete *it;
    
    _atts.clear();

      // copy the new dimensions
    for ( const_att_itr it = atts.begin(); it != atts.end(); ++it )
        _atts.push_back( (*it)->clone() );
}

void              
NcVarBase::setAttributes( const NcAttBase* const atts[], int numatts )
{
      // clear the current dimensions
    for ( att_itr it = _atts.begin(); it != _atts.end(); ++it )
        delete *it;
    
    _atts.clear();

      // copy the new dimensions
    for ( int i=0; i<numatts; i++ )
        _atts.push_back( atts[i]->clone() );
}


void
NcVarBase::addAttribute( const NcAttBase& att )
{
    _atts.push_back( att.clone() );
}

const int* 
NcVarBase::dimids() const
{
    static int dimIds[NC_MAX_DIMS];
    int i = 0;
    for ( const_dim_itr it = _dims.begin(); it != _dims.end(); ++it )
        dimIds[i++] = (*it).id();
    return dimIds;
}


//======================================================================


NcDimension&
NcDimension::operator=( const NcDimension& d ) 
{ 
    _name = d._name; 
    _size = d._size; 
    _id = d._id;
    _isUnlimDim = d._isUnlimDim;
    return *this; 
}

} // end namespace ncfile
