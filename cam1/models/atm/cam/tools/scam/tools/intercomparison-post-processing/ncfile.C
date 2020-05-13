#ifndef lint
static char rcsid[] = "$Id: ncfile.C,v 1.1.6.1 2004/05/18 20:51:36 jet Exp $";
#endif /* lint */

#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>

#include "ncfile.H"
#include "field.H"
 
//
//: Constructor: sets the filename and specifies error handling
//  Inputs: 
//     file - name of the netcdf file;
//     NCerrorHandling - allow NetCDF to handle errors (default false)
//
NCFile::NCFile( const char* file, bool NCerrorHandling )
{
    handleErrs = NCerrorHandling;
    filename = new char[ strlen( file ) + 1 ];
    strcpy( filename, file );
}


NCFile::~NCFile()
{
    delete[] filename;
    errStatus = nc_close( NCID );
}

size_t
NCFile::GetDimSize( const char* dimName )
{
    size_t    size;
    int     dimID;
    
    errStatus = nc_inq_dimid( NCID, dimName, &dimID );
    if ( errStatus != NC_NOERR ) {
        HandleError( errStatus );
        return 0;
    }
    nc_inq_dimlen( NCID, dimID, &size);

    return size;
}

bool
NCFile::GetGlobalAttribute( const char* attName, double* value )
{
    errStatus = nc_get_att_double( NCID, NC_GLOBAL, attName, value );
    if ( errStatus != NC_NOERR ) {
        HandleError( errStatus );
        return false;
    }
    return true;
}

bool
NCFile::GetAttribute(  const char* varname, const char* attName, char* attValue )
{
    int varid;

    errStatus = nc_inq_varid( NCID, varname, &varid );
    if ( errStatus != NC_NOERR ) {
        HandleError( errStatus );
        return false;
    }
    errStatus = nc_get_att_text( NCID, varid, attName, attValue );
    if ( errStatus != NC_NOERR ) {
        HandleError( errStatus );
        return false;
    }
    return true;
}

bool
NCFile::GetAttribute(  const char* varname, const char* attName, float* attValue )
{
    int varid;

    errStatus = nc_inq_varid( NCID, varname, &varid );
    if ( errStatus != NC_NOERR ) {
        HandleError( errStatus );
        return false;
    }
    errStatus = nc_get_att_float( NCID, varid, attName, attValue );
    if ( errStatus != NC_NOERR ) {
        HandleError( errStatus );
        return false;
    }
    return true;
}

void
NCFile::HandleError( int status )
{
    if ( handleErrs ) {
        if ( status != NC_NOERR ) 
            fprintf( stderr, "NetCDF Error:  %s\n", nc_strerror( status ) );
    }
}

bool
NCFile::Open( int mode, bool clobber ) {
//
//  READ mode -----------------------------------------------------
//    
    if ( mode == READ ) {
        errStatus =  nc_open( filename, NC_NOWRITE, &NCID );
        if ( errStatus != NC_NOERR) {
            HandleError( errStatus );
            return false;
        }
        else
            return true;
    }
//
//  WRITE mode -----------------------------------------------------
//
    else if ( mode == WRITE ) {
        errStatus =  nc_open( filename, NC_WRITE, &NCID );
        if ( errStatus != NC_NOERR ) {
            HandleError( errStatus );
            return false;
        }
        else
            return true;
    }
//
//  CREATE mode -----------------------------------------------------
//
    else if ( mode == CREATE ) {
        errStatus = nc_create( filename, ( clobber ? NC_CLOBBER : NC_NOCLOBBER ), &NCID );
        if ( errStatus != NC_NOERR ) {
            HandleError( errStatus );
            return false;
        }
        else 
            return true;
    }
// Unknown mode -
    else {
        fprintf( stderr, "Error: NCFile::Open() - Unknown mode %d\n", mode );
        return false;
    }
}

bool
NCFile::Read( const char* varName, int size, char* outText, int index )
{
    return ( Read( varName, size, ( void* ) outText, index, NC_CHAR ) );
}

bool
NCFile::Read( const char* varName, int size, double* outData, int index )
{
    return ( Read( varName, size, ( void* ) outData, index, NC_DOUBLE ) );
}

bool
NCFile::Read( const char* varName, int size, float* outData, int index )
{
    return ( Read( varName, size, ( void* ) outData, index, NC_FLOAT ) );
}

bool
NCFile::Read( const char* varName, int size, int* outData, int index )
{
    return ( Read( varName, size, ( void* ) outData, index, NC_INT ) );
}

    

bool 
NCFile::Read( const char* varName, float* outData )
{
    // get the variable's id
    int varID;

    errStatus = nc_inq_varid( NCID, varName, &varID );
    if ( errStatus != NC_NOERR ) {
        HandleError( errStatus );
   	fprintf( stderr, "ERROR - NCFile::Read(): Can't find var %s in %s\n",
                 varName, filename );
	return false;
    }
    errStatus =  nc_get_var_float( NCID, varID, outData );
    if ( errStatus != NC_NOERR ) {
        HandleError( errStatus );
   	fprintf( stderr, "ERROR - NCFile::Read(): Can't read var %s in %s\n",
                 varName, filename );
	return false;
    }
    return true;
}

bool 
NCFile::Read( const char* name, Field& f )
{
    int varID, nlev, ntime;

    nlev = GetDimSize( "lev" );
    ntime = GetDimSize( "time" );

    if ( f.GetNlev() > 1 && f.GetNlev() != nlev )
        cerr << "NCFile::Read() - WARNING: field nlev (" <<  f.GetNlev()
             << ") does not match nlev in dataset (" << nlev << "!" << endl;

    float* data = new float[nlev*ntime];

    errStatus = nc_inq_varid( NCID, name, &varID );
    if ( errStatus != NC_NOERR ) {
        HandleError( errStatus );
   	fprintf( stderr, "ERROR - NCFile::Read(): Can't find var %s in %s\n",
                 name, filename );
	return false;
    }
    errStatus =  nc_get_var_float( NCID, varID, data );
    if ( errStatus != NC_NOERR ) {
        HandleError( errStatus );
   	fprintf( stderr, "ERROR - NCFile::Read(): Can't read var %s in %s\n",
                 name, filename );
        delete[] data;
	return false;
    }

    for (int t=0; t<f.GetNtime(); t++ ) 
        for ( int l=0; l<f.GetNlev(); l++ ) 
            f[l][t] = data[t*f.GetNlev() + l];
    delete[] data;
    return true;
}

bool 
NCFile::Read( const char* varName, int* outData )
{
    int varID;
    // get the variable's id
    errStatus = nc_inq_varid( NCID, varName, &varID );
    if ( errStatus != NC_NOERR ) {
        HandleError( errStatus );
   	fprintf( stderr, "ERROR - NCFile::Read():Can't find var %s in %s\n",
                 varName, filename );
	return false;
    }
    errStatus =  nc_get_var_int( NCID, varID, outData );
    if ( errStatus != NC_NOERR ) {
        HandleError( errStatus );
   	fprintf( stderr, "ERROR - NCFile::Read():Can't read var %s in %s\n",
                 varName, filename );
	return false;
    }
    return true;
}

//
//  read: read surface or layer data from netcdf file.
//
bool
NCFile::Read( const char* varName, int size, void* outData, int index, nc_type dataType )
{
    size_t  start[2] = { 0, 0 };
    size_t  count[2] = { 0, 0 };
    int     ndims;
    int     varID;

       // get the variable's id
    errStatus = nc_inq_varid( NCID, varName, &varID );
    if ( errStatus != NC_NOERR ) {
        HandleError( errStatus );
   	fprintf( stderr, "ERROR - NCFile::Read(): Can't find var %s in %s\n",
                 varName, filename );
	return false;
    }
        // get the number of the variable dimensions
    errStatus = nc_inq_varndims( NCID, varID, &ndims );
    if ( errStatus != NC_NOERR ) {
        HandleError( errStatus );
	return false;
    }
    if ( ndims == 0 ) {          // single value in file
            // read single value into outData
        switch( dataType ) {
        case NC_CHAR:
            HandleError( nc_get_var_text( NCID, varID, (char*) outData ) );
            break;
        case NC_DOUBLE:
            HandleError( nc_get_var1_double( NCID, varID,(size_t*)&index , (double*) outData ) );
            break;
        case NC_FLOAT:
            HandleError( nc_get_var1_float( NCID, varID, (size_t*)&index,(float*) outData ) );
            break;
        case NC_INT:
            HandleError( nc_get_var1_int( NCID, varID,(size_t*)
                                          &index, (int*) outData ) );
            break;
        default:
            fprintf( stderr, "ERROR: NCFile::Read() - Unknown Type %d\n", dataType );
            exit( -1 );
        }
        return true;
    }

    start[0] = index;
    if ( ndims >  1 ){
        count[0] = 1;               // get data for single time point
        count[1] = size;
    }
    else
        count[0] = size;
    
        // read the data into outData
    switch( dataType ) {
    case NC_CHAR:
        HandleError( nc_get_vara_text( NCID, varID, start, count,
                                       (char*) outData ) );
        break;
    case NC_DOUBLE:
        HandleError( nc_get_vara_double( NCID, varID, start, count,
                                         (double*) outData ) );
        break;
    case NC_FLOAT:
        HandleError( nc_get_vara_float( NCID, varID, start, count,
                                        (float*) outData ) );
        break;
    case NC_INT:
        HandleError( nc_get_vara_int( NCID, varID, start, count,
                                      (int*) outData ) );
        break;
    default:
        fprintf( stderr, "ERROR: NCFile::Read() - Unknown Type %d\n", dataType );
        exit( -1 );
    }
    return true;
}

bool
NCFile::VarExists( const char* varName )
{
    int varID;
    errStatus = nc_inq_varid( NCID, varName, &varID );
    
    if ( errStatus != NC_NOERR ) 
	return false;
    else
        return true; 
}
