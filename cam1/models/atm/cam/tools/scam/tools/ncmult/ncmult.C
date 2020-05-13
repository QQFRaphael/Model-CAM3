#include <stdlib.h>
#include <sys/param.h>
#include <stdio.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>
#include <netcdf.h>

void
HandleNCerr( int status )
{
    if ( status != NC_NOERR ) {
        cerr << "NETCDF ERROR: " << nc_strerror( status ) << endl;
        exit( 1 );
    }
}

void
PrintUsage()
{
            cerr << "\nusage: ncmult [-o outputfile] [-v varname1,mult1[,varmame2,mult2[...]]] inputfile\n" 
                << endl;
}

int main( int argc, char* argv[] )
{
    int ncid;
    int stat;
    int ndims;
    int nlev, ntime, nlat, nlon;
    char inputfile[MAXPATHLEN];
    char outputfile[MAXPATHLEN];
    char var_name[NC_MAX_VARS][NC_MAX_NAME];
    int  var_id;
    int  dim_id[NC_MAX_DIMS];
    size_t  dim_len[NC_MAX_DIMS];
    int  var_dimids[NC_MAX_DIMS];
    int  var_ndims;
    int  var_natts;
    size_t  var_size;
    nc_type var_type;
    bool  overwrite = true;

    int nvars = 0;
    float mult[NC_MAX_VARS];

    char    c;
    char*   ptr;
    char validOptions[]="v:o:";

    
    if  ( argc == 1 ) {
        PrintUsage();
        exit( 1 );
    }

    while ( (c = getopt( argc, argv, validOptions ) ) != -1 ) {
        switch (c) {
        case 'v':   
            ptr = strtok( optarg,"," );
            do {
                strcpy( var_name[nvars], ptr );
                if ( (ptr = strtok( NULL, "," ) ) == NULL ) {
                    PrintUsage();
                    exit( 1 );
                }
                mult[nvars++] = atof( ptr ); 
            }
            while ( ptr = strtok( NULL, "," ) );
            break;
        case 'o':
            strcpy( outputfile, optarg );
            overwrite = false;
            break;
        default:
            PrintUsage();
            exit( 1 );
        }
    }

    if ( optind != argc -1 ) {
        PrintUsage();
        exit( 1 );
    }
    
    strcpy( inputfile, argv[optind] );
    
    if ( ! overwrite ) {
          // copy the inputfile to the outputfile

        ifstream in( inputfile );
        if ( ! in ) {
            cerr << "ncmult: ERROR: Can't open input file \""
                 << inputfile << "\""<< endl;
            exit( 1 );
        }
        ofstream out( outputfile );
        if ( ! out ) {
            cerr << "ncmult: ERROR: Can't create output file \""
                 << outputfile << "\""<< endl;
            exit( 1 );
        }
        char b;
        while( in.get( b ) ) {
            out.put( b );
        }
    }
    else
        strcpy( outputfile, inputfile );

    if ( nc_open( outputfile, NC_WRITE, &ncid ) != NC_NOERR ) {
        cerr << "ncmult: ERROR: unable to open \""
             << outputfile << "\"" << endl;
        exit( 1 );
    }

    HandleNCerr( nc_inq_ndims( ncid, &ndims ) );
    
    for ( int d=0; d<ndims; d++ ) 
        HandleNCerr( nc_inq_dimlen( ncid, d, &dim_len[d] ) );
    
    for ( int n=0; n<nvars; n++ ) {
        var_size = 1;
        if ( nc_inq_varid( ncid, var_name[n], &var_id ) != NC_NOERR ) {
            cerr << "ncmult: WARNING: Variable \"" << var_name[n] <<
                "\" not found - skipping" << endl;
        }
        HandleNCerr( nc_inq_var( ncid, var_id, 0, &var_type,
                                 &var_ndims, var_dimids, &var_natts ) );
        
          // calculate the size of the variable
        for ( int i=0; i<var_ndims; i++ ) 
            var_size *= dim_len[var_dimids[i]];
        
        float* floatarray;
        double* doublearray;
        int* intarray;

        switch ( var_type ) {
        case NC_INT:
            intarray = new int[var_size];
            HandleNCerr( nc_get_var_int( ncid, var_id, intarray ) );
            for ( int s=0; s<var_size; s++ )
                intarray[s] *= mult[n];
            HandleNCerr( nc_put_var_int( ncid, var_id, intarray ) );
            delete[] intarray;
            break;            
        case NC_FLOAT:
            floatarray = new float[var_size];
            HandleNCerr( nc_get_var_float( ncid, var_id, floatarray ) );
            for ( int s=0; s<var_size; s++ )
                floatarray[s] *= mult[n];
            HandleNCerr( nc_put_var_float( ncid, var_id, floatarray ) );
            delete[] floatarray;
            break;            
        case NC_DOUBLE:
            doublearray = new double[var_size];
            HandleNCerr( nc_get_var_double( ncid, var_id, doublearray ) );
            for ( int s=0; s<var_size; s++ )
                doublearray[s] *= mult[n];
            HandleNCerr( nc_put_var_double( ncid, var_id, doublearray ) );
            delete[] doublearray;
            break;            
        default:
            cerr << "Illegal variable type (not one of NC_INT, NC_DOUBLE, NC_FLOAT)"
                 << endl;
            exit( 1 );
        }
    }

    HandleNCerr( nc_close( ncid ) );
    return 0;
}
