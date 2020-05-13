#include <netcdf.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>


#define DEFAULT_LPRINT 5

struct dimension
{
    int    id;
    int    size;
    char   name[NC_MAX_NAME];
};

typedef struct dimension dimension;

class  Variable
{
public:
    Variable() { }
    int GetDimSize( const char* dim_name ) {
        for ( int i=0; i<ndims; i++ ) {
            if ( !strcmp( dims[dim_ids[i]].name, dim_name ) )
                return dims[dim_ids[i]].size;
        }
        return 1;               // didn't find it, assume 1
    }
    int         id;
    int         ndims;
    int         natts;
    nc_type     type;
    int         dim_ids[NC_MAX_VAR_DIMS];
    dimension*  dims;
    char        name[NC_MAX_NAME];
    double*     data;
    double*     ave;
};


void
PrintUsage()
{

}

void
HandleNCerr( int status )
{
    if ( status != NC_NOERR ) {
        fprintf( stderr, "NETCDF ERROR: %s\n ", nc_strerror( status ));
        exit( 1 );
    }
}


//==================================================================
//==================================================================

int main( int argc, char* argv[] )
{
    int      ncid;
    int      ave_ncid;
    char*     input_file;
    int      nvars, ndims;
    int      n, i;

    dimension *dims;
    Variable*  vars;
    int    timeDimid;


    input_file = strdup(  argv[1] );

    HandleNCerr( nc_open( input_file, NC_NOWRITE, &ncid ));
    
      //  get number of dimensions, number of variables
    
    HandleNCerr( nc_inq_nvars( ncid, &nvars ) );
    HandleNCerr( nc_inq_ndims( ncid, &ndims ) );
    
    dims = new dimension[ndims];
    vars = new Variable[nvars];

      // get dimensions, variables info

    for ( i=0; i< ndims; i++ ) {
        HandleNCerr( nc_inq_dim( ncid, i, dims[i].name, (unsigned int*)&dims[i].size ) );
        dims[i].id = i;
    }
    
      // get the time dimension id
    int tsecVarid;

    HandleNCerr( nc_inq_dimid( ncid, "time", &timeDimid ) );
    int ntime = dims[timeDimid].size;
    int* time = new int[ntime];

    HandleNCerr( nc_inq_varid( ncid, "tsec", &tsecVarid ) );
    HandleNCerr( nc_get_var_int( ncid, tsecVarid, time ) );
    
    int timestep = time[1] - time[0];
    
    int steps_per_day = 86400/timestep;

      // initialize variables
    
    for ( n=0; n<nvars; n++ ) {
        HandleNCerr( nc_inq_var( ncid, n, vars[n].name, &vars[n].type,
                                &vars[n].ndims, vars[n].dim_ids, &vars[n].natts ) );
        vars[n].id = n;
        vars[n].dims = dims;
    } 

    
    /* 
     *  Create new netcdf files
     */
    
    HandleNCerr( nc_create( "diurnal_ave.nc", NC_WRITE, &ave_ncid ) );

    int timeDimSize = dims[timeDimid].size;

    dims[timeDimid].size = steps_per_day;
    for ( i=0;i<ndims;i++ ) 
        HandleNCerr( nc_def_dim( ave_ncid, dims[i].name, dims[i].size, NULL ) );

    dims[timeDimid].size = timeDimSize;
    for ( n=0; n<nvars; n++ ) {
        HandleNCerr( nc_def_var( ave_ncid, vars[n].name, vars[n].type,
                                 vars[n].ndims, vars[n].dim_ids, NULL ) );
          // copy the units attribute
        nc_copy_att( ncid, vars[n].id, "units", ave_ncid, vars[n].id );
    }
    HandleNCerr( nc_enddef( ave_ncid ) );


    
      /*
       * read one variable at a time
       */
    for ( n=0; n<nvars; n++ ) {        
          
          // allocate memory for the variable's data
        int ntime = vars[n].GetDimSize("time");
        int nlev = vars[n].GetDimSize("lev");

        vars[n].ave = new double[nlev*steps_per_day];
        memset( vars[n].ave, 0, sizeof(double)*nlev*steps_per_day );
        vars[n].data = new double[ntime*nlev];

          // read in the entire variable
        HandleNCerr( nc_get_var_double( ncid, vars[n].id, vars[n].data ) );

        int ntime_even = ntime - ntime%steps_per_day;
          // compute the average
        for ( int l=0; l<nlev; l++ ) {
            if ( ntime_even >= steps_per_day ) {
                for (int t=0; t<ntime_even; t++ ) 
                    vars[n].ave[(t%steps_per_day)*nlev+l] += vars[n].data[t*nlev+l];
                for( int t=0; t<steps_per_day; t++ )
                    vars[n].ave[t*nlev+l] /= ntime_even/steps_per_day;
            }
            else 
                vars[n].ave[l] = vars[n].data[l];
        }
            
        HandleNCerr( nc_put_var_double( ave_ncid, vars[n].id, 
                                       vars[n].ave ) );
        delete vars[n].ave;
        delete vars[n].data;
    }
    nc_close( ncid );
    nc_close( ave_ncid );
    return 0;
}




    










