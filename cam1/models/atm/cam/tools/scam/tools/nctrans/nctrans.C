#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/param.h>          // MAXPATHLEN
#include <multiset.h>

struct dimension
{
    int    id;
    uint   size;
    char   name[NC_MAX_NAME];
};

typedef struct dimension dimension;

class  Variable
{
public:
    Variable() { id=-1; ndims=0; natts=0; mult=1; }
    int GetDimSize( const char* dim_name ) {
        for ( int i=0; i<ndims; i++ ) {
            if ( !strcmp( dims[dim_ids[i]].name, dim_name ) )
                return dims[i].size;
        }
        return 0;               // didn't find it
    }
    int        id;
    int        ndims;
    int        natts;
    nc_type    type;
    char       mult_units[NC_MAX_NAME];
    float      mult;
    int        dim_ids[NC_MAX_VAR_DIMS];
    dimension* dims;
    bool       is_modified;
    char       name[NC_MAX_NAME];
    static int nvars;
};

int Variable::nvars=0;

int 
GetVarId( Variable* vars, const char* name ) 
{
    for ( int n=0; n<vars[0].nvars; n++ ) 
        if ( !strcmp( vars[n].name, name ) )
            return vars[n].id;
    return -1;                  // didn't find it
}

void
ApplyMult( Variable& var, void* vardata )
{
    int i;
    int size = 1;
    float*  fp;
    double* dp;
    int*    ip;
    short*  sp;

    if ( ! var.is_modified )
        return;
    for ( i=0; i<var.ndims; i++ )
        size *= var.dims[i].size;

    switch ( var.type ) {
    case NC_FLOAT:
        fp=(float*)vardata;
        for( i=0; i<size; i++ )
           *fp++ *=  var.mult;
        break;
    case NC_INT:
        ip=(int*)vardata;
        for( i=0; i<size; i++ )
           *ip++ *=  var.mult;
        break;
    case NC_SHORT:
        sp=(short*)vardata;
        for( i=0; i<size; i++ )
           *sp++ *=  var.mult;
        break;
    case NC_DOUBLE:
        dp=(double*)vardata;
        for( i=0; i<size; i++ )
           *dp++ *=  var.mult;
        break;
    }
}


void
handle_err( int status )
{
    if ( status != NC_NOERR ) {
        cerr << "NETCDF ERROR: " << nc_strerror( status ) << endl;
        exit( 1 );
    }
}

void
PrintUsage();

void
ParseVaropts( char* varopts, Variable* vars )
{
    int mod_varid;
    char* ptr;
    int pos = 0;
    for ( ptr = strtok( varopts, "," ); ptr != NULL; pos++, ptr = strtok( NULL, "," ) ) {
        if ( pos%3 == 0 ) {
            if ( ( mod_varid = GetVarId( vars,ptr ) ) < 0 ) {
                cerr << "Error: Couldn't find variable " << ptr 
                     << " in input dataset." << endl;
                exit ( -1 );
            }
            vars[mod_varid].is_modified = true;
        }
        if ( pos%3 == 1 )
            strcpy( vars[mod_varid].mult_units, ptr );
        if ( pos%3 == 2 ) 
            vars[mod_varid].mult = atof( ptr );
    }
}

double
GetTimeMult( Variable& timevar, int ncid, char* new_units )
{
    double new_seconds, old_seconds;
    char old_units[NC_MAX_NAME];

    if ( !strcmp( new_units, "seconds"  ) )
        new_seconds=1;
    else if ( !strcmp( new_units, "hours"  ) )
        new_seconds=3600;
    else if ( !strcmp( new_units, "days"  ) )
        new_seconds=86400;
    else {
        cerr << "Error: Don't know how to convert time to " <<
            new_units << endl;
        exit ( -1 );
    }
    
    handle_err( nc_get_att_text( ncid, timevar.id,
                                 "units", old_units ) );
    if ( !strcmp( old_units, "seconds"  ) )
        old_seconds=1;
    else if ( !strcmp( old_units, "s"  ) )
        old_seconds=1;
    else if ( !strcmp( old_units, "hours"  ) )
        old_seconds=3600;
    else if ( !strcmp( old_units, "days"  ) )
        old_seconds=86400;
    else {
        cerr << "Error: Don't know how to convert time to " <<
            old_units << endl;
        exit ( -1 );
    }
    
    return old_seconds/new_seconds;
}
    


int 
main( int argc, char* argv[] )
{
    int      ncid, out_ncid;
    int      lev_dimid, time_dimid, lon_dimid, lat_dimid;
    char     infile[MAXPATHLEN], outfile[MAXPATHLEN];
    int      nvars, ndims;
    int      n, i;
    void*    vardata;
    char     c;
    char     varopts[1024];
    bool     mod_all_vars = false;
    char     outtime[32] = "";
    double   time_mult;
    int      time_varid;

    dimension *dims;
    Variable*  vars;

    
    char opts[] = "av:o:t:";

      // parse command line
    if ( argc < 2 ) 
        PrintUsage();
    while ( (c = getopt( argc, argv, opts ) ) != -1 ) {
        switch (c) {
        case 'a':  
            mod_all_vars = true;
            break;
        case 'v':  
            if ( ! optarg )
                PrintUsage();
            strcpy( varopts, optarg );
            break;
        case 't':  
            if ( ! optarg )
                PrintUsage();
            strcpy( outtime, optarg );
            break;
        case 'o':               // name for saving history output
            if ( ! optarg )
                PrintUsage();
            strcpy( outfile, optarg );
            break;
        case '?':
            PrintUsage();
            break;
        }
    }

    strcpy( infile, argv[optind] );
    handle_err( nc_open( infile, NC_NOWRITE, &ncid ) );
    
      /*
       *   get number of dimensions, number of variables
       *  (assume all datasets are same in this respect)
       */
    
    handle_err( nc_inq_nvars( ncid, &nvars ) );
    handle_err( nc_inq_ndims( ncid, &ndims ) );
    
    Variable::nvars =  nvars;

    dims = new dimension[ndims];
    vars = new Variable[nvars];


      /*
       *  get dimensions, variables info
       */
    for ( i=0; i< ndims; i++ ) {
        handle_err( nc_inq_dim( ncid, i, dims[i].name, &dims[i].size ) );
        dims[i].id = i;
    }
    
      /*
       * initialize variables
       */
    
    for ( n=0; n<nvars; n++ ) {
        handle_err( nc_inq_var( ncid, n, vars[n].name, &vars[n].type,
                                &vars[n].ndims, vars[n].dim_ids, &vars[n].natts ) );
        vars[n].id = n;
        vars[n].dims = dims;
    }
    
    handle_err( nc_inq_dimid( ncid, "lon",  &lon_dimid ) );
    handle_err( nc_inq_dimid( ncid, "lat",  &lat_dimid ) );
    handle_err( nc_inq_dimid( ncid, "lev",  &lev_dimid ) );
    handle_err( nc_inq_dimid( ncid, "time", &time_dimid ) );
    
    /* 
     *  output netcdf file
     */
    
    handle_err( nc_create( outfile, NC_WRITE, &out_ncid ) );
    int id;
    for ( i=0;i<ndims;i++ ) {
        handle_err( nc_def_dim( out_ncid, dims[i].name, dims[i].size, &id ) );
    }
      //
      // rename "tsec" variable to "time", and change its type 
      // to NC_DOUBLE
      // 
    if ( outtime[0] ){
        if ( ( time_varid = GetVarId( vars, "time" ) ) < 0 ) 
        {
            if ( ( time_varid = GetVarId( vars, "tsec" ) ) < 0 )
            {
                cerr << "Error: Couldn't find time variable in "
                     << infile;
                exit ( - 1 );
            }
            
            strcpy( vars[time_varid].name, "time" );
            vars[time_varid].type = NC_DOUBLE;
        }        
    }    

    for ( n=0; n<nvars; n++ ) {
        handle_err( nc_def_var( out_ncid, vars[n].name, vars[n].type,
                                vars[n].ndims, vars[n].dim_ids, &id ) );
    }
    
      // copy the attributes
    char attname[NC_MAX_NAME];
    int  attnum;
    for ( n=0; n<nvars; n++ ) {
        for ( attnum=0; attnum < vars[n].natts; attnum++ ) {
            handle_err( nc_inq_attname( ncid, vars[n].id, attnum, attname ) );
            nc_copy_att( ncid, vars[n].id, attname, out_ncid, vars[n].id );
        }
    }
//--------------------------------------------------------------------
    if ( outtime[0] ) {
        handle_err( nc_put_att_text( out_ncid, time_varid, "units", strlen( outtime ) + 1, outtime ) ); 
        char time_longname[256];
        sprintf(time_longname,"Time in %s", outtime );
        handle_err( nc_put_att_text( out_ncid, time_varid, "long_name", strlen( time_longname ) + 1, time_longname ) ); 
    }
//--------------------------------------------------------------------      //
      // parse the command line variables list
      //
    ParseVaropts( varopts, vars );

    for ( n=0; n<nvars; n++ ) {
        if ( vars[n].is_modified ) {
            handle_err( nc_put_att_text( out_ncid, n, "units", strlen( vars[n].mult_units ) + 1, vars[n].mult_units ) );
        }
    }
//--------------------------------------------------------------------
    if ( mod_all_vars == true ) {
        float one = 1;
        char units[NC_MAX_NAME];

        for ( n=0; n<nvars; n++ ) {
            if ( nc_get_att_text( ncid, n, "plot_units", units )
                 == NC_NOERR ) {
                vars[n].is_modified = true;
                handle_err( nc_put_att_text( out_ncid, n, "units", strlen( units ) + 1, units ) );
                handle_err( nc_get_att_float( ncid, n, "plot_multiplier", &vars[n].mult ) );
                handle_err( nc_put_att_float( out_ncid, n, "plot_multiplier", NC_FLOAT, 1, &one ) );
            }
        }
    }
//--------------------------------------------------------------------

    handle_err( nc_enddef( out_ncid ) );
    
      /*
       * create space for data
       */
    int max_data_size = dims[time_dimid].size * dims[lev_dimid].size * dims[lat_dimid].size * dims[lon_dimid].size;
    vardata = new double[ max_data_size ];
    
//--------------------------------------------------------------------
      /*
       * read in variables one at a time and convert units
       */

    for ( n=0; n<nvars; n++ ) {
        switch( vars[n].type ) {
        case NC_FLOAT :
            handle_err( nc_get_var_float( ncid, vars[n].id, (float*)vardata ) );
            ApplyMult( vars[n], vardata );
            handle_err( nc_put_var_float( out_ncid, vars[n].id, (float*)vardata ) );
            break;
        case NC_DOUBLE :
            handle_err( nc_get_var_double( ncid, vars[n].id, (double*)vardata ) );
            ApplyMult( vars[n], vardata );
            handle_err( nc_put_var_double( out_ncid, vars[n].id, (double*)vardata ) );
            break;
        case NC_INT :
            handle_err( nc_get_var_int( ncid, vars[n].id, (int*)vardata ) );
            ApplyMult( vars[n], vardata );
            handle_err( nc_put_var_int( out_ncid, vars[n].id, (int*)vardata ) );
            break;
        case NC_SHORT :
            handle_err( nc_get_var_short( ncid, vars[n].id, (short*)vardata ) );
            ApplyMult( vars[n], vardata );
            handle_err( nc_put_var_short( out_ncid, vars[n].id, (short*)vardata ) );
            break;
        case NC_CHAR :
            handle_err( nc_get_var_text( ncid, vars[n].id, (char*)vardata ) );
            handle_err( nc_put_var_text( out_ncid, vars[n].id, (char*)vardata ) );
            break;
        default:
            cerr << "Error: can't handle data type of variable " << vars[n].name << endl;
            exit ( -1 );
        }
    }    
    

      /*
       * convert time dimension
       */
    if ( outtime[0] ) {
        int ntime = dims[time_dimid].size;
        time_mult = GetTimeMult( vars[time_varid], ncid, outtime );

        handle_err( nc_get_var_double( ncid, time_varid, (double*)vardata ) );
        double* t = (double*)vardata;
        for ( int i=0; i<ntime; i++ ) {
            *t *= time_mult;
            ++t;
        }
        handle_err( nc_put_var_double( out_ncid, time_varid, (double*)vardata ) );
    }
    
    nc_close( ncid );
    nc_close( out_ncid );
    return 0;
}

void
PrintUsage()
{
    cerr <<  "usage: nctrans [-v varname,units,mult[,varname2,units,mult]] [-t timeunits] [-a] [-o <outputfile>] inputfile\n" << endl 
         <<  "   -v: variables to rescale - name, new units, multiplication factor "
         << endl
         << "   -t: convert time units; can be either seconds, hours, or days"
         << "   -a: rescale all variables according to attributes plot_units, plot_multiplier contained in input dataset"
         << endl
         << "  -o: name of output file"
         << endl
         << "Example:"
         << endl
         << "  nctrans -v lev,mb,.01,divq,gm/kg/day,8.64e7 -t days -o foo.nc arm0795.nc"
         << endl
         << "This will convert the \"lev\" variable from Pascals to Millibars,the \"divq\" variable from kg/kg/sec to gm/kg/day, and the time/tsec variable from seconds to days, putting the output in \"foo.nc\", using \"arm0795.nc\" as the input."
         << endl;

    
    
    exit( -1 );
}






