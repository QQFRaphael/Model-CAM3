#include <netcdf.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include <multimap.h>

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
    Variable() { ave_sum = 0; ave_squared_sum = 0; }
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
    double**    data;
    double*     ave;
    double*     sdev;
    double*     kurt;
    double*     skew;
    double      ave_sum;
    double      ave_squared_sum;
    double*     dev_sums;
    double*     dev_squared_sums;
};


void
PrintUsage()
{
            cerr << "\nusage: sdev [-l[n]] files\n" 
                << "-l[n]: print out highest (n) L1 and L2 normalizations" 
                << endl;
}

void
HandleNCerr( int status )
{
    if ( status != NC_NOERR ) {
        cerr << "NETCDF ERROR: " << nc_strerror( status ) << endl;
        exit( 1 );
    }
}


/*-----------------------------------------------------------------*/
void
CompStats( Variable* var, int num_inputfiles, bool calculate_L )
{
    int  n;
    double ep, sum, dev, variance, kurt, skew, p;
    
    int nlev;

    nlev = var->GetDimSize( "lev" );

    n = num_inputfiles;
    ep = 0;                     // for "corrected two-pass algorithm to minimize round-off
    dev = 0.0;

    if ( num_inputfiles < 2 ) {
        cerr << "num_inputfiles must be at least 2 in comp_stats" 
             << endl;
        exit( 1 );
    }
    
      // loop over time and levels, computing variance 
      // between input files 
      // (adapted from numerical recipes in C)
    int ntime = var->GetDimSize("time");
    for (int t=0; t<ntime; t++ ) {
        for ( int l=0; l<nlev; l++ ) {
            sum = 0.0;
            variance = 0.0;
            kurt = 0.0;
            skew = 0.0;
            p = 0;
            for ( int f=0; f<num_inputfiles; f++ ) 
                sum += var->data[f][t*nlev+l];
            
            var->ave[t*nlev+l] = sum / num_inputfiles;
            if ( calculate_L ) {
                var->ave_sum += fabs( var->ave[t*nlev+l] ); // for L1 calc
                var->ave_squared_sum += var->ave[t*nlev+l] * var->ave[t*nlev+l]; // for L2 calc
            }
            for ( int f=0; f<num_inputfiles; f++ ) {
                dev = var->data[f][t*nlev+l] - var->ave[t*nlev+l];
                var->dev_sums[f] += fabs( dev ); // for L1 calc
                variance += ( p = dev * dev );
                skew += ( p *= sum );
                kurt += ( p *= sum );
                if ( calculate_L ) 
                    var->dev_squared_sums[f] += ( dev * dev ); // for L2 calc
            }
            
            variance = ( variance - ep * ep / num_inputfiles ) / ( num_inputfiles - 1 );
            var->sdev[t*nlev+l] = sqrt( variance );
            if ( variance ) {   // skew and kurtosis only exist if variance != 0
                var->skew[t*nlev+l] = skew / (num_inputfiles*(variance*var->sdev[t*nlev+l]));
                var->kurt[t*nlev+l] = kurt / (num_inputfiles*variance*variance) - 3.0;
            }             
        }
    }
}

//==================================================================
//==================================================================

int main( int argc, char* argv[] )
{
    int*     ncid;
    int      ave_ncid, sdev_ncid,kurt_ncid,skew_ncid;
    int      num_inputfiles=0;
    char*    input_files[1000];
    int      nvars, ndims;
    int      n, i;
    int      tmp_int[100000];
    double   tmp_double[100000];
    
    bool     calculate_L = false;

    dimension *dims;
    Variable*  vars;
    char   tmp[MAX_NC_NAME];
    char    c;
    char validOptions[]="l:?";
    int num_L_to_print = DEFAULT_LPRINT;
    
    if  ( argc == 1 ) {
        PrintUsage();
        exit( 1 );
    }

    while ( (c = getopt( argc, argv, validOptions ) ) != -1 ) {
        switch (c) {
        case 'l':   
            calculate_L = true;
            if ( optarg != NULL )
                num_L_to_print = atoi( optarg );
            if ( num_L_to_print == 0 )
                num_L_to_print = DEFAULT_LPRINT;
            break;
        case '?':
        default:
            PrintUsage();
            exit( 1 );
        }
    }
    
            
    for ( ; optind < argc; optind++ ) {
        input_files[num_inputfiles] = new char[NC_MAX_NAME];
        strcpy( input_files[num_inputfiles], argv[optind] );
        num_inputfiles++;
    }

    ncid = new int[num_inputfiles];
    //    for ( i=0; i<num_inputfiles; i++ ) {
    //        HandleNCerr( nc_open( input_files[i], NC_NOWRITE, &ncid[i] ) );
    //    }
    HandleNCerr( nc_open( input_files[0], NC_NOWRITE, &ncid[0] ) );
    
      //  get number of dimensions, number of variables
      //  (assume all datasets are same in this respect)
    
    HandleNCerr( nc_inq_nvars( ncid[0], &nvars ) );
    HandleNCerr( nc_inq_ndims( ncid[0], &ndims ) );
    
    dims = new dimension[ndims];
    vars = new Variable[nvars];

      // get dimensions, variables info

    for ( i=0; i< ndims; i++ ) {
        HandleNCerr( nc_inq_dim( ncid[0], i, dims[i].name, (unsigned int*)&dims[i].size ) );
        dims[i].id = i;
    }


      // initialize variables
    
    for ( n=0; n<nvars; n++ ) {
        HandleNCerr( nc_inq_var( ncid[0], n, vars[n].name, &vars[n].type,
                                &vars[n].ndims, vars[n].dim_ids, &vars[n].natts ) );
        vars[n].id = n;
        vars[n].dims = dims;
    } 

    
    /* 
     *  Create new netcdf files
     */
    
    HandleNCerr( nc_create( "ave.nc", NC_WRITE, &ave_ncid ) );
    HandleNCerr( nc_create( "sdev.nc", NC_WRITE, &sdev_ncid ) );
    HandleNCerr( nc_create( "kurt.nc", NC_WRITE, &kurt_ncid ) );
    HandleNCerr( nc_create( "skew.nc", NC_WRITE, &skew_ncid ) );

    for ( i=0;i<ndims;i++ ) {
        HandleNCerr( nc_def_dim( ave_ncid, dims[i].name, dims[i].size, NULL ) );
        HandleNCerr( nc_def_dim( sdev_ncid, dims[i].name, dims[i].size, NULL ) );
        HandleNCerr( nc_def_dim( kurt_ncid, dims[i].name, dims[i].size, NULL ) );
        HandleNCerr( nc_def_dim( skew_ncid, dims[i].name, dims[i].size, NULL ) );
    }
    for ( n=0; n<nvars; n++ ) {
            HandleNCerr( nc_def_var( ave_ncid, vars[n].name, vars[n].type,
                                    vars[n].ndims, vars[n].dim_ids, NULL ) );
            HandleNCerr( nc_def_var( sdev_ncid, vars[n].name, vars[n].type,
                                    vars[n].ndims, vars[n].dim_ids, NULL ) );
            HandleNCerr( nc_def_var( kurt_ncid, vars[n].name, vars[n].type,
                                    vars[n].ndims, vars[n].dim_ids, NULL ) );
            HandleNCerr( nc_def_var( skew_ncid, vars[n].name, vars[n].type,
                                    vars[n].ndims, vars[n].dim_ids, NULL ) );
    }
      // copy the units attribute
     for ( n=0; n<nvars; n++ ) {
         nc_copy_att( ncid[0], vars[n].id, "units", ave_ncid, vars[n].id );
         nc_copy_att( ncid[0], vars[n].id, "units", sdev_ncid, vars[n].id );
         nc_copy_att( ncid[0], vars[n].id, "units", kurt_ncid, vars[n].id );
         nc_copy_att( ncid[0], vars[n].id, "units", skew_ncid, vars[n].id );
    }
   

    HandleNCerr( nc_enddef( ave_ncid ) );
    HandleNCerr( nc_enddef( sdev_ncid ) );
    HandleNCerr( nc_enddef( kurt_ncid ) );
    HandleNCerr( nc_enddef( skew_ncid ) );

    
      /*
       * read one variable at a time
       */
    for ( n=0; n<nvars; n++ ) {
      //  lat,lon,lev,time,basedate variables
      //   just copy them to output files unchanged
        if ( nc_get_att_text( ncid[0], vars[n].id, "num_levels", tmp ) != NC_NOERR ) {
//  	  if ( strcmp( vars[n].name, "pres" ) ) // pres should be treated like other fields
//  	    {
//                  if ( vars[n].type == NC_FLOAT ||  vars[n].type == NC_DOUBLE ) {
//                      HandleNCerr( nc_get_var_double( ncid[0], vars[n].id, tmp_double ) );
//                      HandleNCerr( nc_put_var_double( ave_ncid, vars[n].id, tmp_double ) );
//                      HandleNCerr( nc_put_var_double( sdev_ncid, vars[n].id, tmp_double ) );
//                      HandleNCerr( nc_put_var_double( kurt_ncid, vars[n].id, tmp_double ) );
//                      HandleNCerr( nc_put_var_double( skew_ncid, vars[n].id, tmp_double ) );
//                  }
//                  else if ( vars[n].type == NC_INT ) {
//                      HandleNCerr( nc_get_var_int( ncid[0], vars[n].id, tmp_int ) );
//                      HandleNCerr( nc_put_var_int( ave_ncid, vars[n].id, tmp_int ) );
//                      HandleNCerr( nc_put_var_int( sdev_ncid, vars[n].id, tmp_int ) );
//                      HandleNCerr( nc_put_var_int( kurt_ncid, vars[n].id, tmp_int ) );
//                      HandleNCerr( nc_put_var_int( skew_ncid, vars[n].id, tmp_int ) );
//                  }
//  		continue;           /* ignore text variables */
//  	    }
        }
        
          // 
          // field variables
          //
        
	if ( !(vars[n].type == NC_FLOAT ||  vars[n].type == NC_DOUBLE) )
	  continue;           /* ignore text variables */
        cout << "working on \"" <<  vars[n].name << "\"" << endl;
        
          // allocate memory for the variable's data
        int ntime = vars[n].GetDimSize("time");
        int nlev = vars[n].GetDimSize("lev");

        vars[n].dev_sums = new double[num_inputfiles];
        vars[n].dev_squared_sums = new double[num_inputfiles];
        vars[n].ave = new double[ntime*nlev];
        vars[n].sdev = new double[ntime*nlev];
        vars[n].kurt = new double[ntime*nlev];
        vars[n].skew = new double[ntime*nlev];
        memset( vars[n].dev_sums, 0, sizeof(double) * num_inputfiles );
        memset( vars[n].dev_squared_sums, 0, sizeof(double) * num_inputfiles );
        memset( vars[n].ave, 0, sizeof(double) * ntime*nlev );
        memset( vars[n].sdev, 0, sizeof(double) * ntime*nlev );
        memset( vars[n].kurt, 0, sizeof(double) * ntime*nlev );
        memset( vars[n].skew, 0, sizeof(double) * ntime*nlev );

        vars[n].data = new double*[num_inputfiles];
        for ( int f=0; f<num_inputfiles; f++ ) {
            vars[n].data[f] = new double[ntime*nlev];
            if ( vars[n].data[f] == NULL ) {
                cerr << " memory allocation failed\n" ;
                exit (1);
            }
        }
          // read in the entire variable
        for ( int f=0; f<num_inputfiles; f++ ) {
	  if (f>0)
	  HandleNCerr( nc_open( input_files[f], NC_NOWRITE, &ncid[f] ) );
	  HandleNCerr( nc_get_var_double( ncid[f], vars[n].id,
					  vars[n].data[f]) );
	  if (f>0)
	  HandleNCerr( nc_close(ncid[f] ) );
        }
          /*
           * compute statistics
           */
        CompStats( &vars[n], num_inputfiles, calculate_L );
        
          /*
           * write statistical data back to the new files
           */
        HandleNCerr( nc_put_var_double( ave_ncid, vars[n].id, 
                                       vars[n].ave ) );
        HandleNCerr( nc_put_var_double( sdev_ncid, vars[n].id, 
                                       vars[n].sdev ) );
        HandleNCerr( nc_put_var_double( kurt_ncid, vars[n].id, 
                                       vars[n].kurt ) );
        HandleNCerr( nc_put_var_double( skew_ncid, vars[n].id, 
                                       vars[n].skew ) );
    
        if ( calculate_L ) {
              /*
               * compute L1 and L2 values for the ensemble members
               */
            multimap<double, const char*> L1_list;
            multimap<double, const char*> L2_list;
            double L1, L2;
            for ( i=0; i<num_inputfiles; i++ ) {
                  // compute L1
                L1 = vars[n].dev_sums[i]/vars[n].ave_sum;
                L1_list.insert( pair<double, const char*>( L1, input_files[i]) );
                
                  // compute L2
                L2 = sqrt( vars[n].dev_squared_sums[i] ) / 
                    sqrt( vars[n].ave_squared_sum );
                L2_list.insert( pair<double, const char*>( L2, input_files[i] ) );
            }
            
            cout << "\n" << vars[n].name << " - Lowest L1 values:  \n" << endl;
            int x=0;
            for ( multimap<double, const char*>::iterator it = L1_list.begin();
                  it != L1_list.end(); ++it ) {
                if ( x < num_L_to_print )
                    cout << (*it).second << ", " <<  (*it).first << endl;
                x++;
            }
            cout << "\n" << vars[n].name << " - Highest L1 values:  \n" << endl;
            x=0;
            for ( multimap<double, const char*>::iterator it = L1_list.begin();
                  it != L1_list.end(); ++it ) {
                if ( L1_list.size() - x <= num_L_to_print )
                    cout << (*it).second << ", " <<  (*it).first << endl;
                x++;
            }
            cout << "\n" << vars[n].name << " - Lowest L2 values:  \n" << endl;
            x=0;
            for ( multimap<double, const char*>::iterator it = L2_list.begin();
                  it != L2_list.end(); ++it ) {
                if ( x < num_L_to_print ) 
                    cout << (*it).second << ", " << (*it).first << endl;
                x++;
            }
            cout << "\n" << vars[n].name << " - Highest L2 values:  \n" << endl;
            x=0;
            for ( multimap<double, const char*>::iterator it = L2_list.begin();
                  it != L2_list.end(); ++it ) {
                if ( L2_list.size() - x <= num_L_to_print ) 
                    cout << (*it).second << ", " << (*it).first << endl;
                x++;
            }
            cout << "\n --------------------------\n " << endl;
        }

          // free the memory used for this variable
        delete[] vars[n].dev_sums;
        delete[] vars[n].dev_squared_sums;
        delete[] vars[n].ave;
        delete[] vars[n].sdev;
        delete[] vars[n].kurt;
        delete[] vars[n].skew;
        for ( int f=0; f<num_inputfiles; f++ ) 
            delete[] vars[n].data[f];
        delete[] vars[n].data;

    }
    
    //    for ( i=0; i<num_inputfiles; i++ ) 
    //	  nc_close( ncid[i] );
    HandleNCerr( nc_close(ncid[0] ) );
    nc_close( ave_ncid );
    nc_close( sdev_ncid );
    nc_close( kurt_ncid );
    nc_close( skew_ncid );

    return 0;
}




    










