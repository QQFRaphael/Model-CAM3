#include <netcdf.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include <iostream.h>

#define MAX_FILES 1000
#define NO_VAL -99999.0
#define DEFAULT_NBINS 90
#define DEFAULT_SS 1

void
PrintUsage()
{
    cerr << "\nusage: pdf -v <variable> [-l <level>] [-h] [-t <timeslice>] [-b <bin-size>] [-s <stencil diameter>] [-o <outputfile>] <input sccm history files>" << endl
         << "   -v : name of variable to extract from the input files" << endl
         << "   -l : index of level to process (if not specified, all levels will be processed)" << endl
         << "   -h : show histogram on stdout" << endl
         << "   -t : index of timeslice to process (if not specified, all timeslices will be processed)" << endl
         << "   -b : bin size (default divides range into" << DEFAULT_NBINS << " equal bins)" << endl
         << "   -s : stencil size in level and time dimension (default is " << DEFAULT_SS << ")" << endl
         << "   -o : name of output file (if not specified, output will be in pdf.<pid>.nc)" << endl
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

int main( int argc, char* argv[] )
{
    int      ncid[MAX_FILES];
    int      lev_dimID,lat_dimID,lon_dimID,time_dimID;
    int      nlev, ntime;
    int      levidx = -1, timeidx = -1;
    int      time_varID, lev_varID;
    int      lat_varID, lon_varID;
    int      varID;
    bool     showHistogram = false;
    int      var_ndims;
    int*     var_dimids;
    float*   vardata;
    char     varname[256] = "";
    int      nfiles=0;
    char*    input_files[1000];
    char     outputfile[256] = "";
    char     c;
    int      nbins = DEFAULT_NBINS;        // number of bins, default = 10
    float    binSize = 0;
    float*   binvals;
    int**    binCount;
    int      ss = DEFAULT_SS;                // stencil size, default = 1
    float    minval=NO_VAL, maxval = NO_VAL;
    char     validOptions[]="v:l:ht:b:s:o:?";

    if ( argc == 1 ) {
        PrintUsage();
        exit( 1 );
    }

    while ( (c = getopt( argc, argv, validOptions ) ) != -1 ) {
        switch (c) {
        case 'v':   
            strcpy( varname, optarg );
            break;
        case 'b':   
            binSize = atof( optarg );
            break;
        case 's':   
            ss = atoi( optarg );
            break;
        case 'h':   
            showHistogram = true;
            break;
        case 'l':   
            levidx = atoi( optarg ) - 1; // convert to C indexing
            break;
        case 't':   
            timeidx = atoi( optarg ) - 1; // convert to C indexing
            break;           
        case 'o':
            strncpy( outputfile, optarg, 256 );
            break;
        case '?':
        default:
            PrintUsage();
            exit( 1 );
        }
    }
      //
      // make sure that varname has been provided
      //

    if ( !varname[0] ) {
        PrintUsage();
        exit( 1 );
    }
    
    if ( optind == argc ) {
        PrintUsage();
        exit( 1 );
    }


    for ( ; optind < argc; optind++ ) {
        input_files[nfiles] = new char[NC_MAX_NAME];
        strcpy( input_files[nfiles], argv[optind] );
        nfiles++;
    }

    HandleNCerr( nc_open( input_files[0], NC_NOWRITE, &ncid[0] ) );
    HandleNCerr(nc_inq_dimid( ncid[0], "lev",  &lev_dimID ));
    HandleNCerr(nc_inq_dimid( ncid[0], "time", &time_dimID ));
    HandleNCerr(nc_inq_dimid( ncid[0], "lon", &lon_dimID ));
    HandleNCerr(nc_inq_dimid( ncid[0], "lat", &lat_dimID ));

    HandleNCerr(nc_inq_dimlen( ncid[0], lev_dimID, (size_t*)&nlev ));
    HandleNCerr(nc_inq_dimlen( ncid[0], time_dimID, (size_t*)&ntime ));

    HandleNCerr(nc_inq_varid( ncid[0], "lev", &lev_varID ));
    HandleNCerr(nc_inq_varid( ncid[0], "lat", &lat_varID ));
    HandleNCerr(nc_inq_varid( ncid[0], "lon", &lon_varID ));
    HandleNCerr(nc_inq_varid( ncid[0], "time", &time_varID ));
    HandleNCerr(nc_inq_varid( ncid[0], varname, &varID ));
    HandleNCerr(nc_inq_varndims( ncid[0], varID, &var_ndims ));
    var_dimids = new int[var_ndims];
    HandleNCerr(nc_inq_vardimid( ncid[0], varID, var_dimids ));

    float  *lev,*lon, *lat, *time;
    lev = new float[nlev];
    time = new float[ntime];
    lat = new float[1];
    lon = new float[1];

      // get the time, lev, lat, lon variables
    HandleNCerr(nc_get_var_float( ncid[0], time_varID, time ));
    HandleNCerr(nc_get_var_float( ncid[0], lev_varID, lev ));
    HandleNCerr(nc_get_var_float( ncid[0], lon_varID, lon ));
    HandleNCerr(nc_get_var_float( ncid[0], lat_varID, lat ));
    HandleNCerr(nc_close(ncid[0]));

    int stencilSize = 0;
    int ts = 0;                 // stencil size in time dimension
    int ls = 0;                 // stencil size in level dimension




    if ( timeidx == -1 ) { // no time index specified, do all times
        ts = ntime;
        timeidx = 0;
    }
    else {
        ts = ss;
        timeidx = timeidx -ts/2;
    }
        
    if ( nlev == 1 )  { // only one level
        levidx = 0;
        ls = 1;   
    }
    else {
        if ( levidx == -1 ) { // no level index specified, do all levels
            ls = nlev;
            levidx = 0;
        }
        else {
            ls = ss;
            levidx = levidx - ls/2;
        }
    }

    stencilSize = ts * ls;

    vardata = new float[nfiles*stencilSize];

    unsigned int start[4];
    unsigned int count[4];

    for ( int i=0; i< var_ndims; i++ ) {
        if( var_dimids[i] == lev_dimID ) {
            start[i] = levidx;
            count[i] = ls;
        }
        else if ( var_dimids[i] == time_dimID ) {
            start[i] = timeidx;
            count[i] = ts;
        }
        else if ( var_dimids[i] == lon_dimID ) {
            start[i] = 0;
            count[i] = 1;
        }
        else if ( var_dimids[i] == lat_dimID ) {
            start[i] = 0;
            count[i] = 1;
        }
        else {
            fprintf(stderr,"pdf: ERROR: Unknown dimension in dataset\n" );
            exit( -1 );
        }
    }
    for ( int n=0; n<nfiles; n++ ) { 
        HandleNCerr( nc_open( input_files[n], NC_NOWRITE, &ncid[n] ) );
        HandleNCerr( nc_get_vara_float( ncid[n], varID, start, count, &vardata[n*stencilSize] ));
        HandleNCerr( nc_close( ncid[n] ) );
    }
      // find the min and max to get the proper binning
    
    for ( int n=0; n<nfiles; n++ ) {
        for ( int s=0; s<stencilSize; s++ ) {
            if ( minval == NO_VAL )
                minval = vardata[n*stencilSize + s];
            if ( maxval == NO_VAL )
                maxval = vardata[n*stencilSize + s];
            
            if ( vardata[n*stencilSize + s] < minval )
                minval = vardata[n*stencilSize + s];
            if ( vardata[n*stencilSize + s] > maxval )
                maxval = vardata[n*stencilSize + s];
        }
    }
    
      //
      // set up the bins
      //

    if ( binSize != 0 )       // binSize provided
        nbins = (int) ceil( (maxval-minval)/binSize );
    else
        binSize = (maxval-minval)/nbins;

    cout << "\nsetting up bins; minval = " << minval << ", maxval =" << maxval 
         << ", nbins = " << nbins << endl;

    binCount = new int*[stencilSize];

    for ( int s=0; s<stencilSize; s++ ) {
        binCount[s] = new int[nbins+1];
    }
    binvals = new float[nbins+1];

    for ( int s=0; s<stencilSize; s++ ) {
        for ( int nb=0; nb <= nbins; nb++ ) {
            binvals[nb] = minval + nb * binSize;
            binCount[s][nb] = 0;
        }
    }

      // put the data into the bins

    for ( int nf=0;nf<nfiles; nf++ ) {
        for ( int nb=0; nb<=nbins; nb++ ){
            for ( int t=0; t<ts; t++ ) {
                for ( int l=0; l<ls; l++ ) {
                    if ( vardata[nf*stencilSize + t*ls +l] >= binvals[nb]  &&
                         vardata[nf*stencilSize + t*ls +l] <= binvals[nb+1] ) 
                        binCount[t*ls +l][nb]++;
                }
            }
        }
    }

    if ( showHistogram == true ) {
        for ( int l=0; l<ls; l++ ) {
            for ( int t=0; t<ts; t++ ) {
                cout << "level: " << lev[l+levidx] << endl;
                cout << "time:  " << time[t+timeidx] << endl;
                for ( int nb=0; nb<nbins; nb++ ) {
                    cout.width(10);
                    cout << binvals[nb] << ": " ;
                    for ( int i=0; i<binCount[t*ls+l][nb]; i++ )
                        cout << "*";
                    cout << endl;
                }
                cout << "\n=====================================================\n";
            }
        }
        cout << endl;
    }

    int outncid, timeDimID, levDimID, latDimID, lonDimID,  binDimID;
    int latVarID, lonVarID, levVarID, dataVarID, timeVarID, binVarID;

    if ( !outputfile[0] ) {
        sprintf( outputfile, "pdf.%d.nc", (int)getpid() );
        cout << "\n *** netcdf output is in " << outputfile << " ***" << endl;
    }

    HandleNCerr( nc_create( outputfile, NC_CLOBBER, &outncid ) );
    HandleNCerr( nc_def_dim( outncid, "time", (size_t) ts,  &timeDimID ) );
    HandleNCerr( nc_def_dim( outncid, "lev", (size_t) ls, &levDimID ) );
    HandleNCerr( nc_def_dim( outncid, "bin", (size_t)nbins, &binDimID ) );
    HandleNCerr( nc_def_dim( outncid, "lat", 1, &latDimID ) );
    HandleNCerr( nc_def_dim( outncid, "lon", 1, &lonDimID ) );
    
    int dimIDs[] = { timeDimID, levDimID, binDimID, latDimID, lonDimID };

    HandleNCerr( nc_def_var( outncid, "time", NC_FLOAT, 1, &dimIDs[0], &timeVarID ) );
    HandleNCerr( nc_def_var( outncid, "lev",  NC_FLOAT, 1, &dimIDs[1], &levVarID ) );
    HandleNCerr( nc_def_var( outncid, "bin",  NC_FLOAT, 1, &dimIDs[2], &binVarID ) );
    HandleNCerr( nc_def_var( outncid, "lat",  NC_FLOAT, 1, &dimIDs[3], &latVarID ) );
    HandleNCerr( nc_def_var( outncid, "lon",  NC_FLOAT, 1, &dimIDs[4], &lonVarID ) );
    HandleNCerr( nc_def_var( outncid, varname,  NC_INT, 3, dimIDs , &dataVarID ) );
    HandleNCerr( nc_enddef( outncid ) );
    
    HandleNCerr( nc_put_var_float ( outncid, latVarID, lat ) );
    HandleNCerr( nc_put_var_float ( outncid, binVarID, binvals ) );
    HandleNCerr( nc_put_var_float ( outncid, lonVarID, lon ) );
    HandleNCerr( nc_put_var_float ( outncid, levVarID, &lev[levidx] ) );
    HandleNCerr( nc_put_var_float ( outncid, timeVarID, &time[timeidx] ) );

      //
      // need to change the shape of the array to pass to netcdf
      //

    int* tmp = new int[nbins*stencilSize];
    for ( int t=0; t<ts; t++ ) 
        for ( int l=0; l<ls; l++ ) 
            for ( int nb=0; nb<nbins; nb++ ) 
                tmp[(t*ls+l)*nbins+nb] = binCount[t*ls+l][nb];

    HandleNCerr( nc_put_var_int( outncid, dataVarID, tmp ) );

    HandleNCerr( nc_close ( outncid ) );

    return 0;
}    


