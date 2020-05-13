/*
** Read in monthly 1-degree mixed layer depth data (Monterey and Levitus,
** circa 1997) and put it in netcdf format.  Include enough metadata so that
** a tool such as Ferret can produce reasonable plots.  Input files are in
** ASCII. Naming convention is mld.pd.xxx, where the xxx is the 3-digit month
** index (e.g. January is 001).  These represent the "fixed density criterion", 
** or method 2, from the Monterey and Levitus atlas). I got them from Steve Worley 
** in SCD from the web page http://dss.ucar.edu/datasets/ds285.0/data/woa94.
** Land values are coded as -99.9.  This value is used as _FillValue on output.
**
** Usage is "definemld1x1 -m output_1x1_mld_file [-v]".
**
** Jim Rosinski, Oct. 2001
*/

#include <unistd.h>      // getopt
#include <stdio.h>
#include <string.h>      // string functions
#include <sys/types.h>   // size_t
#include <netcdf.h>

#include "definemld1x1.h"

void usage_exit (char *);

int main (int argc, char *argv[])
{
  char *mldfile   = NULL;  // Mixed layer depths file: input
  char filename[11];       // Input filename (mld.pd.xxx). Files from Steve Worley (SCD)
  char string[80];         // attribute string

  int c;                   // character from arg list
  int i, j;                // longitude, latitude indices
  int n;                   // month index
  int ret;                 // return code fro netcdf function calls
  const int ndm[NMO] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};  // Days per month
  const char *timeatt = "days since 0000-01-01";   // attribute for time variable

  const double fillvalue = -99.9;                  // fill value (as on mld.pd.xxx datasets)

  // netcdf ids and lengths

  int dimids[3] = {-1, -1, -1};  // dimension ids for output netcdf variables
  int ncid = -1;                 // input/output initial conditions file
  int londimid = -1;             // longitude dimension id
  int latdimid = -1;             // latitude dimension id
  int timedimid = -1;            // time dimension id
  int lonid = -1;                // longitude variable id
  int latid = -1;                // latitude variable id
  int timeid = -1;               // time variable id
  int mldid = -1;                // mixed layer depth id
  int mldannid = -1;             // mixed layer depth annual average id
  int accum[NLAT][NLON];         // accumulation counter

  size_t start[3] = {0, 0, 0};        // only 1st value will change
  size_t count[3] = {1, NLAT, NLON};  // tells netcdf to grab a single time of lon x lat data

  double mld[NMO][NLAT][NLON];     // Mixed layer depths read from ASCII file
  double mldann[NLAT][NLON];       // Annual average of MLD
  double lon[NLON];                // longitude (degrees)
  double lat[NLAT];                // latitude (degrees)
  double time[NMO];                // time (days since year start)
  double sum;                      // day sum
  
  Boolean verbose = false;         // added printout

  FILE *fp = NULL;                 // input file pointer

  // Parse command line
 
  while ((c = getopt (argc, argv, "m:v")) != -1) {
 
    switch(c) {                                                                 
    case 'm':
      mldfile = optarg;
      break;
    case 'v':
      verbose = true;
      break;
    default:
      printf ("Unknown argument: %c\n", c);
      usage_exit (argv[0]);
    }
  }

  if ( ! mldfile) {
    printf ("-m mldfile not specified\n");
    usage_exit (argv[0]);
  }

  // Open output file

  if ((ret = nc_create (mldfile, NC_WRITE, &ncid)) != NC_NOERR)
    err_exit ("Error creating mld file %s\n%s\n", mldfile, nc_strerror (ret));

  // Define output dimensions and variables

  wrap_nc_def_dim (ncid, "lon", NLON, &londimid);
  wrap_nc_def_dim (ncid, "lat", NLAT, &latdimid);
  wrap_nc_def_dim (ncid, "time", NC_UNLIMITED, &timedimid);

  wrap_nc_def_var (ncid, "lon", NC_DOUBLE, 1, &londimid, &lonid);
  ret = nc_put_att_text (ncid, lonid, "units", strlen ("degrees_east"), "degrees_east");

  wrap_nc_def_var (ncid, "lat", NC_DOUBLE, 1, &latdimid, &latid);
  ret = nc_put_att_text (ncid, latid, "units", strlen ("degrees_north"), "degrees_north");

  wrap_nc_def_var (ncid, "time", NC_DOUBLE, 1, &timedimid, &timeid);
  ret = nc_put_att_text (ncid, timeid, "units", strlen (timeatt), timeatt);

  dimids[0] = timedimid;
  dimids[1] = latdimid;
  dimids[2] = londimid;

  wrap_nc_def_var (ncid, "MLDMONTHLY", NC_DOUBLE, 3, dimids, &mldid);
  (void) nc_put_att_double (ncid, mldid, "_FillValue", NC_DOUBLE, 1, &fillvalue);
  strcpy (string, "monthly mixed layer depths");
  (void) nc_put_att_text (ncid, mldid, "long_name", strlen (string), string);
  (void) nc_put_att_text (ncid, mldid, "units", strlen ("m"), "m");

  dimids[0] = latdimid;
  dimids[1] = londimid;

  wrap_nc_def_var (ncid, "MLDANN", NC_DOUBLE, 2, dimids, &mldannid);
  (void) nc_put_att_double (ncid, mldannid, "_FillValue", NC_DOUBLE, 1, &fillvalue);
  strcpy (string, "yearly average mixed layer depths");
  (void) nc_put_att_text (ncid, mldid, "long_name", strlen (string), string);
  (void) nc_put_att_text (ncid, mldannid, "units", strlen ("m"), "m");

  wrap_nc_enddef (ncid);

  // Write coordinate data

  for (i = 0; i < NLON; i++)
    lon[i] = 0.5 + i*(360./NLON);

  for (j = 0; j < NLAT; j++)
    lat[j] = -89.5 + j*(180./NLAT);

  // Define time for each month as midpoint of the month

  sum = 0.;
  for (n = 0; j < NMO; n++) {
    time[n] = sum + 0.5*ndm[n];
    sum += ndm[n];
  }

  wrap_nc_put_var_double (ncid, lonid, lon);
  wrap_nc_put_var_double (ncid, latid, lat);

  // Initialize variable for annual average

  for (j = 0; j < NLAT; j++)
    for (i = 0; i < NLON; i++) {
      mldann[j][i] = 0.;
      accum[j][i] = 0;
    }

  // loop over months, read and write mixed layer depths

  for (n = 0; n < NMO; n++) {
    sprintf (filename, "mld.pd.%3.3d", n+1);
    if ((fp = fopen (filename, "r")) == NULL)
      err_exit ("fopen failure on file %s:%s\n", filename, strerror (errno));

    if (verbose)
      printf ("Successfully opened file %s\n", filename);

    for (j = 0; j < NLAT; j++) {
      for (i = 0; i < NLON; i++) {
	if (fscanf (fp, "%lf", &mld[n][j][i]) != 1)
	  err_exit ("fscanf failure n=%d j=%d i=%d\n", n, j, i);
	if (mld[n][j][i] != fillvalue) {
	  ++accum[j][i];
	  mldann[j][i] += mld[n][j][i];
	}
      }
    }

    // Write the monthly data to the netcdf file

    start[0] = n;
    wrap_nc_put_vara_double (ncid, timeid, &start[0], &count[0], &time[n]);
    wrap_nc_put_vara_double (ncid, mldid, start, count, (double *) mld[n]);

    if (fclose (fp) != 0)
      err_exit ("fclose failure:%s\n", strerror (errno));
  }

  // Normalize the annual average data and write it to the netcdf file

  for (j = 0; j < NLAT; j++)
    for (i = 0; i < NLON; i++)
      if (accum[j][i] > 0) {
	mldann[j][i] /= accum[j][i];
	if (accum[j][i] < 12)
	  printf ("accum[%d][%d] = %d\n", j, i, accum[j][i]);
      } else {
	mldann[j][i] = fillvalue;
      }

  start[0] = 0;
  start[1] = 0;
  count[0] = NLAT;
  count[1] = NLON;
  wrap_nc_put_vara_double (ncid, mldannid, start, count, (double *) mldann);

  wrap_nc_close (ncid);

  return 0;
}

void usage_exit (char *argv0)
{
  err_exit ("Usage: %s -m MLDFILE [-v]\n", argv0);
}
