/*
** This code adds Q flux information computed from a control run to a Slab Ocean
** Model boundary dataset.  That boundary dataset must already contain mixed layer
** depths (MLDANN) and SST fields.  The specifics are as follows:
**
** sstfile: Boundary monthly SST information
**          input: SST_cpl field used for open water calculation of QFLUX field.
**                 ice_cov field defines where sea ice is for QFLUX calculation.
**                 MLDANN is annually averaged mixed layer depths, also required 
**                 for Q flux calculations.
**          output: QFLUX is net heat flux computed from a control run:
**                  QFLUX = FNET - (H2OTERM - ICETERM)
**                 
** firstfile through lastfile: input: FNET = FSNS - FLNS - LHFLX - SHFLX.
*/

#include <math.h>        /* M_PI */
#include <unistd.h>      // getopt
#include <stdlib.h>      // atoi
#include <ctype.h>       // isdigit
#include <stdio.h>
#include <string.h>      // string functions
#include <sys/types.h>   // size_t
#include <math.h>
#include <netcdf.h>

#include "defineqflux.h"

void usage_exit (char *);
int date (int, int);
char *bldfn (char *, int, int);
void incdate (int *, int *);
void AddToString (char *, char *, char *);
int check_hname (char *, int *);
int getprvnxtindx (int, int, int *, int *);
double daysof (int, int);
void timefactors (int, int, int, int, double *, double *);

const int daysmo[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
int date_sst[12];        // dates from SST file
int datesec_sst[12];     // seconds of date from SST file

int main (int argc, char *argv[])
{
  char *sstfile   = NULL;  // SST boundary file: QFLUX added on output
  char *firstfile = NULL;  // 1st file in sequence from control run
  char *lastfile  = NULL;  // last file in sequence from control run
  char *fname     = NULL;  // next file in monthly sequence
  char tmp[5];             // temporary for copying
  char *cmdline   = NULL;  // cmd line added to history attribute of output file(s)
  char *history   = NULL;  // history attribute
  char caseid[33];         // control case name
  char oldcaseid[33];      // control case name (prv months)
  char *prefix    = "";    // prefix to filenames (e.g. <prefix>xxxx-xx.nc)

  int c;                   // character from arg list
  int i, j;                // longitude, latitude indices
  int n;                   // month index
  int yr;                  // year value
  int mo;                  // month value (1-12)
  int first_yr, last_yr;   // 1st and last years from control run
  int first_mo, last_mo;   // 1st and last months from control run
  int ret;                 // return code fro netcdf function calls
  int npls, nmin;          // next and prv month index
  int moaccum[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};  // accumulation counter for each month
  int nacs[12];            // accumulation counter (a la CAM)
  int prefix_len1;         // prefix length
  int prefix_len2;         // prefix length

  size_t hislen;           // length of history attribute
  size_t caselen;          // length of case attribute on on control run

  // netcdf ids and lengths

  int dimids[3] = {-1, -1, -1};  // lon x lat x time dimension ids
  int xydimids[2] = {-1, -1};    // lon x lat dimension ids
  int ncid_sst = -1;             // input/output SST file
  int ncid_cntl = -1;            // input control run file
  int lonid = -1;                // longitude 
  int latid = -1;                // latitude
  int nlonid = -1;               // number of lons per lat
  int sstid = -1;                // sea surface temperature
  int sicthkid = -1;             // sea ice thickness
  int siccntid = -1;             // sea ice thickness
  int icefracid = -1;            // ice fraction
  int fsnsid = -1;               // net solar flux at surface
  int flnsid = -1;               // net longwave flux at surface
  int lhflxid = -1;              // surface latent heat flux
  int shflxid = -1;              // surface sensible heat flux
  int mldid = -1;                // mixed layer depth
  int qfluxid = -1;              // net atm-ocean heat flux (CSIM sign convention)
  int qfluxyavgid = -1;          // yearly average QFLUX
  int qfluxprediddleid = -1;     // QFLUX before time diddling
  int diag1id = -1;              // diagnostic should be roundoff to qfluxyavg
  int diag2id = -1;              // diagnostic should be roundoff to monthly qflux
  int h2otermid = -1;            // delta-sst term
  int icetermid = -1;            // ice growth/decay term
  int fnetid = -1;               // net atm flux term
  int fnetyavgid = -1;           // yearly average FNET
  int timedimid_sst = -1;        // time dimension id (sst)
  int dateid_sst = -1;           // date id (YYYYMMDD) (sst)
  int datesecid_sst = -1;        // date id (YYYYMMDD) (sst)
  int nstep;                     // current model timestep

  size_t timelen = -1;           // length of time dimension variable
  size_t start[3] = {0, 0, 0};   // only 1st value will change
  size_t startxy[2] = {0, 0};    // always start at x=0, y=0
  size_t count[3] = {1, PLAT, PLON};  // tells netcdf to grab a single time of lon x lat data
  size_t countxy[2] = {PLAT, PLON};   // tells netcdf to grab lon x lat data

  double sst[12][PLAT][PLON];    // SST read from dataset
  double fsns[PLAT][PLON];       // FSNS from control run
  double flns[PLAT][PLON];       // FLNS from control run
  double lhflx[PLAT][PLON];      // LHFLX from control run
  double shflx[PLAT][PLON];      // SHFLX from control run
  double fnet[12][PLAT][PLON];   // fsns - flns - lhflx - shflx
  double fnetyavg[PLAT][PLON];   // yearly average FNET
  double qflux[12][PLAT][PLON];  // net flux
  double qfluxyavg[PLAT][PLON];  // yearly average QFLUX (dataset)
  double diag1[PLAT][PLON];      // diagnostic should be roundoff to qfluxyavg
  double diag2[12][PLAT][PLON];  // diagnostic should be roundoff to qfluxyavg
  double modelqflux;             // QFLUX as time interpolated in model
  double modelqfluxyavg[PLAT][PLON]; // yearly average QFLUX (model)
  double h2oterm[12][PLAT][PLON]; // intermediate in QFLUX computation
  double iceterm[12][PLAT][PLON]; // intermediate in QFLUX computation
  double sicvol[12][PLAT][PLON];  // sea ice volume
  double sicthk[12][PLAT][PLON];  // sea ice thickness (from SST file or cntl run)
  double siccnt[12][PLAT][PLON];  // sea ice concentration (fraction)
  double mldann[PLAT][PLON];      // annual mean mixed layer depths (input)
  double flat[PLAT];              // model latitude (lat coord variable from IC file
  double dayfact;                 // time interpolation factor
  double fact1, fact2;            // time interpolation factors (model)
  double delsst;                  // time change in SST
  double delsicvol;               // time change in sea ice volume
  double siccntnpls, siccntnmin;  // bounded ice concentrations at + and - time levels
  double sstnpls, sstnmin;        // bounded SSTs at + and - time levels
  double timesum;                 // sum of time weights
  double diddle[PLAT][PLON];      // correction to QFLUX
  
  const int dtime = 1200;         // length of model timestep (for diagnostic calc)
  const int ntspy = (86400*365)/dtime;  // number of timesteps per year
  const double volmax = 5.0;      // max ice volume
  const double tsice = -1.7999;   // min temperature of sea water (C)

  //JR 10/30/02: changed constants to match shared constants

  const double rhoocn = 1.026e3;  // mass density of water (kg/m3) (approx mixed layer)
  const double cpocn = 3.996e3;   // specific heat of sea water (joules/kg/k)
  const double latice = 3.337e5;  // latent heat of fusion ice (J/kg)
  const double rhoice = 0.917e3;  // density of ice (kg/m^3)
  const double li = latice*rhoice; // heat of fusion of sea-ice (J/m3)

  int nlon[PLAT];                 // number of longitudes per latitude

  Boolean diffaccum;              // check if number of months matches
  Boolean verbose = false;        // added printout
  Boolean preserve_mavg = true;   // Default: Do apply C. Bitz time diddling algorithm to
                                  // preserve monthly means in time interpolation
  Boolean preserve_yavg = false;  // Default: Do not apply Briegleb time diddling algorithm to
                                  // preserve yearly mean in time interpolation

  // Parse command line, saving contents for adding metadata to output file(s).

  cmdline = (char *) malloc (strlen (argv[0]) + 1);
  strcpy (cmdline, argv[0]);

  while ((c = getopt (argc, argv, "bcf:l:s:v")) != -1) {
 
    switch(c) {                                                                 
    case 'b':
      preserve_yavg = true;
      preserve_mavg = false;
      realloc (cmdline, strlen (cmdline) + 4);
      strcat (cmdline, " -b");
      break;
    case 'c':
      preserve_mavg = true;
      preserve_yavg = false;
      realloc (cmdline, strlen (cmdline) + 4);
      strcat (cmdline, " -b");
      break;
    case 'f':
      firstfile = optarg;
      AddToString (cmdline, " -f ", firstfile);
      break;
    case 'l':
      lastfile = optarg;
      AddToString (cmdline, " -l ", lastfile);
      break;
    case 's':
      sstfile = optarg;
      AddToString (cmdline, " -s ", sstfile);
      break;
    case 'v':
      verbose = true;
      break;
    default:
      printf ("Unknown argument: %c\n", c);
      usage_exit (argv[0]);
    }
  }

  if ( ! firstfile)
    printf ("-f firstfile not specified\n");

  if ( ! lastfile)
    printf ("-l lastfile not specified\n");

  if ( ! sstfile)
    printf ("-s sstfile not specified\n");

  if ( ! firstfile || ! lastfile || ! sstfile)
    usage_exit (argv[0]);

  printf ("%s\n", cmdline);

  // Open input files

  if ((ret = nc_open (sstfile, NC_WRITE, &ncid_sst)) != NC_NOERR)
    err_exit ("Error opening SST file %s\n%s\n", sstfile, nc_strerror (ret));
  
  // Ensure consistent latitude and longitude dimension sizes

  check_dims_consistent (ncid_sst, sstfile);

  wrap_nc_inq_varid (ncid_sst, "lat", &latid);
  wrap_nc_get_var_double (ncid_sst, latid, flat);

  // Check for reduced grid: Currently can only handle full grid

  for (j = 0; j < PLAT; j++)
    nlon[j] = PLON;

  if (nc_inq_varid (ncid_sst, "nlon", &nlonid) == NC_NOERR)  // reduced grid
    wrap_nc_get_var_int (ncid_sst, nlonid, nlon);

  for (j = 0; j < PLAT; j++)
    if (nlon[j] != PLON)
      err_exit ("Cannot currently handle SST on reduced grid\n");

  // Ensure consistent grid characteristics of all input files

  check_grid_consistent (ncid_sst, sstfile, nlon);

  // MLD file

  wrap_nc_inq_varid (ncid_sst, "MLDANN", &mldid);
  wrap_nc_get_vara_double (ncid_sst, mldid, startxy, countxy, (double *) mldann);

  // SST file: must be 12 monthly mean values

  wrap_nc_inq_unlimdim (ncid_sst, &timedimid_sst);
  wrap_nc_inq_dimlen (ncid_sst, timedimid_sst, &timelen);
  
  if (timelen != 12)
    err_exit ("Got %d for length of time dimension on SST file, wanted 12\n", timelen);

  // Time variable id's

  wrap_nc_inq_varid (ncid_sst, "date", &dateid_sst);
  wrap_nc_inq_varid (ncid_sst, "datesec", &datesecid_sst);

  wrap_nc_get_var_int (ncid_sst, dateid_sst, date_sst);
  wrap_nc_get_var_int (ncid_sst, datesecid_sst, datesec_sst);

  // SST and ice fraction id's

  wrap_nc_inq_varid (ncid_sst, "SST_cpl", &sstid);
  wrap_nc_inq_varid (ncid_sst, "ice_cov", &icefracid);

  // loop over months, compute sea ice thickness, volume

  for (n = 0; n < 12; n++) {
    start[0] = n;
    wrap_nc_get_vara_double (ncid_sst, sstid, start, count, (double *) sst[n]);
    wrap_nc_get_vara_double (ncid_sst, icefracid, start, count, (double *) siccnt[n]);

    if (verbose)
      printf ("Read sst and ice conc file %s sst=%g conc=%g\n",
	      sstfile, sst[n][0][0], siccnt[n][PLAT-1][0]);
	      
    mksith (siccnt[n], sicthk[n], flat, nlon, n);

    // sea ice volume

    for (j = 0; j < PLAT; j++) {
      for (i = 0; i < nlon[j]; i++) {
	sicvol[n][j][i] = 0.;                
	if (siccnt[n][j][i] > 0.) {
	  sicvol[n][j][i] = siccnt[n][j][i] * sicthk[n][j][i];
	  
	  if (sicvol[n][j][i] < 0.001 ) 
	    sicvol[n][j][i] = 0.;
	  if (sicvol[n][j][i] > volmax) 
	    sicvol[n][j][i] = volmax;
	}
      }
    }
  }

  // Control run files.
  // Ensure that file name follows assumed naming convention

  if ( check_hname (firstfile, &prefix_len1) != 0)
    err_exit ("First file %s naming format unknown\n", firstfile);

  if ( check_hname (lastfile, &prefix_len2) != 0)
    err_exit ("Last file %s naming format unknown\n", lastfile);

  if (prefix_len1 != prefix_len2 || strncmp (firstfile, lastfile, prefix_len1) != 0)
    err_exit ("%s and %s do not have compatible prefixes\n", firstfile, lastfile);
  
  if (prefix_len1 > 0) {
    prefix = (char *) malloc (prefix_len1 + 1);
    strncpy (prefix, firstfile, prefix_len1);
    prefix[prefix_len1] = '\0';
  }

  strncpy (tmp, &firstfile[prefix_len1], 4);
  tmp[4] = '\0';
  first_yr = atoi (tmp);
    
  strncpy (tmp, &lastfile[prefix_len2], 4);
  tmp[4] = '\0';
  last_yr = atoi (tmp);
    
  strncpy (tmp, &firstfile[prefix_len1+5], 2);
  tmp[2] = '\0';
  first_mo = atoi (tmp);
  
  strncpy (tmp, &lastfile[prefix_len2+5], 2);
  tmp[2] = '\0';
  last_mo = atoi (tmp);
  
  // accumulate net flux for control run

  for (n = 0; n < 12; n++)
    for (j = 0; j < PLAT; j++)
      for (i = 0; i < nlon[j]; i++) {
	fnet[n][j][i] = 0.;
      }

  // Always grab the 1st timeslice from monthly averaged file
  // i.e. assume each file has only a single timeslice on it

  start[0] = 0;  

  // Loop over control run files, computing Q fluxes.

  yr = first_yr;
  mo = first_mo;
  memset (caseid, 0, sizeof caseid);

  for (; date (yr, mo) <= date (last_yr, last_mo); incdate (&yr, &mo)) {
    fname = bldfn (prefix, yr, mo);
    if ((ret = nc_open (fname, NC_NOWRITE, &ncid_cntl)) != NC_NOERR)
      err_exit ("Error opening control case file %s\n%s", fname, nc_strerror (ret));

    // Ensure case names do not change

    if (nc_inq_attlen (ncid_cntl, NC_GLOBAL, "case", &caselen) != NC_NOERR)
      err_exit ("Error inquiring case name of control run\n");

    if (caselen > (sizeof caseid) - 1)
      err_exit ("Length of case name on %s too long\n", fname);

    wrap_nc_get_att_text (ncid_cntl, NC_GLOBAL, "case", caseid);

    if (yr != first_yr || mo != first_mo)
      if (strcmp (oldcaseid, caseid))
	err_exit ("Case names do not match: %s vs. %s\n", oldcaseid, caseid);
    
    strcpy (oldcaseid, caseid);        // Save case id for comparison next iteration

    check_dims_consistent (ncid_cntl, fname);

    // Read the data needed to compute fluxes
    // JR 12/10/02 changed field names (added "OI") to use only non-land
    // JR portion of fluxes from control run.

    wrap_nc_inq_varid (ncid_cntl, "FSNSOI", &fsnsid);
    wrap_nc_inq_varid (ncid_cntl, "FLNSOI", &flnsid);
    wrap_nc_inq_varid (ncid_cntl, "LHFLXOI", &lhflxid);
    wrap_nc_inq_varid (ncid_cntl, "SHFLXOI", &shflxid);

    wrap_nc_get_vara_double (ncid_cntl, fsnsid, start, count, (double *) fsns);
    wrap_nc_get_vara_double (ncid_cntl, flnsid, start, count, (double *) flns);
    wrap_nc_get_vara_double (ncid_cntl, lhflxid, start, count, (double *) lhflx);
    wrap_nc_get_vara_double (ncid_cntl, shflxid, start, count, (double *) shflx);

    if (verbose)
      printf ("Read fluxes file %s fsns=%g flns=%g lhflx=%g shflx=%g\n",
	      fname, fsns[9][0], flns[9][0], lhflx[9][0], shflx[9][0]);
	      
    ++moaccum[mo-1];     // increment per-month counter for fields just read in

    for (j = 0; j < PLAT; j++)
      for (i = 0; i < nlon[j]; i++)
	fnet[mo-1][j][i] += fsns[j][i] - flns[j][i] - lhflx[j][i] - shflx[j][i];

    wrap_nc_close (ncid_cntl);
  }

  // compute average flux

  diffaccum = false;
  for (n = 0; n < 12; n++) {
    if (moaccum[n] == 0)
      err_exit ("Error: counter for month %d is zero\n", n);
      
    npls = (n + 1) % 12;

    if (moaccum[n] != moaccum[npls])
      diffaccum = true;

    for (j = 0; j < PLAT; j++)
      for (i = 0; i < nlon[j]; i++) {
	fnet[n][j][i] /= moaccum[n];
      }
  }

  if (diffaccum) {
    printf ("NOTE: Accumulation counts do not match for all months of control run:\n");
    for (n = 0; n < 12; n++)
      printf ("Month %d=%d\n", n+1, moaccum[n]);
  }

  // Compute QFLUX

  for (n = 0; n < 12; n++) {
    npls = (n + 1) % 12;
    nmin = n == 0 ? 11 : n - 1;

    //JR: 11/6/02 changed dayfact definition for consistency per BPB

    dayfact = 2.*daysmo[n]*86400.;

    if (verbose)
      printf ("Using SST mid-months %d and %d to compute QFLUX mid-month %d\n dayfact=%g\n",
	      nmin, npls, n, dayfact/86400.);

    for (j = 0; j < PLAT; j++)
      for (i = 0; i < nlon[j]; i++) {

	// h2oterm will go away over points completely covered by ice.
	// BPB: 10/25/02 changed delsst to be more accurate
	// JR: 11/4/02 changed ice conc. and sst terms to bound properly for time diddling.
	//     Matches prescribed ice calc. in model.

	siccntnpls = MIN (MAX (siccnt[npls][j][i], 0.), 1.);
	siccntnmin = MIN (MAX (siccnt[nmin][j][i], 0.), 1.);
	sstnpls    = MAX (sst[npls][j][i], tsice);
	sstnmin    = MAX (sst[nmin][j][i], tsice);
	delsst    = (1. - siccntnpls)*sstnpls -  (1. - siccntnmin)*sstnmin;
	delsicvol   = sicvol[npls][j][i] - sicvol[nmin][j][i];
	h2oterm[n][j][i] = (rhoocn * cpocn * mldann[j][i] * delsst) / dayfact;
	iceterm[n][j][i] = (li * delsicvol) / dayfact;
	qflux[n][j][i] = fnet[n][j][i] - (h2oterm[n][j][i] - iceterm[n][j][i]);
      }
  }

  // Modify history attribute

  wrap_nc_redef (ncid_sst);
  if (nc_inq_attlen (ncid_sst, NC_GLOBAL, "history", &hislen) == NC_NOERR) {
    history = (char *) malloc (hislen + 1);
    wrap_nc_get_att_text (ncid_sst, NC_GLOBAL, "history", history);
    //    printf ("history before new stuff:\n%s\n", history);
    hislen = strlen (history) + strlen (cmdline) + 2;    // 2 is for newline and \0
    realloc (history, hislen);
    //    printf ("history after realloc:\n%s\n", history);
    strcat (history, "\n");
    strcat (history, cmdline);
    //    printf ("history after new stuff:\n%s\n", history);
  } else {
    history = cmdline;
    hislen = strlen (history) + 1;
  }
  wrap_nc_put_att_text (ncid_sst, NC_GLOBAL, "history", hislen, history);
      
  // Overwrite or define control_case attribute

  wrap_nc_put_att_text (ncid_sst, NC_GLOBAL, "control_case", caselen, caseid);
  wrap_nc_enddef (ncid_sst);

  // Add QFLUX to SST file

  if ((ret = nc_inq_varid (ncid_sst, "QFLUX", &qfluxid)) != NC_NOERR) {
    wrap_nc_inq_dimid (ncid_sst, "lon", &lonid);
    wrap_nc_inq_dimid (ncid_sst, "lat", &latid);

    dimids[0] = timedimid_sst;
    dimids[1] = latid;
    dimids[2] = lonid;

    xydimids[0] = latid;
    xydimids[1] = lonid;

    wrap_nc_redef (ncid_sst);
    wrap_nc_def_var (ncid_sst, "QFLUXPREDIDDLE", NC_DOUBLE, 3, dimids, &qfluxprediddleid);
    wrap_nc_def_var (ncid_sst, "QFLUX",       NC_DOUBLE, 3, dimids, &qfluxid);
    wrap_nc_def_var (ncid_sst, "H2OTERM",     NC_DOUBLE, 3, dimids, &h2otermid);
    wrap_nc_def_var (ncid_sst, "ICETERM",     NC_DOUBLE, 3, dimids, &icetermid);
    wrap_nc_def_var (ncid_sst, "FNET",        NC_DOUBLE, 3, dimids, &fnetid);
    wrap_nc_def_var (ncid_sst, "SICTHK",      NC_DOUBLE, 3, dimids, &sicthkid);
    wrap_nc_def_var (ncid_sst, "SICCNT",      NC_DOUBLE, 3, dimids, &siccntid);
    wrap_nc_def_var (ncid_sst, "DIAG2",       NC_DOUBLE, 3, dimids, &diag2id);

    wrap_nc_def_var (ncid_sst, "QFLUXYAVG",      NC_DOUBLE, 2, xydimids, &qfluxyavgid);
    wrap_nc_def_var (ncid_sst, "FNETYAVG",    NC_DOUBLE, 2, xydimids, &fnetyavgid);
    wrap_nc_def_var (ncid_sst, "DIAG1",       NC_DOUBLE, 2, xydimids, &diag1id);

    wrap_nc_enddef (ncid_sst);
  }

  for (j = 0; j < PLAT; j++)
    for (i = 0; i < PLON; i++) {
      fnetyavg[j][i] = 0.;
      qfluxyavg[j][i]   = 0.;
    }

  // Write out monthly averaged data and compute yearly averages

  timesum = 0.;

  for (n = 0; n < 12; n++) {
    start[0] = n;
    wrap_nc_put_vara_double (ncid_sst, qfluxprediddleid, start, count, (double *) qflux[n]);
    wrap_nc_put_vara_double (ncid_sst, h2otermid, start, count, (double *) h2oterm[n]);
    wrap_nc_put_vara_double (ncid_sst, icetermid, start, count, (double *) iceterm[n]);
    wrap_nc_put_vara_double (ncid_sst, fnetid, start, count, (double *) fnet[n]);
    wrap_nc_put_vara_double (ncid_sst, sicthkid, start, count, (double *) sicthk[n]);
    wrap_nc_put_vara_double (ncid_sst, siccntid, start, count, (double *) siccnt[n]);

    for (j = 0; j < PLAT; j++)
      for (i = 0; i < PLON; i++) {
	fnetyavg[j][i] += fnet[n][j][i]*daysmo[n];
	qfluxyavg[j][i]   += qflux[n][j][i]*daysmo[n];
      }

    timesum += daysmo[n];
  }

  for (j = 0; j < PLAT; j++)
    for (i = 0; i < PLON; i++) {
      fnetyavg[j][i] /= timesum;
      qfluxyavg[j][i]   /= timesum;
    }
  
  // Write out yearly averages

  wrap_nc_put_vara_double (ncid_sst, fnetyavgid, startxy, countxy, (double *) fnetyavg);
  wrap_nc_put_vara_double (ncid_sst, qfluxyavgid, startxy, countxy, (double *) qfluxyavg);

  // Apply time diddling to QFLUX 

  for (j = 0; j < PLAT; j++)
    for (i = 0; i < PLON; i++)
      modelqfluxyavg[j][i] = 0.;

  // Loop over 1 year of nsteps to compute yearly averaged QFLUX as would be seen by model

  for (nstep = 0; nstep < ntspy; nstep++) {
    (void) getprvnxtindx (nstep, dtime, &nmin, &npls);
    (void) timefactors (nstep, dtime, nmin, npls, &fact1, &fact2);

    if (nstep == 0)
      printf (" fact1,fact2=%20.16f %20.16f\n", fact1, fact2);

    for (j = 0; j < PLAT; j++)
      for (i = 0; i < PLON; i++) {
	modelqflux = qflux[nmin][j][i]*fact1 + qflux[npls][j][i]*fact2;
	modelqfluxyavg[j][i] += modelqflux;
      }
  }

  // Compute Briegleb time diddling factor

  for (j = 0; j < PLAT; j++)
    for (i = 0; i < PLON; i++) {
      modelqfluxyavg[j][i] /= ntspy;
      diddle[j][i] = qfluxyavg[j][i] - modelqfluxyavg[j][i];
    }

  // Apply time diddling factor and write output QFLUX and its negative (QFLUX)
  // Use either Briegleb's procedure for preserving the yearly average, or
  // Taylor's procedure for preserving monthly averages.

  if (preserve_yavg) {
    for (n = 0; n < 12; n++) {
      
      for (j = 0; j < PLAT; j++)
	for (i = 0; i < PLON; i++) {
	  qflux[n][j][i] += diddle[j][i];
	}

      start[0] = n;
      wrap_nc_put_vara_double (ncid_sst, qfluxid, start, count, (double *) qflux[n]);
    }

  } else if (preserve_mavg) {

    timediddle_mavg (qflux);

    for (n = 0; n < 12; n++) {
      start[0] = n;
      wrap_nc_put_vara_double (ncid_sst, qfluxid, start, count, (double *) qflux[n]);
    }
  }

  // Loop over 1 year of nsteps again and compute a monthly averaged QFLUX that should
  // match the pre-time-diddled value.

  for (j = 0; j < PLAT; j++)
    for (i = 0; i < PLON; i++)
      diag1[j][i] = 0.;

  for (n = 0; n < 12; n++) {
    nacs[n] = 0;
    for (j = 0; j < PLAT; j++)
      for (i = 0; i < PLON; i++)
	diag2[n][j][i] = 0.;
  }

  for (nstep = 0; nstep < ntspy; nstep++) {
    n = getprvnxtindx (nstep, dtime, &nmin, &npls);
    ++nacs[n];
    (void) timefactors (nstep, dtime, nmin, npls, &fact1, &fact2);
    for (j = 0; j < PLAT; j++)
      for (i = 0; i < PLON; i++) {
	modelqflux = qflux[nmin][j][i]*fact1 + qflux[npls][j][i]*fact2;
	diag1[j][i]    += modelqflux;
	diag2[n][j][i] += modelqflux;
      }
  }

  for (j = 0; j < PLAT; j++)
    for (i = 0; i < PLON; i++) 
      diag1[j][i] /= ntspy;

  wrap_nc_put_vara_double (ncid_sst, diag1id, startxy, countxy, (double *) diag1);

  for (n = 0; n < 12; n++) {

    for (j = 0; j < PLAT; j++)
      for (i = 0; i < PLON; i++)
	diag2[n][j][i] /= nacs[n];

    start[0] = n;
    wrap_nc_put_vara_double (ncid_sst, diag2id, start, count, (double *) diag2[n]);
  }

  wrap_nc_close (ncid_sst);
  return 0;
}

void usage_exit (char *argv0)
{
  err_exit ("Usage: %s -f firstfile_fromcontrolrun -l lastfile_fromcontrolrun"
	    " -s SSTFILE [-v]\n", argv0);
}

/*
** date: convert input year and month into a single int (e.g.: 199502)
*/

int date (int yr, int mo)
{
  return yr*100 + mo;
}

/*
** bldfn:convert input year and month into a string of the form yyyy-mo.nc
*/

char *bldfn (char *prefix, int yr, int mo)
{
  char *newfile;

  // Create character string yyyy-mo.nc

  newfile = (char *) malloc (strlen (prefix) + 11);
  sprintf (newfile, "%s%4.4d-%2.2d.nc", prefix, yr, mo);
  return newfile;
}
  
/*
** incdate: increment year and month by a single month
*/

void incdate (int *yr, int *mo)
{
  if (*mo == 12) {
    *yr += 1;
    *mo = 1;
  } else {
    *mo += 1;
  }
}

void AddToString (char *cmdline,
		  char *arg,
		  char *setting)
{
  realloc (cmdline, strlen (cmdline) + strlen (arg) + strlen (setting) + 1);
  strcat (cmdline, arg);
  strcat (cmdline, setting);
}

// check_hname: check that input character string follows
// <anything>xxxx-xx.nc naming convention, where x is a digit. 

int check_hname (char *name,          
		 int *prefix_len)
{
  int i;
  int len;

  len = strlen (name);
  if (len < 10)
    return -1;

  i = len - 10;
  if (isdigit (name[i]) && isdigit (name[i+1]) && isdigit (name[i+2]) && 
      isdigit (name[i+3]) && name[i+4] == '-' && isdigit (name[i+5]) && 
      isdigit (name[i+6]) && strcmp (&name[i+7], ".nc") == 0) {
    *prefix_len = i;
    return 0;
  }
  return -1;
}

// getprvnxtindx: find month indices bracketing a given model timestep
// return value is current month index (zero-based)

int getprvnxtindx (int nstep,   // current "model" timestep
		   int dtime,   // timestep in seconds
		   int *prv,    // previous index
		   int *nxt)    // next index
{
  int n;
  int ntspdy;                 // number of timesteps per day
  int curmoindx = -1;         // current month index
  int daysum;
  double daysintegrated;
  double timefile;

  if (86400 % dtime != 0) {
    printf ("getprvindx: dtime must divide evenly into 86400\n");
    exit (1);
  }
    
  ntspdy         = 86400 / dtime;
  nstep          = nstep % (365*ntspdy);  // modulate nstep to this year
  daysintegrated = (nstep*dtime)/86400.;

  // Determine current month index

  daysum = 0;
  for (n = 0; n < 12; n++) {
    daysum += daysmo[n];
    if (daysum > daysintegrated) {
      curmoindx = n;
      break;
    }
  }

  if (curmoindx < 0) {
    printf ("getprvnxtindx: cannot find curmoindx for nstep %d\n", nstep);
    exit (1);
  }

  for (n = 0; n < 12; n++) {
    timefile = daysof (date_sst[n], datesec_sst[n]);
    if (timefile > daysintegrated) {
      if (n == 0) {
	*prv = 11;
	*nxt = 0;
      } else {
	*prv = n - 1;
	*nxt = n;
      }
      return curmoindx;
    }
  }

  // If get to here, must be beyond date of last time sample (typically mid-December)

  if (daysintegrated < 365.) {
    *prv = 11;
    *nxt = 0;
  } else {
    printf ("getprvnxtindx: cannot find indices for nstep %d\n", nstep);
    exit (1);
  }
  return curmoindx;
}

double daysof (int ymd, int datesec)
{
  int moindx;
  int n;

  double total;
  
  if ((moindx = ((ymd/100) % 100) - 1) < 0 || moindx > 11)
    err_exit ("modindx1=%d\n", moindx);

  total = 0.;
  for (n = 0; n < moindx; n++) 
    total += daysmo[n];

  total += ymd % 100 - 1;
  total += datesec/86400.;
  return total;
}

void timefactors (int nstep, 
		  int dtime, 
		  int nmin, 
		  int npls, 
		  double *fact1, 
		  double *fact2)
{
  double daysintegrated;          // model integration time
  double days1, days2;            // valid times of monthly values on dataset

  daysintegrated = (nstep*dtime)/86400.;
  days1 = daysof (date_sst[nmin], datesec_sst[nmin]);
  days2 = daysof (date_sst[npls], datesec_sst[npls]);

  // handle wrapped situation

  if (nmin == 11 && npls == 0) {
    if (daysintegrated < days1)
      days1 -= 365.;
    else if (daysintegrated > days2)
      days2 += 365.;
  } else if (nmin >= npls) {
    printf ("timefactors: bad nmin, npls=%d %d\n", nmin, npls);
    exit (1);
  }
  
  *fact1 = (days2 - daysintegrated) / (days2 - days1);
  *fact2 = (daysintegrated - days1) / (days2 - days1);
  
  if (*fact1 < 0. || *fact1 > 1. || *fact2 < 0. || *fact2 > 1. || 
      fabs (*fact1 + *fact2 - 1.) > 1.e-6) {
    printf ("timefactors: bad fact1,fact2= %f %f\n", *fact1, *fact2);
    exit (1);
  }
}



