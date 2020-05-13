/*------------------------------------------------------------------------*
 * File: dataset.cpp 
 * $Author: jet $
 * $Id: dataset.cpp,v 1.1.6.1 2004/05/18 20:51:51 jet Exp $ *
 *------------------------------------------------------------------------*/

#ifndef lint
static char rcsid[] = "$Id: dataset.cpp,v 1.1.6.1 2004/05/18 20:51:51 jet Exp $";
#endif /* lint */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "realtype.h"
#include "dataset.h"
#include "runtype.h"
#include "msgdlg.h"
#include "manager.h"
#include "MainWndImpl.h"

#include <iostream>

#ifndef _RWSTD_NO_NAMESPACE
  using namespace std;
#endif

// the reqVars must follow the order layed out in runtype.h
const char* reqVars[][MAX_DATASET_VARS] = {
      // analysis 
    { "lon", "lat", "lev", "time", "date", "datesec","hyam", 
      "hyai", "hybi", "hybm", "PHIS", "U", "V", "T", "PS", "OMEGA" },
      // model
    { "lon", "lat", "lev", "time", "date", "datesec", "hyam", "hyai", 
      "hybi", "hybm", "PHIS", "U", "V", "T", "PS", "TS1", "TS2", "TS3", "TS4", "OMEGA" },
      // iop
    { "lon","lat","lev","phis","q","divT","divq","Ps","omega" },
      // lsminit
    { "mcdate", "mcsec", "LANDVECIXY","LANDVECJXY","T_VEG_INI","T_GRND_INI","H2OCAN_INI","H2OSNO_INI","SNOWDP_INI","SNOWAGE_INI","SNLSNO_INI","T_SOISNO_INI","T_LAKE_INI","H2OSOI_LIQ_INI","H2OSOI_ICE_INI","ZSNO_INI","DZSNO_INI","ZISNO_INI" },
      // lsmsurf
    {"LONGXY","LATIXY","LANDMASK","LANDFRAC","SOIL_COLOR","PCT_SAND","PCT_CLAY","PCT_WETLAND","PCT_LAKE","PCT_GLACIER","PCT_URBAN","PFT","PCT_PFT","MONTHLY_LAI","MONTHLY_SAI","MONTHLY_HEIGHT_TOP","MONTHLY_HEIGHT_BOT"},
      // ozone
    { "lat", "lev", "time", "date","datesec","OZONE" }, 
      // press
    { "hyam", "hyai", "hybi", "hybm" },
      // sst
    { "date", "datesec", "SST_cpl","ice_cov" },
      // absems
    { "ah2onw","eh2onw","ah2ow","ln_ah2ow","cn_ah2ow","ln_eh2ow","cn_eh2ow"},
      // aeroptics
    {   "wvl", "wvn", "wvl_sw", "wvn_sw",
        "wvl_ccm", "wvl_min_ccm", "wvl_max_ccm", "rh",
        "ext_suso", "ssa_suso", "asm_suso", "ext_ccm_suso",
        "ssa_ccm_suso", "asm_ccm_suso", "ext_waso", "ssa_waso",
        "asm_waso", "ext_ccm_waso", "ssa_ccm_waso", "asm_ccm_waso",
        "ext_ssam", "ssa_ssam", "asm_ssam", "ext_ccm_ssam",
        "ssa_ccm_ssam", "asm_ccm_ssam", "ext_sscm", "ssa_sscm",
        "asm_sscm", "ext_ccm_sscm", "ssa_ccm_sscm", "asm_ccm_sscm",
        "ext_inso", "ssa_inso", "asm_inso", "ext_ccm_inso",
        "ssa_ccm_inso", "asm_ccm_inso", "ext_soot", "ssa_soot",
        "asm_soot", "ext_ccm_soot", "ssa_ccm_soot", "asm_ccm_soot",
        "ext_minm", "ssa_minm", "asm_minm", "ext_ccm_minm",
        "ssa_ccm_minm", "asm_ccm_minm", "ext_miam", "ssa_miam",
        "asm_miam", "ext_ccm_miam", "ssa_ccm_miam", "asm_ccm_miam",
        "ext_micm", "ssa_micm", "asm_micm", "ext_ccm_micm",
        "ssa_ccm_micm", "asm_ccm_micm", "ext_mitr", "ssa_mitr",
        "asm_mitr", "ext_ccm_mitr", "ssa_ccm_mitr", "asm_ccm_mitr",
        "ext_dst1","ssa_dst1" },

      // aermass
    {   "PS","SO4","SSLT","OCPHO","OCPHI","BCPHO",
	"BCPHI","DSTQ01","DSTQ02","DSTQ03","DSTQ04" },
      // LSMPFT
    {},
      // user
    { "bdate", "lon", "lat", "lev", "tsec" },
      // saved initial conditions
    { "lon", "lat", "lev", "bdate", "tsec", "modelfile", "analysisfile", "iopfile", 
      "ozonfile",  "sstfile", "lsminifile", "switches", "runtype" }
};

const char* impIopVars[]  = { "Ptend", "u", "v","Tg","Tsair" };

const char* datasetDesc[] = { "Global Analysis", 
			      "Global Model", 
			      "IOP", 
			      "LSM initial conditions",
			      "LSM Surfdat", 
			      "Ozone", 
			      "Pressure Levels", 
			      "Sea-Surface Temp", 
			      "Absortivity Emissivity",
			      "Aersol Optics",
			      "Aersol Mass",
			      "CLM PFT",
			      "User",
			      "Saved Initial Conditions" }; 
// these are the names given to the dataset types in scam defaults files
const char* datasetType[] = { "analysisfile", 
			      "modelfile", 
			      "iopfile", 
			      "lsminifile",
			      "lsmsurffile", 
			      "ozonfile",
			      "pressfile", 
			      "sstfile",
			      "absemsfile",
			      "aeropticsfile",
			      "aermassfile",
			      "lsmpftfile",
 			      "userfile",
			      "sicfile" }; 

bool
Dataset::CheckVars( string& missingVars )
{
    bool  status = TRUE;
    
    for ( int n=0; reqVars[filetype][n] != 0; n++) {
        if( hasVariable( reqVars[filetype][n] ) )
            continue;
        
          // missing divq/divt presents special case!
          // they can be substituted by 3d versions
        if ( reqVars[filetype][n] == "divq" ) 
            if ( hasVariable( "divq3d" ) )
                continue;

        if ( reqVars[filetype][n] == "divT" )
            if ( hasVariable( "divT3d" ) )
                continue;
        
        if( !missingVars.empty() )
            missingVars += ", ";            
        missingVars += reqVars[filetype][n];
        status = FALSE;
    }
    
      //
      // check for "important" IOP variables 
      //
    if ( filetype == IOP ) {
        for ( size_t n=0; n < sizeof( impIopVars )/sizeof( string ); n++) {
            if( !hasVariable( impIopVars[n] ) ) {
                if( !missingVars.empty() )
                    missingVars += ", ";            
                missingVars += impIopVars[n];
            }
        }
    }
    return status;
}

bool
Dataset::CheckLat( real_t lat, int latidx )
{
    return ( FindIdx( "lat", lat ) == latidx );
}

bool
Dataset::CheckLon( real_t lon, int lonidx )
{
    return ( FindIdx( "lon", lon ) == lonidx );
}

int
Dataset::FindIdx( string varname, real_t value )
{
    int idx;
    real_t  next;
    
      // find the nearest index 

    try {
        NcVariable<real_t> var = variable( varname );
        var.read();
        idx = 0;
        for ( unsigned i=0; i < var.size(); i++ ) {
            next = var[i];
            if ( fabs( value - next ) < fabs( value - var[idx] ) ) 
                idx = i;
        }
        return idx;
    } catch ( NcErr& e ) {
        cerr << "ERROR: "__FILE__":" << __LINE__ 
             << " Dataset::FindIdx(): " << e << endl;
        exit( -1 );
    }
}

string 
Dataset::TypeDescription( int type )
{
    return datasetDesc[type];
}

string 
Dataset::TypeName( int type )
{
    return datasetType[type];
}




