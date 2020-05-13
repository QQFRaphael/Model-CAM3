/*------------------------------------------------------------------------*
 * File: ncarg.cpp 
 * $Author: jet $
 * $Id: ncarg.cpp,v 1.1.6.1 2004/05/18 20:52:00 jet Exp $ *------------------------------------------------------------------------*/

#ifndef lint
static char rcsid[] = "$Id: ncarg.cpp,v 1.1.6.1 2004/05/18 20:52:00 jet Exp $";
#endif /* lint */


#ifndef NO_NCARG
#include "realtype.h"
#include <sstream>  
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#ifndef _RWSTD_NO_NAMESPACE
  using namespace std;
#endif


#include <ncarg/ncargC.h>
#include <ncarg/gks.h>

#include <netcdf.h>
#include "ncarg.h"
#include "timeconvert.h"
#include "max.h"
#include "msgdlg.h"

#define PREF_NUM_CONTOURS 10
#define MAX_NUM_CONTOURS (2*PREF_NUM_CONTOURS)
#define MAX_LABELS 11
#define MAP_SIZE  9600000
#define AREA_SIZE 400000
#define GROUP_SIZE 20
#define LBL_FONTSIZE .015

extern "C" {
    
struct rgb {
    double red, green, blue;
};
static struct rgb std_rgb[] = { { 1.0, 1.0, 1.0 },  /* white */
                                { 0.0, 0.0, 1.0 },  /* blue */
                                { 0.0, 0.4, 0.0 },  /* green */
                                { 1.0, 1.0, 0.5 },  /* yellow*/
                                { 1.0, 0.2, 0.2 },  /* red */
                                { 0.0, 0.0, 0.0 } }; /* black */
    
static char input_file[MAX_PATH_LEN];
static char field_name[256];
static char field_long_name[256];
static char field_units[256];
static float field_multiplier;
static int  ncid;
static int  start_index;
static int  stop_index;
static int   ntime;
static int basedate;
static float* hours;             /* time data (in hours) from input file */
static int   time_format;
static float steps_per_hour;
static bool GksIsOpen = false;

bool
findBoundingArrayIndices(float value, float *array, int array_size
			     ,int *lowIdx,int *hiIdx)
{
  int i;
  *lowIdx = *hiIdx = 0;

  if(value < array[0])
    return(false);
  if(value > array[array_size - 1])
    return(false);

  for(i = 0;i<array_size;i++) {
      if(value >= array[i] && value < array[i+1]) {
	  *lowIdx = i;
	  *hiIdx = i + 1;
	  return(true);
      }
  }
  return(false);
}

/*-----------------------------------------------------------------------
 * Interpolate
 *
 *  Interpolates a value between points described by x1,y1 and x2,y2.
 *
 * INPUTS:
 *
 *        x1, y1              Point 1.
 *        x2, y2              Point 2.
 *        x                   The value between x1 and x2.
 *
 * OUTPUTS:
 *
 *        Returns the interpolated value between y1 and y2.
 */
float Interpolate(float x1,float y1,float x2,float y2,float x)
{
  return (((y2 - y1)/(x2 - x1)) * (x - x1) + y1); 
}



static int
colram(float *xcs,
       float *ycs,
       int *ncs,
       int *iai,
       int *iag,
       int *nai)
{
    Gint ifll;
    int i;
    Gpoint_list point_list;

    /*
     * Assume the polygon will be filled until we find otherwise.
     */
    ifll = 1;

    /*
     * If any of the area identifiers is negative, don't fill the polygon.
     */
    for (i = 0; i < *nai; i++) {
        if (iai[i] < 0)
            ifll = 0;
    }

    /*
     * Otherwise, fill the polygon in the color implied by its area
     * identifier relative to edge group 3 (the contour-line group).
     */
    if (ifll != 0) {
        ifll = 0;
        for (i = 0; i < *nai; i++) {
            if (iag[i] == 3)
                ifll = iai[i];
        }
        if (ifll > 0 && ifll <= MAX_NUM_CONTOURS) {
            gset_fill_colr_ind(ifll + 1);

            point_list.points = (Gpoint *) malloc(sizeof(Gpoint) * (*ncs));
            point_list.num_points = *ncs;
            for (i = 0; i < *ncs; i++) {
                point_list.points[i].x = xcs[i];
                point_list.points[i].y = ycs[i];
            }
            gfill_area(&point_list);
        }
    }

    return 1;
}


void OpenGKS() {
    if ( GksIsOpen )            // only can open if currently not open
        return;
    Gasfs list_asf[] = { { GASF_INDIV },{ GASF_INDIV },{ GASF_INDIV },
                         { GASF_INDIV },{ GASF_INDIV },{ GASF_INDIV },
                         { GASF_INDIV },{ GASF_INDIV },{ GASF_INDIV },
                         { GASF_INDIV },{ GASF_INDIV } };
    gopen_gks("stdout", 0);
    /*
     * Set some GKS parameters
     */
    gset_clip_ind(GIND_NO_CLIP);       /* turn off clipping              */
    gset_asfs(list_asf);               /* all source flags == individual */
    gset_fill_int_style(GSTYLE_SOLID); /* force solid filling            */
    GksIsOpen = true;
}

void
CloseGKS() 
{
    if ( GksIsOpen )
        gclose_gks();
    GksIsOpen = false;
}

void
ActivateOutputDevice( int deviceID) {
    gactivate_ws(deviceID);
}

void
ClearOutputDevice( int deviceID) {
    gclear_ws(deviceID, (Gctrl_flag) 0);
}

void
DeactivateOutputDevice( int deviceID) {
    gdeactivate_ws(deviceID);
}


/*-----------------------------------------------------------------------
 * OpenOutputDevice
 *
 * INPUTS:
 *
 * OUTPUTS:
 */
int
OpenOutputDevice( int deviceID, int type )
{
    char* connectionID = NULL;
    Gint wkstype;

    wkstype = type;

    gopen_ws( deviceID, connectionID, wkstype);

    return deviceID;
}

/*-----------------------------------------------------------------------
 * CloseOutputDevice
 *
 * INPUTS: none
 *
 * OUTPUTS: none
 */
void
CloseOutputDevice( int id)
{
        gclose_ws(id);
}

/*-----------------------------------------------------------------------
 * ShowFieldLinePlot
 *
 * INPUTS: none
 *
 * OUTPUTS: none
 */
static void
ShowFieldLinePlot(float *Data)
{
    float *local_data;
    float *time_labels;  /*  label values */
    int i;
    int numTimeSlices = stop_index - start_index;


    local_data = (float *) malloc(sizeof(float) * numTimeSlices );
    time_labels = (float *) malloc(sizeof(float) * numTimeSlices );
    
    for ( i = 0; i < numTimeSlices; i++ ) 
        local_data[i] = Data[i+start_index];


      /* 
       * Bottom Label ( time units )
       */
    c_agsetc("LABEL/NAME.", "B");     
    c_agseti("LINE/NUMBER.", -100);
    c_agsetf("LINE/CHARACTER.", .020);
    c_agsetf("AXIS/BOTTOM/NUMERIC/WIDTH.", .02);
    switch( time_format ) {
    case TimeConverter::STEPS:
        for ( i = 0; i < numTimeSlices; i++ )
            time_labels[i] = (int)( hours[i+start_index] * steps_per_hour);
         c_agsetc("LINE/TEXT.", "Steps$");
         break;
    case TimeConverter::HOURS:
        for ( i = 0; i < numTimeSlices; i++ )
            time_labels[i] = hours[i+start_index];
        c_agsetc("LINE/TEXT.", "Hours$");
        break;
    case TimeConverter::DAYS:
    case TimeConverter::DATE:
        // NCAR graphics can't handle dates properly, so just display days
        for ( i = 0; i < numTimeSlices; i++ )
            time_labels[i] = hours[i+start_index] / 24.0;
        c_agsetc("LINE/TEXT.", "Days$");
        break;
    default:
        cerr <<  "ERROR: "__FILE__":" << __LINE__ 
             << " : ShowFieldLinePlot() - bad time format: " << time_format <<endl;
        exit( -1 );
    }

      /*
       * Left Label  ( field units )
       */
    c_agsetc("LABEL/NAME.", "L");       
    c_agseti("LINE/NUMBER.", 100);
    c_agsetf("LINE/CHARACTER.", .02);
    c_agsetf("AXIS/LEFT/NUMERIC/WIDTH.", .02);
    c_agsetc("LINE/TEXT.", field_units);

      /* 
       * Top Label ( field long name )
       */
    c_agsetc("LABEL/NAME.", "T");       
    c_agseti("LINE/NUMBER.", 100);
    c_agsetf("LINE/CHARACTER.", .04);
    c_agsetc("LINE/TEXT.", field_long_name);
    
      /*
       * Filename Label ( at bottom )
       */ 
    {     
        float vl, vr, vb, vt;
        float wl, wr, wb, wt;
        int lf;
        c_getset(&vl, &vr, &vb, &vt, &wl, &wr, &wb, &wt, &lf); /* save original values */
        c_set(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1); /* set viewport to whole window */
        c_pchiqu( .5, .02, input_file, LBL_FONTSIZE, 0.0, 0.0); /* print the label */
        c_set(vl, vr, vb, vt,wl, wr, wb, wt, lf); /* reset to the original values */
    }

      /* next 3 calls take the place of a call to ezxy */
      /* this is necessary because ezxy includes a call to c_frame() */
      /* which causes the plot to wait for a mouse click before returning */
    c_agstup( time_labels, 1, 0, numTimeSlices, 1,
              local_data, 1, 0, numTimeSlices, 1 );
    c_agback();
    c_agcurv( time_labels, 1, local_data, 1, numTimeSlices, 1 );

    
    free((char *)time_labels);
    free((char *)local_data);
}

/*-----------------------------------------------------------------------
 * ShowFieldAverage
 *
 * INPUTS: none
 *
 * OUTPUTS: none
 */
static void
ShowFieldAverage(float *Data, float* lev_data,
                 int nlev)
{
    float *ydra;
    int i;
    
    ydra = (float *) malloc(sizeof(float) * nlev);
    for (i = 0; i < nlev; i++) 
        ydra[i] = lev_data[i];
    
      /* 
       * Top Label (field_long_name)
       */
    c_agsetc("LABEL/NAME.", "T");       
    c_agseti("LINE/NUMBER.", 100);
    c_agsetf("LINE/CHARACTER.", .04);
    c_agsetc("LINE/TEXT.", field_long_name );

      /* 
       * Bottom Label (field_units)
       */
    c_agsetc("LABEL/NAME.", "B");  
    c_agseti("LINE/NUMBER.", -100);
    c_agsetf("LINE/CHARACTER.", .02);
    c_agsetf("AXIS/BOTTOM/NUMERIC/WIDTH.", .02);
    c_agsetc("LINE/TEXT.", field_units );

      /* 
       * Left Label (pressure)
       */
    c_agsetc("LABEL/NAME.", "L");     
    c_agseti("LINE/NUMBER.", 100);
    c_agsetf("LINE/CHARACTER.", .02);
    c_agsetf("AXIS/LEFT/NUMERIC/WIDTH.", .02);
    c_agsetc("LINE/TEXT.", "Pressure (mb)$");
    
      /*
       * Filename Label ( at bottom )
       */ 
    {     
        float vl, vr, vb, vt;
        float wl, wr, wb, wt;
        int lf;
        c_getset(&vl, &vr, &vb, &vt, &wl, &wr, &wb, &wt, &lf); /* save original values */
        c_set(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1); /* set viewport to whole window */
        c_pchiqu( .5, .02, input_file, LBL_FONTSIZE, 0.0, 0.0); /* print the label */
        c_set(vl, vr, vb, vt,wl, wr, wb, wt, lf); /* reset to the original values */
    }

    /* next 3 calls take the place of a call to ezxy */
    /* this is necessary because ezxy includes a call to c_frame() */
    /* which causes the plot to wait for a mouse click before returning */
    c_agstup( Data, 1, 0, nlev, 1, lev_data, 1, 0, nlev, 1 );
    c_agback();
    c_agcurv( Data, 1, lev_data, 1, nlev, 1 ) ;

    free( ydra);
}

/*-----------------------------------------------------------------------
 * ShowFieldContour
 *
 * INPUTS: none
 *
 * OUTPUTS: none
 */
static int
ShowFieldContour(float *Data, float *lev_data, long lev, int wksid )
{
    int*   iamap;
    int  iarea[GROUP_SIZE];
    int  igrp[GROUP_SIZE];
    int   lfin[MAX_NUM_CONTOURS+1];
    float xwrk[AREA_SIZE];
    float ywrk[AREA_SIZE];
    float *local_data;
    int i,j;
    int ntime_to_display;
    float zmin, zmax;
    char label[80];
    char*  llbs[MAX_NUM_CONTOURS+2];
    float clev1, clev2, cont_interval;
    float* time_labels;
    char time_label[256];
    float* press_levels,press,deltaPressure;
    int *iwrk;
    int int_wksize;
    float *rwrk;
    int real_wksize;
    int hiIdx,lowIdx; /*High index and Low index*/
    /*    float valLow = 0;
	  float valHigh = 0; */
    float valLow;
    float valHigh;
    Gcolr_rep  colr;
    int ncl;                    /* number of contour levels chosen by conpack */
    float clev;                 /* "nice" contour levels chosen by conpack */
    double colidx;
    int lo_idx,hi_idx;
    double delta; 
    int mapsize;
    int nstd;
    double color_interval;
    valLow=0;    
    valHigh=0;    
    nstd = ( sizeof( std_rgb )/ sizeof( struct rgb ) );

    ntime_to_display = stop_index - start_index + 1;
    local_data = (float *) malloc(sizeof(float) * lev * ntime_to_display);
      /*
       * interpolate the data from model pressure levels
       * onto evenly spaced pressure levels
       */
    press_levels = (float *)malloc(sizeof(float) * lev);
    deltaPressure = (lev_data[lev-1] - lev_data[0]) / (lev - 1);
    for(i = lev-1,press = lev_data[lev-1];i >= 0; i--,press -= deltaPressure)
	press_levels[i] = press;
    
    for (i = 0; i < lev; i++) {
        for (j = 0; j < ntime_to_display; j++) {
            if (i == lev -1 || i == 0) {
                *(local_data + i*ntime_to_display + j) 
                    = (float)*(Data + i*ntime + j + start_index);
            }
            else {
                
                if( !findBoundingArrayIndices(press_levels[i],lev_data,lev,&lowIdx,&hiIdx)) {
                    cerr << "ERROR: "__FILE__":" <<__LINE__
                         << ": ShowFieldContour() - Couldn't find bounds of data array\n";
                    exit( -1 );
                }
                valLow = (float)*(Data + lowIdx*ntime + j + start_index);
                valHigh = (float)*(Data + hiIdx*ntime + j + start_index);
                *(local_data + i*ntime_to_display + j) 
                    = Interpolate(lev_data[lowIdx],valLow
                                  ,lev_data[hiIdx],valHigh,press_levels[i]);
            }
        }
    }
    real_wksize = int_wksize = ntime_to_display*lev;
    if(real_wksize < 2500)
        real_wksize = 2500;
    
    if(int_wksize < 1000)
        int_wksize = 1000;
    
    
    iwrk = (int *)malloc(int_wksize*sizeof(int));
    rwrk = (float *)malloc(real_wksize*sizeof(float));
    
    if(iwrk == NULL || rwrk == NULL) {
        cerr << "ERROR: "__FILE__":" <<__LINE__
             << "ShowFieldContour() - Can't malloc memory for iwrk or rwrk.\n";
	exit(-1);
    }
    
    /*
     * Set some attributes of the contour plot
     */
    c_cpseti("NOF - NUMERIC OMISSION FLAGS", 0);   /* !trim trailing zeroes  */
    /*
     * ask for PREF_NUM_CONTOURS contours levels
     * we could actually get as many as PREF_NUM_CONTOURS * 2
     */
    c_cpseti("CLS - CONTOUR LEVEL SELECTOR", +PREF_NUM_CONTOURS); 
    /*
     * Set our own borders so that the plot isn't gridded (we use all the
     * window)
     */
    c_cpseti("SET - DO-SET-CALL FLAG", 0);
    /* set viewport */
    time_labels = (float *)malloc(sizeof(float) * ntime);
    switch( time_format ) {
    case TimeConverter::STEPS:
        for ( i = 0; i < ntime; i++ ) 
            time_labels[i] = (int)( hours[i] * steps_per_hour);
        c_labmod("(I5)","(I5)",0,0,15,15,0,0,1);
        strcpy(time_label,"Steps");
        break;
    case TimeConverter::HOURS:
        for ( i = 0; i < ntime; i++ ) { 
            time_labels[i] = hours[i];
        }
        c_labmod("(F6.1)","(I5)",0,0,15,15,0,0,1);
        strcpy(time_label,"Hours");
        break;
    case TimeConverter::DAYS:
    case TimeConverter::DATE:
          // NCAR graphics can't handle dates properly, so just display days
        for ( i = 0; i < ntime; i++ ) 
            time_labels[i] = hours[i] / 24.0;
        c_labmod("(F6.1)","(I5)",0,0,15,15,0,0,1);
        strcpy(time_label,"Days");
        break;
    default:
        cerr << "ERROR: "__FILE__":" <<__LINE__
             << ": ShowFieldContour() - bad time format\n";
        exit( -1 );
    }
    c_set(0.10, 0.83, 0.15, 0.9,// left right bottom top
          time_labels[start_index], time_labels[stop_index], 1000., 0., 1);
      /*
       * Set the ranges of the x (time) and y axis (pressure 1000-0 mb) 
       */
    c_cpsetr("YC1 - Y COORDINATE AT INDEX 1", 0.); 
    c_cpsetr("YCN - Y COORDINATE AT INDEX N", 1000.);
    c_cpsetr("XC1 - X COORDINATE AT INDEX 1",  time_labels[start_index]);
    c_cpsetr("XCM - X COORDINATE AT INDEX M", time_labels[stop_index]);
    
    c_cpsetr("SPV - SPECIAL VALUE", NC_FILL_FLOAT );
      /*
       * Initialize the area map. Guess a map size - needs to be large enough
       *  to display the necessary contours (requires 10 bytes for each vertex),
       *  but there is no good way to predict the number ahead of time. 
       */
    mapsize = (stop_index - start_index) * 10000;
    iamap = (int*)malloc(mapsize * sizeof(int));
    if ( iamap == NULL ) {
        ShowMsg( __FILE__, __LINE__, "ERROR: Memory allocation for contouring failed! \n Try plotting a shorter time interval" );
        return -1;
    }
    c_arinam(iamap, mapsize);
      /*
       * initialize ConPack
       */
    c_cprect(local_data, ntime_to_display, ntime_to_display, lev,
             rwrk, real_wksize, iwrk, int_wksize);
      /*
       * Force conpack to choose "nice" contour levels so we can figure out
       *   a color scheme based on the number chosen
       */
    c_cppkcl( local_data, rwrk, iwrk );
      /*
       * find out how many contours were used and the contour interval
       * so we can set up the colors properly
       */
    c_cpgeti("NCL - NUMBER OF CONTOUR LEVELS", &ncl );
    if ( ncl == 0 ) {
        stringstream msg;
	//        msg << "Field " << field_name << " has constant value " << valLow << " - unable to plot" << ends;
        msg << "Field " << field_name << " has constant value " << valLow << " - unable to plot";
        if (wksid == OUTPUT_FILE_ID)
            cerr << msg.str() << endl;
        else
            ShowMsg( __FILE__, __LINE__,msg.str().c_str());
        return -1;
    }
    
      /*   generate a color table
       *   make the color table use some reasonable standard values
       *   and linearly interpolate between them to complete it
       *   we'll use white, blue, green, yellow, red, black as the std.
       */
    
      /*
       * set the first two indexed colors to black and white 
       */
    colr.rgb.red = 0;
    colr.rgb.green = 0;
    colr.rgb.blue = 0;
    gset_colr_rep( wksid, 0, &colr );
    colr.rgb.red = 1;
    colr.rgb.green = 1;
    colr.rgb.blue = 1;
    gset_colr_rep( wksid, 1, &colr );
    
    colidx = 0;
    color_interval  = (double)(nstd-1)/ (double)ncl; 
    
    for ( i=0; i<=ncl; i++ ) {
        colidx = i*color_interval;
        
        lo_idx = (int) colidx;

        if (i == ncl)
           hi_idx = 0;
        else 
           hi_idx = 1;

        delta = colidx - lo_idx;
        colr.rgb.red = std_rgb[lo_idx].red*(1-delta) +
            std_rgb[lo_idx+hi_idx].red*(delta);
        colr.rgb.green = std_rgb[lo_idx].green*(1-delta) +
            std_rgb[lo_idx+hi_idx].green*(delta);
        colr.rgb.blue = std_rgb[lo_idx].blue*(1-delta) +
            std_rgb[lo_idx+hi_idx].blue*(delta);
        gset_colr_rep( wksid, i+2, &colr );
    }
    
      /*
       * add contour lines to the area map 
       */
    c_cpclam(local_data, rwrk, iwrk, iamap);
    
      /*
       * Color the map. 
       */
    c_arscam(iamap, xwrk, ywrk, AREA_SIZE, iarea, igrp, GROUP_SIZE, colram);
    
    free( iamap );
      /*
       * Put the border and tick marks on the plot
       */
    if (ntime_to_display - 1 <= MAX_LABELS)
        c_gridal(ntime_to_display - 1,0,10,0,1,1,5,0.,0.);
    else 
        c_gridal(MAX_LABELS,0,10,0,1,1,5,0.,0.);
    
      /*
       * Draw a label bar for the plot, relating colors to values.
       */
    c_cpgetr("ZMN", &zmin);
    c_cpgetr("ZMX", &zmax);
    c_cpseti("PAI - PARAMETER ARRAY INDEX", 1 );
    c_cpgetr("CLV - CONTOUR LEVEL VALUES", &clev1 );
    c_cpseti("PAI - PARAMETER ARRAY INDEX", 2 );
    c_cpgetr("CLV - CONTOUR LEVEL VALUES", &clev2 );
      /*
       * print the actual values of the contour levels into the label bar
       */
    cont_interval = clev2 - clev1;
    llbs[0] = (char *) malloc(sizeof(char) * 10);
    sprintf(llbs[0], "%.4g", clev1-cont_interval );
    for (i = 0; i < ncl; i++) {
        lfin[i] = i+2;           /* set color index array */
        c_cpseti("PAI - PARAMETER ARRAY INDEX", i+1 ); 
        c_cpgetr("CLV - CONTOUR LEVEL VALUES", &clev );
        llbs[i+1] = (char *) malloc(sizeof(char) * 10);
        sprintf(llbs[i+1], "%.4g", clev ); 
    }
    llbs[ncl+1] = (char *) malloc(sizeof(char) * 10);
    sprintf(llbs[ncl+1], "%.4g", clev+cont_interval );
    lfin[ncl] = ncl+2;
        
    c_lbseti("CBL - COLOR OF BOX LINES", 1);
    c_lblbar(1,                 // vertical
             0.85,              // left edge
             1.0,               // right edge
             0.15,              // bottom edge
             0.9,               // top edge
             ncl+1,             // number of boxes
             .25,               // horizontal fill fraction
             1.0,               // vertical fill fraction
             lfin,              // fill color indices
             0,                 // routine used to fill 
             llbs,              // labels
             ncl+2,             // # of labels
             1                  // put labels to the right
        );

      /* LABELS */
    {

        float vl, vr, vb, vt;
        float wl, wr, wb, wt;
        int lf;
        c_getset(&vl, &vr, &vb, &vt, &wl, &wr, &wb, &wt, &lf); /* save original values */
        c_set(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1); /* set viewport */

          /* Bottom label - time */
        c_pchiqu( 0.5, 0.06, time_label, LBL_FONTSIZE, 0.0, 0.0);
          /*  Left Label (pressure) */
        c_pchiqu( 0.02 , 0.5, "Pressure (mb)", LBL_FONTSIZE, 90.0, 0.0);
          /* Title Label */
        sprintf(label, "%s (%s)",field_long_name, field_units); 
        c_pchiqu(0.5, 0.95, label, LBL_FONTSIZE * 1.3, 0.0, 0.0); 
          /* Filename Label  */
        c_pchiqu(0.5, 0.01, input_file, LBL_FONTSIZE, 0.0, 0.0); 

        c_set(vl, vr, vb, vt,wl, wr, wb, wt, lf); /* reset to the original values */ 
    }

    free((char *)local_data);
    for ( i=0; i<ncl+2; i++ )
        free( llbs[i] );
    free((char *)time_labels);
    free((char *)press_levels);
    return 0;
}

static void
rotateField(float *in,
            float *out,
            int x,
            int y)
{
    int i,j;

    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            *(out + x*j + i) = *(in + y*i + j);
        }
    }
}

static int
OpenInputFile()
{
  if ((nc_open(input_file, NC_NOWRITE, &ncid )) != NC_NOERR) {
      cerr << "ERROR: "__FILE__", line " << __LINE__ 
           << " - OpenInputFile():Can not open input file " << input_file << endl;
      return(false);
  }
  return(true);
}

static int
CloseInputFile()
{
  nc_close(ncid);
  return(true);
}

/*-----------------------------------------------------------------------
 * CreatePlot
 *      Performs the actual Plot, based on the current values of the form.
 *
 * INPUTS:
 *
 *      w               Unused
 *      client_data     Unused
 *      call_data       Unused
 *
 * OUTPUTS: none
 */
int
CreatePlot( const char* _input_file, const char* _field_name,
            const char* _field_long_name, const char* _field_units,
            real_t _field_multiplier,
            int start_step, int stop_step, int _time_format,
            real_t _steps_per_hour, int _basedate, bool do_average,
            int deviceID )
{
    int     plev_dimID, time_dimID, ilev_dimID;
    int     nlev;
    int     ilev;
    int     varID;
    float   *field_data;
    int     levID;
    int     timeID;
    float   *lev_data;
    float   *tmp_data;
    int     ndims;
    int     dimids[MAX_VAR_DIMS];
    int     i;
    int     status = 0;

    time_format = _time_format;
    basedate = _basedate;
    
    strcpy( input_file, _input_file );
    strcpy( field_name, _field_name );
    strcpy( field_long_name, _field_long_name );
    strcpy( field_units, _field_units );

    
    field_multiplier = _field_multiplier;
    OpenInputFile();
     
    /*
     * Get dimension IDs to use in accessing data
     */
    nc_inq_dimid( ncid, "lev", &plev_dimID );
    nc_inq_dimid( ncid, "ilev", &ilev_dimID );
    nc_inq_dimid( ncid, "time", &time_dimID );

    nc_inq_dimlen( ncid, plev_dimID, (size_t*)&nlev );
    nc_inq_dimlen( ncid, ilev_dimID, (size_t*)&ilev );
    nc_inq_dimlen( ncid, time_dimID, (size_t*)&ntime );

    nc_inq_varid( ncid, field_name, &varID );
    nc_inq_varndims( ncid, varID, &ndims );

    nc_inq_vardimid( ncid, varID, dimids );

    if ( ndims == 4 && dimids[1] == ilev_dimID ) {
      nlev=ilev;
      nc_inq_varid( ncid, "ilev", &levID );
    }else
      nc_inq_varid( ncid, "lev", &levID );
      
    lev_data = (float *) malloc( sizeof(float) * nlev );
    nc_get_var_float( ncid, levID, lev_data );

    hours = (float*) malloc( sizeof( float ) * ntime );
    nc_inq_varid( ncid, "time", &timeID );
    nc_get_var_float( ncid, timeID, hours );

    steps_per_hour = _steps_per_hour;

    stop_index = (int) (((stop_step/steps_per_hour) - hours[0] )/(hours[1] - hours[0])) ;
    start_index = (int) (((start_step/steps_per_hour) - hours[0] ) / (hours[1] - hours[0])) ;

    if ( ndims == 4 ) {
        field_data = (float *) malloc(sizeof(float) * nlev * ntime);

        if ( dimids[0] == time_dimID ) {
	    tmp_data = (float *) malloc(sizeof(float) * nlev * ntime);
	    nc_get_var_float( ncid, varID, tmp_data );
	    rotateField( tmp_data, field_data, ntime, nlev );
	    free( (char *) tmp_data );
        }
        else
            nc_get_var_float(ncid, varID, field_data);
        
	CloseInputFile();
        /*
         * apply the field_multiplier
         */
        for ( i=0; i<nlev*ntime; i++ )
            field_data[i] *= field_multiplier;

        if (do_average) {
            float *line_data;
            int i,j;

            line_data = (float *)malloc(sizeof(float) * nlev);
            for (i = 0; i < nlev; i++) {
                *(line_data + i) = 0.0;
                for (j = start_index; j <= stop_index; j++) {
                    *(line_data + i) += *(field_data + j + i*ntime);
                }
                *(line_data + i) /= ( stop_index - start_index + 1 );
            }
            
            ShowFieldAverage(line_data, lev_data, nlev);
            free(line_data);
        } else {
            status = ShowFieldContour( field_data, lev_data, nlev, deviceID );
        }
        free((char*) field_data);
        free((char*) lev_data);
        free((char*) hours );
    }
    else if (ndims == 3) {
        field_data = (float *) malloc(sizeof(float) * ntime);
        nc_get_var_float(ncid, varID, field_data);
        /*
         * apply the field_multiplier
         */
        for ( i=0; i<ntime; i++ )
            field_data[i] *= field_multiplier;
        
	CloseInputFile();
        ShowFieldLinePlot(field_data);
        free((char*) field_data);
    }
    else { 
        ShowMsg( __FILE__, __LINE__,"ERROR: field %s has wrong number of dimensions",  field_name );
        return 0;
    }

    if ( deviceID == OUTPUT_FILE_ID ) /* output is to a file */
        c_ngpict(deviceID, 1); /* Action '1' == Update and clear (send END signal)*/
    else
        c_ngpict(deviceID, 0); /* Action '0' == Update  */
    return status;
}

} // extern "C"

#endif  // NO_NCARG
