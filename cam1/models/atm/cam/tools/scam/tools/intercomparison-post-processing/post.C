#include <map.h>
#include "field.H"

#define KEY first
#define FLD second

#define FIELD(f) fields.find(f)->FLD
#define APPEND(x,y) fields.insert( fldpair( x, y) )

void
PrintUsage()
{
    cerr << "usage: post <histfile> <simulation> <\"brief description\">" << endl
         << "   histfile: an SCCM history file" << endl
         << "   simulation: the simulation name (used for output)" <<endl
         << "   description: a brief description of the run enclosed by quotes" << endl;
}

//-----------------------------------------------------------
int 
main( int argc, char* argv[] )
{
    char histfile[256];
    char desc[256];
    char sim[256];
    float* time;
    float* lev;
    int ntime;
    int nlev;

    if ( argc != 4 ) {
        PrintUsage();
        exit( 1 );
    }

    strcpy( histfile, argv[1] );
    strcpy( sim, argv[2] );
    strcpy( desc, argv[3] );

    NCFile h( histfile );

    if( ! h.Open() ) {
        cerr << "post: ERROR: can't open \"" << histfile << "\"" << endl;
        exit( 1 );
    }
    
    ntime = h.GetDimSize( "time" );
    nlev = h.GetDimSize( "lev" );
    time = new float[ ntime ];
    lev = new float[ nlev ];

    h.Read( "time", time );
    float timeInterval = time[1] - time[0];
    for ( int t=0;t<ntime;t++ )
        time[t] -= timeInterval/2;         // we want the time at the midpoint (everything is in hours)
    
    h.Read( "lev", lev );


      // use the STL map class to build a field list
      // note that we're not providing a functor to 
      // compare items for sorting purposes, and the
      // default (less()) leaves the fields in the order
      // in which they are added.

    typedef map<const char*,Field>  fldlist;
    typedef map<const char*,Field>::iterator flditr;
    typedef pair<const char*,Field> fldpair;

    fldlist fields;
    flditr it;


      //---------------------------------------------
      // 2-D fields
      //---------------------------------------------
    
    APPEND( "ZMID",  Field( "HEIGHT (km)", "%7.3f", nlev, ntime, 1, .001 ) );
    APPEND( "Tk",     Field( "TEMPERATURE (K)", "%7.2f", nlev, ntime, 2 ) );
    APPEND( "Q",     Field( "WATER VAPOR MIXING RATIO (g/kg)", "%7.3f", nlev, ntime, 3, 1000.0 ) );
    APPEND( "RELHUM",Field( "RELATIVE HUMIDITY", "%6.3f", nlev, ntime, 4 ) );
    APPEND( "CLWMR", Field( "CLOUD WATER MIXING RATIO (g/kg)", "%7.4f", nlev, ntime, 5 ) );
    APPEND( "CIWMR", Field( "CLOUD ICE MIXING RATIO (g/kg)", "%7.4f", nlev, ntime, 6 ) );
    APPEND( "RMR",   Field( "RAIN MIXING RATIO (N/A)", "%7.4f", nlev, ntime, 7 ) );
    APPEND( "SMR",   Field( "SNOW MIXING RATIO (N/A)", "%7.4f", nlev, ntime, 8 ) );
    APPEND( "GMR",   Field( "GRAUPEL MIXING RATIO (N/A)", "%7.4f", nlev, ntime, 9 ) );
    APPEND( "CLOUD", Field( "CLOUD FRACTION", "%6.3f", nlev, ntime, 10 ) );
    APPEND( "U",     Field( "HORIZONTAL WIND VELOCITY IN X-DIRECTION (m/s)", "%7.2f", nlev, ntime, 11 ) );
    APPEND( "V",     Field( "HORIZONTAL WIND VELOCITY IN Y-DIRECTION (m/s)", "%7.2f", nlev, ntime, 12 ) );
    APPEND( "Q1MQR", Field( "APPARENT HEAT SOURCE (K/day)", "%7.2f", nlev, ntime, 13, 86400.0 ) );
    APPEND( "Q2",    Field( "APPARENT MOISTURE SINK (K/day)", "%7.2f", nlev, ntime, 14, 86400.0 ) );
    APPEND( "Q1CC",  Field( "CONVECTIVE Q1C (N/A)", "%7.2f", nlev, ntime, 15 ) );
    APPEND( "Q1CS",  Field( "CONVECTIVE Q1C (N/A)", "%7.2f", nlev, ntime, 16 ) );
    APPEND( "Q2C",   Field( "CONVECTIVE Q2 (N/A)", "%7.2f", nlev, ntime, 17 ) );
    APPEND( "Q2S",   Field( "CONVECTIVE Q2 (N/A)", "%7.2f", nlev, ntime, 18 ) );
    APPEND( "QR",    Field( "RADIATIVE HEATING RATE (K/day)", "%7.2f", nlev, ntime, 19, 86400.0 ) );
    APPEND( "QRS",   Field( "SOLAR RADIATIVE HEATING RATE (K/day)", "%7.2f", nlev, ntime, 20, 86400.0 ) ); 
    APPEND( "QRL",   Field( "INFRARED HEATING RATE (K/day)", "%7.2f", nlev, ntime, 21, 86400.0 ) ); 
    APPEND( "QCLR",  Field( "CLEAR RADIATIVE HEATING RATE (N/A)", "%7.2f", nlev, ntime, 22 ) ); 
    APPEND( "QCLD",  Field( "CLOUDY RADIATIVE HEATING RATE (N/A)", "%7.2f", nlev, ntime, 23 ) ); 
    APPEND( "CMFMC", Field( "CLOUD MASS FLUX (mb/s)", "%8.5f", nlev, ntime, 24, 0.0980616 ) );
    APPEND( "CMFU",  Field( "UPDRAFT CLOUD MASS FLUX (N/A)", "%7.5f", nlev, ntime, 25 ) ); 
    APPEND( "CMFD",  Field( "DOWNDRAFT DRAFT CLOUD MASS FLUX (N/A)", "%7.5f", nlev, ntime, 26 ) ); 
    APPEND( "CUF",   Field( "FRACTIONAL AREA OF UPDRAFT CORES (N/A)", "%6.3f", nlev, ntime, 27 ) ); 
    APPEND( "CDF",   Field( "FRACTIONAL AREA OF DOWNDRAFT CORES (N/A)", "%6.3f", nlev, ntime, 28 ) ); 
    APPEND( "CUSA",  Field( "AVERAGE CORE UPDRAFT SPEED (N/A)", "%7.3f", nlev, ntime, 29 ) ); 
    APPEND( "CDSA",  Field( "AVERAGE CORE DOWNDRAFT SPEED (N/A)", "%7.3f", nlev, ntime, 30 ) ); 


      //=========================
      // 1-D fields
      //=========================
    
      // 
      // GROUP 1
      // 

    APPEND( "TS",   Field( "SURFACE SKIN TEMPERATURE (K)", "%7.2f", 1, ntime, 2, 1.0, 1 ) );
    APPEND( "NSS",  Field( "NEAR-SURFACE DRY STATIC ENERGY (kJ/kg)", "%7.2f", 1, ntime, 3,  0.001, 1 ) );
    APPEND( "NSQ",  Field( "NEAR-SURFACE WATER VAPOR MIXING RATIO (g/kg)", "%6.2f", 1, ntime, 4, 1000.0, 1 ) );
    APPEND( "NSH",  Field( "NEAR-SURFACE MOIST STATIC ENERGY (kJ/kg)", "%7.2f", 1, ntime, 5, 0.001, 1 ) );
    APPEND( "NSU",  Field( "NEAR-SURFACE HORIZONTAL WIND VELOCITY IN X-DIRECTION (m/s)", "%7.2f", 1, ntime, 6, 1.0, 1 ) );
    APPEND( "NSV",  Field( "NEAR-SURFACE HORIZONTAL WIND VELOCITY IN Y-DIRECTION (m/s)", "%7.2f", 1, ntime, 7, 1.0, 1 ) );
    APPEND( "SHFLX",Field( "SURFACE TURBULENT FLUX OF SENSIBLE HEAT (W/m^2)", "%6.1f", 1, ntime, 8, 1.0, 1 ) );
    APPEND( "LHFLX",Field( "SURFACE TURBULENT FLUX OF LATENT HEAT (W/m^2)", "%6.1f", 1, ntime, 9, 1.0, 1 ) );
    APPEND( "TAUX", Field( "SURFACE TURBULENT FLUX OF HORIZONTAL MOMENTUM COMPONENT IN X-DIRECTION (N/m^2)", "%8.4f", 1, ntime, 10, 1.0, 1 ) );
    APPEND( "TAUY", Field( "SURFACE TURBULENT FLUX OF HORIZONTAL MOMENTUM COMPONENT IN Y-DIRECTION (N/m^2)", "%8.4f", 1, ntime, 11, 1.0, 1 ) );
    
      // 
      // group 2
      // 
    
    APPEND( "FSDS",  Field( "SURFACE DOWNWELLING SOLAR RADIATIVE FLUX (W/m^2)", "%7.1f", 1, ntime, 2, 1.0, 2 ) );
    APPEND( "FSUS",   Field( "SURFACE UPWELLING SOLAR RADIATIVE FLUX (W/m^2)", "%7.1f", 1, ntime, 3, 1.0, 2 ) );
    APPEND( "FLDS",  Field( "SURFACE DOWNWELLING INFRARED RADIATIVE FLUX (W/m^2)", "%6.1f", 1, ntime, 4, 1.0, 2 ) );
    APPEND( "FLUS",   Field( "SURFACE UPWELLING INFRARED RADIATIVE FLUX (W/m^2)", "%6.1f", 1, ntime, 5, 1.0, 2 ) );
    APPEND( "SOLIN", Field( "TOA DOWNWELLING SOLAR RADIATIVE FLUX (W/m^2)", "%7.1f", 1, ntime, 6, 1.0, 2 ) );
    APPEND( "FSUT", Field( "TOA UPWELLING SOLAR RADIATIVE FLUX (W/m^2)", "%6.1f", 1, ntime, 7, 1.0, 2 ) );
    APPEND( "FLNT",  Field( "TOA UPWELLING INFRARED RADIATIVE FLUX (W/m^2)", "%6.1f", 1, ntime, 8, 1.0, 2 ) );
    APPEND( "CLDTOT",Field( "CLOUD AMOUNT", "%6.3f", 1, ntime, 9, 1.0, 2 ) );
    APPEND( "ACLD",  Field( "COLD CLOUD TOP AREA (N/A)", "%6.3f", 1, ntime, 10, 1.0, 2 ) ); 
    APPEND( "PW",    Field( "PRECIPITABLE WATER (kg/m^2)", "%6.2f", 1, ntime, 11, 1.0, 2 ) );
    APPEND( "TOTLWP",Field( "CLOUD LIQUID WATER PATH (kg/m^2)", "%10.3e", 1, ntime, 12, 0.001, 2 ) );
    APPEND( "TOTIWP",Field( "CLOUD ICE PATH (kg/m^2)", "%10.3e", 1, ntime, 13, 0.001, 2 ) );
    
      // 
      // group3
      // 
    
    APPEND( "RP",   Field( "VERTICALLY INTEGRATED RAIN (N/A)", "%10.3e", 1, ntime, 2, 1.0, 3 ) );
    APPEND( "SP",   Field( "VERTICALLY INTEGRATED SNOW (N/A)", "%10.3e", 1, ntime, 3, 1.0, 3 ) );
    APPEND( "GP",   Field( "VERTICALLY INTEGRATED GRAUPEL (N/A)", "%10.3e", 1, ntime, 4, 1.0, 3 ) );
    APPEND( "PREC", Field( "SURFACE RAINFALL RATE (mm/day)", "%7.2f", 1, ntime, 5, 8.64e+07, 3 ) );
    APPEND( "PRECC",Field( "CONVECTIVE SURFACE RAINFALL RATE (mm/day)", "%7.2f", 1, ntime, 6, 8.64e+07, 3 ) );
    APPEND( "PRECL",Field( "STRATIFORM SURFACE RAINFALL RATE (mm/day)", "%7.2f", 1, ntime, 7, 8.64e+07, 3 ) );
    APPEND( "AR",   Field( "RAIN FRACTIONAL AREA (N/A)", "%6.3f", 1, ntime, 8, 1.0, 3 ) );
    APPEND( "AC",   Field( "CONVECTIVE FRACTIONAL AREA (N/A)", "%6.3f", 1, ntime, 9, 1.0, 3 ) );
    APPEND( "AS",   Field( "STRATIFORM FRACTIONAL AREA (N/A)", "%6.3f", 1, ntime, 10, 1.0, 3 ) );
    
    
      //===========================================
      // read the fields in from the netcdf file
      //  (fill in with missing values if n/a)
      //===========================================
    
    for ( it = fields.begin(); it != fields.end(); ++it ) {
        if ( ! h.Read( it->KEY, it->FLD ) ) {
            cout << "\nWarning! Missing field: "<< it->KEY << "\n" << endl;
            it->FLD.SetDataIsMissing( true );
        }
    }

    
      //========================================
      // initialize the derived fields
      //========================================
    
    const float Cp = 1.00464e3; // specific heat cap. of dry air (J/kg/K)
    const float g = 9.80616;    // acceleration due to gravity (m/s^2)
    const float L = 2.5104e6;   // latent heat of vaporization
    
    for ( int t=0; t<ntime; t++ )
        FIELD("NSS")[0][t] = ( Cp * FIELD("T")[nlev-1][t] ) + ( g * FIELD("ZMID")[nlev-1][t] ); 
    FIELD("NSS").SetDataIsMissing( false );

    for ( int t=0; t<ntime; t++ )
        FIELD("NSQ")[0][t] = FIELD("Q")[nlev-1][t];
    FIELD("NSQ").SetDataIsMissing( false );

    for ( int t=0; t<ntime; t++ )
        FIELD("NSH")[0][t] = FIELD("NSS")[0][t] + ( L * FIELD("NSQ")[0][t] );
    FIELD("NSH").SetDataIsMissing( false );
    
    for ( int t=0; t<ntime; t++ )
        FIELD("NSU")[0][t] = FIELD("U")[nlev-1][t];
    FIELD("NSU").SetDataIsMissing( false );

    for ( int t=0; t<ntime; t++ )
        FIELD("NSV")[0][t] = FIELD("V")[nlev-1][t];
    FIELD("NSV").SetDataIsMissing( false );

    Field flns( 1, ntime );
    h.Read( "FLNS", flns ); 
    for ( int t=0; t<ntime; t++ )
        FIELD("FLDS")[0][t] = FIELD("FLUS")[0][t] - flns[0][t];
    FIELD("FLDS").SetDataIsMissing( false );

    Field fsnt( 1, ntime );
    h.Read( "FSNT", fsnt ); 
    for ( int t=0; t<ntime; t++ )
        FIELD("FSUT")[0][t] = FIELD("SOLIN")[0][t] - fsnt[0][t];
    FIELD("FSUT").SetDataIsMissing( false );
    
    Field fsns( 1, ntime );
    h.Read( "FSNS", fsns );
    for ( int t=0; t<ntime; t++ )
        FIELD("FSUS")[0][t] = FIELD("FSDS")[0][t] - fsns[0][t];
    FIELD("FSUS").SetDataIsMissing( false );

      //============================
      // write out the files
      //============================

    char filename[256];
    char fstring[256];

      //----------------------------
      // 2D fields
      //----------------------------

    for ( it = fields.begin(); it != fields.end(); ++it ) {
        if ( it->FLD.GetNlev() < 2 ) {
            continue;
        }
          // if the data is missing, just create an empty file
          //  with "not_submitted" appended to the name
        if ( it->FLD.DataIsMissing() ) {
            cout << "Creating empty file for "<< it->KEY << endl;
            sprintf( filename, "%s.p%d.hack.not_submitted", sim, it->FLD.GetId() ); 
            ofstream f( filename );
            continue;
        }

        sprintf( filename, "%s.p%d.hack", sim, it->FLD.GetId() ); 
        ofstream f( filename );
        f << "# " << filename << "  " << desc << "  " << it->FLD.GetDescrip() << endl;;
        
        f << " " << nlev << " " << ntime << endl;
        f << " .000000  .000000" << endl;

          // print out the levels
        for ( int i=0; i<nlev; i++ ) {
            sprintf( fstring, " %7.3f", lev[i] );
            f << fstring;
            if ( !((i+1)%10) && i!=nlev-1 )
                f << endl;
        }
        f << endl;

          // print out the time at the midpoint
        for ( int i=0; i<ntime; i++ ) {
            sprintf( fstring, " %6.1f", time[i] );
            f << fstring;
            if ( !((i+1)%10) && i!=ntime-1 )
                f << endl;
        }
        f << endl << endl;
        
        it->FLD.WriteData( f );
        f << endl;
    }

      //---------------------------------------------
      // 1D fields 
      //---------------------------------------------
    
      //
      // group 1
      //
    sprintf( filename, "%s.t1.hack", sim ); 
    
    ofstream g1( filename );
    g1 << "# "<< filename  << " " << desc << endl;
    for ( int t=0; t<ntime; t++ ) {
        sprintf( fstring, "%6.1f", time[t] ); // time at midpoint
        g1 << fstring;
        for ( it = fields.begin(); it != fields.end(); ++it ) {
            if ( it->FLD.GetGroup() == 1 ) 
                it->FLD.WriteData( g1, t );
        }
        g1 << endl;
    }
    g1 << endl;
    
      //
      // group 2
      //
    sprintf( filename, "%s.t2.hack", sim ); 
    
    ofstream g2( filename );
    g2 << "# " << filename  << " " << desc << endl;

    for ( int t=0; t<ntime; t++ ) {
        sprintf( fstring, "%6.1f", time[t] ); // time at midpoint
        g2 << fstring;
        for ( it = fields.begin(); it != fields.end(); ++it ) {
            if ( it->FLD.GetGroup() == 2 ) 
                it->FLD.WriteData( g2, t );
        }
        g2 << endl;
    }
    g2 << endl;
    
      //
      // group 3
      //
    sprintf( filename, "%s.t3.hack", sim ); 
    
    ofstream g3( filename );
    g3 << "# " << filename  << " " << desc << endl;

    for ( int t=0; t<ntime; t++ ) {
        sprintf( fstring, "%6.1f", time[t] ); // midpoint
        g3 << fstring ;
        for ( it = fields.begin(); it != fields.end(); ++it ) {
            if ( it->FLD.GetGroup() == 3 ) 
                it->FLD.WriteData( g3, t );
        }
        g3 << endl;
    }
    g3 << endl;
    
    return 0;
}







