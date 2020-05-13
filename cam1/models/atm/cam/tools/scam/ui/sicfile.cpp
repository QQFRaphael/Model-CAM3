/*------------------------------------------------------------------------*
 * File: sicfile.cpp 
 * $Author: jet $
 * $Id: sicfile.cpp,v 1.1.6.1 2004/05/18 20:52:02 jet Exp $ *
 *------------------------------------------------------------------------*/
#include "field.h"
#include "dataset.h"
#include "ncfile.h"
#include "manager.h"

#include "sicfile.h"

#define STR_DIMLEN 256             // length of string dimension

using namespace ncfile;

Sicfile::Sicfile( const string& name, OpenMode mode, bool clobber ) :
    NcFile( name, mode, clobber )
{
    NcVariable<real_t> lat, lon, lev;
    NcIntVar bdate, tsec, runtype, sw; 
    NcCharVar swdesc, dataset[NUM_INIT_FILES];
    NcDimension dims[4];
    int rtype;
    if ( mode == CREATE ) {
        NcDimension timeDim( "time", 1 );
        NcDimension levDim( "lev", (size_t)MANAGER.NumLevs() );
        NcDimension latDim( "lat", 1 );
        NcDimension lonDim( "lon", 1 );
        NcDimension switchDim( "switch", (size_t)NUM_SWITCHES );
        NcDimension stringDim( "string", (size_t)STR_DIMLEN );
        
        addDimension( timeDim );
        addDimension( levDim );
        addDimension( latDim );
        addDimension( lonDim );
        addDimension( switchDim );
        addDimension( stringDim );
        
        addAttribute( NcCharAtt( "Title", "SCAM Saved Initial Conditions File" ) );
        addAttribute( NcIntAtt( "time_step_length", MANAGER.StepLen() ) );
        
          //
          // add the standard variables
          //

          // latitude variable
        dims[0] = latDim;
        lat = NcVariable<real_t>( "lat", dims, 1 );
        lat.addAttribute( NcCharAtt( "long_name", "latitude" ) );
        lat.addAttribute( NcCharAtt( "units", "degrees north" ) );
        addVariable( lat );

          // longitude variable
        dims[0] = lonDim;
        lon = NcVariable<real_t>( "lon", dims, 1 );
        lon.addAttribute( NcCharAtt( "long_name", "longitude" ) );
        lon.addAttribute( NcCharAtt( "units", "degrees east" ) );
        addVariable( lon );

          // level variable
        dims[0] = levDim;
        lev = NcVariable<real_t>( "lev", dims, 1 );
        lev.addAttribute( NcCharAtt( "long_name", "baseline vertical pressure coordinate profile" ) );
        lev.addAttribute( NcCharAtt( "units", "hPa" ) );
        addVariable( lev );

          // switches variable
        dims[0] = switchDim;
        sw = NcIntVar( "switches", dims, 1 );
        sw.addAttribute( NcCharAtt( "long_name", "model logical switches" ) );
        addVariable( sw );

          // switch descriptions variable
        dims[0] = switchDim;
        dims[1] = stringDim;
        swdesc = NcCharVar( "switch_desc", dims, 2 );
        swdesc.addAttribute( NcCharAtt( "long_name", "logical switches descriptions" ) );
        addVariable( swdesc );

          // basedate variable
        bdate = NcIntVar( "bdate", dims, 0 );
        bdate.addAttribute( NcCharAtt( "units", "yymmdd" ) );
        addVariable( bdate );

          // basesecs variable
        tsec = NcIntVar( "tsec", dims, 0 );
        tsec.addAttribute( NcCharAtt( "units", "seconds" ) );
        tsec.addAttribute( NcCharAtt( "long_name", "seconds since basedate" ) );
        addVariable( tsec );

          // run type variable
        runtype = NcIntVar( "runtype", dims, 0 );
        runtype.addAttribute( NcCharAtt( "long_name", "scam run type" ) );
        addVariable( runtype );

          // dataset names
        dims[0] = stringDim;
        for ( int i=0; i<NUM_INIT_FILES; i++ ) {
            dataset[i] = NcCharVar( Dataset::TypeName(i), dims, 1 );
            dataset[i].addAttribute( NcCharAtt( "long_name", Dataset::TypeDescription(i).c_str() ));
            addVariable( dataset[i] );
        }

          //
          // add the modifiable fields
          //
        NcDimension levDims[4] = { timeDim, levDim, latDim, lonDim };
        NcDimension srfDims[3] = { timeDim, latDim, lonDim };
        
        vector<NcVariable<real_t>*> fldVars;  // holds the netCDF field variables

        const FieldList& fl =  MANAGER.GetFieldList();

        for ( FieldListConstItr it = fl.begin(); it != fl.end(); ++it ) {
            Field *f = (*it).second;
            if ( ! f->IsModifiable() ) 
                continue;       // only add modifiable fields
              // define the variable
            NcVariable<real_t>* fvar;
            if ( f->IsMultiLevel() ) 
                fvar = new NcVariable<real_t>( f->Name(), levDims, 4 );
            else
                fvar = new NcVariable<real_t>( f->Name(), srfDims, 3 );
              // define the variable attributes 
            NcCharAtt longName( string("long_name"), f->LongName().c_str());
            NcCharAtt units( string("units"), f->Units().c_str() );
            NcIntAtt numLevs( string("num_levels"), f->NumLevs() );
            NcAttBase* atts[3] = { &longName, &units, &numLevs };
            fvar->setAttributes( atts, 3 );
            
              // add the variable to the file
            addVariable( *fvar ); 
              // get the field's data from the model and
              // copy it to the FieldVar for writing below

            MANAGER.ReadField( *f );
            for ( int i=0; i<f->NumLevs(); i++ )
                (*fvar)[i] = (*f)[i];
              // store it for later when we can write it
              // (can't write variables while in define mode)
            fldVars.push_back( fvar );
        }
        
        endDefineMode();
        
          //************************************************
          // write the variables to the file
          //************************************************
        
        lat = MANAGER.Lat();
        lat.write();
        lon = MANAGER.Lon();
        lon.write();
        lev = MANAGER.BaseLevels();
        lev.write();
        bdate = MANAGER.BaseDate();
        bdate.write();
        tsec = MANAGER.BaseSecs();
        tsec.write();
        if (MANAGER.RunType()>SIC) 
	  rtype=MANAGER.RunType()-SIC;
	else
	  rtype=MANAGER.RunType();
        runtype = rtype;

        runtype.write();

          // dataset names
        for ( int i=0; i<NUM_INIT_FILES; i++ ) {
              // initialize its value
            strcpy( &dataset[i][0], MANAGER.GetDataset( i ).name().c_str() );
            dataset[i].write();
        }        
          // model switches        
        for ( int i=0; i<NUM_SWITCHES; i++ ) {
            sw[i] = MANAGER.SwitchState(i);
            strcpy( &swdesc[i*STR_DIMLEN], MANAGER.SwitchDesc(i).c_str() );
        }
        sw.write();    
        swdesc.write();
 
          // write the modifiable fields' data

        for( vector<NcVariable<real_t>*>::iterator it = fldVars.begin(); it!= fldVars.end(); ++it ) {
            (*it)->write();
              // delete the NcVariable<real_t> after it has been written
            delete *it;         
        }
        
                    
    } // if ( mode == CREATE )

    if ( mode == NcFile::READ ) {
        lat = variable("lat");
        lat.read();
        MANAGER.SetLat( lat );
        
        lon = variable("lon");
        lon.read();
        MANAGER.SetLon( lon );
        
        bdate = variable("bdate");
        bdate.read();
        MANAGER.SetBaseDate( bdate );
        
        tsec = variable("tsec");
        tsec.read();
        MANAGER.SetBaseSecs( tsec[0] );
        
        sw = variable("switches");
        sw.read();            
        swdesc = variable("switch_desc");
        swdesc.read();
        int strdimlen = dimension("string").size();
        for ( int i=0; i<NUM_SWITCHES; i++ ) {
            MANAGER.SetSwitch( sw[i], i );
            MANAGER.SetSwitchDesc( &swdesc[i*strdimlen], i );
        }
        
        runtype = variable("runtype");
        runtype.read();
        MANAGER.SetRunType( runtype );
        MANAGER.SetDataset( name, MODEL );
        
        for ( int i=0; i<SIC; i++ ) {
            dataset[i] = variable( Dataset::TypeName( i ) );
            dataset[i].read();
            MANAGER.SetDataset( dataset[i].data(), i );
        }
        MANAGER.SetDataset( name, SIC );
    }
}

Sicfile::~Sicfile()
{
}


