/*------------------------------------------------------------------------*
 * FILE: history.cpp 
 * $Author: jet $
 * $Id: history.cpp,v 1.1.6.1 2004/05/18 20:51:54 jet Exp $ *
 *------------------------------------------------------------------------*/
#include "field.h"
#include "ncfile.h"
#include "manager.h"

#include "history.h"


using namespace ncfile;

#define STR_DIMLEN 256             // length of string dimension

History::History( const string& name, const string& case_name, OpenMode mode, bool clobber ) :
    NcFile( name, mode, clobber )
{
    if ( mode == CREATE ) {
        NcDimension timeDim( "time", NC_UNLIMITED, true );
        NcDimension levDim( "lev", (size_t)MANAGER.NumLevs() );
        NcDimension ilevDim( "ilev", (size_t)MANAGER.NumLevs()+1 );
        NcDimension latDim( "lat", 1 );
        NcDimension lonDim( "lon", 1 );
        NcDimension switchDim( "switch", (size_t)NUM_SWITCHES );
        NcDimension stringDim( "string", (size_t)STR_DIMLEN );
        
        addDimension( timeDim );
        addDimension( levDim );
        addDimension( ilevDim );
        addDimension( latDim );
        addDimension( lonDim );
        addDimension( switchDim );
        addDimension( stringDim );
        
        addAttribute( NcCharAtt( "Title", "SCAM History File" ) );
        addAttribute( NcIntAtt( "time_step_length", MANAGER.StepLen() ) );
        addAttribute( NcCharAtt( "case", case_name.c_str() ));
        
          // -----------------------------------------------------
          // add the standard variables
          // -----------------------------------------------------
        NcDimension dims[4];
      
          // latitude variable
        dims[0] = latDim;
        latVar = NcVariable<real_t>( "lat", dims, 1 );
        latVar.addAttribute( NcCharAtt("long_name", "latitude") );
        latVar.addAttribute( NcCharAtt( "units", "degrees north" ));
        addVariable( latVar );

          // longitude variable
        dims[0] = lonDim;
        lonVar = NcVariable<real_t>( "lon", dims, 1 );
        lonVar.addAttribute( NcCharAtt( "long_name","longitude" ) );
        lonVar.addAttribute( NcCharAtt( "units", "degrees east" ) );
        addVariable( lonVar );

          // level variable
        dims[0] = levDim;
        levVar = NcVariable<real_t>( "lev", dims, 1 );
        levVar.addAttribute( NcCharAtt( "long_name","baseline vertical pressure coordinate profile" ));
        levVar.addAttribute( NcCharAtt( "units", "hPa" ) ); 
        addVariable( levVar );

          // interface level variable
        dims[0] = ilevDim;
        ilevVar = NcVariable<real_t>( "ilev", dims, 1 );
        ilevVar.addAttribute( NcCharAtt( "long_name","baseline vertical interface pressure coordinate profile" ));
        ilevVar.addAttribute( NcCharAtt( "units", "hPa" ) ); 
        addVariable( ilevVar );

          // time variable
        dims[0] = timeDim;
        timeVar = NcDoubleVar( "time", dims, 1 );
        timeVar.addAttribute( NcCharAtt( "long_name", "elapsed model time" ) );
        timeVar.addAttribute( NcCharAtt( "units", "hours" ) );
        addVariable( timeVar );

          // tsec variable
        dims[0] = timeDim;
        tsecVar = NcIntVar( "tsec", dims, 1 );
        tsecVar.addAttribute( NcCharAtt( "long_name", "seconds since basedate" ) );
        tsecVar.addAttribute( NcCharAtt( "units", "seconds" ) );
        addVariable( tsecVar );

          // switch descriptions variable
        dims[0] = switchDim;
        dims[1] = stringDim;
        swDescVar = NcCharVar( "switch_desc", dims, 2 );
        swDescVar.addAttribute( NcCharAtt("long_name", "logical switches descriptions" ) );
        addVariable( swDescVar );
        
          // switches variable
        dims[0] = switchDim;
        switchVar = NcIntVar( "switches", dims, 1 );
        switchVar.addAttribute( NcCharAtt( "long_name", "model logical switches" ) );
        addVariable( switchVar );
        

          // ipres variable
        dims[0] =  timeDim, dims[1] = ilevDim, dims[2] = latDim, dims[3] = lonDim;
        ipresVar = NcVariable<real_t>( "ipres", dims, 4 );
        ipresVar.addAttribute( NcCharAtt( "long_name","model interface pressure levels" ) );
        ipresVar.addAttribute( NcCharAtt( "units", "hPa" ) );
        addVariable( ipresVar );


          // pres variable
        dims[0] =  timeDim, dims[1] = levDim, dims[2] = latDim, dims[3] = lonDim;
        presVar = NcVariable<real_t>( "pres", dims, 4 );
        presVar.addAttribute( NcCharAtt( "long_name","model pressure levels" ) );
        presVar.addAttribute( NcCharAtt( "units", "hPa" ) );
        addVariable( presVar );


          // basedate variable (no dimensions)
        bdateVar = NcIntVar( "bdate", dims, 0 );
        bdateVar.addAttribute( NcCharAtt( "units", "yymmdd" ) );
        addVariable( bdateVar );

          // 
          // add the field variables 
          //           //
        NcDimension levDims[4] = { timeDim, levDim, latDim, lonDim };
        NcDimension ilevDims[4] = { timeDim, ilevDim, latDim, lonDim };
        NcDimension srfDims[3] = { timeDim, latDim, lonDim };

        for ( FieldListConstItr it = MANAGER.GetFieldList().begin(); 
              it !=  MANAGER.GetFieldList().end(); ++it ) {
            Field* f = (*it).second;
            delete f->_histVar;
            f->_histVar = 0;
            if (!f->IsSaved() ) 
                continue;

              // set up dimensions and attributes for the fields
            if ( f->IsMultiLevel() )
	      if ( f->NumLevs() == MANAGER.NumLevs()+1)
                f->_histVar = new NcVariable<real_t>( f->Name(), ilevDims, 4 );
	      else
                f->_histVar = new NcVariable<real_t>( f->Name(), levDims, 4 );
            else
                f->_histVar = new NcVariable<real_t>( f->Name(), srfDims, 3 );
        
            NcCharAtt longName( string("long_name"), f->LongName().c_str() );
            NcCharAtt units( string("units"), f->Units().c_str() );
            NcCharAtt plotUnits( string("plot_units"), f->PlotUnits().c_str() );
            NcIntAtt numLevs( string("num_levels"), f->NumLevs() );
            NcFloatAtt plotMult( string("plot_multiplier"), f->PlotMult() );
            
            NcAttBase* atts[5] = { &longName, &units, &plotUnits, &numLevs, &plotMult };
            f->_histVar->setAttributes( atts, 5 );
              // add the field variable to the file
            addVariable( *f->_histVar ); 
        }
          // we're done with the metadata specification
        endDefineMode();
        
          //
          // write the initial conditions variables
          //
        lonVar = MANAGER.Lon();
        lonVar.write();

        latVar = MANAGER.Lat();
        latVar.write();

        levVar = MANAGER.BaseLevels();
        levVar.write();

        ilevVar = MANAGER.BaseILevels();
        ilevVar.write();

        bdateVar = MANAGER.BaseDate();
        bdateVar.write();

        for ( size_t i=0; i<switchDim.size(); i++ ) {
            switchVar[i] = MANAGER.SwitchState(i);
            size_t j=0;
            for ( ; j<MANAGER.SwitchDesc(i).size(); j++ )
                swDescVar[i*stringDim.size()+j] = MANAGER.SwitchDesc(i)[j];
            swDescVar[i*stringDim.size()+j] = 0; // null terminate the string
        }
        swDescVar.write();
        switchVar.write();
    }
    
    if ( mode == READ ) {
        timeVar = variable( "time" );
        tsecVar = variable( "tsec" );
        presVar = variable( "pres" );
	ipresVar = variable( "ipres" );
    }
}

//
// get a list of the variables (fields) that are plottable.
//
FieldList
History::PlottableFields()
{
    FieldList fieldList;

      // clear the fieldList
    for ( FieldListItr it = fieldList.begin(); it != fieldList.end(); ++it )
        delete (*it).second;
    fieldList.clear();

    for ( int id=0; id<numVars(); id++ ) {
        NcVariable<real_t> v = variable( id );
        if ( v.hasAttribute( "num_levels" ) ) 
            fieldList.insert( fieldpair( v.name(), new Field( v ) ) );
    }
    return fieldList;
}

void History::WriteTime( double hours, int seconds, int index ) 
{
    timeVar = hours;
    tsecVar = seconds;
    timeVar.writeRecord(index );
    tsecVar.writeRecord(index );
}

void History::WriteLevs( const vector<real_t>& levels, int index )
{
    presVar = levels;
    presVar.writeRecord( index );
}
void History::WriteILevs( const vector<real_t>& levels, int index )
{
    ipresVar = levels;
    ipresVar.writeRecord( index );
}

History::~History() 
{
}
