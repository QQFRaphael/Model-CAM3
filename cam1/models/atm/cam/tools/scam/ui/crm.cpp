
#include <iostream.h>
#include <qstring.h>
#include <string.h>
#include <fstream.h>
#include <stdio.h>
// #include <stdlib.h>
// #include <ctype.h>

#include "ncfile.h"
#include "timeconvert.h"
#include "crm.h"
#include "defaults.h"

using namespace ncfile;

// initialize static member 
crm* crm::_instance = 0;

crm::crm() {
}

crm::~crm() {
}

string
crm::Setup( string& startupfile, const string& infilename )
{
  struct crm_input_data inputfields;
  string filename;

  Defaults d( startupfile, Defaults::READ );
  string userdatadir = d.GetStringDefault( "userdatadir" );

  cout << "Reading CRM input from file " << infilename << endl;
  Read_Input_Data( infilename, inputfields);
  
  Write_Pressure_NC_File(userdatadir + PRESSFILE,inputfields);
  
  Write_Ozone_NC_File(userdatadir + OZONEFILE,inputfields);
  
  Write_SCAM_NC_File(userdatadir + SCAMFILE,inputfields);
  
  Write_SCAM_Quickstart_File(userdatadir + SCAMQSFILE,inputfields);

  return  userdatadir + SCAMQSFILE ;
}

void
crm::PostProcess( const string& ncoutfile)
{
  struct crm_output_data output;

  try {
    Read_Output_Data( ncoutfile, output);
  }
  catch ( NcErr ) {
    cerr << "Couldn't read CRM output file -- text output file not made." << endl;
  }
  Print_Output_Data( CRMTextFile, output);
}

void  
crm::Read_Input_Data( const string& infilename, struct crm_input_data& input) 
{
  char buffer[1024] = {0};
  int lcount = -1;
  real_t tcount,julday; 

  ifstream infile(infilename.data());
  if (!infile) {
    cerr << "Error: cannot open file " << infilename.data() << endl;
    exit(1); 
  }
  
  infile.getline(buffer,1024); // Comment line
  infile.getline(buffer,1024); // Comment line
  infile.getline(buffer,1024); // Comment line
  infile >> julday; infile.getline(buffer,1024);
  infile >> input.lat; infile.getline(buffer,1024);
  infile.getline(buffer,1024); // Comment line
  infile >> tcount;
  while ((int)(tcount - 1) == ++lcount) {
    infile >> input.pmid[lcount];
    input.pmid[lcount] = input.pmid[lcount] * 100.;
    infile >> input.T[lcount];
    infile >> input.q[lcount];
    infile >> input.o3mmr[lcount];
    infile >> input.cld[lcount];
    infile >> input.clwp[lcount]; 
    infile.getline(buffer,1024);
    infile >> tcount;
  }
  cout << "Found " << lcount << " levels." << endl;
  input.plev = lcount;
  input.Ps = tcount * 100.;
  infile.getline(buffer,1024);
  
  infile >> input.Tsair	;    infile.getline(buffer,1024);
  infile >> input.Tg	;    infile.getline(buffer,1024);
  infile >> input.oro	;    infile.getline(buffer,1024);
  infile.getline(buffer,1024); // Obsolete input line
  infile >> input.snowh	;    infile.getline(buffer,1024);
  infile >> input.asdir	;    infile.getline(buffer,1024);
  infile >> input.asdif	;    infile.getline(buffer,1024);
  infile >> input.aldir	;    infile.getline(buffer,1024);
  infile >> input.aldif	;    infile.getline(buffer,1024);
  infile.getline(buffer,1024); // Obsolete input line
  infile >> input.co2vmr  ;    infile.getline(buffer,1024);
  infile >> input.n2ovmr  ;    infile.getline(buffer,1024);
  infile >> input.ch4vmr  ;    infile.getline(buffer,1024);
  infile >> input.f11vmr  ;    infile.getline(buffer,1024);
  infile >> input.f12vmr  ;    infile.getline(buffer,1024);
  infile >> input.tauvis  ;    infile.getline(buffer,1024);
  infile >> input.scon    ;    infile.getline(buffer,1024);
  infile >> input.iyear_AD;    infile.getline(buffer,1024);
  infile >> input.lon     ;    infile.getline(buffer,1024);

  // Fill in other values
  input.phis = 0.;
  TimeConverter::SecondsToDate( (int)(julday * 86400), input.iyear_AD*10000+100 ,
				input.bdate, input.tsec );
  for (lcount = 0; lcount < input.plev ; ++lcount ) {
    input.hyam[lcount] = 0.0;
    input.hybm[lcount] = input.pmid[lcount] / input.Ps;
  }
  input.pint[lcount] = input.Ps;
  input.hybi[lcount] = 1.0;
  input.hyai[lcount] = 0.;
  for (--lcount ; lcount > 0 ; --lcount ) {
    input.hyai[lcount] = 0.0;
    input.pint[lcount] = 0.5 * (input.pmid[lcount]+input.pmid[lcount-1]);
    input.hybi[lcount] = 0.5 * (input.hybm[lcount]+input.hybm[lcount-1]);
  }
  input.pint[lcount]=0.5*input.pmid[lcount];
  input.hybi[lcount]=0.5*input.hybm[lcount];
  input.hyai[lcount]=0.;
}
    
void
crm::Write_Pressure_NC_File( const string& pressfile, struct crm_input_data& input) 
{
  int lcount;

  NcFile pf = NcFile(pressfile,NcFile::CREATE,true);

  NcDimension levDim( "lev", input.plev );
  NcDimension ilevDim( "ilev", input.plev+1 );
  pf.addDimension( levDim );
  pf.addDimension( ilevDim );

  NcDimension dims[4];
  dims[0] = levDim;
  NcVariable<real_t> levVar = NcVariable<real_t>( "lev", dims, 1 );
  NcVariable<real_t> hyamVar = NcVariable<real_t>( "hyam", dims, 1 );
  NcVariable<real_t> hybmVar = NcVariable<real_t>( "hybm", dims, 1 );
  dims[0] = ilevDim;
  NcVariable<real_t> hyaiVar = NcVariable<real_t>( "hyai", dims, 1 );
  NcVariable<real_t> hybiVar = NcVariable<real_t>( "hybi", dims, 1 );
  pf.addVariable( levVar );
  pf.addVariable( hyamVar );
  pf.addVariable( hybmVar );
  pf.addVariable( hyaiVar );
  pf.addVariable( hybiVar );
  pf.endDefineMode();

  for (lcount = 0; lcount < input.plev ; ++lcount ) {
    levVar[lcount] = input.pmid[lcount];
    hyamVar[lcount] = input.hyam[lcount];
    hybmVar[lcount] = input.hybm[lcount];
    hyaiVar[lcount] = input.hyai[lcount];
    hybiVar[lcount] = input.hybi[lcount];
  }
  hyaiVar[lcount] = input.hyai[lcount];
  hybiVar[lcount] = input.hybi[lcount];
  levVar.write();
  hyamVar.write();
  hybmVar.write();
  hyaiVar.write();
  hybiVar.write();

  //  NcDimension timeDim( "time", NC_UNLIMITED, true );
  //  NcDimension latDim( "lat", 1 );
  //  NcDimension lonDim( "lon", 1 );
//   pf.addDimension( timeDim );
//   pf.addDimension( latDim );
//   pf.addDimension( lonDim );
  
}

void
crm::Write_Ozone_NC_File( const string& ozonefile, struct crm_input_data& input) 
{
  NcFile of = NcFile(ozonefile,NcFile::CREATE,true);
  
  NcDimension levDim( "lev", input.plev );
  NcDimension timeDim( "time", 12);
  NcDimension latDim( "lat", 2 );
  NcDimension lonDim( "lon", 1 );
  of.addDimension( lonDim );
  of.addDimension( latDim );
  of.addDimension( levDim );
  of.addDimension( timeDim );
  
  NcDimension dims[4];
  dims[0] = lonDim;
  NcVariable<real_t> lonVar = NcVariable<real_t>( "lon", dims, 1 );
  dims[0] = latDim;
  NcVariable<real_t> latVar = NcVariable<real_t>( "lat", dims, 1 );
  dims[0] = levDim;
  NcVariable<real_t> levVar = NcVariable<real_t>( "lev", dims, 1 );
  dims[0] = timeDim;
  NcVariable<int> dateVar = NcVariable<int>( "date", dims, 1 );
  NcVariable<int> datesecVar = NcVariable<int>( "datesec", dims, 1 );
  NcVariable<real_t> timeVar = NcVariable<real_t>( "time", dims, 1 );
  dims[0] = timeDim;
  dims[1] = latDim;
  dims[2] = levDim;
  dims[3] = lonDim;
  NcVariable<real_t> ozoneVar = NcVariable<real_t>( "OZONE", dims, 4 );
  
  of.addVariable( lonVar );
  of.addVariable( latVar );
  of.addVariable( levVar );
  of.addVariable( dateVar );
  of.addVariable( datesecVar );
  of.addVariable( timeVar );
  of.addVariable( ozoneVar );
  of.endDefineMode();

  lonVar = 0.;
  latVar[0] = -45.0;
  latVar[1] = 45.0;
  for (int levcount = 0; levcount < input.plev ; ++levcount )
    levVar[levcount] = input.pmid[levcount] / 100.0000;
  for (int mcount = 0; mcount < 12 ; ++mcount ) {
    for (int latcount = 0 ; latcount < 2 ; ++latcount ) 
      for ( int levcount = 0; levcount < input.plev ; ++levcount )
	ozoneVar[mcount*2*input.plev+latcount*input.plev+levcount] = 
	  input.o3mmr[levcount] * 28.9644 / 48.0000;
    dateVar[mcount] = 19900101+mcount*100;
    datesecVar[mcount] = 0;
    timeVar[mcount] = 1 + mcount*31;
  }
  
  latVar.write();
  lonVar.write();
  levVar.write();
  ozoneVar.write();
  dateVar.write();
  datesecVar.write();
  timeVar.write();
  
}

void
crm::Write_SCAM_NC_File( const string& scamfile, struct crm_input_data& input) 
{
  NcFile of = NcFile(scamfile,NcFile::CREATE,true);
  
  NcDimension lonDim( "lon", 1 );
  NcDimension latDim( "lat", 1 );
  NcDimension levDim( "lev", input.plev );
  NcDimension timeDim( "time", 2);
  of.addDimension( lonDim );
  of.addDimension( latDim );
  of.addDimension( levDim );
  of.addDimension( timeDim );
  
  NcDimension dims[4];
  dims[0] = lonDim;
  NcVariable<real_t> lonVar = NcVariable<real_t>( "lon", dims, 1 );
  dims[0] = latDim;
  NcVariable<real_t> latVar = NcVariable<real_t>( "lat", dims, 1 );
  dims[0] = levDim;
  NcVariable<real_t> levVar = NcVariable<real_t>( "lev", dims, 1 );
  dims[0] = lonDim;
  dims[1] = latDim;
  NcVariable<real_t> phisVar = NcVariable<real_t>( "phis", dims, 2 );
  NcVariable<real_t> oroVar = NcVariable<real_t>( "landfrac", dims, 2 );
  NcVariable<real_t> snowhiceVar = NcVariable<real_t>( "snowhice", dims, 2 );
  NcVariable<real_t> co2vmrVar = NcVariable<real_t>( "co2vmr", dims, 2 );
  NcVariable<real_t> n2ovmrVar = NcVariable<real_t>( "n2ovmr", dims, 2 );
  NcVariable<real_t> ch4vmrVar = NcVariable<real_t>( "ch4vmr", dims, 2 );
  NcVariable<real_t> f11vmrVar = NcVariable<real_t>( "f11vmr", dims, 2 );
  NcVariable<real_t> f12vmrVar = NcVariable<real_t>( "f12vmr", dims, 2 );
  NcVariable<real_t> sconVar = NcVariable<real_t>( "scon", dims, 2 );
  NcVariable<real_t> tauvisVar = NcVariable<real_t>( "tauvis", dims, 2 );
  NcVariable<real_t> aldirVar = NcVariable<real_t>( "aldir", dims, 2 );
  NcVariable<real_t> aldifVar = NcVariable<real_t>( "aldif", dims, 2 );
  NcVariable<real_t> asdirVar = NcVariable<real_t>( "asdir", dims, 2 );
  NcVariable<real_t> asdifVar = NcVariable<real_t>( "asdif", dims, 2 );
  dims[0] = timeDim;
  NcVariable<int> tsecVar = NcVariable<int>( "tsec", dims, 1 );
  NcVariable<int> bdateVar = NcVariable<int>( "bdate",dims,0 );
  dims[0] = timeDim;
  dims[1] = levDim;
  dims[2] = latDim;
  dims[3] = lonDim;
  NcVariable<real_t> TVar = NcVariable<real_t>( "T", dims, 4 );
  NcVariable<real_t> qVar = NcVariable<real_t>( "q", dims, 4 );
  NcVariable<real_t> divTVar = NcVariable<real_t>( "divT", dims, 4 );
  NcVariable<real_t> divqVar = NcVariable<real_t>( "divq", dims, 4 );
  NcVariable<real_t> cldVar = NcVariable<real_t>( "cld", dims, 4 );
  NcVariable<real_t> clwpVar = NcVariable<real_t>( "clwp", dims, 4 );
  dims[0] = timeDim;
  dims[1] = latDim;
  dims[2] = lonDim;
  NcVariable<real_t> PsVar = NcVariable<real_t>( "Ps", dims, 3 );
  NcVariable<real_t> omegaVar = NcVariable<real_t>( "omega", dims, 3 );
  NcVariable<real_t> TgVar = NcVariable<real_t>( "Tg", dims, 3 );
  NcVariable<real_t> TsairVar = NcVariable<real_t>( "Tsair", dims, 3 );
  
  of.addVariable( lonVar );
  of.addVariable( latVar  );
  of.addVariable( levVar  );
  of.addVariable( phisVar  );
  of.addVariable( oroVar  );
  of.addVariable( snowhiceVar );
  of.addVariable( co2vmrVar );
  of.addVariable( n2ovmrVar );
  of.addVariable( ch4vmrVar );
  of.addVariable( f11vmrVar );
  of.addVariable( f12vmrVar );
  of.addVariable( sconVar  );
  of.addVariable( tauvisVar );
  of.addVariable( aldirVar  );
  of.addVariable( aldifVar  );
  of.addVariable( asdirVar  );
  of.addVariable( asdifVar  );
  of.addVariable( tsecVar  );
  of.addVariable( bdateVar  );
  of.addVariable( TVar  );
  of.addVariable( qVar  );
  of.addVariable( divTVar  );
  of.addVariable( divqVar  );
  of.addVariable( cldVar  );
  of.addVariable( clwpVar  );
  of.addVariable( PsVar  );
  of.addVariable( omegaVar  );
  of.addVariable( TgVar  );
  of.addVariable( TsairVar  );
  of.endDefineMode();

  lonVar   = input.lon;
  latVar   = input.lat;
  phisVar   = input.phis;
  oroVar   = input.oro;
  snowhiceVar= input.snowh;
  co2vmrVar  = input.co2vmr;
  n2ovmrVar  = input.n2ovmr;
  ch4vmrVar  = input.ch4vmr;
  f11vmrVar  = input.f11vmr;
  f12vmrVar  = input.f12vmr;
  sconVar   = input.scon * 1000.;
  tauvisVar  = input.tauvis;
  aldirVar   = input.aldir;
  aldifVar   = input.aldif;
  asdirVar   = input.asdir;
  asdifVar   = input.asdif;
  bdateVar   = input.bdate;
  for (int timecount = 0 ; timecount < 2 ; ++timecount) {
    tsecVar[timecount]   = input.tsec+timecount*100000;
    PsVar[timecount]   = input.Ps;
    omegaVar[timecount]   = 0.;
    TgVar[timecount]   = input.Tg;
    TsairVar[timecount]   = input.Tsair;
    for (int levcount = 0 ; levcount < input.plev ; ++levcount) {
      TVar[timecount*input.plev+levcount]   = input.T[levcount];
      qVar[timecount*input.plev+levcount]   = input.q[levcount];
      divTVar[timecount*input.plev+levcount]   = 0.;
      divqVar[timecount*input.plev+levcount]   = 0.;
      cldVar[timecount*input.plev+levcount]   = input.cld[levcount];
      clwpVar[timecount*input.plev+levcount]   = input.clwp[levcount];
    }
  }
  for (int levcount = 0; levcount < input.plev ; ++levcount )
    levVar[levcount] = input.pmid[levcount];

  lonVar.write();
  latVar.write();
  levVar.write();
  phisVar.write();
  oroVar.write();
  snowhiceVar.write();
  co2vmrVar.write();
  n2ovmrVar.write();
  ch4vmrVar.write();
  f11vmrVar.write();
  f12vmrVar.write();
  sconVar.write();
  tauvisVar.write();
  aldirVar.write();
  aldifVar.write();
  asdirVar.write();
  asdifVar.write();
  tsecVar.write();
  bdateVar.write();
  TVar.write();
  qVar.write();
  divTVar.write();
  divqVar.write();
  cldVar.write();
  clwpVar.write();
  PsVar.write();
  omegaVar.write();
  TgVar.write();
  TsairVar.write();
}

void
crm::Write_SCAM_Quickstart_File( const string& scamqsfile, struct crm_input_data& input) 
{
  ofstream qsf(scamqsfile.data());
  if (!qsf) {
    cerr << "Error: cannot open file " << scamqsfile << endl;
    exit(1); 
  }
  qsf << "# SCCM Defaults File #" << endl;
  qsf << "runtype=\"2\"" << endl;
  qsf << "basedate=\"" << input.bdate << "\"" << endl;
  qsf << "basesecs=\"" << input.tsec << "\"" << endl;
  qsf << "iopstartoffset=\"0\"" << endl;
  qsf << "globaldatadir=\"./data/global/\"" << endl;
  qsf << "iopdatadir=\"./data/iop/\"" << endl;
  qsf << "boundarydatadir=\"./data/boundary/\"" << endl;
  qsf << "userdatadir=\"./userdata/\"" << endl;
  qsf << "histfile=\"crm.nc\"" << endl;
  qsf << "lat=\"" << input.lat << "\"" << endl;
  qsf << "lon=\"" << input.lon << "\"" << endl;
  qsf << "steplen=\"1800\"" << endl;
  qsf << "endstep=\"1\"" << endl;
  qsf << "savefreq=\"1\"" << endl;
  qsf << "timedisplayformat=\"0\"" << endl;
  qsf << "showsettings=\"0\"" << endl;
  qsf << "analysisfile=\"./data/global/cami_0000-09-01_64x128_L26_c030918.nc\"" << endl;
  qsf << "modelfile=\"./data/global/cami_0000-09-01_64x128_L26_c030918.nc\"" << endl;
  qsf << "userfile=\"\"" << endl;
  qsf << "iopfile=\"./userdata/iop_crm.nc\"" << endl;
  qsf << "lsminifile=\"./data/boundary/clmi_0000-09-01_64x128_T42_c020514.landvec.nc\"" << endl;
  qsf << "ozonfile=\"./userdata/ozone_crm.nc\"" << endl;
  qsf << "pressfile=\"./userdata/press_crm.nc\"" << endl;
  qsf << "absemsfile=\"./data/boundary/abs_ems_factors_fastvx.052001.nc\"" << endl;
  qsf << "aeropticsfile=\"./data/boundary/AerosolOptics_c040105.nc\"" << endl;
  qsf << "aermassfile=\"./data/boundary/AerosolMass_V_64x128_clim_c031022.nc\"" << endl;
  qsf << "lsmsurffile=\"./data/boundary/clms_64x128_c020514.nc\"" << endl;
  qsf << "lsmpftfile=\"./data/boundary/pft-physiology\"" << endl;
  qsf << "sstfile=\"./data/boundary/sst_HadOIBl_bc_64x128_clim_c020411.nc\"" << endl;
  qsf << "sicfile=\"\"" << endl;
  qsf << "savefields=\"T Q CLOUD QRL QRS SRFRAD FUL FDL FSUL FSDL FUS FDS FSUS FSDS \"" << endl;
  qsf << "switch_desc1=\"Logical Switch 1\"" << endl;
  qsf << "switch1=\"0\"" << endl;
  qsf << "switch_desc2=\"Logical Switch 2\"" << endl;
  qsf << "switch2=\"0\"" << endl;
  qsf << "switch_desc3=\"Logical Switch 3\"" << endl;
  qsf << "switch3=\"0\"" << endl;
  qsf << "switch_desc4=\"Logical Switch 4\"" << endl;
  qsf << "switch4=\"0\"" << endl;
  qsf << "switch_desc5=\"Logical Switch 5\"" << endl;
  qsf << "switch5=\"0\"" << endl;
  qsf << "switch_desc6=\"Logical Switch 6\"" << endl;
  qsf << "switch6=\"0\"" << endl;
  qsf << "switch_desc7=\"Logical Switch 7\"" << endl;
  qsf << "switch7=\"0\"" << endl;
  qsf << "switch_desc8=\"Logical Switch 8\"" << endl;
  qsf << "switch8=\"0\"" << endl;
  qsf << "switch_desc9=\"Logical Switch 9\"" << endl;
  qsf << "switch9=\"0\"" << endl;
  qsf << "switch_desc10=\"Logical Switch 10\"" << endl;
  qsf << "switch10=\"0\"" << endl;
  qsf << "switch_desc11=\"Logical Switch 11\"" << endl;
  qsf << "switch11=\"0\"" << endl;
  qsf << "switch_desc12=\"Logical Switch 12\"" << endl;
  qsf << "switch12=\"0\"" << endl;
  qsf << "switch_desc13=\"Logical Switch 13\"" << endl;
  qsf << "switch13=\"0\"" << endl;
  qsf << "switch_desc14=\"Logical Switch 14\"" << endl;
  qsf << "switch14=\"0\"" << endl;
  qsf << "switch_desc15=\"Logical Switch 15\"" << endl;
  qsf << "switch15=\"0\"" << endl;
  qsf << "switch_desc16=\"Logical Switch 16\"" << endl;
  qsf << "switch16=\"0\"" << endl;
  qsf << "switch_desc17=\"Logical Switch 17\"" << endl;
  qsf << "switch17=\"0\"" << endl;
  qsf << "switch_desc18=\"Logical Switch 18\"" << endl;
  qsf << "switch18=\"0\"" << endl;
  qsf << "switch_desc19=\"Logical Switch 19\"" << endl;
  qsf << "switch19=\"0\"" << endl;
  qsf << "switch_desc20=\"CRM Switch\"" << endl;
  qsf << "switch20=\"1\"" << endl;
}

//  Set the CRM text output file
void
crm::SaveCRMTextFile( const QString& filename )
{
  CRMTextFile = filename;
}

bool
crm::IsCRMTextFile()
{
  return !CRMTextFile.isEmpty();
}

void
crm::Read_Output_Data( const string& ncoutfile, struct crm_output_data& output)
{
  NcFile of = NcFile(ncoutfile.data(),NcFile::READ);
  NcDimension levDim = of.dimension( "lev" );
  NcDimension ilevDim = of.dimension( "ilev" );
  output.lev = levDim.size();
  output.ilev = ilevDim.size();
  NcVariable<real_t> qrsVar = of.variable( "QRS" );
  NcVariable<real_t> qrlVar = of.variable( "QRL" );
  NcVariable<real_t> fulVar = of.variable( "FUL" );
  NcVariable<real_t> fdlVar = of.variable( "FDL" );
  NcVariable<real_t> fusVar = of.variable( "FUS" );
  NcVariable<real_t> fdsVar = of.variable( "FDS" );
  NcVariable<real_t> levVar = of.variable( "lev" );
  NcVariable<real_t> ilevVar = of.variable( "ilev" );
  qrsVar.read();
  qrlVar.read();
  fulVar.read();
  fdlVar.read();
  fusVar.read();
  fdsVar.read();
  levVar.read();
  ilevVar.read();
  for (int lcount = 0 ; lcount < output.lev ; ++lcount) {
    output.pmid[lcount] = levVar[lcount];
    output.qrs[lcount] = qrsVar[lcount];
    output.qrl[lcount] = qrlVar[lcount];
  }
  for (int lcount = 0 ; lcount < output.ilev ; ++lcount) {
    output.pint[lcount] = ilevVar[lcount];
    output.ful[lcount] = fulVar[lcount];
    output.fdl[lcount] = fdlVar[lcount];
    output.fus[lcount] = fusVar[lcount];
    output.fds[lcount] = fdsVar[lcount];
  }
}

void
crm::Print_Output_Data( const QString& txtfile, struct crm_output_data& output)
{
  ofstream of(txtfile);
  if (!of) {
    cerr << "Error: cannot open file " << txtfile << endl;
    exit(1); 
  }
  of << "Begin CCM3 Column Radiation Model" << endl;
  of << "CCM3 CRM Results:" << endl;
  of << "Conventions:" << endl;
  of << "Shortwave fluxes are positive downward" << endl;
  of << "Longwave fluxes are positive upward" << endl;
  of << "Net Radiative fluxes are positive downward (into the system)'" << endl;
  of << "Fluxes defined to be zero are not reported (e.g., LWdown flx TOA)' " << endl;
  of << "Abbreviations, Acronyms and Definitions:" << endl;
  of << "LW   = Longwave" << endl;
  of << "LWCF = Longwave Cloud Forcing" << endl;
  of << "NCF  = Net Cloud Forcing = SWCF+LWCF" << endl;
  of << "NIR  = Near Infrared (0.7 < lambda < 5.0 microns)" << endl;
  of << "N7   = NOAA7 satellite NIR instrument weighted flux" << endl;
  of << "NRF  = Net Radiative Flux: sum of SW and LW fluxes" << endl;
  of << "SW   = Shortwave" << endl;
  of << "SWCF = Shortwave Cloud Forcing" << endl;
  of << "TOA  = Top of Atmosphere" << endl;
  of << "Vis  = Visible (0.2 < lambda < 0.7 microns)" << endl;
  of << "atm  = Atmosphere" << endl;
  of << "clr  = Clear sky (diagnostic computation with no clouds)" << endl;
  of << "ctr  = Center" << endl;
  of << "dff  = Diffuse flux" << endl;
  of << "drc  = Direct flux" << endl;
  of << "dwn  = Downwelling" << endl;
  of << "frc  = Fraction" << endl;
  of << "lqd  = Liquid" << endl;
  of << "mpc  = Mass path column" << endl;
  of << "net  = Net flux = downwelling minus upwelling flux" << endl;
  of << "spc  = Spectral" << endl;
  of << "sfc  = Surface level" << endl;
  of << "vmr  = Volume mixing ratio" << endl;
  of << "wvl  = Wavelength" << endl;
  of << "um   = Microns" << endl;
  of << "up   = Upwelling" << endl;
  of << " " << endl;
  of << "Heating rates:" << endl;
  of << "Level   Pressure       SW          LW          Net" << endl;
  of << "   #       mb        K day-1     K day-1     K day-1" << endl;
  for (int lcount = 0 ; lcount < output.lev ; ++lcount) {
    of.width(4);
    of << lcount + 1;
    of.width(12);
    of.precision(3);
    of.flags( ios::fixed );
    of << output.pmid[lcount];
    of.width(12);
    of.precision(5);
    of << output.qrs[lcount]*86400;
    of.width(12);
    of << output.qrl[lcount]*86400;
    of.width(12);
    of << (output.qrs[lcount] + output.qrl[lcount]) *86400;
    of << endl;
  }

  of << " " << endl;
  of << "Fluxes:" << endl;
  of << "Level   Pressure      FUL          FDL      FUS       FDS " << endl;
  of << "   #       mb        W/m**2     W/m**2     W/m**2   W/m**2 " << endl;
  for (int lcount = 0 ; lcount < output.ilev ; ++lcount) {
    of.width(4);
    of << lcount + 1;
    of.width(12);
    of.precision(3);
    of.flags( ios::fixed );
    of << output.pint[lcount];
    of.width(12);
    of.precision(5);
    of << output.ful[lcount];
    of.width(12);
    of << output.fdl[lcount];
    of.width(12);
    of << output.fus[lcount];
    of.width(12);
    of << output.fds[lcount];
    of << endl;
  }

  of << "End CAM2 CRM" << endl;
  
}
