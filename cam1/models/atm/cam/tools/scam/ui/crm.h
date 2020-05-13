#ifndef crm_h
#define crm_h

#include <max.h>
#include <qstring.h>
#define CRM crm::Instance()
#define OZONEFILE "/ozone_crm.nc"
#define PRESSFILE "/press_crm.nc"
#define SCAMFILE "/iop_crm.nc"
#define SCAMQSFILE "/crm_quickstart"

class crm 
{
public:
  ~crm();
  string   Setup( string& startupfile, const string& infilename );
  void   PostProcess( const string& ncoutfile );
  void   SaveCRMTextFile( const QString& filename );
  bool   IsCRMTextFile();
  static crm& Instance();
  
private:
  static crm* _instance;
  QString CRMTextFile; //CRM text output file
  
  struct crm_input_data {
    int plev;
    real_t lat;
    real_t pmid[MAX_LEVELS];
    real_t pint[MAX_LEVELS];
    real_t hybi[MAX_LEVELS];
    real_t hyai[MAX_LEVELS];
    real_t hybm[MAX_LEVELS];
    real_t hyam[MAX_LEVELS];
    real_t T[MAX_LEVELS];
    real_t q[MAX_LEVELS];
    real_t o3mmr[MAX_LEVELS];
    real_t cld[MAX_LEVELS];
    real_t clwp[MAX_LEVELS];
    real_t Ps;
    real_t Tsair;
    real_t Tg;
    real_t oro;
    real_t snowh;
    real_t asdir;
    real_t asdif;
    real_t aldir;
    real_t aldif;
    real_t co2vmr;
    real_t n2ovmr;
    real_t ch4vmr;
    real_t f11vmr;
    real_t f12vmr;
    real_t tauvis;
    real_t scon;
    int iyear_AD;
    real_t lon;
    real_t phis;
    int tsec;
    int bdate;
  };

  
  struct crm_output_data {
    int lev;
    int ilev;
    real_t pmid[MAX_LEVELS];
    real_t pint[MAX_LEVELS];
    real_t qrl[MAX_LEVELS];
    real_t qrs[MAX_LEVELS];
    real_t ful[MAX_LEVELS];
    real_t fdl[MAX_LEVELS];
    real_t fus[MAX_LEVELS];
    real_t fds[MAX_LEVELS];
  };

  crm();
  void Read_Input_Data( const string& infilename, struct crm_input_data& input);
  void Write_Pressure_NC_File( const string& pressfile, struct crm_input_data& input);
  void Write_Ozone_NC_File( const string& ozonefile, struct crm_input_data& input);
  void Write_SCAM_NC_File( const string& scamfile, struct crm_input_data& input);
  void Write_SCAM_Quickstart_File( const string& scamqsfile, struct crm_input_data& input);
  void Read_Output_Data( const string& ncoutfile, struct crm_output_data& output);
  void Print_Output_Data( const QString& txtfile, struct crm_output_data& output);
};

inline crm& crm::Instance() { 
    if ( _instance == NULL )  _instance = new crm();
    return *_instance; }
#endif
