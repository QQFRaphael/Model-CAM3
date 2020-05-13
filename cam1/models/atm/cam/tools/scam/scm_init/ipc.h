#ifndef _IPC_H
#define _IPC_H

/* INTERPROCESS COMMUNICATION */
/* used by scam to facilitate communication */
/* between the GUI and the model */


#include "realtype.h"
#include "max.h"       


struct field_info {
    char   name[MAX_FIELDNAME_LEN];
    char   longname[MAX_FIELDNAME_LEN];
    char   units[64];
    char   std_units[64];
    real_t  mult;
    real_t  min;
    real_t  max;
    int    is_shown;
    int    is_modifiable;
    int    is_averaged;
    int    size;
};

struct field_data {
    int id;
    real_t data[MAX_LEVELS];
}; 
  
struct level_data {
    int nlev;
    int ilev;
    real_t levels[MAX_LEVELS];
    real_t ilevels[MAX_LEVELS];
};

struct init_data {
    real_t lat;
    real_t lon;
    int    sw[NUM_SWITCHES];
    int    bd;
    int    bs;
    int    rt;
    int    sl;
    int    ie;
    int    rst;
    char   ana[MAX_PATH_LEN];
    char   ini[MAX_PATH_LEN];
    char   iop[MAX_PATH_LEN];
    char   lsmini[MAX_PATH_LEN];
    char   ozo[MAX_PATH_LEN];
    char   prs[MAX_PATH_LEN];
    char   sic[MAX_PATH_LEN];
    char   sst[MAX_PATH_LEN];
    char   abs[MAX_PATH_LEN];
    char   optics[MAX_PATH_LEN];
    char   mass[MAX_PATH_LEN];
    char   lsmpft[MAX_PATH_LEN];
    char   lsmsurf[MAX_PATH_LEN];
    char   usr[MAX_PATH_LEN];
};

/*-----------------------------------------------------*/    
/* structs and parameters used only by rpc version     */
/*-----------------------------------------------------*/    

/*
 * the '%' in front of the following #define makes it get
 *  placed unchanged into the rpc header file by rpcgen
 *   (rpcgen defines RPC_HDR when compiling the header)
 */
#ifdef RPC_HDR 
%#define RPC_SVC_FG /* make server run in the foreground */
#endif
struct dataset_var_real {
    int      found;
    real_t    data[MAX_TIME_DIM];
};

struct dataset_var_int {
    int      found;
    int      data[MAX_TIME_DIM];
};

struct dataset_req { 
    char      dataset[MAX_PATH_LEN];
    char      varname[64]; 
};
/*--------------------------------------*/    


/*-----------------------------------------*/    
/* parameters and structs for fifo version */
/*-----------------------------------------*/    

#define INIT_MODEL_REQ_ID        1
#define FIELD_INFO_REQ_ID        2
#define GET_FIELD_REQ_ID         3
#define GET_LEVELS_REQ_ID        4
#define SET_FIELD_REQ_ID         5
#define RESET_FIELD_REQ_ID       6
#define STEP_REQ_ID              7
/* Increased the following to allow lots of levels for CRM */
#define REQ_BUFSIZE           4096
#define REPLY_BUFSIZE         4096

/* real_t rpad was added below to get some compilers (IRIX) to */
/* align the srv_req.buf on on real_t boundary.  If no real_t  */
/* is present then it at best alligns things on int boundaries */

struct srv_req {
    real_t rpad;    /* padding to align buf on 8 byte boundary */
    int id;                    /* id of request */
    int pad;                    /* padding to align buf on 8 byte boundary */
    char buf[REQ_BUFSIZE];
};

struct srv_reply {
    real_t rpad;    /* padding to align buf on 8 byte boundary */
    int status;
    int pad;                    /* padding to align buf on 8 byte boundary */
    char buf[REPLY_BUFSIZE];
};

/*--------------------------------------*/    


typedef struct field_info field_info;
typedef struct field_data field_data;
typedef struct level_data level_data;
typedef struct init_data  init_data;
typedef struct dataset_var_real dataset_var_real;
typedef struct dataset_var_int dataset_var_int;
typedef struct dataset_req  dataset_req;
typedef struct srv_req    srv_req;
typedef struct srv_reply  srv_reply;    

#endif /* _IPC_H */








