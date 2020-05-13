#include <sys/types.h>   // size_t
#include <netcdf.h>

// Grid constants

#define PLON 96
#define PLAT 48

// macros

#define land(i,j) (ioro[j][i] == 1)

// typedefs 

typedef enum {false = 0, true = 1} Boolean;


// Function prototypes

void err_exit (const char *, ...);
void usage_exit (char *);

void wrap_nc_inq_varid (int, char *, int *);
void wrap_nc_get_var_double (int, int, double *);
void wrap_nc_get_var_int (int, int, int *);
void wrap_nc_get_vara_double (int, int, size_t *, size_t *, double *);
void wrap_nc_put_vara_double (int, int, size_t *, size_t *, double *);
void wrap_nc_inq_unlimdim (int, int *);
void wrap_nc_inq_dimlen (int, int, size_t *);
void wrap_nc_inq_dimid (int, char *, int *);
void wrap_nc_def_var (int, char *, nc_type, int, int *, int *);
void wrap_nc_enddef (int);
void wrap_nc_redef (int);
void wrap_nc_close (int);
