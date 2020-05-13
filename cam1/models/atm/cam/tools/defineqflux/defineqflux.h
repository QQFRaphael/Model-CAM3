#include <sys/types.h>   // size_t
#include <netcdf.h>

// Grid constants

#define PLON 128
#define PLAT 64

// Macros

#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

// typedefs 

typedef enum {false = 0, true = 1} Boolean;

// Function prototypes

void err_exit (const char *, ...);
void mksith (double [PLAT][PLON], double [PLAT][PLON], double *, int *, int);
void check_dims_consistent (int, char *);
void check_grid_consistent (const int ncid, const char *, const int *);
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

int gepp (double **, double **, const int, double[], int[]);
int backsolve (double **, const int, const double[], double[]);
void printeq (double **, const int, const double[]);
void print_matrix (double **, const int);
