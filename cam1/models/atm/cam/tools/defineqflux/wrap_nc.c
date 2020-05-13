#include <sys/types.h>     // size_t
#include <netcdf.h>

#include "definesomic.h"

int ret;

void wrap_nc_inq_varid (int ncid, char *name, int *varid)
{
  if ((ret = nc_inq_varid (ncid, name, varid)) != NC_NOERR)
    err_exit ("nc_inq_varid %s failure: %s\n", name, nc_strerror (ret));
}

void wrap_nc_get_att_text (int ncid, int varid, const char *name, char *tp)
{
  if ((ret = nc_get_att_text (ncid, varid, name, tp)) != NC_NOERR)
    err_exit ("nc_get_att_text %s %s failure: %s\n", name, tp, nc_strerror (ret));
}

void wrap_nc_put_att_text (int ncid, int varid, const char *name, size_t len, const char *tp)
{
  if ((ret = nc_put_att_text (ncid, varid, name, len, tp)) != NC_NOERR)
    err_exit ("nc_put_att_text %s %s failure: %s\n", name, tp, nc_strerror (ret));
}

void wrap_nc_get_var_double (int ncid, int varid, double *values)
{
  if ((ret = nc_get_var_double (ncid, varid, values)) != NC_NOERR)
    err_exit ("nc_get_var_double %d failure: %s\n", varid, nc_strerror (ret));
}

void wrap_nc_get_var_int (int ncid, int varid, int *values)
{
  if ((ret = nc_get_var_int (ncid, varid, values)) != NC_NOERR)
    err_exit ("nc_get_var_int %d failure: %s\n", varid, nc_strerror (ret));
}

void wrap_nc_get_vara_double (int ncid, int varid, size_t *start, size_t *count, double *values)
{
  ret = nc_get_vara_double (ncid, varid, start, count, values);
  if (ret != NC_NOERR)
    err_exit ("nc_get_vara_double %d failure: %s\n", varid, nc_strerror (ret));
}

void wrap_nc_put_vara_double (int ncid, int varid, size_t *start, size_t *count, double *values)
{
  ret = nc_put_vara_double (ncid, varid, start, count, values);
  if (ret != NC_NOERR)
    err_exit ("nc_put_vara_double %d failure: %s\n", varid, nc_strerror (ret));
}

void wrap_nc_inq_unlimdim (int ncid, int *dimid)
{
  if ((ret = nc_inq_unlimdim (ncid, dimid)) != NC_NOERR)
    err_exit ("nc_inq_unlimdim %d failure: %s\n", dimid, nc_strerror (ret));
}

void wrap_nc_inq_dimlen (int ncid, int dimid, size_t *length)
{
  if ((ret = nc_inq_dimlen (ncid, dimid, length)) != NC_NOERR)
    err_exit ("nc_inq_dimlen %d failure: %s\n", dimid, nc_strerror (ret));
}

void wrap_nc_inq_dimid (int ncid, char *name, int *dimid)
{
  if ((ret = nc_inq_dimid (ncid, name, dimid)) != NC_NOERR)
    err_exit ("nc_inq_dimid %s failure: %s\n", name, nc_strerror (ret));
}

void wrap_nc_def_var (int ncid, char *name, nc_type xtype, int ndims, int *dimids, int *varid)
{
  int ret;

  if ((ret = nc_def_var (ncid, name, xtype, ndims, dimids, varid)) != NC_NOERR)
    err_exit ("nc_def_var %s failure: %s\n", name, nc_strerror (ret));
}

void wrap_nc_enddef (int ncid)
{
  if ((ret = nc_enddef (ncid)) != NC_NOERR)
    err_exit ("nc_enddef failure: %s\n", nc_strerror (ret));
}

void wrap_nc_redef (int ncid)
{
  if ((ret = nc_redef (ncid)) != NC_NOERR)
    err_exit ("nc_redef failure: %s\n", nc_strerror (ret));
}

void wrap_nc_close (int ncid)
{
  if ((ret = nc_close (ncid)) != NC_NOERR)
    err_exit ("nc_close failure: %s\n", nc_strerror (ret));
}

