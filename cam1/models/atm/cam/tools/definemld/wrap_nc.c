#include <sys/types.h>     // size_t
#include <netcdf.h>

void wrap_nc_inq_varid (int ncid, char *name, int *varid)
{
  int ret;

  if ((ret = nc_inq_varid (ncid, name, varid)) != NC_NOERR)
    err_exit ("nc_inq_varid %s failure: %s\n", name, nc_strerror (ret));
}

void wrap_nc_get_var_double (int ncid, int varid, double *values)
{
  int ret;

  if ((ret = nc_get_var_double (ncid, varid, values)) != NC_NOERR)
    err_exit ("nc_get_var_double %d failure: %s\n", varid, nc_strerror (ret));
}

void wrap_nc_get_var_int (int ncid, int varid, int *values)
{
  int ret;

  if ((ret = nc_get_var_int (ncid, varid, values)) != NC_NOERR)
    err_exit ("nc_get_var_int %d failure: %s\n", varid, nc_strerror (ret));
}

void wrap_nc_get_vara_double (int ncid, int varid, size_t *start, size_t *count, double *values)
{
  int ret;

  ret = nc_get_vara_double (ncid, varid, start, count, values);
  if (ret != NC_NOERR)
    err_exit ("nc_get_vara_double %d failure: %s\n", varid, nc_strerror (ret));
}

void wrap_nc_put_vara_double (int ncid, int varid, size_t *start, size_t *count, double *values)
{
  int ret;

  ret = nc_put_vara_double (ncid, varid, start, count, values);
  if (ret != NC_NOERR)
    err_exit ("nc_put_vara_double %d failure: %s\n", varid, nc_strerror (ret));
}

void wrap_nc_put_var_double (int ncid, int varid, double *values)
{
  int ret;

  ret = nc_put_var_double (ncid, varid, values);
  if (ret != NC_NOERR)
    err_exit ("nc_put_var_double %d failure: %s\n", varid, nc_strerror (ret));
}

void wrap_nc_inq_unlimdim (int ncid, int *dimid)
{
  int ret;

  if ((ret = nc_inq_unlimdim (ncid, dimid)) != NC_NOERR)
    err_exit ("nc_inq_unlimdim %d failure: %s\n", dimid, nc_strerror (ret));
}

void wrap_nc_inq_dimlen (int ncid, int dimid, size_t *length)
{
  int ret;

  if ((ret = nc_inq_dimlen (ncid, dimid, length)) != NC_NOERR)
    err_exit ("nc_inq_dimlen %d failure: %s\n", dimid, nc_strerror (ret));
}

void wrap_nc_inq_dimid (int ncid, char *name, int *dimid)
{
  int ret;

  if ((ret = nc_inq_dimid (ncid, name, dimid)) != NC_NOERR)
    err_exit ("nc_inq_dimid %s failure: %s\n", name, nc_strerror (ret));
}

void wrap_nc_def_dim (int ncid, char *name, size_t len, int *dimidp)
{
  int ret;

  if ((ret = nc_def_dim (ncid, name, len, dimidp)) != NC_NOERR)
    err_exit ("nc_def_dim %s failure: %s\n", name, nc_strerror (ret));
}

void wrap_nc_def_var (int ncid, char *name, nc_type xtype, int ndims, int *dimids, int *varid)
{
  int ret;

  if ((ret = nc_def_var (ncid, name, xtype, ndims, dimids, varid)) != NC_NOERR)
    err_exit ("nc_def_var %s failure: %s\n", name, nc_strerror (ret));
}

void wrap_nc_enddef (int ncid)
{
  int ret;

  if ((ret = nc_enddef (ncid)) != NC_NOERR)
    err_exit ("nc_enddef failure: %s\n", nc_strerror (ret));
}

void wrap_nc_redef (int ncid)
{
  int ret;

  if ((ret = nc_redef (ncid)) != NC_NOERR)
    err_exit ("nc_redef failure: %s\n", nc_strerror (ret));
}

void wrap_nc_close (int ncid)
{
  int ret;

  if ((ret = nc_close (ncid)) != NC_NOERR)
    err_exit ("nc_close failure: %s\n", nc_strerror (ret));
}

