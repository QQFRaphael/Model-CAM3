#include <sys/types.h>     // size_t

#include "defineqflux.h"

void check_dims_consistent (int ncid,
			    char *filename)
{
  int latid;
  int lonid;
  size_t latsiz;
  size_t lonsiz;

  wrap_nc_inq_dimid (ncid, "lat", &latid);
  wrap_nc_inq_dimlen (ncid, latid, &latsiz);

  if (latsiz != PLAT)
    err_exit ("Compile-time PLAT=%d inconsistent with dimension lat=%d from file %s\n",
	      PLAT, latsiz, filename);

  wrap_nc_inq_dimid (ncid, "lon", &lonid);
  wrap_nc_inq_dimlen (ncid, lonid, &lonsiz);

  if (lonsiz != PLON)
    err_exit ("Compile-time PLON=%d inconsistent with dimension lon=%d from file %s\n",
	      PLON, lonsiz, filename);
}
  
void check_grid_consistent (const int ncid,
			    const char *filename,
			    const int nlonin[PLAT])
{
  int j;
  int nlonid;
  int nlon[PLAT];

  if (nc_inq_varid (ncid, "nlon", &nlonid) == NC_NOERR)    // reduced grid
    wrap_nc_get_var_int (ncid, nlonid, nlon);
  else                                                    // full grid
    for (j = 0; j < PLAT; j++)
      nlon[j] = PLON;

  for (j = 0; j < PLAT; j++)
    if (nlon[j] != nlonin[j])
      err_exit ("file %s nlon[%d]=%d inconsistent with nlonin[%d]=%d\n",
		filename, j, nlon[j], j, nlonin[j]);
}
  
