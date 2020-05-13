#include <misc.h>

module rgrid

  use pmgrid, only: plat
  use pspect, only: pmmax, pmax
  use infnan, only: bigint

  implicit none

  integer :: nlon(plat)        = bigint ! num longitudes per latitude
  integer :: beglatpair(pmmax) = bigint
#if ( defined SCAM )
  integer :: nmmax(plat)     = bigint
#else
  integer :: nmmax(plat/2)     = bigint
#endif
  integer :: wnummax(plat)     = bigint ! cutoff Fourier wavenumber

  logical :: fullgrid                   ! true => no grid reduction towards poles
end module rgrid
