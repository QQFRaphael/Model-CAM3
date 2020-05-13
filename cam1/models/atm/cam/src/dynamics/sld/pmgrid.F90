#include <misc.h>
#include <params.h>

module pmgrid

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize grid point resolution parameters
! 
! Author: 
! 
!-----------------------------------------------------------------------
  implicit none

  private

  public plon, plev, plat, plnlv, plevp
  public masterproc, iam, begirow, endirow, beglat, endlat
  public numlats, dyngrid_set

! Grid point resolution parameters

  integer, parameter :: plon       = PLON  ! number of longitudes
  integer, parameter :: plev       = PLEV  ! number of vertical levels
  integer, parameter :: plat       = PLAT  ! number of latitudes
  integer, parameter :: plevp      = plev + 1 ! plev + 1
  integer, parameter :: plnlv      = plon*plev     ! Length of multilevel field slice

  integer iam
  integer begirow    ! beg. index for lat pairs owned by a proc
  integer endirow    ! end. index for lat pairs owned by a proc
  integer beglat     ! beg. index for lats owned by a given proc
  integer endlat     ! end. index for lats owned by a given proc
  integer numlats    ! number of lats owned by a given proc
  logical masterproc 
  logical :: dyngrid_set = .false. ! flag indicates dynamics grid has been set
#if ( ! defined SPMD )
  parameter (iam        = 0)
  parameter (begirow    = 1)
  parameter (endirow    = plat/2)
  parameter (beglat     = 1)
  parameter (endlat     = plat)
  parameter (numlats    = plat)
  parameter (masterproc = .true.)
#endif
end module pmgrid


