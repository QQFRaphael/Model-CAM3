#include <misc.h>
#include <params.h>

subroutine herzin(pkdim   ,pf      ,f       ,fst     ,fsb     , &
                  sig     ,dsig    ,sigdp   ,kdp     ,fdp     , &
                  nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interpolate field on vertical slice to vertical departure point using
! Hermite cubic interpolation.
! 
! Method: 
! 
! Author: 
! Original version:  J. Olson
! Standardized:      J. Rosinski, June 1992
! Reviewed:          D. Williamson, P. Rasch, August 1992
! Reviewed:          D. Williamson, P. Rasch, March 1996
!
!-----------------------------------------------------------------------
!
! $Id: herzin.F90,v 1.1.44.2 2004/09/23 17:29:23 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plev
!-----------------------------------------------------------------------
   implicit none
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: pkdim             ! vertical dimension
   integer, intent(in) :: pf                ! dimension (number of fields)
!
   real(r8), intent(in) :: f    (plon,pkdim,pf) ! fields
   real(r8), intent(in) :: fst  (plon,pkdim,pf) ! z-derivatives at top edge of interval
   real(r8), intent(in) :: fsb  (plon,pkdim,pf) ! z-derivatives at bot edge of interval
   real(r8), intent(in) :: sig  (pkdim)         ! vertical grid coordinates
   real(r8), intent(in) :: dsig (pkdim)         ! intervals between vertical grid pts.
   real(r8), intent(in) :: sigdp(plon,plev)     ! vertical coord. of departure point
!
   integer, intent(in) :: kdp  (plon,plev)  ! vertical index  of departure point
   integer, intent(in) :: nlon
!
! Output arguments
!
   real(r8), intent(out) :: fdp(plon,plev,pf)   ! z-interpolants
!
!-----------------------------------------------------------------------
!
!  pkdim   Vertical dimension of vertical slice arrays.
!  pf      Number of fields being interpolated.
!  f       Vertical slice of data to be interpolated.
!  fst     z-derivatives at the top edge of each interval contained in f
!  fsb     z-derivatives at the bot edge of each interval contained in f
!  sig     Sigma values corresponding to the vertical grid
!  dsig    Increment in sigma value for each interval in vertical grid.
!  sigdp   Sigma value at the trajectory midpoint or endpoint for each
!          gridpoint in a vertical slice from the global grid.
!  kdp     Vertical index for each gridpoint.  This index points into a
!          vertical slice array whose vertical grid is given by sig.
!          E.g.,   sig(kdp(i,j)) .le. sigdp(i,j) .lt. sig(kdp(i,j)+1) .
!  fdp     Value of field at the trajectory midpoints or endpoints.
!
!---------------------------Local variables-----------------------------
!
   integer i,k,m             ! indices
!
   real(r8) dzk                  ! vert interval containing the dep. pt.
   real(r8) zt                   ! |
   real(r8) zb                   ! |
   real(r8) ht (plon)            ! | -- interpolation coefficients
   real(r8) hb (plon)            ! |
   real(r8) dht(plon)            ! |
   real(r8) dhb(plon)            ! |
!
!-----------------------------------------------------------------------
!
   do k=1,plev
      do i=1,nlon
         dzk = dsig(kdp(i,k))
         zt = ( sig(kdp(i,k)+1) - sigdp(i,k) )/dzk
         zb = 1. - zt
         ht (i) = ( 3.0 - 2.0*zt )*zt**2
         hb (i) = ( 3.0 - 2.0*zb )*zb**2
         dht(i) = -dzk*( zt - 1. )*zt**2
         dhb(i) =  dzk*( zb - 1. )*zb**2
      end do
!
! Loop over fields.
!
      do m=1,pf
         do i=1,nlon
            fdp(i,k,m) = f(i,kdp(i,k)  ,m)* ht(i) +  &
               fst(i,kdp(i,k),m)*dht(i) +  &
               f(i,kdp(i,k)+1,m)* hb(i) +  &
               fsb(i,kdp(i,k),m)*dhb(i)
         end do
      end do
   end do
!
   return
end subroutine herzin
