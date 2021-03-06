#include <misc.h>
#include <params.h>
   subroutine tstep(lm      ,zdt     ,ztdtsq  )
!-----------------------------------------------------------------------
!
! Solution of the vertically coupled system of equations arising
! from the semi-impicit equations for each spectral element along
! two dimensional wavenumber n.  The inverse matrix depends
! only on two dimensional wavenumber and the reference atmosphere.
! It is precomputed and stored for use during the forecast. The routine
! overwrites the d,T and lnps coefficients with the new values.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, D. Williamson, August 1992
! Reviewed:          B. Boville, D. Williamson, April 1996
!
!-----------------------------------------------------------------------
!
! $Id: tstep.F90,v 1.2.28.5 2004/03/30 22:10:59 rosinski Exp $
! $Author: rosinski $
!
!-----------------------------------------------------------------------
      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid
      use pspect
      use comspe
      use commap

      implicit none

!-----------------------------------------------------------------------
#include <comhyb.h>
!-----------------------------------------------------------------------
!
! Input arguments
!
      integer, intent(in) :: lm            ! local Fourier wavenumber index

      real(r8), intent(in) :: zdt             ! timestep, dt (seconds)
      real(r8), intent(in) :: ztdtsq(pnmax)   ! dt*(n(n+1)/a^2 where n is 2-d wavenumber
!
!---------------------------Local workspace-----------------------------
!
      real(r8) z(2*pnmax,plev) ! workspace for computation of spectral array d
      real(r8) hhref           ! href/2 (reference hydrostatic matrix / 2)
      real(r8) hbps            ! bps/2 (ref. coeff. for lnps term in div. eq. / 2)
      real(r8) ztemp           ! temporary workspace

      integer m            ! global wavenumber index
      integer n,j          ! 2-d wavenumber index
      integer k,kk         ! level indices
      integer lmr,lmc      ! real and imaginary spectral indices
      integer ir,ii        ! real and imaginary spectral indices
      integer nn           ! real and imaginary spectral indices
!
!-----------------------------------------------------------------------
!
! Complete rhs of helmholtz eq.
!
      m  = locm(lm,iam)
      lmr = lnstart(lm)
      lmc = 2*lmr
      do k=1,plev
!
! Coefficients for diagonal terms
!
         hhref = 0.5*href(k,k)
         hbps = 0.5*bps(k)
!
! Loop along total wavenumber index (in spectral space)
! Add lnps and diagonal (vertical space) T terms to d(t-1)
!
         do n=1,nlen(m)
            ir = lmc + 2*n - 1
            ii = ir + 1
            d(ir,k) = d(ir,k) + ztdtsq(n+m-1)*(hhref*t(ir,k) + hbps*alps(ir))
            d(ii,k) = d(ii,k) + ztdtsq(n+m-1)*(hhref*t(ii,k) + hbps*alps(ii))
         end do
         if (k.lt.plev) then
            do kk=k+1,plev
!
! Add off-diagonal (vertical space) T terms to d(t-1)
!
               hhref = 0.5*href(kk,k)
               do n=1,nlen(m)
                  ir = lmc + 2*n - 1
                  ii = ir + 1
                  d(ir,k) = d(ir,k) + ztdtsq(n+m-1)*hhref*t(ir,kk)
                  d(ii,k) = d(ii,k) + ztdtsq(n+m-1)*hhref*t(ii,kk)
               end do
            end do
         end if
      end do                    ! k=1,plev (calculation level)
!
! Solution of helmholtz equation
! First: initialize temporary space for solution
!
      do k=1,plev
         do j=1,2*pnmax
            z(j,k) = 0.
         end do
      end do
!
! Multiply right hand side by inverse matrix
!
      if (nlen(m) >= plev) then
         do kk=1,plev
            do k=1,plev
               do n=1,nlen(m)
                  ir = lmc + 2*n - 1
                  ii = ir + 1
                  z(2*n-1,k) = z(2*n-1,k) + bm1(kk,k,m+n-1)*d(ir,kk)
                  z(2*n  ,k) = z(2*n  ,k) + bm1(kk,k,m+n-1)*d(ii,kk)
               end do
            end do                  ! inner loop over levels
         end do                    ! outer loop over levels
      else
         do n=1,nlen(m)
            ir = lmc + 2*n - 1
            ii = ir + 1
            do kk=1,plev
               do k=1,plev
                  z(2*n-1,k) = z(2*n-1,k) + bm1(kk,k,m+n-1)*d(ir,kk)
                  z(2*n  ,k) = z(2*n  ,k) + bm1(kk,k,m+n-1)*d(ii,kk)
               end do
            end do                  ! inner loop over levels
         end do                    ! outer loop over levels
      endif
!
! Move solution for divergence to d
!
      do k=1,plev
         do n=1,nlen(m)
            ir = lmc + 2*n - 1
            ii = ir + 1
            d(ir,k) = z(2*n-1,k)
            d(ii,k) = z(2*n  ,k)
         end do
      end do
!
! Complete ln(pstar) and T forecasts
! Add semi-implicit part to surface pressure (vector multiply)
!
      do k=1,plev
         ztemp = zdt*hypd(k)/hypi(plevp)
         do n=1,nlen(m)
            ir = lmc + 2*n - 1
            ii = ir + 1
            alps(ir) = alps(ir) - ztemp*d(ir,k)
            alps(ii) = alps(ii) - ztemp*d(ii,k)
         end do
      end do
!
! Add semi-implicit part to temperature (matrix multiply)
!
      if (2*nlen(m) >= plev) then
         do kk=1,plev
            do k=1,plev
               do nn = lmc+1, lmc+2*nlen(m)
                  t(nn,k) = t(nn,k) - zdt*tau(kk,k)*d(nn,kk)
               end do
            end do
         end do
      else
         do n=1,nlen(m)
            ir = lmc + 2*n - 1
            ii = ir + 1
            do kk=1,plev
               do k=1,plev
                  t(ir,k) = t(ir,k) - zdt*tau(kk,k)*d(ir,kk)
                  t(ii,k) = t(ii,k) - zdt*tau(kk,k)*d(ii,kk)
               end do
            end do
         end do
      endif
!
      return
   end subroutine tstep

