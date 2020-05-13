#include <misc.h>

module pfixer

!----------------------------------------------------------------------- 
!
! BOP
!
! !MODULE: metdata
!
! !DESCRIPTION
! Corrects (or fixes) mass fluxes and edge pressures to be consistent
! with offline surface pressure at next model time step.
!
! !USES
  use shr_kind_mod, only: r8 => shr_kind_r8
  use abortutils, only: endrun

! !PUBLIC MEMBER FUNCTIONS
  public :: adjust_press

! !REVISION HISTORY:
!   04.01.30  Stacy Walters    Creation
!   04.02.15  F Vitt  Fixed bug in edge pressure corrections
!   04.08.27  F Vitt  Added ability to handle 2D decomposition
!
! EOP
!----------------------------------------------------------------------- 

  private

contains

!-----------------------------------------------------------------------
!       ... adjust mass fluxes and pressures for lin-rood transport
!-----------------------------------------------------------------------

  subroutine adjust_press( im, jm, km, jfirst, jlast, kfirst, klast, &
       pmet1, pmet2, mfx, mfy, pexy, ifirstxy, ilastxy, jfirstxy, jlastxy )
!!$       pmet1, pmet2, padj, mfx, mfy, pexy, ifirstxy, ilastxy, jfirstxy, jlastxy )

#if defined( SPMD )
    use pmgrid,   only : iam, npr_y
    use mod_comm, only : mp_send3d, mp_recv3d
    use mod_comm,     only : mp_sendirr, mp_recvirr
    use spmd_dyn, only : ijk_yz_to_xy, ijk_xy_to_yz, comm_z
    use parutilitiesmodule, only: parcollective2d, SUMOP
#endif
    use time_manager, only : get_nstep
    use pmgrid,       only : twod_decomp

    use history,        only: outfld

    implicit none

    !-----------------------------------------------------------------------
    !       ... dummy arguments
    !-----------------------------------------------------------------------
    integer, intent(in)     :: im                                    ! east-west dimension
    integer, intent(in)     :: jm                                    ! south-north dimension
    integer, intent(in)     :: km                                    ! vertical dimension
    integer, intent(in)     :: jfirst                                ! first latitude
    integer, intent(in)     :: jlast                                 ! last latitude
    integer, intent(in)     :: ifirstxy                              ! first longitude
    integer, intent(in)     :: ilastxy                               ! last longitude
    integer, intent(in)     :: jfirstxy                              ! first latitude
    integer, intent(in)     :: jlastxy                               ! last latitude
    integer, intent(in)     :: kfirst                                ! first level
    integer, intent(in)     :: klast                                 ! last level
    real(r8), intent(in)    :: pmet1(im,jfirst:jlast)                ! surface pressure at t(n) (Pa)
    real(r8), intent(in)    :: pmet2(im,jfirst:jlast)                ! surface pressure at t(n+1) from met fields (Pa)
!!$    real(r8), intent(inout) :: padj(im,jfirst:jlast)                 ! adjusted surface pressure at t(n+1) (Pa)
    real(r8), intent(inout) :: mfx(im,jfirst:jlast,kfirst:klast)     ! zonal mass flux
    real(r8), intent(inout) :: mfy(im,jfirst:jlast+1,kfirst:klast)   ! meridional mass flux
    real(r8), intent(inout) :: pexy(ifirstxy:ilastxy,km+1,jfirstxy:jlastxy)

    !-----------------------------------------------------------------------
    !       ... local variables
    !-----------------------------------------------------------------------
    integer, parameter :: nstep0 = -1

    integer  :: i, j, k, km1
    integer  :: nstep
#if defined( SPMD )
    integer  :: dest, src
#endif
    integer  :: ndx(2)
    real(r8) :: dpi(im,jfirst:jlast,kfirst:klast)
    real(r8) :: dpixy(ifirstxy:ilastxy,jfirstxy:jlastxy,km)
    real(r8) :: dpi_in(im,jfirst:jlast,kfirst:klast)
    real(r8) :: dpi_inxy(ifirstxy:ilastxy,jfirstxy:jlastxy,1:km)
    real(r8) :: dpi_c(im,jfirst:jlast,kfirst:klast)
    real(r8) :: dps(im,jfirst:jlast)
    real(r8) :: dps_in(im,jfirst:jlast)
    real(r8) :: dps_inxy(ifirstxy:ilastxy,jfirstxy:jlastxy)
    real(r8) :: ps_diffxy(ifirstxy:ilastxy,jfirstxy:jlastxy)

    real(r8) :: dmfx(im,jfirst:jlast,kfirst:klast)     ! original zonnal mass flux
    real(r8) :: dmfy(im,jfirst:jlast+1,kfirst:klast)   ! original meridional mass flux 
    real(r8) :: emfx(im,jfirst:jlast)     ! zonal mass flux error
    real(r8) :: emfy(im,jfirst:jlast+1)   ! meridional mass flux error
 
    real(r8) :: tmp2d(im,kfirst:klast)
    logical :: debug = .false.
    logical :: method1 = .true.

#if defined( SPMD )
    ! Send one latitude of mfy to the south
    if( mod(iam,npr_y) /= 0 ) then
       dest = iam-1
    else
       dest = -1
    end if
    if( mod(iam+1,npr_y) /= 0 ) then
       src  = iam+1
    else
       src = -1
    end if
    call mp_send3d( dest, src, im, jm, km, &
         1, im, jfirst, jlast+1, kfirst, klast, &
         1, im, jfirst, jfirst, kfirst, klast, mfy)
#endif

    do j = jfirst,jlast
       dps(:,j) = pmet2(:,j) - pmet1(:,j)
    end do

#if defined( SPMD )
    call mp_recv3d( src, im, jm, km, &
         1, im, jfirst, jlast+1, kfirst, klast, &
         1, im, jlast+1, jlast+1, kfirst, klast, mfy)
#endif

    nstep = get_nstep()
!-----------------------------------------------------------------------
!       ... store incoming mass fluxes
!-----------------------------------------------------------------------
    if (debug) then
       do k = kfirst,klast
          do j = jfirst,jlast
             dmfx(:,j,k) = mfx(:,j,k)
             dmfy(:,j,k) = mfy(:,j,k)
          end do
          if( jlast /= jm ) then
             dmfy(:,jlast+1,k) = mfy(:,jlast+1,k)
          end if
       end do
    endif

!-----------------------------------------------------------------------
!       ... incoming mass flux divergence
!-----------------------------------------------------------------------
    call calc_divergence( im, jm, km, jfirst, jlast, &
         kfirst, klast, mfx, mfy, dpi_in )

    if (debug) then
       do j = jfirst, jlast
          do k = kfirst, klast
             do i = 1, im
                tmp2d(i,k) = dpi_in(i,j,k)
             enddo
          enddo
          call outfld('DPI_IN  ', tmp2d, im   ,j   )
       enddo
    endif

!-----------------------------------------------------------------------
!       ... surface pressure from mass flux divergence
!-----------------------------------------------------------------------
!  Two different methods to compute change in ps give differnt 
!  results if 2D decomp is used (round off error).  Method 1 gives
!  identical 1D vs 2D decomposition results.
!-----------------------------------------------------------------------
    if (method1) then 

       ! xfer dpi_in to dpi_inxy
       if (twod_decomp .eq. 1) then
#if defined (SPMD)
          call mp_sendirr( dpi_in,   ijk_yz_to_xy%SendDesc, ijk_yz_to_xy%RecvDesc, dpi_inxy)
          call mp_recvirr( dpi_inxy, ijk_yz_to_xy%RecvDesc )
#endif
       else            
          dpi_inxy(:,:,:) = dpi_in(:,:,:)
       endif

       ! vertical sum
       do j = jfirstxy,jlastxy
          do i = ifirstxy,ilastxy
             dps_inxy(i,j) = sum( dpi_inxy(i,j,1:km) )
          end do
       end do

       ! xfer dps_inxy to dps_in
       ! Embed in 3D array since transpose machinery cannot handle 2D arrays
       if (twod_decomp .eq. 1) then
#if defined (SPMD)
          do k = 1,km
             do j = jfirstxy,jlastxy
                do i = ifirstxy,ilastxy
                   dpixy(i,j,k) = dps_inxy(i,j)
                enddo
             enddo
          enddo

          call mp_sendirr( dpixy, ijk_xy_to_yz%SendDesc, ijk_xy_to_yz%RecvDesc, dpi )
          call mp_recvirr( dpi,   ijk_xy_to_yz%RecvDesc )

          do j = jfirst,jlast
             do i = 1,im
                dps_in(i,j) = dpi(i,j,kfirst)
             enddo
          enddo
#endif
       else            
          dps_in(:,:) = dps_inxy(:,:)
       endif

    else ! method1

       ! this method does not give identical results as the above method
       ! when two dimensional decomposition is used

       do j = jfirst,jlast
          do i = 1,im
             dps_in(i,j) = sum( dpi_in(i,j,kfirst:klast) )
          end do
       end do

#if ( defined SPMD )
       if (twod_decomp .eq. 1) then
          call parcollective2d( comm_z, SUMOP, im, jlast-jfirst+1,  dps_in )
       endif
#endif

    endif ! method1

    if (debug) then
       do j = jfirst, jlast
          call outfld('DPS1    ',   dps(1,j), im   ,j   )
          call outfld('DPS2    ',dps_in(1,j), im   ,j   )
       enddo
    endif

!-----------------------------------------------------------------------
!       ... modify (fix) mass fluxes
!-----------------------------------------------------------------------
    call do_press_fix_llnl( im, jm, km, jfirst, jlast, &
         kfirst, klast, dps, dps_in, mfx, mfy )

!-----------------------------------------------------------------------
!       ... modified mass flux divergence
!-----------------------------------------------------------------------
    call calc_divergence( im, jm, km, jfirst, jlast, &
         kfirst, klast, mfx, mfy, dpi_c )
      
!-----------------------------------------------------------------------
!       ... differential mass flux divergence
!-----------------------------------------------------------------------
    do k = kfirst,klast
       do j = jfirst,jlast
          dpi(:,j,k) = dpi_c(:,j,k) - dpi_in(:,j,k)
       end do
    end do

    if (twod_decomp .eq. 1) then
#if defined (SPMD)
       call mp_sendirr( dpi,   ijk_yz_to_xy%SendDesc,ijk_yz_to_xy%RecvDesc, dpixy)
       call mp_recvirr( dpixy, ijk_yz_to_xy%RecvDesc )
#endif
    else            
       dpixy(:,:,:) = dpi(:,:,:)
    endif


!-----------------------------------------------------------------------
!       ... modify pe
!-----------------------------------------------------------------------

    if (debug) then
       write(*,*) ' '
       write(*,*) 'adjust_press: max pe diff %  @ nstep,ifirstxy,ilastxy,jfirstxy,jlastxy = ',&
            nstep,ifirstxy,ilastxy,jfirstxy,jlastxy
    endif

    do k = 1+1,km+1
       km1 = k - 1

       if (debug) then
          do j = jfirstxy,jlastxy
             do i = ifirstxy,ilastxy
                ps_diffxy(i,j) = sum( dpixy(i,j,1:km1) )/ pexy(i,k,j )
             end do
          end do
       endif

       if( nstep > nstep0 ) then
          do j = jfirstxy,jlastxy
             do i = ifirstxy,ilastxy
                pexy(i,k,j) = pexy(i,k,j) + sum( dpixy(i,j,1:km1) ) 
             end do
          end do
       end if
       if (debug) then

          ndx(:)       = maxloc( abs( ps_diffxy(:,:) ) )

          ndx(1)       = ndx(1) + ifirstxy - 1
          ndx(2)       = ndx(2) + jfirstxy - 1

          write(*,'("pfixer press change error (% error,press adjmnt,new pe)",1x,3i5,1p,3g15.7)') &
               k,ndx(:),100._r8*abs( ps_diffxy(ndx(1),ndx(2)) ), &
               dpixy(ndx(1),ndx(2),km1),pexy(ndx(1),k,ndx(2))

       endif
    end do

    if (debug) then 
       write(*,*) ' '
       write(*,*) 'adjust_press: max  mass flux error  @ nstep,jfirst,jlast,kfirst,klast = ',&
            nstep,jfirst,jlast,kfirst,klast

       do k = kfirst,klast

          do j=jfirst,jlast
             do i=1,im
                emfx(i,j) = ( mfx(i,j,k)-dmfx(i,j,k) ) 
             enddo
          enddo

          ndx(:)       = maxloc( abs( emfx(:,:) ) )
          ndx(2)       = ndx(2) + jfirst - 1

          write(*,'("pfixer max x flux error (diff,fixed,orig) ",1x,3i5,1p,3g15.7)') &
               k,ndx(:), emfx( ndx(1),ndx(2) ) , &
               mfx(ndx(1),ndx(2),k),  dmfx(ndx(1),ndx(2),k)

          do j=jfirst,jlast+1
             do i=1,im
                emfy(i,j) = ( mfy(i,j,k)-dmfy(i,j,k) ) 
             enddo
          enddo

          ndx(:)       = maxloc( abs( emfy(:,:) ) )
          ndx(2)       = ndx(2) + jfirst - 1

          write(*,'("pfixer max y flux error (diff,fixed,orig) ",1x,3i5,1p,3g15.7)') &
               k,ndx(:), emfy( ndx(1),ndx(2) ) , &
               mfy(ndx(1),ndx(2),k),  dmfy(ndx(1),ndx(2),k)

       enddo
    endif

  end subroutine adjust_press

!-----------------------------------------------------------------------
!       ... calculate horizontal mass flux divergence
!-----------------------------------------------------------------------
  subroutine calc_divergence( im, jm, km, jfirst, jlast, &
       kfirst, klast, mfx, mfy, dpi )

    use dynamics_vars, only : rcap, acosp

    implicit none

    !-----------------------------------------------------------------------
    !       ... dummy arguments
    !-----------------------------------------------------------------------
    integer, intent(in)     :: im                                         ! east-west dimension
    integer, intent(in)     :: jm                                         ! south-north dimension
    integer, intent(in)     :: km                                         ! vertical dimension
    integer, intent(in)     :: jfirst                                     ! first latitude
    integer, intent(in)     :: jlast                                      ! last latitude
    integer, intent(in)     :: kfirst                                     ! first level
    integer, intent(in)     :: klast                                      ! last level
    real(r8), intent(in)    :: mfx(im,jfirst:jlast,kfirst:klast)          ! zonal mass flux
    real(r8), intent(in)    :: mfy(im,jfirst:jlast+1,kfirst:klast)        ! meridional mass flux
    real(r8), intent(inout) :: dpi(im,jfirst:jlast,kfirst:klast)          ! horizontal mass flux divergence

!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
    integer  :: i, j, k, js2g0, jn2g0
    real(r8) :: sum1

    js2g0 = max( 2,jfirst )
    jn2g0 = min( jm-1,jlast )

#ifdef SPMD
!$omp parallel do private( j, k )
#endif
    do k = kfirst,klast
!-----------------------------------------------------------------------
!       ... north-south component
!-----------------------------------------------------------------------
       do j = js2g0,jn2g0
          dpi(:,j,k) = (mfy(:,j,k) - mfy(:,j+1,k)) * acosp(j)
       end do
!-----------------------------------------------------------------------
!       ... east-west component
!-----------------------------------------------------------------------
       do j = js2g0,jn2g0
          dpi(:im-1,j,k) = dpi(:im-1,j,k) + mfx(:im-1,j,k) - mfx(2:im,j,k)
          dpi(im,j,k)    = dpi(im,j,k) + mfx(im,j,k) - mfx(1,j,k)
       end do
!-----------------------------------------------------------------------
!       ... poles
!-----------------------------------------------------------------------
       if( jfirst == 1 ) then
          sum1 = -sum( mfy(:,2,k) )*rcap
          dpi(:,1,k) = sum1
       end if
       if( jlast == jm ) then
          sum1 = sum( mfy(:,jm,k) ) * rcap
          dpi(:,jm,k) = sum1
       end if
    end do
#ifdef SPMD
!$omp end parallel do
#endif

  end subroutine calc_divergence

!-----------------------------------------------------------------------
!       ... fix the mass fluxes to match the met field pressure tendency
!     See: http://asd.llnl.gov/pfix/index.html
!-----------------------------------------------------------------------
  subroutine do_press_fix_llnl( im, jm, km, jfirst, jlast, &
       kfirst, klast, dps, dps_ctm, mfx, mfy )

    use pmgrid,        only : plon, plat, plev, plevp, masterproc
    use dynamics_vars, only : cosp, rcap
    use commap,        only : gw => w
#ifdef SPMD
    use mpishorthand,  only : mpicom, mpi_double_precision, mpi_success
    use spmd_dyn,      only : npes, compute_gsfactors
#endif

    implicit none

!-----------------------------------------------------------------------
!       ... dummy arguments
!-----------------------------------------------------------------------
    integer, intent(in)     :: im                                         ! east-west dimension
    integer, intent(in)     :: jm                                         ! south-north dimension
    integer, intent(in)     :: km                                         ! vertical dimension
    integer, intent(in)     :: jfirst                                     ! first latitude
    integer, intent(in)     :: jlast                                      ! last latitude
    integer, intent(in)     :: kfirst                                     ! first level
    integer, intent(in)     :: klast                                      ! last level
    real(r8), intent(in)    :: dps(im,jfirst:jlast)                       ! surface pressure change from met fields
    real(r8), intent(in)    :: dps_ctm(im,jfirst:jlast)                   ! vert. sum of dpi from original mass fluxes
    real(r8), intent(inout) :: mfx(im,jfirst:jlast,kfirst:klast)          ! zonal mass flux
    real(r8), intent(inout) :: mfy(im,jfirst:jlast+1,kfirst:klast)        ! meridional mass flux

!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
    integer     :: i, j, jglob, k, astat, ierr
    integer     :: jn2g0, js2g0, jn2g1
    integer     :: cnt
#ifdef SPMD
    integer     :: numrecv(0:npes-1)
    integer     :: displs(0:npes-1)
#endif
    real(r8)    :: dpress_g                       ! global pressure error
    real(r8)    :: fxmean, factor
    real(r8)    :: ddps(plon,jfirst:jlast)        ! surface pressure change error
    real(r8)    :: dpresslat(plat)
    real(r8)    :: mmfd(plat)
    real(r8)    :: mmf(plat+1)
    real(r8)    :: fxintegral(plon+1)
    real(r8)    :: xcolmass_fix(plon,jfirst:jlast)

#include <comhyb.h>

    js2g0 = max( 2,jfirst )
    jn2g0 = min( jm-1,jlast )
    jn2g1 = min( jm-1,jlast+1 )

    do j = jfirst,jlast
       ddps(:,j) = dps(:,j) - dps_ctm(:,j)
    end do
    factor = .5_r8/im
    do j = jfirst,jlast
       dpresslat(j) = sum( ddps(:,j) ) * gw(j) * factor
    end do

#ifdef SPMD
    call compute_gsfactors( 1, cnt, numrecv, displs )
    cnt = jlast - jfirst + 1
    call mpi_allgatherv( dpresslat(jfirst:jlast), cnt, mpi_double_precision, &
         dpresslat, numrecv, displs, mpi_double_precision, mpicom, ierr )
    if( ierr /= mpi_success ) then
       write(*,*) 'do_press_fix_llnl: mpi_allgatherv failed; error code = ',ierr
       call endrun
    end if
#endif
    dpress_g = sum( dpresslat(:) )
    if( masterproc ) then
       write(*,*) 'do_press_fix_llnl: dpress_g = ',dpress_g
    end if

!-----------------------------------------------------------------------
!     calculate mean meridional flux divergence (df/dy).
!     note that mmfd is actually the zonal mean pressure change,
!     which is related to df/dy by geometrical factors.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     	... handle non-pole regions.
!-----------------------------------------------------------------------
    factor = 1./im
    do j = jfirst,jlast
       mmfd(j) = dpress_g - sum( ddps(:,j) ) * factor
    end do

#ifdef SPMD
    cnt = jlast - jfirst + 1
    call mpi_allgatherv( mmfd(jfirst:jlast), cnt, mpi_double_precision, &
         mmfd, numrecv, displs, mpi_double_precision, mpicom, ierr )
    if( ierr /= mpi_success ) then
       write(*,*) 'do_press_fix_llnl: mpi_allgatherv failed; error code = ',ierr
       call endrun
    end if
#endif

!-----------------------------------------------------------------------
!     calculate mean meridional fluxes (cosp*fy).
!     nb: this calculation is being done for global lats, i.e., (1,plat)
!-----------------------------------------------------------------------
    mmf(2) = mmfd(1) / (rcap*plon)
    do j = 2,plat-1
       mmf(j+1) = mmf(j) + mmfd(j) * cosp(j)
    end do

!-----------------------------------------------------------------------
!     fix latitude bands.
!     note that we do not need to worry about geometry here because
!     all boxes in a latitude band are identical.
!     note also that fxintegral(plon+1) should equal fxintegral(1),
!     i.e., zero.
!-----------------------------------------------------------------------
#ifdef SPMD
!$omp parallel do private( i, j, k, fxmean, fxintegral )
#endif
    do j = js2g0,jn2g0
       fxintegral(1) = 0._r8
       do i = 1,im
          fxintegral(i+1)  = fxintegral(i) - (ddps(i,j) - dpress_g) - mmfd(j)
       end do
       fxintegral(1)       = fxintegral(im+1)
       fxmean              = sum( fxintegral(:im) ) * factor
       xcolmass_fix(:im,j) = fxintegral(:im) - fxmean
    end do
#ifdef SPMD
!$omp end parallel do
#endif

!-----------------------------------------------------------------------
!     	... distribute colmass_fix in vertical
!-----------------------------------------------------------------------
#ifdef SPMD
!$omp parallel do private( j, k )
#endif
    do k = kfirst,klast
       do j = js2g0,jn2g0
          mfx(:,j,k) = mfx(:,j,k) + xcolmass_fix(:,j) * hybd(k)
       end do
       do j = js2g0,jn2g1
          mfy(:,j,k) = mfy(:,j,k) + mmf(j) * hybd(k)
       end do
       if( jlast == jm ) then
          mfy(:,jm,k) = mfy(:,jm,k) + mmf(jm) * hybd(k)
       end if
    end do
#ifdef SPMD
!$omp end parallel do
#endif

  end subroutine do_press_fix_llnl

end module pfixer
