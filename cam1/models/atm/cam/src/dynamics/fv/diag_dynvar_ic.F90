#include <misc.h>
#include <params.h>

  subroutine diag_dynvar_ic(phis, ps, t3, u3s, v3s, q3)
!
!----------------------------------------------------------------------- 
! 
! Purpose: record state variables to IC file
!
!-----------------------------------------------------------------------
!
    use shr_kind_mod , only: r8 => shr_kind_r8
    use pmgrid
    use history      , only: outfld, write_inithist
    use constituents , only: ppcnst, cnst_name
    use dynamics_vars, only: ng_d, ng_s

    implicit none
!
!-----------------------------------------------------------------------
!
! Arguments
!
    real(r8), intent(in) :: phis(plon,beglat:endlat)  ! surface geopotential (grav*zs)
    real(r8), intent(in) :: ps  (plon,beglat:endlat)  ! Surface pressure (pa) 
    real(r8), intent(in) :: t3  (plon,beglev:endlev,beglat:endlat)  ! Temperature (K)
    real(r8), intent(in) :: u3s (plon,beglat-ng_d:endlat+ng_s,beglev:endlev)  ! u wind velocities, staggered grid
    real(r8), intent(in) :: v3s (plon,beglat-ng_s:endlat+ng_d,beglev:endlev)  ! v wind velocities, staggered grid
    real(r8), intent(in) :: q3  (plon,beglat-ng_d:endlat+ng_d,beglev:endlev,ppcnst)  ! Tracers
!
!---------------------------Local workspace-----------------------------
!
    integer i, j, k, m   ! indices
    real(r8) tmp(plon,beglev:endlev)
!
!-----------------------------------------------------------------------
!
    if( write_inithist() ) then

!$OMP PARALLEL DO PRIVATE (I, J, K, M, TMP)
       do j = beglat, endlat

          call outfld ('PS&IC      ', ps  (1,j) , plon, j)
          call outfld ('T&IC       ', t3(1,1,j) , plon, j)

          do k = beglev, endlev
             do i = 1, plon
                tmp(i,k) = u3s(i,j,k)
             enddo
          enddo
          call outfld ('US&IC      ', tmp       , plon, j)

          do k = beglev, endlev
             do i = 1, plon
                tmp(i,k) = v3s(i,j,k)
             enddo
          enddo
          call outfld ('VS&IC      ', tmp       , plon, j)

          do m = 1, ppcnst
             do k = beglev, endlev
                do i = 1, plon
                   tmp(i,k) = q3(i,j,k,m)
                enddo
             enddo
             call outfld(trim(cnst_name(m))//'&IC' , tmp  , plon, j)
          end do

       enddo

    end if

    return
  end subroutine diag_dynvar_ic
