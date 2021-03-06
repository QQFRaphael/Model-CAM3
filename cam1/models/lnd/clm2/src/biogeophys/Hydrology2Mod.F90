#include <misc.h>
#include <preproc.h>

module Hydrology2Mod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: Hydrology2Mod
!
! !DESCRIPTION:
! Calculation of soil/snow hydrology.
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: Hydrology2        ! Calcultes soil/snow hydrology
!
! !REVISION HISTORY:
! 2/28/02 Peter Thornton: Migrated to new data structures.
! 7/12/03 Forrest Hoffman ,Mariana Vertenstein : Migrated to vector code
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Hydrology2
!
! !INTERFACE:
  subroutine Hydrology2(lbc, ubc, num_nolakec, filter_nolakec, &
                        num_soilc, filter_soilc, num_snowc, filter_snowc, &
                        num_nosnowc, filter_nosnowc)
!
! !DESCRIPTION:
! This is the main subroutine to execute the calculation of soil/snow
! hydrology
! Calling sequence is:
!  Hydrology2:                 surface hydrology driver
!    -> SnowWater:             change of snow mass and snow water onto soil
!    -> SurfaceRunoff:         surface runoff
!    -> Infiltration:          infiltration into surface soil layer
!    -> SoilWater:             soil water movement between layers
!          -> Tridiagonal      tridiagonal matrix solution
!    -> Drainage:              subsurface runoff
!    -> SnowCompaction:        compaction of snow layers
!    -> CombineSnowLayers:     combine snow layers that are thinner than minimum
!    -> DivideSnowLayers:      subdivide snow layers that are thicker than maximum
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use clm_varcon      , only : denh2o, denice, istice, istwet, istsoil, spval
    use clm_varpar      , only : nlevsoi, nlevsno
    use SnowHydrologyMod, only : SnowCompaction, CombineSnowLayers, DivideSnowLayers, &
                                 SnowWater, BuildSnowFilter
    use SoilHydrologyMod, only : Infiltration, SoilWater, Drainage, SurfaceRunoff
#if (defined COUP_CAM)
    use time_manager    , only : get_step_size, get_nstep, is_perpetual
#else
    use time_manager    , only : get_step_size, get_nstep
#endif
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc                    ! column bounds
    integer, intent(in) :: num_nolakec                 ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(ubc-lbc+1)   ! column filter for non-lake points
    integer, intent(in) :: num_soilc                   ! number of column soil points in column filter
    integer, intent(in) :: filter_soilc(ubc-lbc+1)     ! column filter for soil points
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: cgridcell(:)     ! column's gridcell
    integer , pointer :: clandunit(:)     ! column's landunit
    integer , pointer :: ityplun(:)       ! landunit type
    integer , pointer :: snl(:)           ! number of snow layers
    real(r8), pointer :: h2ocan(:)        ! canopy water (mm H2O)
    real(r8), pointer :: h2osno(:)        ! snow water (mm H2O)
    real(r8), pointer :: watsat(:,:)      ! volumetric soil water at saturation (porosity)
    real(r8), pointer :: sucsat(:,:)      ! minimum soil suction (mm)
    real(r8), pointer :: bsw(:,:)         ! Clapp and Hornberger "b"
    real(r8), pointer :: z(:,:)           ! layer depth  (m)
    real(r8), pointer :: forc_rain(:)     ! rain rate [mm/s]
    real(r8), pointer :: forc_snow(:)     ! snow rate [mm/s]
    real(r8), pointer :: begwb(:)         ! water mass begining of the time step
    real(r8), pointer :: qflx_evap_tot(:) ! qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
!
! local pointers to implicit inout arguments
!
    real(r8), pointer :: dz(:,:)          ! layer thickness depth (m)
    real(r8), pointer :: zi(:,:)          ! interface depth (m)
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: endwb(:)         ! water mass end of the time step
    real(r8), pointer :: snowage(:)       ! non dimensional snow age [-]
    real(r8), pointer :: wf(:)            ! soil water as frac. of whc for top 0.5 m
    real(r8), pointer :: snowice(:)       ! average snow ice lens
    real(r8), pointer :: snowliq(:)       ! average snow liquid water
    real(r8), pointer :: t_snow(:)        ! vertically averaged snow temperature
    real(r8), pointer :: t_grnd(:)        ! ground temperature (Kelvin)
    real(r8), pointer :: t_soisno(:,:)    ! soil temperature (Kelvin)
    real(r8), pointer :: h2osoi_ice(:,:)  ! ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)  ! liquid water (kg/m2)
    real(r8), pointer :: h2osoi_vol(:,:)  ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
    real(r8), pointer :: qflx_drain(:)    ! sub-surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_surf(:)     ! surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_infl(:)     ! infiltration (mm H2O /s)
    real(r8), pointer :: qflx_qrgwl(:)    ! qflx_surf at glaciers, wetlands, lakes
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: g,l,c,j,fc                 ! indices
    integer  :: num_snowc                  ! number of column snow points
    integer  :: filter_snowc(ubc-lbc+1)    ! column filter for snow points
    integer  :: num_nosnowc                ! number of column non-snow points
    integer  :: filter_nosnowc(ubc-lbc+1)  ! column filter for non-snow points
    integer  :: nstep                      ! time step number
    real(r8) :: dtime                      ! land model time step (sec)
    real(r8) :: zwice(lbc:ubc)             ! the sum of ice mass of soil (kg/m2)
    real(r8) :: vol_liq(lbc:ubc,1:nlevsoi) ! partial volume of liquid water in layer
    real(r8) :: s(lbc:ubc,1:nlevsoi)       ! wetness of soil (including ice)
    real(r8) :: zwt(lbc:ubc)               ! water table depth
    real(r8) :: fcov(lbc:ubc)              ! fractional area with water table at surface
    real(r8) :: dwat(lbc:ubc,1:nlevsoi)    ! change in soil water
    real(r8) :: hk(lbc:ubc,1:nlevsoi)      ! hydraulic conductivity (mm h2o/s)
    real(r8) :: dhkdw(lbc:ubc,1:nlevsoi)   ! d(hk)/d(vol_liq)
#if (defined DGVM)
    real(r8) :: watdry                     ! temporary
    real(r8) :: rwat(lbc:ubc)              ! soil water wgted by depth to maximum depth of 0.5 m
    real(r8) :: swat(lbc:ubc)              ! same as rwat but at saturation
    real(r8) :: rz(lbc:ubc)                ! thickness of soil layers contributing to rwat (m)
    real(r8) :: tsw                        ! volumetric soil water to 0.5 m
    real(r8) :: stsw                       ! volumetric soil water to 0.5 m at saturation
#endif
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (gridcell-level)

    forc_rain => clm3%g%a2lf%forc_rain
    forc_snow => clm3%g%a2lf%forc_snow

    ! Assign local pointers to derived subtypes components (landunit-level)

    ityplun => clm3%g%l%itype

    ! Assign local pointers to derived subtypes components (column-level)

    cgridcell     => clm3%g%l%c%gridcell
    clandunit     => clm3%g%l%c%landunit
    snl           => clm3%g%l%c%cps%snl
    snowage       => clm3%g%l%c%cps%snowage
    t_snow        => clm3%g%l%c%ces%t_snow
    t_grnd        => clm3%g%l%c%ces%t_grnd
    h2ocan        => clm3%g%l%c%cws%pws_a%h2ocan
    h2osno        => clm3%g%l%c%cws%h2osno
    wf            => clm3%g%l%c%cps%wf
    snowice       => clm3%g%l%c%cws%snowice
    snowliq       => clm3%g%l%c%cws%snowliq
    watsat        => clm3%g%l%c%cps%watsat
    sucsat        => clm3%g%l%c%cps%sucsat
    bsw           => clm3%g%l%c%cps%bsw
    z             => clm3%g%l%c%cps%z
    dz            => clm3%g%l%c%cps%dz
    zi            => clm3%g%l%c%cps%zi
    t_soisno      => clm3%g%l%c%ces%t_soisno
    h2osoi_ice    => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq    => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_vol    => clm3%g%l%c%cws%h2osoi_vol
    qflx_evap_tot => clm3%g%l%c%cwf%pwf_a%qflx_evap_tot
    qflx_drain    => clm3%g%l%c%cwf%qflx_drain
    qflx_surf     => clm3%g%l%c%cwf%qflx_surf
    qflx_infl     => clm3%g%l%c%cwf%qflx_infl
    qflx_qrgwl    => clm3%g%l%c%cwf%qflx_qrgwl
    endwb         => clm3%g%l%c%cwbal%endwb
    begwb         => clm3%g%l%c%cwbal%begwb

    ! Determine time step and step size

    nstep = get_nstep()
    dtime = get_step_size()

    ! Determine initial snow/no-snow filters (will be modified possibly by
    ! routines CombineSnowLayers and DivideSnowLayers below

    call BuildSnowFilter(lbc, ubc, num_nolakec, filter_nolakec, &
         num_snowc, filter_snowc, num_nosnowc, filter_nosnowc)

    ! Determine the change of snow mass and the snow water onto soil

    call SnowWater(lbc, ubc, num_snowc, filter_snowc, num_nosnowc, filter_nosnowc)

    ! Determine soil hydrology

    call SurfaceRunoff(lbc, ubc, num_soilc, filter_soilc, &
         zwice, vol_liq, s, zwt, fcov)

    call Infiltration(lbc, ubc,  num_soilc, filter_soilc)

    call SoilWater(lbc, ubc, num_soilc, filter_soilc, &
         vol_liq, dwat, hk, dhkdw)

    call Drainage(lbc, ubc, num_soilc, filter_soilc, &
         zwice, vol_liq, s, zwt, fcov, hk, dhkdw, dwat)

#if (!defined COUP_CAM)

    ! Natural compaction and metamorphosis.

    call SnowCompaction(lbc, ubc, num_snowc, filter_snowc)

    ! Combine thin snow elements

    call CombineSnowLayers(lbc, ubc, num_snowc, filter_snowc)

    ! Divide thick snow elements

    call DivideSnowLayers(lbc, ubc, num_snowc, filter_snowc)

#else

    if (.not. is_perpetual()) then

       ! Natural compaction and metamorphosis.

       call SnowCompaction(lbc, ubc, num_snowc, filter_snowc)

       ! Combine thin snow elements

       call CombineSnowLayers(lbc, ubc, num_snowc, filter_snowc)

       ! Divide thick snow elements

       call DivideSnowLayers(lbc, ubc, num_snowc, filter_snowc)

    else

       do fc = 1, num_snowc
          c = filter_snowc(fc)
          h2osno(c) = 0._r8
       end do
       do j = -nlevsno+1,0
          do fc = 1, num_snowc
             c = filter_snowc(fc)
             if (j >= snl(c)+1) then
                h2osno(c) = h2osno(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
             end if
          end do
       end do

    end if

#endif

    ! Set empty snow layers to zero

!dir$ concurrent
!cdir nodep
    do fc = 1, num_snowc
       c = filter_snowc(fc)
       if (snl(c) > -nlevsno) then
          snowage(c) = 0._r8
       end if
    end do
    do j = -nlevsno+1,0
!dir$ concurrent
!cdir nodep
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j <= snl(c) .and. snl(c) > -nlevsno) then
             h2osoi_ice(c,j) = 0._r8
             h2osoi_liq(c,j) = 0._r8
             t_soisno(c,j) = 0._r8
             dz(c,j) = 0._r8
             z(c,j) = 0._r8
             zi(c,j-1) = 0._r8
          end if
       end do
    end do

    ! Build new snow filter

    call BuildSnowFilter(lbc, ubc, num_nolakec, filter_nolakec, &
         num_snowc, filter_snowc, num_nosnowc, filter_nosnowc)

    ! Vertically average t_soisno and sum of h2osoi_liq and h2osoi_ice
    ! over all snow layers for history output

!dir$ concurrent
!cdir nodep
    do fc = 1, num_snowc
       c = filter_snowc(fc)
       t_snow(c)  = 0._r8
       snowice(c) = 0._r8
       snowliq(c) = 0._r8
    end do
!dir$ concurrent
!cdir nodep
    do fc = 1, num_nosnowc
       c = filter_nosnowc(fc)
       t_snow(c)  = spval
       snowice(c) = spval
       snowliq(c) = spval
    end do

    do j = -nlevsno+1, 0
!dir$ concurrent
!cdir nodep
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             t_snow(c)  = t_snow(c) + t_soisno(c,j)
             snowice(c) = snowice(c) + h2osoi_ice(c,j)
             snowliq(c) = snowliq(c) + h2osoi_liq(c,j)
          end if
       end do
    end do

    ! Determine ground temperature, ending water balance and volumetric soil water

!dir$ concurrent
!cdir nodep
    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       if (snl(c) < 0) t_snow(c) = t_snow(c)/abs(snl(c))
       t_grnd(c) = t_soisno(c,snl(c)+1)
       endwb(c) = h2ocan(c) + h2osno(c)
    end do

    do j = 1, nlevsoi
!dir$ concurrent
!cdir nodep
       do fc = 1, num_nolakec
          c = filter_nolakec(fc)
          endwb(c) = endwb(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
          h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)
       end do
    end do

    ! Determine wetland and land ice hydrology (must be placed here
    ! since need snow updated from CombineSnowLayers)

!dir$ concurrent
!cdir nodep
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       l = clandunit(c)
       g = cgridcell(c)
       if (ityplun(l)==istwet .or. ityplun(l)==istice) then
          qflx_drain(c) = 0._r8
          qflx_surf(c) = 0._r8
          qflx_infl(c) = 0._r8
          qflx_qrgwl(c) = forc_rain(g) + forc_snow(g) - qflx_evap_tot(c) - (endwb(c)-begwb(c))/dtime
       end if
    end do

#if (defined DGVM)

    ! Available soil water up to a depth of 0.5 m.
    ! Potentially available soil water (=whc) up to a depth of 0.5 m.
    ! Water content as fraction of whc up to a depth of 0.5 m.

!dir$ concurrent
!cdir nodep
    do c = lbc,ubc
       l = clandunit(c)
       if (ityplun(l) == istsoil) then
          rwat(c) = 0._r8
          swat(c) = 0._r8
          rz(c)   = 0._r8
       end if
    end do

    do j = 1, nlevsoi
!dir$ concurrent
!cdir nodep
       do c = lbc,ubc
          l = clandunit(c)
          if (ityplun(l) == istsoil) then
             if (z(c,j)+0.5*dz(c,j) <= 0.5) then
                watdry = watsat(c,j) * (316230./sucsat(c,j)) ** (-1./bsw(c,j))
                rwat(c) = rwat(c) + (h2osoi_vol(c,j)-watdry) * dz(c,j)
                swat(c) = swat(c) + (watsat(c,j)    -watdry) * dz(c,j)
                rz(c) = rz(c) + dz(c,j)
             end if
          end if
       end do
    end do

!dir$ concurrent
!cdir nodep
    do c = lbc,ubc
       l = clandunit(c)
       if (ityplun(l) == istsoil) then
          if (rz(c) /= 0._r8) then
             tsw  = rwat(c)/rz(c)
             stsw = swat(c)/rz(c)
          else
             watdry = watsat(c,1) * (316230./sucsat(c,1)) ** (-1./bsw(c,1))
             tsw = h2osoi_vol(c,1) - watdry
             stsw = watsat(c,1) - watdry
          end if
          wf(c) = tsw/stsw
       else
          wf(c) = 1.0
       end if
    end do

#endif

  end subroutine Hydrology2

end module Hydrology2Mod
