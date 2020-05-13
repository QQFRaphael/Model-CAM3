#include <misc.h>
#include <params.h>

module diagnostics

!---------------------------------------------------------------------------------
! Module to compute a variety of diagnostics quantities for history files
!---------------------------------------------------------------------------------

use shr_kind_mod,  only: r8 => shr_kind_r8
use infnan,        only: nan
use physics_types, only: physics_state
use ppgrid,        only: pcols, pver, pvermx, begchunk, endchunk
use history,       only: outfld, write_inithist
use constituents,  only: pcnst, ppcnst, cnst_name, cnst_longname, cnst_cam_outfld
use chemistry,     only: chem_is
use abortutils,    only: endrun

implicit none
private
save

! Public interfaces

public :: &
   diag_defaultopts,   &! set default values of namelist variables
   diag_setopts,       &! get namelist input
   diag_init,          &! initialization
   diag_allocate,      &! allocate memory for module variables
   diag_deallocate,    &! deallocate memory for module variables
   diag_conv_tend_ini, &! initialize convective tendency calcs
   diag_phys_writeout, &! output diagnostics of the dynamics
   diag_conv,          &! output diagnostics of convective processes
   diag_surf,          &! output diagnostics of the surface
   diag_physvar_ic

! Private data

real(r8), allocatable :: &
   dtcond(:,:,:),  &! temperature tendency due to convection
   dqcond(:,:,:,:)  ! constituent tendencies due to convection

character(len=8) :: diag_cnst_conv_tend = 'q_only' ! output constituent tendencies due to convection
                                                   ! 'none', 'q_only' or 'all'

character(len=8), public :: dcconnam(ppcnst)       ! names of convection tendencies

contains

!===============================================================================

subroutine diag_defaultopts(diag_cnst_conv_tend_out)
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
!-----------------------------------------------------------------------

   character(len=*), intent(out), optional :: diag_cnst_conv_tend_out
!-----------------------------------------------------------------------
   if ( present(diag_cnst_conv_tend_out) ) then
      if (chem_is('waccm_ghg') .or. chem_is('waccm_mozart')) then
         diag_cnst_conv_tend_out = 'none'
      else
         diag_cnst_conv_tend_out = diag_cnst_conv_tend
      end if
   endif
end subroutine diag_defaultopts

!================================================================================================

subroutine diag_setopts(diag_cnst_conv_tend_in)
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
!-----------------------------------------------------------------------

   character(len=*), intent(in), optional :: diag_cnst_conv_tend_in
!-----------------------------------------------------------------------
   if ( present(diag_cnst_conv_tend_in) ) then
      diag_cnst_conv_tend = diag_cnst_conv_tend_in
   endif
end subroutine diag_setopts

!================================================================================================

subroutine diag_init

! Declare the history fields for which this module contains outfld calls.

   use history,        only: addfld, add_default, phys_decomp
   use comsrf,         only: plevmx, tsnam
   use constituents,   only:  cnst_need_pdeldry

   integer :: k, m

! outfld calls in diag_phys_writeout

   call addfld ('NSTEP   ','timestep',1,'A','Model timestep',phys_decomp)
   call addfld ('PHIS    ','m2/s2   ',1,    'I','Surface geopotential',phys_decomp)
   call addfld ('PS      ','Pa      ',1,    'A','Surface pressure',phys_decomp)
   call addfld ('T       ','K       ',pver, 'A','Temperature',phys_decomp)
   call addfld ('U       ','m/s     ',pver, 'A','Zonal wind',phys_decomp)
   call addfld ('V       ','m/s     ',pver, 'A','Meridional wind',phys_decomp)
   call addfld (cnst_name(1),'kg/kg ',pver, 'A',cnst_longname(1),phys_decomp)

   call addfld ('Z3      ','m       ',pver, 'A','Geopotential Height (above sea level)',phys_decomp)
   call addfld ('Z700    ','m       ',1,    'A','Geopotential Z at 700 mbar pressure surface',phys_decomp)
   call addfld ('Z500    ','m       ',1,    'A','Geopotential Z at 500 mbar pressure surface',phys_decomp)
   call addfld ('Z300    ','m       ',1,    'A','Geopotential Z at 300 mbar pressure surface',phys_decomp)
   call addfld ('Z050    ','m       ',1,    'A','Geopotential Z at 50 mbar pressure surface',phys_decomp)

   call addfld ('ZZ      ','m2      ',pver, 'A','Eddy height variance' ,phys_decomp)
   call addfld ('VZ      ','m2/s    ',pver, 'A','Meridional transport of geopotential energy',phys_decomp)
   call addfld ('VT      ','K m/s   ',pver, 'A','Meridional heat transport',phys_decomp)
   call addfld ('VU      ','m2/s2   ',pver, 'A','Meridional flux of zonal momentum' ,phys_decomp)
   call addfld ('VV      ','m2/s2   ',pver, 'A','Meridional velocity squared' ,phys_decomp)
   call addfld ('VQ      ','m/skg/kg',pver, 'A','Meridional water transport',phys_decomp)

   call addfld ('UU      ','m2/s2   ',pver, 'A','Zonal velocity squared' ,phys_decomp)
   call addfld ('WSPEED  ','m/s     ',pver, 'X','Horizontal total wind speed' ,phys_decomp)

   call addfld ('OMEGA   ','Pa/s    ',pver, 'A','Vertical velocity (pressure)',phys_decomp)
   call addfld ('OMEGAT  ','K Pa/s  ',pver, 'A','Vertical heat flux' ,phys_decomp)
   call addfld ('OMEGAU  ','m Pa/s2 ',pver, 'A','Vertical flux of zonal momentum' ,phys_decomp)
   call addfld ('OMEGA850','Pa/s    ',1,    'A','Vertical velocity at 850 mbar pressure surface',phys_decomp)
   call addfld ('OMEGA500','Pa/s    ',1,    'A','Vertical velocity at 500 mbar pressure surface',phys_decomp)

   call addfld ('MQ      ','kg/m2   ',pver, 'A','Water vapor mass in layer',phys_decomp)
   call addfld ('TMQ     ','kg/m2   ',1,    'A','Total (vertically integrated) precipitatable water',phys_decomp)
   call addfld ('RELHUM  ','percent ',pver, 'A','Relative humidity',phys_decomp)
   call addfld ('PSL     ','Pa      ',1,    'A','Sea level pressure',phys_decomp)

   call addfld ('T850    ','K       ',1,    'A','Temperature at 850 mbar pressure surface',phys_decomp)
   call addfld ('T300    ','K       ',1,    'A','Temperature at 300 mbar pressure surface',phys_decomp)
   call addfld ('Q850    ','kg/kg   ',1,    'A','Specific Humidity at 850 mbar pressure surface',phys_decomp)
   call addfld ('Q200    ','kg/kg   ',1,    'A','Specific Humidity at 700 mbar pressure surface',phys_decomp)
   call addfld ('U850    ','m/s     ',1,    'A','Zonal wind at 850 mbar pressure surface',phys_decomp)
   call addfld ('U200    ','m/s     ',1,    'A','Zonal wind at 200 mbar pressure surface',phys_decomp)
   call addfld ('V850    ','m/s     ',1,    'A','Meridional wind at 850 mbar pressure surface',phys_decomp)
   call addfld ('V200    ','m/s     ',1,    'A','Meridional wind at 200 mbar pressure surface',phys_decomp)

   call addfld ('TT      ','K2      ',pver, 'A','Eddy temperature variance' ,phys_decomp)

   call addfld ('UBOT    ','m/s     ',1,    'A','Lowest model level zonal wind',phys_decomp)
   call addfld ('VBOT    ','m/s     ',1,    'A','Lowest model level meridional wind',phys_decomp)
   call addfld ('QBOT    ','kg/kg   ',1,    'A','Lowest model level water vapor mixing ratio',phys_decomp)
   call addfld ('ZBOT    ','m       ',1,    'A','Lowest model level height', phys_decomp)

   ! defaults
   call add_default ('PHIS    ', 1, ' ')
   call add_default ('PS      ', 1, ' ')
   call add_default ('T       ', 1, ' ')
   call add_default ('U       ', 1, ' ')
   call add_default ('V       ', 1, ' ')
   call add_default (cnst_name(1), 1, ' ')
   call add_default ('Z3      ', 1, ' ')
   call add_default ('OMEGA   ', 1, ' ')
   call add_default ('RELHUM  ', 1, ' ')

   if ( cnst_need_pdeldry ) then
      call addfld ('PDELDRY ','Pa      ',pver, 'A','Dry pressure difference between levels',phys_decomp)
      call add_default ('PDELDRY ', 1, ' ')
      call addfld ('PSDRY   ','Pa      ',1,    'A','Surface pressure',phys_decomp)
      call add_default ('PSDRY   ', 1, ' ')
   endif

   if (chem_is('cam_default') .or. chem_is('cam_ghg')) then
      call add_default ('VT      ', 1, ' ')
      call add_default ('VU      ', 1, ' ')
      call add_default ('VV      ', 1, ' ')
      call add_default ('VQ      ', 1, ' ')
      call add_default ('UU      ', 1, ' ')
      call add_default ('OMEGAT  ', 1, ' ')
      call add_default ('TMQ     ', 1, ' ')
      call add_default ('PSL     ', 1, ' ')
   end if
   if (chem_is('waccm_ghg') .or. chem_is('waccm_mozart')) then
      call add_default ('PS      ', 2, ' ')
      call add_default ('T       ', 2, ' ')
   end if

! outfld calls in diag_conv

   call addfld ('DTCOND  ','K/s     ',pver, 'A','T tendency - moist processes',phys_decomp)
   do m=1,ppcnst
      dcconnam(m) = 'DC'//cnst_name(m)
      call addfld (dcconnam(m), 'kg/kg/s',pver,'A',trim(cnst_name(m))//' tendency due to moist processes',phys_decomp)
   end do

   ! defaults
   call add_default ('DTCOND  ', 1, ' ')
   if (diag_cnst_conv_tend == 'q_only' .or. diag_cnst_conv_tend == 'all') then
      call add_default (dcconnam(1),  1, ' ')
      if (diag_cnst_conv_tend == 'all') then
         do m = 2, ppcnst
            call add_default (dcconnam(m),  1, ' ')
         end do
      end if
   end if

! outfld calls in diag_surf

   call addfld ('SHFLX   ','W/m2    ',1,    'A','Surface sensible heat flux',phys_decomp)
   call addfld ('LHFLX   ','W/m2    ',1,    'A','Surface latent heat flux',phys_decomp)
   call addfld ('QFLX    ','kg/m2/s ',1,    'A','Surface water flux',phys_decomp)

   call addfld ('TAUX    ','N/m2    ',1,    'A','Zonal surface stress',phys_decomp)
   call addfld ('TAUY    ','N/m2    ',1,    'A','Meridional surface stress',phys_decomp)
   call addfld ('TREFHT  ','K       ',1,    'A','Reference height temperature',phys_decomp)
   call addfld ('TREFHTMN','K       ',1,    'M','Minimum reference height temperature over output period',phys_decomp)
   call addfld ('TREFHTMX','K       ',1,    'X','Maximum reference height temperature over output period',phys_decomp)
   call addfld ('QREFHT  ','kg/kg   ',1,    'A','Reference height humidity',phys_decomp)

   call addfld ('LANDFRAC','fraction',1,    'A','Fraction of sfc area covered by land',phys_decomp)
   call addfld ('ICEFRAC ','fraction',1,    'A','Fraction of sfc area covered by sea-ice',phys_decomp)
   call addfld ('OCNFRAC ','fraction',1,    'A','Fraction of sfc area covered by ocean',phys_decomp)

   call addfld ('TREFMNAV','K       ',1,    'A','Average of TREFHT daily minimum',phys_decomp)
   call addfld ('TREFMXAV','K       ',1,    'A','Average of TREFHT daily maximum',phys_decomp)

   do k=1,plevmx
      call addfld (tsnam(k),'K       ',1,    'A',tsnam(k)//' subsoil temperature',phys_decomp)
   end do
   call addfld ('SICTHK  ','m       ',1,    'A','Sea ice thickness',phys_decomp)
   call addfld ('TSICE   ','K       ',1,    'A','Ice temperature',phys_decomp)

   call addfld ('TS      ','K       ',1,    'A','Surface temperature (radiative)',phys_decomp)
   call addfld ('TSMN    ','K       ',1,    'M','Minimum surface temperature over output period',phys_decomp)
   call addfld ('TSMX    ','K       ',1,    'X','Maximum surface temperature over output period',phys_decomp)
   call addfld ('SNOWHLND','m       ',1,    'A','Water equivalent snow depth',phys_decomp)
   call addfld ('SNOWHICE','m       ',1,    'A','Water equivalent snow depth',phys_decomp)
   call addfld ('TBOT    ','K       ',1,    'A','Lowest model level temperature', phys_decomp)

   call addfld ('ASDIR',   '1',       1,    'A','albedo: shortwave, direct', phys_decomp)
   call addfld ('ASDIF',   '1',       1,    'A','albedo: shortwave, diffuse', phys_decomp)
   call addfld ('ALDIR',   '1',       1,    'A','albedo: longwave, direct', phys_decomp)
   call addfld ('ALDIF',   '1',       1,    'A','albedo: longwave, diffuse', phys_decomp)
   call addfld ('SST',     'K',       1,    'A','sea surface temperature', phys_decomp)

   ! defaults
   call add_default ('SHFLX   ', 1, ' ')
   call add_default ('LHFLX   ', 1, ' ')
   call add_default ('QFLX    ', 1, ' ')
   call add_default ('TAUX    ', 1, ' ')
   call add_default ('TAUY    ', 1, ' ')
   call add_default ('TREFHT  ', 1, ' ')
   call add_default ('LANDFRAC', 1, ' ')
   call add_default ('OCNFRAC ', 1, ' ')
#if ( defined COUP_SOM )
   call add_default ('QREFHT  ', 1, ' ')
#endif
   if (chem_is('cam_default') .or. chem_is('cam_ghg')) then
      call add_default ('ICEFRAC ', 1, ' ')
      call add_default ('TS      ', 1, ' ')
      call add_default ('TSMN    ', 1, ' ')
      call add_default ('TSMX    ', 1, ' ')
      call add_default ('SNOWHLND', 1, ' ')
      call add_default ('SNOWHICE', 1, ' ')
   end if

! outfld calls in diag_physvar_ic

   call addfld ('QCWAT&IC   ','kg/kg   ',pver, 'I','q associated with cloud water'                 ,phys_decomp)
   call addfld ('TCWAT&IC   ','kg/kg   ',pver, 'I','T associated with cloud water'                 ,phys_decomp)
   call addfld ('LCWAT&IC   ','kg/kg   ',pver, 'I','Cloud water (ice + liq'                        ,phys_decomp)
   call addfld ('CLOUD&IC   ','fraction',pver, 'I','Cloud fraction'                                ,phys_decomp)
   call addfld ('PBLH&IC    ','m       ',1,    'I','PBL height'                                    ,phys_decomp)
   call addfld ('TPERT&IC   ','K       ',1,    'I','Perturbation temperature (eddies in PBL)'      ,phys_decomp)
   call addfld ('QPERT&IC   ','kg/kg   ',1,    'I','Perturbation specific humidity (eddies in PBL)',phys_decomp)

   call addfld ('TSICERAD&IC','K       ',1,    'I','Radiatively equivalent ice temperature'        ,phys_decomp)
   call addfld ('TSICE&IC   ','K       ',1,    'I','Ice temperature'                               ,phys_decomp)
   call addfld ('SNOWHICE&IC','m       ',1,    'I','Water equivalent snow depth'                   ,phys_decomp)
   call addfld ('ICEFRAC&IC ','fraction',1,    'I','Fraction of sfc area covered by sea-ice'       ,phys_decomp)
   call addfld ('SICTHK&IC  ','m       ',1,    'I','Sea ice thickness'                             ,phys_decomp)
   call addfld ('TSOCN&IC   ','m       ',1,    'I','Ocean tempertare'                              ,phys_decomp)
   call addfld ('TS&IC      ','K       ',1,    'I','Surface temperature (radiative)'               ,phys_decomp)
   call addfld ('TBOT&IC    ','K       ',1,    'I','Lowest model level temperature'                ,phys_decomp)
   do k = 1,plevmx
      call addfld (trim(tsnam(k))//'&IC','K       ',1,   'I',tsnam(k)//' subsoil temperature'      ,phys_decomp)
   end do

   ! defaults (can't move the add_defaults calls for &IC fields here because mtapes
   !           isn't known yet)


end subroutine diag_init

!===============================================================================

subroutine diag_allocate()

! Allocate memory for module variables.

! Local variables
   character(len=*), parameter :: sub = 'diag_allocate'
   integer :: istat

   allocate(dtcond(pcols,pver,begchunk:endchunk),        &
            dqcond(pcols,pver,ppcnst,begchunk:endchunk), &
            stat=istat)
   if ( istat /= 0 ) then
      call endrun (sub//': ERROR: allocate failed')
   end if
   dtcond = nan
   dqcond = nan

end subroutine diag_allocate

!===============================================================================

subroutine diag_deallocate()

! Deallocate memory for module variables.

! Local variables
   character(len=*), parameter :: sub = 'diag_deallocate'
   integer :: istat

   deallocate(dtcond, dqcond, stat=istat)
   if ( istat /= 0 ) then
      call endrun (sub//': ERROR: deallocate failed')
   end if
end subroutine diag_deallocate
!===============================================================================

subroutine diag_conv_tend_ini(state)

! Initialize convective tendency calcs.

! Argument:

   type(physics_state), intent(in) :: state

! Local variables:

   integer :: i, k, m, lchnk, ncol

   lchnk = state%lchnk
   ncol  = state%ncol

   do k = 1, pver
      do i = 1, ncol
         dtcond(i,k,lchnk) = state%s(i,k)
      end do
   end do

   do m = 1, ppcnst
      do k = 1, pver
         do i = 1, ncol
            dqcond(i,k,m,lchnk) = state%q(i,k,m)
         end do
      end do
   end do
end subroutine diag_conv_tend_ini
!===============================================================================

  subroutine diag_phys_writeout(state)

!----------------------------------------------------------------------- 
! 
! Purpose: record dynamics variables on physics grid
!
!-----------------------------------------------------------------------
    use physconst,     only: gravit, rga, rair
    use wv_saturation, only: aqsat
    use constituents, only:  cnst_need_pdeldry
    use time_manager,    only: get_nstep
#if ( defined COUP_CSM )
    use ccsm_msg, only: psl   ! Store sea-level pressure for CCSM
#endif
#if ( defined SCAM )
use pmgrid,    only: plev, plevp
#include <comfrc.h>
#endif

#include <comctl.h>
!-----------------------------------------------------------------------
!
! Arguments
!
   type(physics_state), intent(inout) :: state
!
!---------------------------Local workspace-----------------------------
!
    real(r8) ftem(pcols,pver) ! temporary workspace
    real(r8) psl_tmp(pcols)   ! Sea Level Pressure
    real(r8) z3(pcols,pver)   ! geo-potential height
    real(r8) p_surf(pcols)    ! data interpolated to a pressure surface
    real(r8) tem2(pcols,pver) ! temporary workspace
    real(r8) timestep(pcols)  ! used for outfld call

    integer k, m, lchnk, ncol, nstep
!
!-----------------------------------------------------------------------
!
    lchnk = state%lchnk
    ncol  = state%ncol

!
! Output NSTEP for debugging
!
   nstep = get_nstep()
   timestep(:ncol) = nstep
   call outfld ('NSTEP   ',timestep, pcols, lchnk)

!This field has the same name as the one that is needed for BFB SCAM
!IOP datasets so don't outfield it here (outfield in tfiltmassfix)
!
    call outfld('T       ',state%t , pcols   ,lchnk   )
    call outfld('PS      ',state%ps, pcols   ,lchnk   )
    call outfld('U       ',state%u , pcols   ,lchnk   )
    call outfld('V       ',state%v , pcols   ,lchnk   )
    do m=1,ppcnst
       if ( cnst_cam_outfld(m) ) then
          call outfld(cnst_name(m),state%q(1,1,m),pcols ,lchnk )
       end if
    end do
    if ( cnst_need_pdeldry .and. ( .not. adiabatic .and. .not. ideal_phys ) ) &
         call outfld('PDELDRY ',state%pdeldry, pcols,   lchnk     )

    if ( cnst_need_pdeldry )  call outfld ('PSDRY',  state%psdry, pcols, lchnk) 

    call outfld('PHIS    ',state%phis,    pcols,   lchnk     )



#if (defined BFB_CAM_SCAM_IOP )
    call outfld('phis    ',state%phis,    pcols,   lchnk     )
#endif

!
! Add height of surface to midpoint height above surface 
!
    do k = 1, pver
       z3(:ncol,k) = state%zm(:ncol,k) + state%phis(:ncol)*rga
    end do
    call outfld('Z3      ',z3,pcols,lchnk)
!           
! Output Z3 on 500mb, 300, 50 and 700 mb surface
!
    call vertinterp(ncol, pcols, pver, state%pmid, 70000._r8, z3, p_surf)
    call outfld('Z700    ', p_surf, pcols, lchnk)
    call vertinterp(ncol, pcols, pver, state%pmid, 50000._r8, z3, p_surf)
    call outfld('Z500    ', p_surf, pcols, lchnk)
    call vertinterp(ncol, pcols, pver, state%pmid,  5000._r8, z3, p_surf)
    call outfld('Z050    ', p_surf, pcols, lchnk)
    call vertinterp(ncol, pcols, pver, state%pmid, 30000._r8, z3, p_surf)
    call outfld('Z300    ', p_surf, pcols, lchnk)
!
! Quadratic height fiels Z3*Z3
!
    ftem(:ncol,:) = z3(:ncol,:)*z3(:ncol,:)
    call outfld('ZZ      ',ftem,pcols,lchnk)

    ftem(:ncol,:) = z3(:ncol,:)*state%v(:ncol,:)*gravit
    call outfld('VZ      ',ftem,  pcols,lchnk)
!
! Meridional advection fields
!
    ftem(:ncol,:) = state%v(:ncol,:)*state%t(:ncol,:)
    call outfld ('VT      ',ftem    ,pcols   ,lchnk     )

    ftem(:ncol,:) = state%v(:ncol,:)*state%q(:ncol,:,1)
    call outfld ('VQ      ',ftem    ,pcols   ,lchnk     )

    ftem(:ncol,:) = state%v(:ncol,:)**2
    call outfld ('VV      ',ftem    ,pcols   ,lchnk     )

    ftem(:ncol,:) = state%v(:ncol,:) * state%u(:ncol,:)
    call outfld ('VU      ',ftem    ,pcols   ,lchnk     )

! zonal advection

    ftem(:ncol,:) = state%u(:ncol,:)**2
    call outfld ('UU      ',ftem    ,pcols   ,lchnk     )

! Wind speed
    ftem(:ncol,:) = sqrt( state%u(:ncol,:)**2 + state%v(:ncol,:)**2)
    call outfld ('WSPEED  ',ftem    ,pcols   ,lchnk     )

! Vertical velocity and advection

#if ( defined SCAM )
    call outfld('OMEGA   ',wfld,    pcols,   lchnk     )
#else
    call outfld('OMEGA   ',state%omega,    pcols,   lchnk     )
#endif

#if (defined BFB_CAM_SCAM_IOP )
    call outfld('omega   ',state%omega,    pcols,   lchnk     )
#endif

    ftem(:ncol,:) = state%omega(:ncol,:)*state%t(:ncol,:)
    call outfld('OMEGAT  ',ftem,    pcols,   lchnk     )
    ftem(:ncol,:) = state%omega(:ncol,:)*state%u(:ncol,:)
    call outfld('OMEGAU  ',ftem,    pcols,   lchnk     )
!
! Output omega at 850 and 500 mb pressure levels
!
    call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%omega, p_surf)
    call outfld('OMEGA850', p_surf, pcols, lchnk)
    call vertinterp(ncol, pcols, pver, state%pmid, 50000._r8, state%omega, p_surf)
    call outfld('OMEGA500', p_surf, pcols, lchnk)
!     
! Mass of q, by layer and vertically integrated
!
    ftem(:ncol,:) = state%q(:ncol,:,1) * state%pdel(:ncol,:) * rga
    call outfld ('MQ      ',ftem    ,pcols   ,lchnk     )

    do k=2,pver
       ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
    end do
    call outfld ('TMQ     ',ftem, pcols   ,lchnk     )
!
! Relative humidity
!
    call aqsat (state%t    ,state%pmid  ,tem2    ,ftem    ,pcols   , &
         ncol ,pver  ,1       ,pver    )
    ftem(:ncol,:) = state%q(:ncol,:,1)/ftem(:ncol,:)*100.
    call outfld ('RELHUM  ',ftem    ,pcols   ,lchnk     )
!
! Sea level pressure
!
    call cpslec (ncol, state%pmid, state%phis, state%ps, state%t,psl_tmp, gravit, rair) 
    call outfld ('PSL     ',psl_tmp  ,pcols, lchnk     )
#if ( defined COUP_CSM )
    psl(:ncol,lchnk) = psl_tmp(:ncol)
#endif
!
! Output T,q,u,v fields on pressure surfaces
!
    call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%t, p_surf)
    call outfld('T850    ', p_surf, pcols, lchnk )
    call vertinterp(ncol, pcols, pver, state%pmid, 30000._r8, state%t, p_surf)
    call outfld('T300    ', p_surf, pcols, lchnk )
    call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%q(1,1,1), p_surf)
    call outfld('Q850    ', p_surf, pcols, lchnk )
    call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, state%q(1,1,1), p_surf)
    call outfld('Q200    ', p_surf, pcols, lchnk )
    call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%u, p_surf)
    call outfld('U850    ', p_surf, pcols, lchnk )
    call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, state%u, p_surf)
    call outfld('U200    ', p_surf, pcols, lchnk )
    call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%v, p_surf)
    call outfld('V850    ', p_surf, pcols, lchnk )
    call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, state%v, p_surf)
    call outfld('V200    ', p_surf, pcols, lchnk )

    ftem(:ncol,:) = state%t(:ncol,:)*state%t(:ncol,:)
    call outfld('TT      ',ftem    ,pcols   ,lchnk   )

#if ( defined COUP_CSM )
!
! Output U, V, T, Q, P and Z at bottom level
!
    call outfld ('UBOT    ', state%u(1,pver)  ,  pcols, lchnk)
    call outfld ('VBOT    ', state%v(1,pver)  ,  pcols, lchnk)
    call outfld ('QBOT    ', state%q(1,pver,1),  pcols, lchnk)
    call outfld ('ZBOT    ', state%zm(1,pver) , pcols, lchnk)
#endif

    return
  end subroutine diag_phys_writeout
!===============================================================================

subroutine diag_conv(state, ztodt)

!----------------------------------------------------------------------- 
! 
! Output diagnostics associated with all convective processes.
!
!-----------------------------------------------------------------------
   use physconst,     only: cpair

! Arguments:

   real(r8),            intent(in) :: ztodt   ! timestep for computing physics tendencies
   type(physics_state), intent(in) :: state

! Local variables:
   
   integer :: i, k, m, lchnk, ncol
   real(r8) :: rtdt

   lchnk = state%lchnk
   ncol  = state%ncol

   rtdt = 1./ztodt

! Total convection tendencies.

   do k = 1, pver
      do i = 1, ncol
         dtcond(i,k,lchnk) = (state%s(i,k) - dtcond(i,k,lchnk))*rtdt / cpair
      end do
   end do
   call outfld('DTCOND  ', dtcond(:,:,lchnk), pcols, lchnk)

   do m = 1, ppcnst
      do k = 1, pver
         do i = 1, ncol
            dqcond(i,k,m,lchnk) = (state%q(i,k,m) - dqcond(i,k,m,lchnk))*rtdt
         end do
      end do
   end do
   do m = 1, ppcnst
      if ( cnst_cam_outfld(m) ) then
         call outfld(dcconnam(m), dqcond(:,:,m,lchnk), pcols, lchnk)
      end if
   end do
end subroutine diag_conv
!===============================================================================

subroutine diag_surf (srfflx, surface, icefrac, ocnfrac, landfrac, &
                      sicthk, snowhland, snowhice, tsice, trefmxav, &
                      trefmnav )

!----------------------------------------------------------------------- 
! 
! Purpose: record surface diagnostics
!
!-----------------------------------------------------------------------

   use comsrf, only: srfflx_state, surface_state, tsnam

#if ( defined COUP_CSM )
    use time_manager, only: is_end_curr_day
#endif
!-----------------------------------------------------------------------
!
! Input arguments
!
    type(srfflx_state),  intent(in) :: srfflx
    type(surface_state), intent(in) :: surface

    real(r8), intent(in) :: icefrac(pcols)   ! ice fraction
    real(r8), intent(in) :: ocnfrac(pcols)   ! ocn fraction
    real(r8), intent(in) :: landfrac(pcols)  ! land fraction
    real(r8), intent(in) :: sicthk(pcols)    ! sea-ice thickness
    real(r8), intent(in) :: snowhland(pcols) ! equivalent liquid water snow depth
    real(r8), intent(in) :: snowhice(pcols)  ! equivalent liquid water snow depth
    real(r8), intent(in) :: tsice(pcols)     ! surface T over seaice

    real(r8), intent(inout) :: trefmnav(pcols) ! daily minimum tref  
    real(r8), intent(inout) :: trefmxav(pcols) ! daily maximum tref
!
!---------------------------Local workspace-----------------------------
!
    integer i,k             ! indexes
    integer :: lchnk        ! chunk identifier
    integer :: ncol         ! longitude dimension
!
!-----------------------------------------------------------------------
!
    lchnk = srfflx%lchnk
    ncol  = srfflx%ncol

    call outfld('SHFLX',    srfflx%shf,       pcols, lchnk)
    call outfld('LHFLX',    srfflx%lhf,       pcols, lchnk)
    call outfld('QFLX',     srfflx%cflx(1,1), pcols, lchnk)

    call outfld('TAUX',     srfflx%wsx,       pcols, lchnk)
    call outfld('TAUY',     srfflx%wsy,       pcols, lchnk)
    call outfld('TREFHT  ', srfflx%tref,      pcols, lchnk)
    call outfld('TREFHTMX', srfflx%tref,      pcols, lchnk)
    call outfld('TREFHTMN', srfflx%tref,      pcols, lchnk)
#if ( defined COUP_CSM )
    call outfld('QREFHT',   srfflx%qref,      pcols, lchnk)
#endif
#if (defined BFB_CAM_SCAM_IOP )
    call outfld('shflx   ',srfflx%shf,   pcols,   lchnk)
    call outfld('lhflx   ',srfflx%lhf,   pcols,   lchnk)
    call outfld('trefht  ',srfflx%tref,  pcols,   lchnk)
#endif
!
! Ouput ocn and ice fractions
!
    call outfld('LANDFRAC', landfrac,         pcols, lchnk)
    call outfld('ICEFRAC',  icefrac,          pcols, lchnk)
    call outfld('OCNFRAC',  ocnfrac,          pcols, lchnk)

#if ( defined COUP_CSM )
!
! Compute daily minimum and maximum of TREF
!
    do i = 1,ncol
       trefmxav(i) = max(srfflx%tref(i),trefmxav(i))
       trefmnav(i) = min(srfflx%tref(i),trefmnav(i))
    end do
    if (is_end_curr_day()) then
       call outfld('TREFMXAV', trefmxav,pcols,   lchnk     )
       call outfld('TREFMNAV', trefmnav,pcols,   lchnk     )
       trefmxav(:ncol) = -1.0e36
       trefmnav(:ncol) =  1.0e36
    endif

#else

    do k=1,pvermx
       call outfld(tsnam(k), surface%tssub(1,k), pcols, lchnk)
    end do
    call outfld('SICTHK',   sicthk,           pcols, lchnk)
    call outfld('TSICE',    tsice,            pcols, lchnk)
#endif

    call outfld('TS',       srfflx%ts,        pcols, lchnk)
    call outfld('TSMN',     srfflx%ts,        pcols, lchnk)
    call outfld('TSMX',     srfflx%ts,        pcols, lchnk)
    call outfld('SNOWHLND', snowhland,        pcols, lchnk)
    call outfld('SNOWHICE', snowhice ,        pcols, lchnk)
    call outfld('TBOT',     surface%tbot,     pcols, lchnk)

    call outfld('ASDIR',    srfflx%asdir,     pcols, lchnk)
    call outfld('ASDIF',    srfflx%asdif,     pcols, lchnk)
    call outfld('ALDIR',    srfflx%aldir,     pcols, lchnk)
    call outfld('ALDIF',    srfflx%aldif,     pcols, lchnk)
    call outfld('SST',      srfflx%sst,       pcols, lchnk)

end subroutine diag_surf

!#######################################################################

   subroutine diag_physvar_ic (lchnk, pbuf)
!
!---------------------------------------------
!
! Purpose: record physics variables on IC file
!
!---------------------------------------------
!
   use phys_buffer, only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
   use buffer     , only: pblht, tpert, qpert
   use comsrf
#if (!defined COUP_CSM)
#if (!defined COUP_SOM)
   use sst_data, only: sst
#endif
#endif
!
! Arguments
!
   integer       , intent(in) :: lchnk  ! chunk identifier
   type(pbuf_fld), intent(in), dimension(pbuf_size_max) :: pbuf
!
!---------------------------Local workspace-----------------------------
!
   integer  :: k                 ! indices
   integer  :: itim, ifld        ! indices
   real(r8), pointer, dimension(:,:) :: cwat_var
!
!-----------------------------------------------------------------------
!
   if( write_inithist() ) then
!
! Associate pointers with physics buffer fields
!
      itim = pbuf_old_tim_idx()
      ifld = pbuf_get_fld_idx('QCWAT')
      cwat_var => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
      call outfld('QCWAT&IC   ',cwat_var, pcols,lchnk)
      ifld = pbuf_get_fld_idx('TCWAT')
      cwat_var => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
      call outfld('TCWAT&IC   ',cwat_var, pcols,lchnk)
      ifld = pbuf_get_fld_idx('LCWAT')
      cwat_var => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
      call outfld('LCWAT&IC   ',cwat_var, pcols,lchnk)
      ifld = pbuf_get_fld_idx('CLD')
      cwat_var => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
      call outfld('CLOUD&IC   ',cwat_var, pcols,lchnk)

      call outfld('PBLH&IC    '     , pblht    (1,  lchnk), pcols, lchnk)
      call outfld('TPERT&IC   '     , tpert    (1,  lchnk), pcols, lchnk)
      call outfld('QPERT&IC   '     , qpert    (1,1,lchnk), pcols, lchnk)

#if ( ! defined COUP_CSM )

      call outfld('TSICERAD&IC'     , tsice_rad(1,  lchnk), pcols, lchnk)
      call outfld('TSICE&IC   '     , tsice    (1,  lchnk), pcols, lchnk)
      call outfld('SNOWHICE&IC'     , snowhice (1,  lchnk), pcols, lchnk)
      call outfld('ICEFRAC&IC '     , icefrac  (1,  lchnk), pcols, lchnk)
      call outfld('SICTHK&IC  '     , sicthk   (1,  lchnk), pcols, lchnk)
#if ( defined COUP_SOM )
      call outfld('TSOCN&IC   '     , tsocn    (1,  lchnk), pcols, lchnk)
#else
      call outfld('TSOCN&IC   '     , sst      (1,  lchnk), pcols, lchnk)
#endif
      call outfld('TS&IC      '     , srfflx_state2d (lchnk)%ts   , pcols, lchnk)
      call outfld('TBOT&IC    '     , surface_state2d(lchnk)%tbot , pcols, lchnk)
      do k=1,pvermx
         call outfld(trim(tsnam(k))//'&IC', surface_state2d(lchnk)%tssub(1,k), pcols, lchnk)
      end do
#endif
   end if

   return
   end subroutine diag_physvar_ic

end module diagnostics
