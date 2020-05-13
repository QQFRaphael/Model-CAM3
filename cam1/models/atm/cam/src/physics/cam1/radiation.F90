#include <params.h>
#include <misc.h>

module radiation

!---------------------------------------------------------------------------------
! Purpose:
!
! Provides the CAM interface to the radiation code
!
! Revision history:
! May  2004, D. B. Coleman,  Initial version of interface module.
! July 2004, B. Eaton,       Use interfaces from new shortwave, longwave, and ozone modules.
!---------------------------------------------------------------------------------

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use ppgrid,        only: pcols, pver, pverp
  use physconst,     only: cpair, cappa
  use abortutils,    only: endrun
  use error_messages, only: handle_err

  implicit none
  private
  save

  public :: radiation_register  ! registers radiation physics buffer fields
  public :: radiation_init  ! calls radini
  public :: radiation_tend  ! moved from radctl.F90

!===============================================================================
contains
!===============================================================================

  subroutine radiation_register
!-----------------------------------------------------------------------
! 
! Register radiation fields in the physics buffer
!
!-----------------------------------------------------------------------

    use phys_buffer,  only: pbuf_times, pbuf_add

    integer idx
    call pbuf_add('QRS' , 'global', 1,pver,1, idx) ! shortwave radiative heating rate 
    call pbuf_add('QRL' , 'global', 1,pver,1, idx) ! longwave  radiative heating rate 


  end subroutine radiation_register

!===============================================================================

  subroutine radiation_init()
!-----------------------------------------------------------------------
!
! Initialize the radiation parameterization, add fields to the history buffer
! 
!-----------------------------------------------------------------------

    use history,       only: addfld, add_default, phys_decomp
    use physconst,     only: gravit, cpair, epsilo, stebol, &
                             pstd, mwdry, mwco2, mwo3
    use radsw,         only: radsw_init
    use radlw,         only: radlw_init
    use radae,         only: radae_init

    call radsw_init(gravit)
    call radlw_init(gravit, stebol)
    call radae_init(gravit, epsilo, stebol, pstd, mwdry, mwco2, mwo3)

    ! Shortwave radiation
    call addfld ('SOLIN   ','W/m2    ',1,    'A','Solar insolation',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('SOLL    ','W/m2    ',1,    'A','Solar downward near infrared direct  to surface',phys_decomp, &
                                                                                                          sampling_seq='rad_lwsw')
    call addfld ('SOLS    ','W/m2    ',1,    'A','Solar downward visible direct  to surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('SOLLD   ','W/m2    ',1,    'A','Solar downward near infrared diffuse to surface',phys_decomp, &
                                                                                                          sampling_seq='rad_lwsw')
    call addfld ('SOLSD   ','W/m2    ',1,    'A','Solar downward visible diffuse to surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('QRS     ','K/s     ',pver, 'A','Solar heating rate',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNS    ','W/m2    ',1,    'A','Net solar flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNT    ','W/m2    ',1,    'A','Net solar flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNTOA  ','W/m2    ',1,    'A','Net solar flux at top of atmosphere',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNTOAC ','W/m2    ',1,    'A','Clearsky net solar flux at top of atmosphere',phys_decomp, &
                                                                                                          sampling_seq='rad_lwsw')
    call addfld ('FSN200  ','W/m2    ',1,    'A','Net shortwave flux at 200 mb',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSN200C ','W/m2    ',1,    'A','Clearsky net shortwave flux at 200 mb',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNTC   ','W/m2    ',1,    'A','Clearsky net solar flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNSC   ','W/m2    ',1,    'A','Clearsky net solar flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSDSC   ','W/m2    ',1,    'A','Clearsky downwelling solar flux at surface',phys_decomp, &
                                                                                                        sampling_seq='rad_lwsw')
    call addfld ('FSDS    ','W/m2    ',1,    'A','Downwelling solar flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FUS     ','W/m2    ',pverp,'I','Shortwave upward flux',phys_decomp)
    call addfld ('FDS     ','W/m2    ',pverp,'I','Shortwave downward flux',phys_decomp)
    call addfld ('FUSC    ','W/m2    ',pverp,'I','Shortwave clear-sky upward flux',phys_decomp)
    call addfld ('FDSC    ','W/m2    ',pverp,'I','Shortwave clear-sky downward flux',phys_decomp)
    call addfld ('FSNIRTOA','W/m2    ',1,    'A','Net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', &
            phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNRTOAC','W/m2    ',1,    'A','Clearsky net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', &
            phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNRTOAS','W/m2    ',1,    'A','Net near-infrared flux (>= 0.7 microns) at top of atmosphere',phys_decomp, &
                                                                                                        sampling_seq='rad_lwsw')
    call addfld ('SWCF    ','W/m2    ',1,    'A','Shortwave cloud forcing',phys_decomp, sampling_seq='rad_lwsw')

    call add_default ('SOLIN   ', 1, ' ')
    call add_default ('QRS     ', 1, ' ')
    call add_default ('FSNS    ', 1, ' ')
    call add_default ('FSNT    ', 1, ' ')
    call add_default ('FSNTOA  ', 1, ' ')
    call add_default ('FSNTOAC ', 1, ' ')
    call add_default ('FSNTC   ', 1, ' ')
    call add_default ('FSNSC   ', 1, ' ')
    call add_default ('FSDSC   ', 1, ' ')
    call add_default ('FSDS    ', 1, ' ')
    call add_default ('SWCF    ', 1, ' ')

    ! aerosol forcing-only calculations
    call addfld ('FSNT_RF ','W/m^2   ',1, 'I','Total column absorbed solar flux (radforce)' ,phys_decomp)
    call addfld ('FSNTC_RF','W/m^2   ',1, 'I','Clear sky total column absorbed solar flux (radforce)' ,phys_decomp)
    call addfld ('FSNS_RF ','W/m^2   ',1, 'I','Surface absorbed solar flux (radforce)' ,phys_decomp)
    call addfld ('FSNSC_RF','W/m^2   ',1, 'I','Clear sky surface absorbed solar flux (radforce)' ,phys_decomp)
    call addfld ('QRS_RF  ','K/s     ',pver, 'I','Solar heating rate (radforce)' ,phys_decomp)

    ! Longwave radiation
    call addfld ('QRL     ','K/s     ',pver, 'A','Longwave heating rate',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLNS    ','W/m2    ',1,    'A','Net longwave flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLNT    ','W/m2    ',1,    'A','Net longwave flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLUT    ','W/m2    ',1,    'A','Upwelling longwave flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLUTC   ','W/m2    ',1,    'A','Clearsky upwelling longwave flux at top of model',phys_decomp, &
                                                                                                        sampling_seq='rad_lwsw')
    call addfld ('FLNTC   ','W/m2    ',1,    'A','Clearsky net longwave flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLN200  ','W/m2    ',1,    'A','Net longwave flux at 200 mb',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLN200C ','W/m2    ',1,    'A','Clearsky net longwave flux at 200 mb',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLNSC   ','W/m2    ',1,    'A','Clearsky net longwave flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('LWCF    ','W/m2    ',1,    'A','Longwave cloud forcing',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FUL     ','W/m2    ',pverp,'I','Longwave upward flux',phys_decomp)
    call addfld ('FDL     ','W/m2    ',pverp,'I','Longwave downward flux',phys_decomp)
    call addfld ('FULC    ','W/m2    ',pverp,'I','Longwave clear-sky upward flux',phys_decomp)
    call addfld ('FDLC    ','W/m2    ',pverp,'I','Longwave clear-sky downward flux',phys_decomp)
    call add_default ('QRL     ', 1, ' ')
    call add_default ('FLNS    ', 1, ' ')
    call add_default ('FLNT    ', 1, ' ')
    call add_default ('FLUT    ', 1, ' ')
    call add_default ('FLUTC   ', 1, ' ')
    call add_default ('FLNTC   ', 1, ' ')
    call add_default ('FLNSC   ', 1, ' ')
    call add_default ('LWCF    ', 1, ' ')

    ! Heating rate needed for d(theta)/dt computation
    call addfld ('HR      ','K/s     ',pver, 'A','Heating rate needed for d(theta)/dt computation',phys_decomp)

  end subroutine radiation_init

!===============================================================================
  
  subroutine radiation_tend(state,ptend,pbuf, &
       surface_state2d, srfflx_state2d, &
       landfrac,landm,icefrac,snowh, &
       fsns,    fsnt, flns,    flnt,  &
       fsds, concld, net_flx)

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Driver for radiation computation.
    ! 
    ! Method: 
    ! Radiation uses cgs units, so conversions must be done from
    ! model fields to radiation fields.
    !
    ! Revision history:
    ! May 2004    D.B. Coleman     Merge of code from radctl.F90 and parts of tphysbc.F90.
    ! 2004-08-09  B. Eaton         Add pointer variables for constituents.
    ! 2004-08-24  B. Eaton         Access O3 and GHG constituents from chem_get_cnst.
    ! 2004-08-30  B. Eaton         Replace chem_get_cnst by rad_constituent_get.
    !-----------------------------------------------------------------------


    use phys_buffer,     only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
    use phys_grid,       only: get_rlat_all_p, get_rlon_all_p
    use param_cldoptics, only: param_cldoptics_calc
    use physics_types,   only: physics_state, physics_ptend
    use time_manager,    only: get_curr_calday
    use comsrf,          only: surface_state, srfflx_state
    use history,         only: outfld
    use cloudsimulator,  only: doisccp, cloudsimulator_run
    use radheat,         only: radheat_tend
    use radsw,           only: radcswmx
    use radlw,           only: radclwmx
    use ppgrid
    use pspect
    use constituents,    only: ppcnst, cnst_get_ind
    use prescribed_aerosols, only: get_aerosol, naer_all, aerosol_diagnostics, &
         aerosol_indirect, get_rf_scales, get_int_scales, radforce, idxVOLC
    use wv_saturation,    only: aqsat
    use rad_constituents, only: rad_constituents_get
    use physconst,        only: cpair, epsilo, stebol
    use aer_optics,       only: idxVIS
    use aerosol_intr,     only: set_aerosol_from_prognostics

#if ( defined SCAM )
#include <max.h>
   use scamMod, only: switch,have_cld,cldobs,have_clwp,clwpobs,have_tg,tground
#endif
                  
#include <comctl.h>
#include <comsol.h>

    !
    ! Arguments
    !

    real(r8), intent(in) :: landfrac(pcols)                ! land fraction
    real(r8), intent(in) :: landm(pcols)                   ! land fraction ramp
    real(r8), intent(in) :: icefrac(pcols)                ! land fraction
    real(r8), intent(in) :: snowh(pcols)                  ! Snow depth (liquid water equivalent)
    real(r8), intent(inout) :: fsns(pcols)                   ! Surface solar absorbed flux
    real(r8), intent(inout) :: fsnt(pcols)                   ! Net column abs solar flux at model top
    real(r8), intent(inout) :: flns(pcols)                   ! Srf longwave cooling (up-down) flux
    real(r8), intent(inout) :: flnt(pcols)                   ! Net outgoing lw flux at model top
    real(r8), intent(out) :: fsds(pcols)                   ! Surface solar down flux
    real(r8), intent(in) ::  concld(pcols,pver)             ! should be pbuf
    real(r8), intent(inout):: net_flx(pcols)

    type(physics_state), intent(in), target :: state
    type(physics_ptend), intent(out)        :: ptend
    type(pbuf_fld),      intent(inout), dimension(pbuf_size_max) :: pbuf
    type(surface_state), intent(inout)      :: surface_state2d
    type(srfflx_state),  intent(in)         :: srfflx_state2d


    !
    ! Local variables
    !
    integer nmxrgn(pcols)                      ! Number of maximally overlapped regions
    real(r8) pmxrgn(pcols,pverp)               ! Maximum values of pressure for each
    !                                             !    maximally overlapped region.
    !                                             !    0->pmxrgn(i,1) is range of pressure for
    !                                             !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
    !                                             !    2nd region, etc
    real(r8) emis(pcols,pver)                  ! Cloud longwave emissivity
    real(r8) :: cicewp(pcols,pver)             ! in-cloud cloud ice water path
    real(r8) :: cliqwp(pcols,pver)             ! in-cloud cloud liquid water path
    real(r8) cltot(pcols)                      ! Diagnostic total cloud cover
    real(r8) cllow(pcols)                      !       "     low  cloud cover
    real(r8) clmed(pcols)                      !       "     mid  cloud cover
    real(r8) clhgh(pcols)                      !       "     hgh  cloud cover
    real(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables

    integer itim, ifld
    real(r8), pointer, dimension(:,:) :: rel     ! liquid effective drop radius (microns)
    real(r8), pointer, dimension(:,:) :: rei     ! ice effective drop size (microns)
    real(r8), pointer, dimension(:,:) :: cld        ! cloud fraction
    real(r8), pointer, dimension(:,:) :: qrs        ! shortwave radiative heating rate 
    real(r8), pointer, dimension(:,:) :: qrl        ! longwave  radiative heating rate 

    integer lchnk, ncol
    real(r8) :: calday                        ! current calendar day
    real(r8) :: clat(pcols)                   ! current latitudes(radians)
    real(r8) :: clon(pcols)                   ! current longitudes(radians)
    real(r8) coszrs(pcols)                     ! Cosine solar zenith angle
   logical  :: conserve_energy = .true.       ! flag to carry (QRS,QRL)*dp across time steps


    !
    ! Local variables from radctl
    !
    integer nspint            ! Num of spctrl intervals across solar spectrum
    integer naer_groups       ! Num of aerosol groups for optical diagnostics
    parameter ( nspint = 19 )
    parameter ( naer_groups = 7 )    ! current groupings are sul, sslt, all carbons, all dust, background, and all aerosols
    real(r8) :: pmxrgnrf(pcols,pverp)             ! temporary copy of pmxrgn
    integer  :: nmxrgnrf(pcols)     ! temporary copy of nmxrgn
    integer i, k              ! index
    integer :: istat
    real(r8) solin(pcols)         ! Solar incident flux
    real(r8) fsntoa(pcols)        ! Net solar flux at TOA
    real(r8) fsntoac(pcols)       ! Clear sky net solar flux at TOA
    real(r8) fsnirt(pcols)        ! Near-IR flux absorbed at toa
    real(r8) fsnrtc(pcols)        ! Clear sky near-IR flux absorbed at toa
    real(r8) fsnirtsq(pcols)      ! Near-IR flux absorbed at toa >= 0.7 microns
    real(r8) fsntc(pcols)         ! Clear sky total column abs solar flux
    real(r8) fsnsc(pcols)         ! Clear sky surface abs solar flux
    real(r8) fsdsc(pcols)         ! Clear sky surface downwelling solar flux
    real(r8) flut(pcols)          ! Upward flux at top of model
    real(r8) lwcf(pcols)          ! longwave cloud forcing
    real(r8) swcf(pcols)          ! shortwave cloud forcing
    real(r8) flutc(pcols)         ! Upward Clear Sky flux at top of model
    real(r8) flntc(pcols)         ! Clear sky lw flux at model top
    real(r8) flnsc(pcols)         ! Clear sky lw flux at srf (up-down)
    real(r8) fln200(pcols)        ! net longwave flux interpolated to 200 mb
    real(r8) fln200c(pcols)       ! net clearsky longwave flux interpolated to 200 mb
    real(r8) fns(pcols,pverp)     ! net shortwave flux
    real(r8) fcns(pcols,pverp)    ! net clear-sky shortwave flux
    real(r8) fsn200(pcols)        ! fns interpolated to 200 mb
    real(r8) fsn200c(pcols)       ! fcns interpolated to 200 mb
    real(r8) fnl(pcols,pverp)     ! net longwave flux
    real(r8) fcnl(pcols,pverp)    ! net clear-sky longwave flux

    real(r8) pbr(pcols,pver)      ! Model mid-level pressures (dynes/cm2)
    real(r8) pnm(pcols,pverp)     ! Model interface pressures (dynes/cm2)
    real(r8) eccf                 ! Earth/sun distance factor
    real(r8) rh(pcols,pver)       ! level relative humidity (fraction)
    real(r8) lwupcgs(pcols)       ! Upward longwave flux in cgs units

    real(r8) esat(pcols,pver)     ! saturation vapor pressure
    real(r8) qsat(pcols,pver)     ! saturation specific humidity

    real(r8) :: frc_day(pcols) ! = 1 for daylight, =0 for night colums
    real(r8) :: aertau(pcols,nspint,naer_groups) ! Aerosol column optical depth
    real(r8) :: aerssa(pcols,nspint,naer_groups) ! Aerosol column averaged single scattering albedo
    real(r8) :: aerasm(pcols,nspint,naer_groups) ! Aerosol column averaged asymmetry parameter
    real(r8) :: aerfwd(pcols,nspint,naer_groups) ! Aerosol column averaged forward scattering

    real(r8) aerosol(pcols, pver, naer_all) ! aerosol mass mixing ratios
    real(r8) scales(naer_all)               ! scaling factors for aerosols

    real(r8), parameter :: cgs2mks = 1.e-3_r8
 
    real(r8), pointer, dimension(:,:) :: n2o   ! nitrous oxide mass mixing ratio
    real(r8), pointer, dimension(:,:) :: ch4   ! methane mass mixing ratio
    real(r8), pointer, dimension(:,:) :: cfc11 ! cfc11 mass mixing ratio
    real(r8), pointer, dimension(:,:) :: cfc12 ! cfc12 mass mixing ratio
    real(r8), pointer, dimension(:,:) :: o3    ! Ozone mass mixing ratio

    character(*), parameter :: name = 'radiation_tend'
!----------------------------------------------------------------------

    call t_startf ('radiation_tend')

    lchnk = state%lchnk
    ncol = state%ncol

    calday = get_curr_calday()

    itim = pbuf_old_tim_idx()
    ifld = pbuf_get_fld_idx('CLD')
    cld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
    ifld = pbuf_get_fld_idx('QRS')
    qrs => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
    ifld = pbuf_get_fld_idx('QRL')
    qrl => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
    ifld = pbuf_get_fld_idx('REL')
    rel  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
    ifld = pbuf_get_fld_idx('REI')
    rei  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
    
#if ( defined SCAM )
    !  For CRM, make cloud equal to input observations:
    if(switch(CRM_SW+1).and.have_cld) then
       do k = 1,pver
          cld(:ncol,k)= cldobs(k)
       enddo
    endif
#endif

    !
    ! Cosine solar zenith angle for current time step
    !
    call get_rlat_all_p(lchnk, ncol, clat)
    call get_rlon_all_p(lchnk, ncol, clon)
    call zenith (calday, clat, clon, coszrs, ncol)


    if (dosw .or. dolw) then

       ! Compute cloud water/ice paths and optical properties for input to radiation
       call t_startf('cldoptics')
       call param_cldoptics_calc(state, cld, landfrac, landm,icefrac, &
            cicewp, cliqwp, emis, rel, rei, pmxrgn, nmxrgn, snowh)
       call t_stopf('cldoptics')

#if ( defined SCAM )
       if(switch(CRM_SW+1).and.have_clwp)then
          do k=1,pver
             ! For CRM, make cloud liquid water path equal to input observations
             cliqwp(:ncol,k) = clwpobs(k)
          end do
       endif
#endif

       ! Get ozone mass mixing ratio.
       call rad_constituents_get('O3', state, o3)

       ! Set chunk dependent radiation input
       call radinp(ncol, state%pmid, state%pint, pbr, pnm, eccf)

       ! Solar radiation computation
       !
       if (dosw) then

          !
          ! calculate heating with aerosols
          !
          call aqsat(state%t, state%pmid, esat, qsat, pcols, &
               ncol, pver, 1, pver)

          ! calculate relative humidity
          rh(1:ncol,1:pver) = state%q(1:ncol,1:pver,1) / qsat(1:ncol,1:pver) * &
               ((1.0_r8 - epsilo) * qsat(1:ncol,1:pver) + epsilo) / &
               ((1.0_r8 - epsilo) * state%q(1:ncol,1:pver,1) + epsilo)

          if (radforce) then

             pmxrgnrf = pmxrgn
             nmxrgnrf = nmxrgn

             call get_rf_scales(scales)

             call get_aerosol(lchnk, state%pint, aerosol, scales)

             ! overwrite with prognostics aerosols
             call set_aerosol_from_prognostics (state, aerosol)

             call aerosol_indirect(ncol,lchnk,landfrac,state%pmid,state%t,state%q,cld,state%zm,rel)

             call t_startf('radcswmx_rf')

             call radcswmx(lchnk   ,ncol ,                            &
                  pnm     ,pbr     ,state%q(1,1,1)   ,rh      ,o3  , &
                  aerosol ,cld     ,cicewp  ,cliqwp  ,rel     , &
                  rei     ,eccf    ,coszrs  ,scon    ,solin   , &
                  srfflx_state2d%asdir   ,srfflx_state2d%asdif   ,srfflx_state2d%aldir   ,srfflx_state2d%aldif   ,nmxrgnrf, &
                  pmxrgnrf,qrs     ,fsnt    ,fsntc   ,fsntoa  , &
                  fsntoac ,fsnirt  ,fsnrtc  ,fsnirtsq,fsns    , &
                  fsnsc   ,fsdsc   ,fsds    ,surface_state2d%sols    ,surface_state2d%soll    , &
                  surface_state2d%solsd   ,surface_state2d%solld   ,frc_day ,                   &
                  aertau  ,aerssa  ,aerasm  ,aerfwd  ,fns     , &
                  fcns    )

             call t_stopf('radcswmx_rf')

             !
             ! Convert units of shortwave fields needed by rest of model from CGS to MKS
             !

             do i = 1, ncol
 
                solin(i) = solin(i)*cgs2mks
                fsnt(i)  = fsnt(i) *cgs2mks
                fsns(i)  = fsns(i) *cgs2mks
                fsntc(i) = fsntc(i)*cgs2mks
                fsnsc(i) = fsnsc(i)*cgs2mks
             end do
             ftem(:ncol,:pver) = qrs(:ncol,:pver)/cpair

             !
             ! Dump shortwave radiation information to history tape buffer (diagnostics)
             !
             call outfld('QRS_RF  ',ftem  ,pcols,lchnk)
             call outfld('FSNT_RF ',fsnt  ,pcols,lchnk)
             call outfld('FSNS_RF ',fsns  ,pcols,lchnk)
             call outfld('FSNTC_RF',fsntc ,pcols,lchnk)
             call outfld('FSNSC_RF',fsnsc ,pcols,lchnk)

          endif ! if (radforce)

          call get_int_scales(scales)

          call get_aerosol(lchnk, state%pint, aerosol, scales)

          ! overwrite with prognostics aerosols
          call set_aerosol_from_prognostics (state, aerosol)

          call aerosol_indirect(ncol,lchnk,landfrac,state%pmid,state%t,state%q,cld,state%zm,rel)

          call t_startf('radcswmx')

          call radcswmx(lchnk   ,ncol    ,                            &
               pnm     ,pbr     ,state%q(1,1,1)   ,rh      ,o3  , &
               aerosol ,cld     ,cicewp  ,cliqwp  ,rel     , &
               rei     ,eccf    ,coszrs  ,scon    ,solin   , &
               srfflx_state2d%asdir   ,srfflx_state2d%asdif   ,srfflx_state2d%aldir   ,srfflx_state2d%aldif   ,nmxrgn  , &
               pmxrgn  ,qrs     ,fsnt    ,fsntc   ,fsntoa  , &
               fsntoac ,fsnirt  ,fsnrtc  ,fsnirtsq,fsns    , &
               fsnsc   ,fsdsc   ,fsds    ,surface_state2d%sols    ,surface_state2d%soll    , &
               surface_state2d%solsd   ,surface_state2d%solld   ,frc_day ,                   &
               aertau  ,aerssa  ,aerasm  ,aerfwd  ,fns     , &
               fcns)

          call t_stopf('radcswmx')

          !  Output net fluxes at 200 mb
          call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fcns, fsn200c)
          call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fns, fsn200)

          !
          ! Convert units of shortwave fields needed by rest of model from CGS to MKS
          !
          do i=1,ncol
 
             solin(i) = solin(i)*cgs2mks
             fsds(i)  = fsds(i)*cgs2mks
             fsnirt(i)= fsnirt(i)*cgs2mks
             fsnrtc(i)= fsnrtc(i)*cgs2mks
             fsnirtsq(i)= fsnirtsq(i)*cgs2mks
             fsnt(i)  = fsnt(i) *cgs2mks
             fsns(i)  = fsns(i) *cgs2mks
             fsntc(i) = fsntc(i)*cgs2mks
             fsnsc(i) = fsnsc(i)*cgs2mks
             fsdsc(i) = fsdsc(i)*cgs2mks
             fsntoa(i)=fsntoa(i)*cgs2mks
             fsntoac(i)=fsntoac(i)*cgs2mks
             fsn200(i)  = fsn200(i)*cgs2mks
             fsn200c(i) = fsn200c(i)*cgs2mks
          end do
          ftem(:ncol,:pver) = qrs(:ncol,:pver)/cpair
          !
          ! Dump shortwave radiation information to history tape buffer (diagnostics)
          !

          call outfld('frc_day ', frc_day, pcols, lchnk)
          call outfld('SULOD_v ', aertau(:,idxVIS,1) ,pcols,lchnk)
          call outfld('SSLTOD_v', aertau(:,idxVIS,2) ,pcols,lchnk)
          call outfld('CAROD_v ', aertau(:,idxVIS,3) ,pcols,lchnk)
          call outfld('DUSTOD_v', aertau(:,idxVIS,4) ,pcols,lchnk)
          call outfld('BGOD_v  ', aertau(:,idxVIS,5) ,pcols,lchnk)
          call outfld('VOLCOD_v', aertau(:,idxVIS,6) ,pcols,lchnk)
          call outfld('AEROD_v ', aertau(:,idxVIS,7) ,pcols,lchnk)
          call outfld('AERSSA_v', aerssa(:,idxVIS,7) ,pcols,lchnk)
          call outfld('AERASM_v', aerasm(:,idxVIS,7) ,pcols,lchnk)
          call outfld('AERFWD_v', aerfwd(:,idxVIS,7) ,pcols,lchnk)
          call aerosol_diagnostics (state, aerosol)

          call outfld('QRS     ',ftem  ,pcols,lchnk)
          call outfld('SOLIN   ',solin ,pcols,lchnk)
          call outfld('FSDS    ',fsds  ,pcols,lchnk)
          call outfld('FSNIRTOA',fsnirt,pcols,lchnk)
          call outfld('FSNRTOAC',fsnrtc,pcols,lchnk)
          call outfld('FSNRTOAS',fsnirtsq,pcols,lchnk)
          call outfld('FSNT    ',fsnt  ,pcols,lchnk)
          call outfld('FSNS    ',fsns  ,pcols,lchnk)
          call outfld('FSNTC   ',fsntc ,pcols,lchnk)
          call outfld('FSNSC   ',fsnsc ,pcols,lchnk)
          call outfld('FSDSC   ',fsdsc ,pcols,lchnk)
          call outfld('FSNTOA  ',fsntoa,pcols,lchnk)
          call outfld('FSNTOAC ',fsntoac,pcols,lchnk)
          call outfld('SOLS    ',surface_state2d%sols  ,pcols,lchnk)
          call outfld('SOLL    ',surface_state2d%soll  ,pcols,lchnk)
          call outfld('SOLSD   ',surface_state2d%solsd ,pcols,lchnk)
          call outfld('SOLLD   ',surface_state2d%solld ,pcols,lchnk)
          call outfld('FSN200  ',fsn200,pcols,lchnk)
          call outfld('FSN200C ',fsn200c,pcols,lchnk)

       end if   ! dosw
       !
       ! Longwave radiation computation
       !
       if (dolw) then
          !
          ! Convert upward longwave flux units to CGS
          !
          do i=1,ncol
             lwupcgs(i) = srfflx_state2d%lwup(i)*1000._r8
 
#if ( defined SCAM )
             if(switch(CRM_SW+1).and.have_tg) lwupcgs(i) = 1000*stebol*tground(1)**4
#endif
          end do

          ! Get gas phase constituents.
          call rad_constituents_get('N2O',   state, n2o)
          call rad_constituents_get('CH4',   state, ch4)
          call rad_constituents_get('CFC11', state, cfc11)
          call rad_constituents_get('CFC12', state, cfc12)

          call t_startf("radclwmx")
          call radclwmx(lchnk, ncol, doabsems, &
             lwupcgs, state%t, state%q(1,1,1), o3, pbr, &
             pnm, state%lnpmid, state%lnpint, n2o, ch4, &
             cfc11, cfc12, cld, emis, pmxrgn, &
             nmxrgn, qrl, flns, flnt, flnsc, &
             flntc, surface_state2d%flwds, flut, flutc, aerosol(1,1,idxVOLC), &
             fnl, fcnl)
          call t_stopf("radclwmx")
          !
          !  Output fluxes at 200 mb
          !
          call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fnl, fln200)
          call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fcnl, fln200c)
          !
          ! Convert units of longwave fields needed by rest of model from CGS to MKS
          !
          do i=1,ncol
             flnt(i)  = flnt(i)*cgs2mks
             flut(i)  = flut(i)*cgs2mks
             flutc(i) = flutc(i)*cgs2mks
             flns(i)  = flns(i)*cgs2mks
             flntc(i) = flntc(i)*cgs2mks
             fln200(i)  = fln200(i)*cgs2mks
             fln200c(i) = fln200c(i)*cgs2mks
             flnsc(i) = flnsc(i)*cgs2mks
             surface_state2d%flwds(i) = surface_state2d%flwds(i)*cgs2mks
             lwcf(i)=flutc(i) - flut(i)
             swcf(i)=fsntoa(i) - fsntoac(i)
          end do
          !
          ! Dump longwave radiation information to history tape buffer (diagnostics)
          !
          call outfld('QRL     ',qrl(:ncol,:)/cpair,ncol,lchnk)
          call outfld('FLNT    ',flnt  ,pcols,lchnk)
          call outfld('FLUT    ',flut  ,pcols,lchnk)
          call outfld('FLUTC   ',flutc ,pcols,lchnk)
          call outfld('FLNTC   ',flntc ,pcols,lchnk)
          call outfld('FLNS    ',flns  ,pcols,lchnk)
          call outfld('FLNSC   ',flnsc ,pcols,lchnk)
          call outfld('LWCF    ',lwcf  ,pcols,lchnk)
          call outfld('SWCF    ',swcf  ,pcols,lchnk)
          call outfld('FLN200  ',fln200,pcols,lchnk)
          call outfld('FLN200C ',fln200c,pcols,lchnk)

       end if  !dolw

       ! Cloud cover diagnostics
       ! radctl can change pmxrgn and nmxrgn so cldsav needs to follow 
       ! radctl.
       !
       call cldsav (lchnk, ncol, cld, state%pmid, cltot, &
            cllow, clmed, clhgh, nmxrgn, pmxrgn)
       !
       ! Dump cloud field information to history tape buffer (diagnostics)
       !
       call outfld('CLDTOT  ',cltot  ,pcols,lchnk)
       call outfld('CLDLOW  ',cllow  ,pcols,lchnk)
       call outfld('CLDMED  ',clmed  ,pcols,lchnk)
       call outfld('CLDHGH  ',clhgh  ,pcols,lchnk)

       call outfld('CLOUD   ',cld    ,pcols,lchnk) !This should be written out by stratiform package
       if (doisccp) then
          call cloudsimulator_run(state, srfflx_state2d%ts, concld, cld, cliqwp, &
               cicewp, rel, rei, emis, coszrs  )
       end if
    else   !  if (dosw .or. dolw) then

       ! convert radiative heating rates from Q*dp to Q for energy conservation
       if (conserve_energy) then
!DIR$ CONCURRENT
          do k =1 , pver
!DIR$ CONCURRENT
             do i = 1, ncol
                qrs(i,k) = qrs(i,k)/state%pdel(i,k)
                qrl(i,k) = qrl(i,k)/state%pdel(i,k)
             end do
          end do
       end if

    end if   !  if (dosw .or. dolw) then

    ! Compute net radiative heating tendency
    call radheat_tend(state, ptend, qrl, qrs, fsns, &
                      fsnt, flns, flnt, srfflx_state2d%asdir, net_flx)

    ! Compute heating rate for dtheta/dt 
    do k=1,pver
       do i=1,ncol
          ftem(i,k) = (qrs(i,k) + qrl(i,k))/cpair * (1.e5/state%pmid(i,k))**cappa
       end do
    end do
    call outfld('HR      ',ftem    ,pcols   ,lchnk   )

    ! convert radiative heating rates to Q*dp for energy conservation
    if (conserve_energy) then
!DIR$ CONCURRENT
       do k =1 , pver
!DIR$ CONCURRENT
          do i = 1, ncol
             qrs(i,k) = qrs(i,k)*state%pdel(i,k)
             qrl(i,k) = qrl(i,k)*state%pdel(i,k)
          end do
       end do
    end if

    call t_stopf ('radiation_tend')


  end subroutine radiation_tend

!===============================================================================

subroutine radinp(ncol, pmid, pint, pmidrd, pintrd, eccf)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Set latitude and time dependent arrays for input to solar
! and longwave radiation.
! Convert model pressures to cgs.
! 
! Author: CCM1, CMS Contact J. Kiehl
!-----------------------------------------------------------------------
   use shr_orb_mod
   use time_manager, only: get_curr_calday

#include <comsol.h>
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: ncol                 ! number of atmospheric columns

   real(r8), intent(in) :: pmid(pcols,pver)    ! Pressure at model mid-levels (pascals)
   real(r8), intent(in) :: pint(pcols,pverp)   ! Pressure at model interfaces (pascals)
!
! Output arguments
!
   real(r8), intent(out) :: pmidrd(pcols,pver)  ! Pressure at mid-levels (dynes/cm*2)
   real(r8), intent(out) :: pintrd(pcols,pverp) ! Pressure at interfaces (dynes/cm*2)
   real(r8), intent(out) :: eccf                ! Earth-sun distance factor

!
!---------------------------Local variables-----------------------------
!
   integer i                ! Longitude loop index
   integer k                ! Vertical loop index

   real(r8) :: calday       ! current calendar day
   real(r8) :: delta        ! Solar declination angle
!-----------------------------------------------------------------------
!
   calday = get_curr_calday()
   call shr_orb_decl (calday  ,eccen     ,mvelpp  ,lambm0  ,obliqr  , &
                      delta   ,eccf)

!
! Convert pressure from pascals to dynes/cm2
!
   do k=1,pver
      do i=1,ncol
         pmidrd(i,k) = pmid(i,k)*10.0
         pintrd(i,k) = pint(i,k)*10.0
      end do
   end do
   do i=1,ncol
      pintrd(i,pverp) = pint(i,pverp)*10.0
   end do

end subroutine radinp

!===============================================================================

end module radiation

