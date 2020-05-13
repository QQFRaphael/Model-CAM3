#include <misc.h>
#include <params.h>

subroutine tphysbc (ztodt,   pblht,   tpert,   qpert,   snowh,   &
                    fsns,    fsnt,    flns,    flnt,    state,   &
                    tend,    pbuf,    fsds,    landm,   &
                    landfrac,ocnfrac, icefrac, surface_state2d, srfflx_state2d )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Evaluate and apply physical processes that are calculated BEFORE 
! coupling to land, sea, and ice models.  
!
! Processes currently included are: 
! dry adjustment, moist convection, stratiform, wet deposition, radiation
!
! Pass surface fields for separate surface flux calculations
! Dump appropriate fields to history file.
! 
! Method: 
!
! Each parameterization should be implemented with this sequence of calls:
!  1)  Call physics interface
!  2)  Check energy
!  3)  Call physics_update
! See Interface to Column Physics and Chemistry Packages 
!   http://www.ccsm.ucar.edu/models/atm-cam/docs/phys-interface/index.html
! 
! Author: CCM1, CMS Contact: J. Truesdale
! 
!-----------------------------------------------------------------------

   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid
   use phys_buffer,     only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
   use stratiform,      only: stratiform_tend
   use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_update, physics_ptend_init
   use diagnostics,     only: diag_conv_tend_ini, diag_phys_writeout, diag_conv
   use history,         only: outfld
   use physconst,       only: cpair
   use constituents,    only: pcnst, pnats, ppcnst, qmin
   use convect_deep,    only: convect_deep_tend, convect_deep_tend_2
   use time_manager,    only: is_first_step, get_nstep
   use convect_shallow, only: convect_shallow_tend
   use check_energy,    only: check_energy_chng, check_energy_fix
   use check_energy,    only: check_tracers_data, check_tracers_init, check_tracers_chng
   use dycore,          only: dycore_is
   use aerosol_intr,    only: aerosol_wet_intr
   use comsrf,          only: surface_state, srfflx_state  ! note these are user defined types, not variables
   use radiation,       only: radiation_tend

   implicit none

#include <comctl.h>
!
! Arguments
!
   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)
   real(r8), intent(in) :: pblht(pcols)                ! Planetary boundary layer height
   real(r8), intent(in) :: tpert(pcols)                ! Thermal temperature excess
   real(r8), intent(inout) :: qpert(pcols,ppcnst)         ! Thermal humidity & constituent excess1
   real(r8), intent(in) :: snowh(pcols)                  ! Snow depth (liquid water equivalent)
   real(r8), intent(inout) :: fsns(pcols)                   ! Surface solar absorbed flux
   real(r8), intent(inout) :: fsnt(pcols)                   ! Net column abs solar flux at model top
   real(r8), intent(inout) :: flns(pcols)                   ! Srf longwave cooling (up-down) flux
   real(r8), intent(inout) :: flnt(pcols)                   ! Net outgoing lw flux at model top
   real(r8), intent(out) :: fsds(pcols)                   ! Surface solar down flux
   real(r8), intent(in) :: landm(pcols)                   ! land fraction ramp
   real(r8), intent(in) :: landfrac(pcols)                ! land fraction
   real(r8), intent(in) :: ocnfrac(pcols)                ! land fraction
   real(r8), intent(in) :: icefrac(pcols)                ! land fraction

   type(physics_state), intent(inout) :: state
   type(physics_tend ), intent(inout) :: tend
   type(pbuf_fld),      intent(inout), dimension(pbuf_size_max) :: pbuf
   type (surface_state), intent(inout) :: surface_state2d
   type (srfflx_state), intent(in)  :: srfflx_state2d

!
!---------------------------Local workspace-----------------------------
!

   type(physics_ptend)   :: ptend                  ! indivdual parameterization tendencies

   integer :: nstep                          ! current timestep number
   integer      lat(pcols)                   ! current latitudes(indices)
   integer      lon(pcols)                   ! current longtitudes(indices)

   real(r8) :: net_flx(pcols)

   real(r8) :: zdu(pcols,pver)               ! detraining mass flux from deep convection
   real(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables
   real(r8) :: cmfmc(pcols,pverp)            ! Convective mass flux--m sub c

   real(r8) cmfcme(pcols,pver)                ! cmf condensation - evaporation
   real(r8) cmfmc2(pcols,pverp)               ! Moist convection cloud mass flux
   real(r8) concld(pcols,pver)             
   real(r8) coszrs(pcols)                     ! Cosine solar zenith angle
   real(r8) dlf(pcols,pver)                   ! Detraining cld H20 from convection
   real(r8) pflx(pcols,pverp)                 ! Conv rain flux thru out btm of lev
   real(r8) rtdt                              ! 1./ztodt

   integer lchnk                              ! chunk identifier
   integer ncol                               ! number of atmospheric columns

   integer  i,k,m                             ! Longitude, level, constituent indices
   integer :: ixcldice, ixcldliq              ! constituent indices for cloud liquid and ice water.
                                           

   real(r8) dellow(pcols)                     ! delta p for bottom three levels of model
   real(r8) tavg(pcols)                       ! mass weighted average temperature for 

! physics buffer fields to compute tendencies for stratiform package
   integer itim, ifld
   real(r8), pointer, dimension(:,:) :: cld        ! cloud fraction


! physics buffer fields for total energy and mass adjustment
   real(r8), pointer, dimension(:  ) :: teout
   real(r8), pointer, dimension(:,:) :: qini
   real(r8), pointer, dimension(:,:) :: tini

   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble

! convective precipitation variables
   real(r8) :: prec_zmc(pcols)                ! total precipitation from ZM convection
   real(r8) :: snow_zmc(pcols)                ! snow from ZM convection
   real(r8) :: prec_cmf(pcols)                ! total precipitation from Hack convection
   real(r8) :: snow_cmf(pcols)                ! snow from Hack convection

! stratiform precipitation variables
   real(r8) :: prec_str(pcols)    ! sfc flux of precip from stratiform (m/s)
   real(r8) :: snow_str(pcols)     ! sfc flux of snow from stratiform   (m/s)
   real(r8) :: prec_pcw(pcols)     ! total precip from prognostic cloud scheme
   real(r8) :: snow_pcw(pcols)     ! snow from prognostic cloud scheme
   real(r8) :: prec_sed(pcols)     ! total precip from cloud sedimentation
   real(r8) :: snow_sed(pcols)     ! snow from cloud ice sedimentation


! energy checking variables
   real(r8) :: zero(pcols)                    ! array of zeros
   real(r8) :: rliq(pcols)                    ! vertical integral of liquid not yet in q(ixcldliq)
   real(r8) :: rliq2(pcols)                   ! vertical integral of liquid from shallow scheme
   real(r8) :: flx_cnd(pcols)
   real(r8) :: flx_heat(pcols)
   type(check_tracers_data):: tracerint             ! energy integrals and cummulative boundary fluxes

!-----------------------------------------------------------------------
    zero = 0.
!
   lchnk = state%lchnk
   ncol  = state%ncol

   rtdt = 1./ztodt

   nstep = get_nstep()


! Associate pointers with physics buffer fields
   itim = pbuf_old_tim_idx()
   ifld = pbuf_get_fld_idx('CLD')
   cld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

   ifld = pbuf_get_fld_idx('TEOUT')
   teout => pbuf(ifld)%fld_ptr(1,1:pcols,1,lchnk,itim)
   ifld = pbuf_get_fld_idx('QINI')
   qini  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld = pbuf_get_fld_idx('TINI')
   tini  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

   ifld = pbuf_get_fld_idx('FRACIS')
   fracis  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1:ppcnst)
!
! Set physics tendencies to 0
   tend %dTdt(:ncol,:pver)  = 0.
   tend %dudt(:ncol,:pver)  = 0.
   tend %dvdt(:ncol,:pver)  = 0.

   call physics_ptend_init (ptend) ! Initialize parameterization tendency structure
!
! Make sure that input tracers are all positive (probably unnecessary)
!
    
   call qneg3('TPHYSBCb',lchnk  ,ncol    ,pcols   ,pver    , &
              1, ppcnst,qmin  ,state%q )

! Save state for convective tendency calculations.
   call diag_conv_tend_ini(state)

   fracis (:ncol,:,1:ppcnst) = 1.
! compute mass integrals of input tracers state
   call check_tracers_init(state, tracerint)

!===================================================
! Global mean total energy fixer
!===================================================
   !*** BAB's FV heating kludge *** save the initial temperature
   tini(:ncol,:pver) = state%t(:ncol,:pver)
   if (dycore_is('LR')) then
      call check_energy_fix(state, ptend, nstep, flx_heat)
      call physics_update(state, tend, ptend, ztodt)
      call check_energy_chng(state, tend, "chkengyfix", nstep, ztodt, zero, zero, zero, flx_heat)
   end if
   qini(:ncol,:pver) = state%q(:ncol,:pver,1)

   call outfld('TEOUT', teout       , pcols, lchnk   )
   call outfld('TEINP', state%te_ini, pcols, lchnk   )
   call outfld('TEFIX', state%te_cur, pcols, lchnk   )
!
!===================================================
! Dry adjustment
! This code block is not a good example of interfacing a parameterization
!===================================================

! Copy state info for input to dadadj
! This is a kludge, so that dadadj does not have to be correctly reformulated in dry static energy

   ptend%s(:ncol,:pver)   = state%t(:ncol,:pver)
   ptend%q(:ncol,:pver,1) = state%q(:ncol,:pver,1)

   call dadadj (lchnk, ncol, state%pmid,  state%pint,  state%pdel,  &
                ptend%s, ptend%q(1,1,1))
   ptend%name  = 'dadadj'
   ptend%ls    = .TRUE.
   ptend%lq(1) = .TRUE.
   ptend%s(:ncol,:)   = (ptend%s(:ncol,:)   - state%t(:ncol,:)  )/ztodt * cpair
   ptend%q(:ncol,:,1) = (ptend%q(:ncol,:,1) - state%q(:ncol,:,1))/ztodt
   call physics_update (state, tend, ptend, ztodt)
!
!===================================================
! Moist convection
!===================================================
!
! Since the PBL doesn't pass constituent perturbations, they
! are zeroed here for input to the moist convection routine
!
   call convect_deep_tend(  prec_zmc,      &
        pblht,    cmfmc,      cmfcme,             &
        tpert,    dlf,        pflx,    zdu,       &
        rliq,    &
        ztodt,    snow_zmc,  &
        state,   ptend, pbuf ) 

   call physics_update(state, tend, ptend, ztodt)

! Check energy integrals, including "reserved liquid"
   flx_cnd(:ncol) = prec_zmc(:ncol) + rliq(:ncol)
   call check_energy_chng(state, tend, "zm_evap", nstep, ztodt, zero, flx_cnd, snow_zmc, zero)


!
! Call Hack (1994) convection scheme to deal with shallow/mid-level convection
!

   call convect_shallow_tend (ztodt   ,&
        qpert     ,   &
        pblht      ,     &
        cmfmc      ,cmfmc2  ,  prec_cmf,   &
        dlf        ,rliq      , rliq2, & 
        snow_cmf   , state, ptend,  pbuf       )
   

   call physics_update (state, tend, ptend, ztodt)

   flx_cnd(:ncol) = prec_cmf(:ncol) + rliq2(:ncol)
   call check_energy_chng(state, tend, "convect_shallow", nstep, ztodt, zero, flx_cnd, snow_cmf, zero)

!===================================================
! Calculate stratiform tendencey (sedimentation, detrain, cloud fraction and microphysics )
!===================================================

   call stratiform_tend(state, ptend, ztodt, &
        icefrac, landfrac, ocnfrac, &
        landm, snowh, & ! sediment
        dlf, & ! detrain
        rliq  , & ! check energy after detrain
        cmfmc,   cmfmc2,  concld,    &
        srfflx_state2d%ts,      srfflx_state2d%sst,        zdu,  &
        prec_str, snow_str, prec_sed, snow_sed, prec_pcw, snow_pcw, & 
        pbuf)

   call physics_update (state, tend, ptend, ztodt)
   call check_energy_chng(state, tend, "cldwat_tend", nstep, ztodt, zero, prec_str, snow_str, zero)

!===================================================
!  Aerosol wet chemistry determines scavenging fractions, and transformations
!
!  Then do convective transport of all trace species except water vapor and
!     cloud liquid and ice (we needed to do the scavenging first
!     to determine the interstitial fraction) 
!===================================================
   call aerosol_wet_intr (state, ptend, ztodt,  pbuf)
   call physics_update (state, tend, ptend, ztodt)

   call convect_deep_tend_2( state,   ptend,  ztodt,  pbuf ) 
   call physics_update (state, tend, ptend, ztodt)

   ! check tracer integrals
   call check_tracers_chng(state, tracerint, "cmfmca", nstep, ztodt, ptend%cflx_srf)

   !===================================================
   ! Moist physical parameteriztions complete: 
   ! send dynamical variables, and derived variables to history file
   !===================================================
   call diag_phys_writeout(state)
   call diag_conv(state, ztodt)

   !===================================================
   ! Radiation computations
   !===================================================
   call radiation_tend(state,ptend,pbuf, &
        surface_state2d, srfflx_state2d, &
        landfrac,landm,icefrac,snowh, &
        fsns,    fsnt, flns,    flnt,  &
        fsds, concld, net_flx)
   ! Set net flux used by spectral dycores
   do i=1,ncol
      tend%flx_net(i) = net_flx(i)
   end do
   call physics_update(state, tend, ptend, ztodt)
   call check_energy_chng(state, tend, "radheat", nstep, ztodt, zero, zero, zero, net_flx)

   ! Save atmospheric fields to force surface models
   call srfxfer (state,surface_state2d,prec_zmc,snow_zmc, &
        prec_cmf,snow_cmf,prec_sed,snow_sed, &
        prec_pcw,snow_pcw,fsns)

end subroutine tphysbc
