#include <misc.h>
#include <params.h>

module aerosol_intr

!---------------------------------------------------------------------------------
! Module to interface the aerosol parameterizations with CAM
! Phil Rasch, Jan 2003
!---------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,      only: pcols, pver
  use physconst,   only: mwdry, mwh2o, gravit
  use constituents,only: ppcnst, cnst_add, cnst_name, advected
  use abortutils, only: endrun

  implicit none

  private          ! Make default type private to the module
  save
!
! Public interfaces
!
  public aerosol_register_cnst                        ! register consituents
  public aerosol_implements_cnst                      ! returns true if consituent is implemented by this package
  public aerosol_init_cnst                            ! initialize mixing ratios if not read from initial file
  public aerosol_defaultopts                          ! get default runtime options
  public aerosol_setopts                              ! set runtime options
  public prognostic_aerosol_initialize                ! initialize (history) variables
  public aerosol_drydep_intr                          ! interface to dry deposition
  public aerosol_wet_intr                             ! interface to wet deposition
  public aerosol_srcsnk_intr                          ! interface to aerosol transformations
  public aerosol_emis_intr                            ! interface to surface emissions
  public aerosol_time_interp                          ! interpolate oxidant and fluxes to current time
  public set_aerosol_from_prognostics

!
! Private data
!

  character(len=16), private, parameter :: def_prognostic_sulfur = 'off'  ! default
! Set this flag to .TRUE. to turn on sulfur
  logical, private :: sulfur = .FALSE.
! Set this flag to .TRUE. to turn on sulfur feedback in 
! set_aerosol_from_prognostics()
  logical, private :: feedback_sulfur = .FALSE.

! Set this flag to .TRUE. to turn on carbon
  logical, private, parameter :: def_carbon = .FALSE.    ! default
  logical, private :: carbon = def_carbon
! Set this flag to .TRUE. to turn on carbon feedback in 
! set_aerosol_from_prognostics()
  logical, private, parameter :: def_feedback_carbon = .FALSE.  ! default
  logical, private :: feedback_carbon = def_feedback_carbon

! Set this flag to .TRUE. to turn on sea salt
  logical, private, parameter :: def_sea_salt = .FALSE.  ! default
  logical, private :: sea_salt = def_sea_salt
! Set this flag to .TRUE. to turn on sea salt feedback in 
! set_aerosol_from_prognostics()
  logical, private, parameter :: def_feedback_sea_salt = .FALSE.  ! default
  logical, private :: feedback_sea_salt = def_feedback_sea_salt



contains

!===============================================================================
  subroutine aerosol_register_cnst
!----------------------------------------------------------------------- 
! 
! Purpose: register aerosols
! 
! Method: 

! Author: P. J. Rasch
! 
!-----------------------------------------------------------------------
    use sulfur_intr, only: sulfur_register_cnst
    use carbon_intr, only: carbon_register_cnst
#if ( defined DUST )
    use dust_intr,   only: dust_register_cnst
#endif
    use seasalt_intr,only: seasalt_register_cnst
    use phys_buffer,  only: pbuf_times, pbuf_add

    implicit none
!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index
    integer idx   ! pbuf idx

!-----------------------------------------------------------------------

! Set names of sulfur species and declare them as history variables
    if ( sulfur ) then
      call sulfur_register_cnst
    endif

    if ( carbon ) then
! Set names of carbon species and declare them as history variables
      call carbon_register_cnst
    endif

#if ( defined DUST )
! Set names of dust species and declare them as history variables
    call dust_register_cnst
#endif

    if ( sea_salt ) then
! Set names of sea salt species and declare them as history variables
      call seasalt_register_cnst
    endif

! Request physics buffer space for fields that don't persist across timesteps.
    call pbuf_add('FRACIS' , 'physpkg', 1,pver,ppcnst, idx)


    return
  end subroutine aerosol_register_cnst


!=======================================================================
  function aerosol_implements_cnst(name)
!----------------------------------------------------------------------- 
! 
! Purpose: return true if specified constituent is implemented by this 
!          package
! 
! Author: T. Henderson
! 
!-----------------------------------------------------------------------
    use sulfur_intr,  only: sulfur_implements_cnst
    use carbon_intr,  only: carbon_implements_cnst
#if ( defined DUST )
    use dust_intr,    only: dust_implements_cnst
#endif
    use seasalt_intr, only: seasalt_implements_cnst
     implicit none
!-----------------------------Arguments---------------------------------

     character(len=*), intent(in) :: name   ! constituent name
     logical :: aerosol_implements_cnst        ! return value
!---------------------------Local workspace-----------------------------
     integer :: m
!-----------------------------------------------------------------------

     aerosol_implements_cnst = .false.

     if ( sulfur ) then
       aerosol_implements_cnst = &
         (aerosol_implements_cnst.OR.sulfur_implements_cnst (name))
     endif

     if ( carbon ) then
       aerosol_implements_cnst = &
         (aerosol_implements_cnst.OR.carbon_implements_cnst (name))
     endif

#if ( defined DUST )
     aerosol_implements_cnst = &
       (aerosol_implements_cnst.OR.dust_implements_cnst (name))
#endif

     if ( sea_salt ) then
       aerosol_implements_cnst = &
         (aerosol_implements_cnst.OR.seasalt_implements_cnst (name))
     endif

  end function aerosol_implements_cnst


!=======================================================================
  subroutine aerosol_init_cnst(name, q)
!----------------------------------------------------------------------- 
! 
! Purpose:  
! Set initial mass mixing ratios.  
!
!-----------------------------------------------------------------------
    use sulfur_intr,  only: sulfur_implements_cnst, sulfur_init_cnst
    use carbon_intr,  only: carbon_implements_cnst, carbon_init_cnst
#if ( defined DUST )
    use dust_intr,    only: dust_implements_cnst, dust_init_cnst
#endif
    use seasalt_intr, only: seasalt_implements_cnst, seasalt_init_cnst
    implicit none
!-----------------------------Arguments---------------------------------

    character(len=*), intent(in) :: name         ! constituent name
    real(r8), intent(out) :: q(:,:,:)            !  mass mixing ratio
!-----------------------------------------------------------------------

    if ( sulfur ) then
      if (sulfur_implements_cnst(name)) then
        call sulfur_init_cnst(name, q)
      endif
    endif

    if ( carbon ) then
      if (carbon_implements_cnst(name)) then
        call carbon_init_cnst(name, q)
      endif
    endif

#if ( defined DUST )
    if (dust_implements_cnst(name)) then
      call dust_init_cnst(name, q)
    endif
#endif

    if ( sea_salt ) then
      if (seasalt_implements_cnst(name)) then
        call seasalt_init_cnst(name, q)
      endif
    endif

  end subroutine aerosol_init_cnst



!=======================================================================
  subroutine aerosol_defaultopts(prognostic_sulfur_out,    &
                                 aero_carbon_out,          &
                                 aero_feedback_carbon_out, &
                                 aero_sea_salt_out,        &
                                 aero_feedback_sea_salt_out)
!----------------------------------------------------------------------- 
!
! Purpose: Return default runtime options
!
! Author: Tom Henderson
!
!-----------------------------------------------------------------------
!
! Arguments:
!
    ! prognostic aerosol sulfur option
    character(len=16), intent(out), optional :: prognostic_sulfur_out
    ! prognostic aerosol carbon option
    logical, intent(out), optional :: aero_carbon_out
    ! prognostic aerosol feedback option for carbon
    logical, intent(out), optional :: aero_feedback_carbon_out
    ! prognostic aerosol sea salt option
    logical, intent(out), optional :: aero_sea_salt_out
    ! prognostic aerosol feedback option for sea salt
    logical, intent(out), optional :: aero_feedback_sea_salt_out
!-----------------------------------------------------------------------
    if ( present(prognostic_sulfur_out) ) then
      prognostic_sulfur_out = def_prognostic_sulfur
    endif

    if ( present(aero_carbon_out) ) then
      aero_carbon_out = def_carbon
    endif
    if ( present(aero_feedback_carbon_out) ) then
      aero_feedback_carbon_out = def_feedback_carbon
    endif

    if ( present(aero_sea_salt_out) ) then
      aero_sea_salt_out = def_sea_salt
    endif
    if ( present(aero_feedback_sea_salt_out) ) then
      aero_feedback_sea_salt_out = def_feedback_sea_salt
    endif
  end subroutine aerosol_defaultopts



!=======================================================================
  subroutine aerosol_setopts(prognostic_sulfur_in,    &
                             aero_carbon_in,          &
                             aero_feedback_carbon_in, &
                             aero_sea_salt_in,        &
                             aero_feedback_sea_salt_in)
!----------------------------------------------------------------------- 
!
! Purpose: Set runtime options and print them
!
! Author: Tom Henderson
!
!-----------------------------------------------------------------------
    use pmgrid, only: masterproc
!-----------------------------------------------------------------------
!
! Arguments:
!
    ! prognostic aerosol sulfur option
    character(len=16), intent(in), optional :: prognostic_sulfur_in
    ! prognostic aerosol carbon option
    logical, intent(in), optional :: aero_carbon_in
    ! prognostic aerosol feedback option for carbon
    logical, intent(in), optional :: aero_feedback_carbon_in
    ! prognostic aerosol sea salt option
    logical, intent(in), optional :: aero_sea_salt_in
    ! prognostic aerosol feedback option for sea salt
    logical, intent(in), optional :: aero_feedback_sea_salt_in
!-----------------------------------------------------------------------
!
! Local variables:
!
    character(len=8) :: onoff
!-----------------------------------------------------------------------
    if ( present(prognostic_sulfur_in) ) then
      if ( prognostic_sulfur_in== "off") then
        sulfur = .FALSE.
        feedback_sulfur = .FALSE.
      elseif( prognostic_sulfur_in == "passive") then
        sulfur = .TRUE.
        feedback_sulfur = .FALSE.
      elseif( prognostic_sulfur_in== "direct") then
        sulfur = .TRUE.
        feedback_sulfur = .TRUE.
      else
        write(6,*) "prognostic_sulfur must be one of off, passive or direct"
        call endrun
      endif
    endif
    if ( present(aero_carbon_in) ) then
      carbon = aero_carbon_in
    endif
    if ( present(aero_feedback_carbon_in) ) then
      feedback_carbon = aero_feedback_carbon_in
    endif
    if ( present(aero_sea_salt_in) ) then
      sea_salt = aero_sea_salt_in
    endif
    if ( present(aero_feedback_sea_salt_in) ) then
      feedback_sea_salt = aero_feedback_sea_salt_in
    endif
    ! error check...
    if (feedback_carbon.and.(.not.carbon)) then
      write(6,*) 'ERROR:  AEROSOL_SETOPTS:  Feedback of prognostic'
      write(6,*) 'carbon aerosols is only allowed when these'
      write(6,*) 'prognostic aerosols are enabled.' 
      call endrun
    endif
    if (feedback_sea_salt.and.(.not.sea_salt)) then
      write(6,*) 'ERROR:  AEROSOL_SETOPTS:  Feedback of prognostic'
      write(6,*) 'sea salt aerosols is only allowed when these'
      write(6,*) 'prognostic aerosols are enabled.' 
      call endrun
    endif
    if (masterproc) then
      ! TBH:  refactor this mess...  
      !----------------------------------
      if (sulfur) then
        onoff = 'on      '
      else
        onoff = 'off     '
      endif
      write(6,*) 'AEROSOL_SETOPTS:  prognostic sulfur aerosols are ', &
        onoff
      if (feedback_sulfur) then
        onoff = 'enabled '
      else
        onoff = 'disabled'
      endif
      write(6,*) 'AEROSOL_SETOPTS:  feedback of prognostic sulfur aerosols is ', &
        onoff
      !----------------------------------
      if (carbon) then
        onoff = 'on      '
      else
        onoff = 'off     '
      endif
      write(6,*) 'AEROSOL_SETOPTS:  prognostic carbon aerosols are ', &
        onoff
      if (feedback_carbon) then
        onoff = 'enabled '
      else
        onoff = 'disabled'
      endif
      write(6,*) 'AEROSOL_SETOPTS:  feedback of prognostic carbon aerosols is ', &
        onoff
      !----------------------------------
      if (sea_salt) then
        onoff = 'on      '
      else
        onoff = 'off     '
      endif
      write(6,*) 'AEROSOL_SETOPTS:  prognostic sea salt aerosols are ', &
        onoff
      if (feedback_sea_salt) then
        onoff = 'enabled '
      else
        onoff = 'disabled'
      endif
      write(6,*) 'AEROSOL_SETOPTS:  feedback of prognostic sea salt aerosols is ', &
        onoff
    endif
  end subroutine aerosol_setopts



!===============================================================================
  subroutine prognostic_aerosol_initialize
!----------------------------------------------------------------------- 
! 
! Purpose: initialize aerosol parameterizations
!          (declare history variables)
! 
! Method: 
! call subroutines
! 
! Author: Phil Rasch
! 
!-----------------------------------------------------------------------
    use sulfur_intr, only: sulfur_initialize
    use carbon_intr, only: carbon_initialize
#if ( defined DUST )
    use dust_intr, only: dust_initialize
#endif
    use seasalt_intr, only: seasalt_initialize

    implicit none

    if ( sulfur ) then
      call sulfur_initialize
    endif

    if ( carbon ) then
      call carbon_initialize
    endif

#if ( defined DUST )
    call dust_initialize
#endif

    return
  end subroutine prognostic_aerosol_initialize

!===============================================================================
  subroutine aerosol_wet_intr (state, ptend, dt,  pbuf)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to wet processing of all aerosols
! 
! Method: 
!  use a modified version of the scavenging parameterization described in
!     Barth et al, 2000, JGR (sulfur cycle paper)
!     Rasch et al, 2001, JGR (INDOEX paper)
! 
! Author: Phil Rasch
! 
!-----------------------------------------------------------------------
    use history,       only: outfld
    use physics_types, only: physics_state, physics_ptend
    use phys_buffer,     only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx

    use phys_grid,     only: get_lat_all_p, get_lon_all_p, get_rlat_all_p
    use sulfur_intr,   only: sulfur_wet_intr
    use carbon_intr,   only: carbon_wet_intr
#if ( defined DUST )
    use dust_intr,     only: dust_wet_intr
#endif
    use wetdep,        only: clddiag
    use time_manager, only: get_curr_date, get_perp_date, get_curr_calday
    use time_manager, only: is_perpetual, get_nstep
    use scyc,       only: scyc_idx1
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
!
    real(r8),            intent(in)  :: dt             ! time step
    type(physics_state), intent(in ) :: state          ! Physics state variables

    type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies
!    real(r8), intent(inout)  :: cflx(pcols,ppcnst)       ! Surface constituent flux (kg/m^2/s)

!    real(r8), intent(inout) :: fracis(pcols,pver,ppcnst)         ! fraction of transported species that are insoluble
   type(pbuf_fld),      intent(inout), dimension(pbuf_size_max) :: pbuf



! local vars
    real(r8) :: calday        ! current calendar day
    integer :: yr, mon, day, ncsec
    integer :: ncdate
    real(r8):: rainmr(pcols,pver)  ! mixing ratio of rain within cloud volume
    real(r8):: cldv(pcols,pver)     ! cloudy volume undergoing wet chem and scavenging
    integer ix
    integer m
    integer  lat(pcols)                  ! latitude indices
    real(r8) clat(pcols)                    ! latitudes
    integer  lon(pcols)                  ! longtitude indices
    integer  nstep
    real(r8) conicw(pcols,pver)          ! convective in-cloud water
    real(r8) cmfdqr(pcols,pver)          ! convective production of rain
    real(r8) cldc(pcols,pver)            ! convective cloud fraction, currently empty


! physics buffer 
    integer itim, ifld
    real(r8), pointer, dimension(:,:) :: cldn       ! cloud fraction
    real(r8), pointer, dimension(:,:) :: cme
    real(r8), pointer, dimension(:,:) :: prain
    real(r8), pointer, dimension(:,:) :: evapr
    real(r8), pointer, dimension(:,:) ::  icwmrdp ! in cloud water mixing ratio, deep convection
    real(r8), pointer, dimension(:,:) ::  rprddp ! rain production, deep convection
    real(r8), pointer, dimension(:,:) ::  icwmrsh ! in cloud water mixing ratio, deep convection
    real(r8), pointer, dimension(:,:) ::  rprdsh ! rain production, deep convection
    real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble


    integer ncol,lchnk

    nstep = get_nstep()

    calday = get_curr_calday()
    if ( is_perpetual() ) then
       call get_perp_date(yr, mon, day, ncsec)
    else
       call get_curr_date(yr, mon, day, ncsec)
    end if
    ncdate = yr*10000 + mon*100 + day

    ncol = state%ncol
    lchnk = state%lchnk

   call get_lat_all_p(lchnk, ncol, lat)
   call get_lon_all_p(lchnk, ncol, lon)
   call get_rlat_all_p(lchnk, ncol, clat)

! Associate pointers with physics buffer fields
   itim = pbuf_old_tim_idx()
   ifld = pbuf_get_fld_idx('CLD')
   cldn => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

   ifld = pbuf_get_fld_idx('QME')
   cme  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld = pbuf_get_fld_idx('PRAIN')
   prain  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld = pbuf_get_fld_idx('NEVAPR')
   evapr  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

   ifld = pbuf_get_fld_idx('FRACIS')
   fracis  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1:ppcnst)



   ifld = pbuf_get_fld_idx('ICWMRDP')
   icwmrdp  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld = pbuf_get_fld_idx('RPRDDP')
   rprddp  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

   ifld = pbuf_get_fld_idx('ICWMRSH')
   icwmrsh  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld = pbuf_get_fld_idx('RPRDSH')
   rprdsh  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

   ! sum deep and shallow convection contributions
   conicw(:ncol,:) = icwmrdp(:ncol,:) + icwmrsh(:ncol,:)
   cmfdqr(:ncol,:) = rprddp(:ncol,:)  + rprdsh(:ncol,:)


!   fields needed for wet scavenging
    call clddiag( state%t, state%pmid, state%pdel, cmfdqr, cldn, cme, evapr, prain, &
         cldv, rainmr, ncol )

    if ( sulfur ) then
!   wet chemistry and scavenging of fields for sulfur cycle
      call sulfur_wet_intr (state, ptend, nstep, dt, lat, clat, lon, cme, prain, &
         evapr, cldv, cldc, cldn, fracis, calday, cmfdqr, conicw, rainmr)
    endif

    if ( carbon ) then
!   wet scavenging for carbon
      call carbon_wet_intr (state, ptend, nstep, dt, lat, clat, cme, prain, &
         evapr, cldv, cldc, cldn, fracis, calday, cmfdqr, conicw, rainmr)
    endif

#if ( defined DUST )
!   wet scavenging for dust
    call dust_wet_intr (state, ptend, nstep, dt, lat, clat, cme, prain, &
       evapr, cldv, cldc, cldn, fracis, calday, cmfdqr, conicw, rainmr)
#endif

    return

  end subroutine aerosol_wet_intr

  subroutine aerosol_drydep_intr (state, ptend, wvflx, dt, &
       fsds, obklen, ts, ustar, prect, snowh, pblh, hflx, landfrac, icefrac, &
       ocnfrac, fv,ram1)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to dry deposition parameterizions and sedimentation of all aerosols
! 
! Method: 
! Use prescribed dry deposition velocities for sulfur and carbon
! Use calculated dry dep velocities from CLM for dust
! 
! Author: Phil Rasch
! 
!-----------------------------------------------------------------------
    use history,       only: outfld
    use physics_types, only: physics_state, physics_ptend
    use phys_grid,     only: get_lat_all_p, get_rlat_all_p
    use sulfur_intr,   only: sulfur_drydep_intr
    use carbon_intr,   only: carbon_drydep_intr
#if ( defined DUST )
    use dust_intr,     only: dust_drydep_intr
#endif
    use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
                            is_perpetual
    use scyc,       only: scyc_idx1
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
    real(r8),            intent(in)  :: dt             ! time step
    type(physics_state), intent(in ) :: state          ! Physics state variables
    real(r8), intent(in) :: fsds(pcols)                 ! longwave down at sfc
    real(r8), intent(in) :: obklen(pcols)                 ! longwave down at sfc
    real(r8), intent(in) :: ustar(pcols)                  ! sfc fric vel
    real(r8), intent(in) :: ts(pcols)                     ! sfc temp
    real(r8), intent(in) :: landfrac(pcols)               ! land fraction
    real(r8), intent(in) :: icefrac(pcols)               ! ice fraction
    real(r8), intent(in) :: ocnfrac(pcols)               ! ocn fraction
    real(r8), intent(in) :: hflx(pcols)                  ! sensible heat flux
    real(r8), intent(in) :: prect(pcols)                 ! prect
    real(r8), intent(in) :: snowh(pcols)                 ! snow depth
    real(r8), intent(in) :: pblh(pcols)                  ! pbl height
!   real(r8), intent(in)  :: wvflx(pcols,ppcnst)         ! water vapor flux
    real(r8), intent(in)  :: wvflx(pcols)         ! water vapor flux
    real(r8), intent(in)  :: fv(pcols)         !  friction velocity
    real(r8), intent(in) :: ram1(pcols)        ! aerodynamic resistance from land model

    type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies

    integer  lat(pcols)                  ! latitude index for S->N storage
    real(r8) clat(pcols)                 ! latitude 
    integer lchnk
    integer ncol
    integer ix
    integer m

    integer :: yr, mon, day, ncsec
    integer :: ncdate

    real(r8) fvtmp(pcols)         ! temporary friction velocity
    real(r8) ram1tmp(pcols)       ! temporary aerodynamic resistance from land model

    if ( is_perpetual() ) then
       call get_perp_date(yr, mon, day, ncsec)
    else
       call get_curr_date(yr, mon, day, ncsec)
    end if

    lchnk = state%lchnk
    ncol = state%ncol
    
    call get_lat_all_p(lchnk, ncol, lat)
    call get_rlat_all_p(lchnk, ncol, clat)

!   note that tendencies are not only in sfc layer (because of sedimentation)
!   and that ptend is updated within each subroutine for different species

    if ( sulfur ) then
      call sulfur_drydep_intr (state, ptend, wvflx, dt, lat, clat, &
         fsds, obklen, ts, ustar, prect, snowh, pblh, hflx, mon, landfrac, icefrac, ocnfrac)
    endif

    if ( carbon ) then
      call carbon_drydep_intr (state, ptend, wvflx, dt, lat, clat, &
         fsds, obklen, ts, ustar, prect, snowh, pblh, hflx, mon, landfrac, icefrac, ocnfrac)
    endif

!    fvtmp = 0.1
!    ram1tmp = 0.1
#if ( defined DUST )
    call dust_drydep_intr (state, ptend, wvflx, dt, lat, clat, &
       fsds, obklen, ts, ustar, prect, snowh, pblh, hflx, mon, landfrac, icefrac, ocnfrac,fv,ram1)
#endif

!    ptend%name  = 'aerosol_dry_intr'
!    ptend%ls    = .FALSE.
!    ptend%lu    = .FALSE.
!    ptend%lv    = .FALSE.
!    ix = scyc_idx1()
!    do m=2,ppcnst
!       if (m >= ix .and. m <= ix+3) ptend%lq(m) = .true.
!    end do

    return
  end subroutine aerosol_drydep_intr

  subroutine aerosol_time_interp

!-----------------------------------------------------------------------
! Purpose: 
! interpolate the emission inventory to correct point in time
! 
! Method: 
! read in data from netcdf files for surrounding time periods
! do linear interpolation
! 
! Author: Phil Rasch
! 
!-----------------------------------------------------------------------

    use sulfur_intr, only: sulfur_time_interp
    use carbon_intr, only: carbon_time_interp


    if ( sulfur ) then
      call sulfur_time_interp
    endif

    if ( carbon ) then
      call carbon_time_interp
    endif

  end subroutine aerosol_time_interp

  subroutine aerosol_emis_intr (state, ptend, cflx, dt)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! return surface fluxes of aerosol species and tendencies in surface layer
! due to surface emissions
! 
! Method: 
! 
! Author: Phil Rasch
! 
!-----------------------------------------------------------------------
    use history,       only: outfld
    use physics_types, only: physics_state, physics_ptend
    use phys_grid,     only: get_lon_all_p, get_lat_all_p, get_rlat_all_p
    use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
                            is_perpetual
    use sulfur_intr,   only: sulfur_emis_intr
!    use scyc, only: doscyc, scyc_idx1, add_volc_emis
!    use sulemis, only: addsulemis
!    use volcemission, only: volcemist

    use carbon_intr, only: carbon_emis_intr

#if ( defined DUST )
    use dust_intr,     only: dust_emis_intr
#endif
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
    real(r8),            intent(in)  :: dt             ! time step
    type(physics_state), intent(in ) :: state          ! Physics state variables

    type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies
    real(r8), intent(inout)  :: cflx(pcols,ppcnst)     ! Surface constituent flux (kg/m^2/s)

    integer  lat(pcols)                  ! latitude index 
    integer  lon(pcols)                  ! longitude index
    real(r8) clat(pcols)                 ! latitude 
    integer lchnk
    integer ncol
    integer i
    integer m
    real astmp(pcols,pver,3)
!
    real(r8) :: calday        ! current calendar day
    integer :: yr, mon, day, ncsec
    integer :: ncdate
    real(r8) :: so2sf(pcols), so4sf(pcols), dmssf(pcols)
    integer :: k

    if ( sulfur ) then
      call sulfur_emis_intr(state, ptend, dt)
    endif

    if ( carbon ) then
      call carbon_emis_intr(state, ptend, dt)
    endif

#if ( defined DUST )
    call dust_emis_intr(state,ptend,dt,cflx)
#endif

    return
  end subroutine aerosol_emis_intr 

  subroutine aerosol_srcsnk_intr (state, ptend, dt, ocnfrc)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! convert aerosols from one category to another by non-wet processes
! gas phase chemistry could go here also
! 

! Author: Phil Rasch
! 
!-----------------------------------------------------------------------
    use history,       only: outfld
    use physics_types, only: physics_state, physics_ptend
    use phys_grid,     only: get_lon_all_p, get_lat_all_p, get_rlat_all_p
    use caer,          only: caer_idx1, caercv
    use seasalt_intr,  only: seasalt_srcsnk
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
    real(r8),            intent(in)  :: dt             ! time step
    type(physics_state), intent(in ) :: state          ! Physics state variables
    real(r8),            intent(in)  :: ocnfrc(pcols)  !

    type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies

    integer i
    integer m

    if ( carbon ) then
      i = caer_idx1()-1
      ! conversion from hydrophobic to hydrophilic form
      ptend%lq(i+1:i+4) = .true.
      do m = 1,2
         call caercv(state%ncol, state%q(:,:,i+m), dt, ptend%q(:,:,i+m), ptend%q(:,:,i+m+2))
      end do
    endif

    if ( sea_salt ) then
      call seasalt_srcsnk(state, dt, state%u(:,pver), ocnfrc, ptend)
    endif

    return
  end subroutine aerosol_srcsnk_intr 

  subroutine set_aerosol_from_prognostics(state, aeromix)

    use prescribed_aerosols, only: naer_all, idxSUL, idxSSLT, idxOCPHO, idxBCPHO, idxOCPHI, idxBCPHI
    use prescribed_aerosols, only: idxDUSTfirst, numDUST,  idxBG
    use constituents, only: ppcnst, cnst_name
    use physics_types, only: physics_state

    implicit none

    type(physics_state), intent(in ) :: state                            ! Physics state variables
    real(r8), intent(inout)          :: aeromix(pcols, pver, naer_all)   ! aerosol mass mixing ratios
    integer m
    integer ncol

    ncol = state%ncol

! TBH:  Commented out the following line at request of Andrew Conley.  9/17/03
!    aeromix(:ncol,:,idxBG) = 0 ! background aerosol to zero
    do m = 1,ppcnst
!       write (6,*) '  set_aerosol_from_prognostics >', trim(cnst_name(m)),'<'
       if (feedback_sulfur) then    
          if (trim(cnst_name(m)) == 'SO4') then
!          write (6,*) ' sulfate found '
             aeromix(:ncol,:,idxSUL) = state%q(:ncol,:,m)
          endif
       endif
       if (feedback_carbon) then    
          if (trim(cnst_name(m)) == 'BCPHO') then
!          write (6,*) ' BCPHO found '
             aeromix(:ncol,:,idxBCPHO) = state%q(:ncol,:,m)
          else if (trim(cnst_name(m)) == 'BCPHI') then
!          write (6,*) ' BCPHI found '
             aeromix(:ncol,:,idxBCPHI) = state%q(:ncol,:,m)
          else if (trim(cnst_name(m)) == 'OCPHO') then
!          write (6,*) ' OCPHO found '
             aeromix(:ncol,:,idxOCPHO) = state%q(:ncol,:,m)
          else if (trim(cnst_name(m)) == 'OCPHI') then
!          write (6,*) ' OCPHI found '
             aeromix(:ncol,:,idxOCPHI) = state%q(:ncol,:,m)
          endif
       endif
       if (feedback_sea_salt) then    
          if (trim(cnst_name(m)) == 'SSLT') then
!          write (6,*) ' SSLT found '
             aeromix(:ncol,:,idxSSLT) = state%q(:ncol,:,m)
          endif
       endif
#if ( defined DUST )
       if (feedback_dust) then    
          if (trim(cnst_name(m)) == 'DSTX01') then
!          write (6,*) ' DSTX01 found '
             aeromix(:ncol,:,idxDUSTfirst) = state%q(:ncol,:,m)
          else if (trim(cnst_name(m)) == 'DSTX02') then
!          write (6,*) ' DSTX02 found '
             aeromix(:ncol,:,idxDUSTfirst+1) = state%q(:ncol,:,m)
          else if (trim(cnst_name(m)) == 'DSTX03') then
!          write (6,*) ' DSTX03 found '
             aeromix(:ncol,:,idxDUSTfirst+2) = state%q(:ncol,:,m)
          else if (trim(cnst_name(m)) == 'DSTX04') then
!          write (6,*) ' DSTX04 found '
             aeromix(:ncol,:,idxDUSTfirst+3) = state%q(:ncol,:,m)
          endif
       endif
#endif
    end do

  end subroutine set_aerosol_from_prognostics

end module aerosol_intr
