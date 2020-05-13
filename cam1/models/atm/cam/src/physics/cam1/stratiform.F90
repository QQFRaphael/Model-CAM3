#undef DEBUG
#include <misc.h>

module stratiform

!---------------------------------------------------------------------------------
! Purpose:
!
! Provides the CAM interface to the prognostic cloud water and ice parameterization
!
! Author: Byron Boville  Sept 04, 2002
!         modified by D.B. Coleman May 2004
!
!---------------------------------------------------------------------------------

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use ppgrid,        only: pcols, pver
  use physconst,     only: gravit, latvap, latice
  use abortutils,    only: endrun

  implicit none
  private
  save

  public :: stratiform_register, stratiform_init_cnst, stratiform_implements_cnst
  public :: stratiform_init
  public :: stratiform_tend

! Private module data

  integer, parameter :: ncnst=2                      ! number of constituents
  character(len=8), dimension(ncnst), parameter :: & ! constituent names
     cnst_names = (/'CLDLIQ', 'CLDICE'/)

  integer :: &
       ixcldice,     &! cloud ice water index
       ixcldliq,     &! cloud liquid water index
       qcwat_idx,    &! qcwat index in physics buffer
       lcwat_idx,    &! lcwat index in physics buffer
       tcwat_idx,    &! tcwat index in physics buffer
       cld_idx        ! cld index in physics buffer

contains


!===============================================================================

  subroutine stratiform_register
!----------------------------------------------------------------------
!
! Register the constituents (cloud liquid and cloud ice) and the fields
! in the physics buffer.
! 
!-----------------------------------------------------------------------
    use constituents, only: cnst_add, advected, nonadvec
    use physconst,    only: mwdry, cpair
    use phys_buffer,  only: pbuf_times, pbuf_add

    implicit none

    logical, parameter :: cldw_adv=.true.  ! true => cloud water is treated as advected tracer

    integer flag
!-----------------------------------------------------------------------

! Register cloud water and determine index (either advected or non-adv).
    if (cldw_adv) then
       flag = advected
    else
       flag = nonadvec
    endif
    call cnst_add(cnst_names(1), flag, mwdry, cpair, 0._r8, ixcldliq, &
         longname='Grid box averaged liquid condensate amount')
    call cnst_add(cnst_names(2), flag, mwdry, cpair, 0._r8, ixcldice, &
         longname='Grid box averaged ice condensate amount')

! Request physics buffer space for fields that persist across timesteps.
    call pbuf_add('QCWAT', 'global', 1,pver,pbuf_times, qcwat_idx)
    call pbuf_add('LCWAT', 'global', 1,pver,pbuf_times, lcwat_idx)
    call pbuf_add('TCWAT', 'global', 1,pver,pbuf_times, tcwat_idx)
    call pbuf_add('CLD',   'global', 1,pver,pbuf_times, cld_idx)

    call pbuf_add('QINI' , 'physpkg', 1,pver,      1, flag)
    call pbuf_add('TINI' , 'physpkg', 1,pver,      1, flag)

    call pbuf_add('QME' , 'physpkg', 1,pver,      1, flag)
    call pbuf_add('PRAIN' , 'physpkg', 1,pver,      1, flag)
    call pbuf_add('NEVAPR' , 'physpkg', 1,pver,      1, flag)


    call pbuf_add('REI' , 'physpkg', 1,pver,      1, flag)
    call pbuf_add('REL' , 'physpkg', 1,pver,      1, flag)


  end subroutine stratiform_register

!===============================================================================

  function stratiform_implements_cnst(name)
!----------------------------------------------------------------------- 
! 
! Purpose: return true if specified constituent is implemented by this package
! 
! Author: B. Eaton
! 
!-----------------------------------------------------------------------
     implicit none
!-----------------------------Arguments---------------------------------

     character(len=*), intent(in) :: name   ! constituent name
     logical :: stratiform_implements_cnst     ! return value
!---------------------------Local workspace-----------------------------
     integer :: m
!-----------------------------------------------------------------------

     stratiform_implements_cnst = .false.
     do m = 1, ncnst
        if (name == cnst_names(m)) then
           stratiform_implements_cnst = .true.
           return
        end if
     end do
  end function stratiform_implements_cnst

!===============================================================================
  subroutine stratiform_init_cnst(name, q)
!----------------------------------------------------------------------- 
! 
! Initialize the cloud water mixing ratios (liquid and ice), if they are
! not read from the initial file
! 
!-----------------------------------------------------------------------
    use pmgrid,        only: plon, plev, plat

    implicit none

! Arguments
    character(len=*), intent(in)  :: name                ! constituent name
    real(r8),         intent(out) :: q(plon,plev,plat)   ! mass mixing ratio
!-----------------------------------------------------------------------

    if ( name == 'CLDLIQ' ) then
       q = 0.0
       return
    else if ( name == 'CLDICE' ) then
       q = 0.0
       return
    end if

  end subroutine stratiform_init_cnst

!===============================================================================
  subroutine stratiform_init()
!----------------------------------------------------------------------
!
! Initialize the cloud water parameterization
! 
!-----------------------------------------------------------------------
    use cldwat,        only: inimc
    use constituents,  only: cnst_get_ind, cnst_name, cnst_longname, sflxnam
    use history,       only: addfld, add_default, phys_decomp
    use physconst,     only: tmelt, rh2o, rhodair

    integer :: m, mm
!-----------------------------------------------------------------------

    ! initialization routine for prognostic cloud water
    call inimc (tmelt, rhodair/1000.0, gravit, rh2o )

    ! register history variables
    do m = 1,ncnst
       call cnst_get_ind(cnst_names(m), mm)
       call addfld (cnst_name(mm), 'kg/kg   ', pver, 'A', cnst_longname(mm), phys_decomp)
       call addfld (sflxnam(mm),   'kg/m2/s ',    1, 'A', trim(cnst_name(mm))//' surface flux', phys_decomp)
       call add_default (cnst_name(mm), 1, ' ')
       call add_default (sflxnam(mm),   1, ' ')
    end do

    call addfld ('FWAUT   ','fraction',pver, 'A','Relative importance of liquid autoconversion' ,phys_decomp)
    call addfld ('FSAUT   ','fraction',pver, 'A','Relative importance of ice autoconversion'    ,phys_decomp)
    call addfld ('FRACW   ','fraction',pver, 'A','Relative  importance of rain accreting liquid',phys_decomp)
    call addfld ('FSACW   ','fraction',pver, 'A','Relative  importance of snow accreting liquid',phys_decomp)
    call addfld ('FSACI   ','fraction',pver, 'A','Relative  importance of snow accreting ice'   ,phys_decomp)
    call addfld ('CME     ','kg/kg/s ',pver, 'A','Rate of cond-evap within the cloud'           ,phys_decomp)
    call addfld ('CMEICE  ','kg/kg/s ',pver, 'A','Rate of cond-evap of ice within the cloud'    ,phys_decomp)
    call addfld ('CMELIQ  ','kg/kg/s ',pver, 'A','Rate of cond-evap of liq within the cloud'    ,phys_decomp)
    call addfld ('ICE2PR  ','kg/kg/s ',pver, 'A','Rate of conversion of ice to precip'          ,phys_decomp)
    call addfld ('LIQ2PR  ','kg/kg/s ',pver, 'A','Rate of conversion of liq to precip'          ,phys_decomp)
    call addfld ('ZMDLF   ','kg/kg/s ',pver, 'A','Detrained liquid water from ZM convection'    ,phys_decomp)
    call addfld ('PRODPREC','kg/kg/s ',pver, 'A','Rate of conversion of condensate to precip'   ,phys_decomp)
    call addfld ('EVAPPREC','kg/kg/s ',pver, 'A','Rate of evaporation of falling precip'        ,phys_decomp)
    call addfld ('EVAPSNOW','kg/kg/s ',pver, 'A','Rate of evaporation of falling snow'          ,phys_decomp)
    call addfld ('HPROGCLD','W/kg'    ,pver, 'A','Heating from prognostic clouds'               ,phys_decomp)
    call addfld ('HCME    ','W/kg'    ,pver, 'A','Heating from cond-evap within the cloud'      ,phys_decomp)
    call addfld ('HEVAP   ','W/kg'    ,pver, 'A','Heating from evaporation of falling precip'   ,phys_decomp)
    call addfld ('HFREEZ  ','W/kg'    ,pver, 'A','Heating rate due to freezing of precip'       ,phys_decomp)
    call addfld ('HMELT   ','W/kg'    ,pver, 'A','Heating from snow melt'                       ,phys_decomp)
    call addfld ('HREPART ','W/kg'    ,pver, 'A','Heating from cloud ice/liquid repartitioning' ,phys_decomp)
    call addfld ('REPARTICE','kg/kg/s',pver, 'A','Cloud ice tendency from cloud ice/liquid repartitioning' ,phys_decomp)
    call addfld ('REPARTLIQ','kg/kg/s',pver, 'A','Cloud liq tendency from cloud ice/liquid repartitioning' ,phys_decomp)
    call addfld ('FICE    ','fraction',pver, 'A','Fractional ice content within cloud'          ,phys_decomp)
    call addfld ('ICWMR   ','kg/kg   ',pver, 'A','Prognostic in-cloud water mixing ratio'       ,phys_decomp)
    call addfld ('ICIMR   ','kg/kg   ',pver, 'A','Prognostic in-cloud ice mixing ratio'         ,phys_decomp)
    call addfld ('PCSNOW  ','m/s'     ,1   , 'A','Snow fall from prognostic clouds'             ,phys_decomp)

    call addfld ('DQSED   ','kg/kg/s' ,pver, 'A','Water vapor tendency from cloud sedimentation',phys_decomp)
    call addfld ('DLSED   ','kg/kg/s' ,pver, 'A','Cloud liquid tendency from sedimentation'     ,phys_decomp)
    call addfld ('DISED   ','kg/kg/s' ,pver, 'A','Cloud ice tendency from sedimentation'        ,phys_decomp)
    call addfld ('HSED    ','W/kg'    ,pver, 'A','Heating from cloud sediment evaporation'      ,phys_decomp)
    call addfld ('SNOWSED ','m/s'     ,1   , 'A','Snow from cloud ice sedimentation'            ,phys_decomp)
    call addfld ('RAINSED ','m/s'     ,1   , 'A','Rain from cloud liquid sedimentation'         ,phys_decomp)
    call addfld ('PRECSED ','m/s'     ,1   , 'A','Precipitation from cloud sedimentation'       ,phys_decomp)

    call add_default ('FICE    ', 1, ' ')

    return
  end subroutine stratiform_init

!===============================================================================

  subroutine stratiform_tend(state, ptend_all, dtime, &
       icefrac, landfrac, ocnfrac, &
       landm,   snowh,    dlf, & 
       rliq ,   &  
       cmfmc,   cmfmc2,   concld,    &
       ts,      sst,      zdu,  &
       prec_str,snow_str, prec_sed, snow_sed, prec_pcw, snow_pcw, & 
       pbuf )

  !----------------------------------------------------------------------- 
  ! 
  ! Purpose: 
  !
  ! Interface to sedimentation, detrain, cloud fraction and 
  !        microphysics subroutines
  !
  ! Author: D.B. Coleman
  ! Date: Apr 2004
  ! 
  !-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid
  use cloud_fraction, only: cldfrc
  use pmgrid,        only: plon, plev, plat
  use physics_types, only: physics_state, physics_ptend, physics_tend
  use physics_types, only: physics_ptend_init, physics_update, physics_tend_init
  use physics_types, only: physics_ptend_sum,  physics_state_copy
  use history,         only: outfld
  use phys_buffer,     only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
  use constituents,    only: cnst_get_ind, ppcnst
  use pkg_cld_sediment, only: cld_sediment_vel, cld_sediment_tend
  use cldwat,        only: pcond, cldwat_fice
  use pkg_cldoptics, only: cldefr


  implicit none

!
! Parameters
!
  real(r8) pnot                  ! reference pressure
  parameter (pnot = 1.e5)


!
! Input arguments
!

  type(physics_state), intent(in)    :: state   ! state variables
  type(physics_ptend), intent(out) :: ptend_all   ! package tendencies
  type(pbuf_fld),      intent(inout), dimension(pbuf_size_max) :: pbuf

  real(r8), intent(in)  :: dtime                ! timestep
  real(r8), intent(in)  :: icefrac (pcols)      ! sea ice fraction (fraction)
  real(r8), intent(in)  :: landfrac(pcols)      ! land fraction (fraction)
  real(r8), intent(in)  :: ocnfrac (pcols)      ! ocean fraction (fraction)
  real(r8), intent(in) :: landm(pcols)          ! land fraction ramped over water
  real(r8), intent(in) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)

  real(r8), intent(in) :: dlf(pcols,pver)       ! detrained water from ZM
  real(r8), intent(in) :: rliq(pcols)                    ! vertical integral of liquid not yet in q(ixcldliq)
  real(r8), intent(in) :: cmfmc(pcols,pverp)    ! convective mass flux--m sub c
  real(r8), intent(in) :: cmfmc2(pcols,pverp)   ! shallow convective mass flux--m sub c

  real(r8), intent(in) :: ts(pcols)             ! surface temperature
  real(r8), intent(in) :: sst(pcols)            ! sea surface temperature
  real(r8), intent(in) :: zdu(pcols,pver)       ! detrainment rate from deep convection

  real(r8), intent(out)  :: prec_str(pcols)  ! [Total] sfc flux of precip from stratiform (m/s) 
  real(r8), intent(out)  :: snow_str(pcols)  ! [Total] sfc flux of snow from stratiform   (m/s)
  real(r8), intent(out)  :: prec_sed(pcols)  ! surface flux of total cloud water from sedimentation
  real(r8), intent(out)  :: snow_sed(pcols) ! surface flux of cloud ice from sedimentation
  real(r8), intent(out)  :: prec_pcw(pcols) ! sfc flux of precip from microphysics(m/s)
  real(r8), intent(out)  :: snow_pcw(pcols) ! sfc flux of snow from microphysics (m/s)

  real(r8), intent(out) :: concld(pcols,pver)    ! convective cloud cover


  !
  ! Local variables
  !

  type(physics_state)  :: state1   ! local copy of the state variable
  type(physics_tend ) :: tend          ! Physics tendencies (empty, needed for physics_update call)
  type(physics_ptend)  :: ptend_loc   ! package tendencies


  integer i,k,m
  integer :: lchnk                  ! chunk identifier
  integer :: ncol                   ! number of atmospheric columns

  ! physics buffer fields
  integer itim, ifld
  real(r8), pointer, dimension(:,:) :: qcwat   ! cloud water old q
  real(r8), pointer, dimension(:,:) :: tcwat   !cloud water old temperature
  real(r8), pointer, dimension(:,:) :: lcwat   ! cloud liquid water old q
  real(r8), pointer, dimension(:,:) ::  cld    ! cloud fraction
   real(r8), pointer, dimension(:,:) :: qme
   real(r8), pointer, dimension(:,:) :: prain
   real(r8), pointer, dimension(:,:) :: nevapr
   real(r8), pointer, dimension(:,:) :: rel     ! liquid effective drop radius (microns)
   real(r8), pointer, dimension(:,:) :: rei     ! ice effective drop size (microns)

  ! local variables for stratiform_sediment
  real(r8) :: rain(pcols)                       ! surface flux of cloud liquid
  real(r8) :: pvliq(pcols,pver+1)               ! vertical velocity of cloud liquid drops (Pa/s)
  real(r8) :: pvice(pcols,pver+1)               ! vertical velocity of cloud ice particles (Pa/s)

  ! local variables for cldfrc
  real(r8)  cldst(pcols,pver)     ! cloud fraction
  real(r8)  clc(pcols)            ! column convective cloud amount
  real(r8) rhdfda(pcols,pver)    ! d_RH/d_cloud_fraction    ====wlin
  real(r8) rhu00(pcols,pver)     ! RH limit, U00             ====wlin
  real(r8) relhum(pcols,pver)         ! RH, output to determine drh/da
  real(r8) rhu002(pcols,pver)         ! same as rhu00 but for perturbed rh 
  real(r8) cld2(pcols,pver)          ! same as cld but for perturbed rh  
  real(r8) concld2(pcols,pver)        ! same as concld but for perturbed rh 
  real(r8) cldst2(pcols,pver)         ! same as cldst but for perturbed rh 
  real(r8) relhum2(pcols,pver)        ! RH after  perturbation            
  real(r8) :: pmid(pcols,pver)      ! midpoint pressures
  real(r8) :: t(pcols,pver)         ! temperature
  real(r8) :: q(pcols,pver)         ! specific humidity
  real(r8) :: omga(pcols,pver)      ! vertical pressure velocity
  real(r8) :: phis(pcols)           ! surface geopotential
  real(r8) :: pdel(pcols,pver)      ! pressure depth of layer
  real(r8) :: ps(pcols)             ! surface pressure

  ! local variables for microphysics
  real(r8) :: rdtime                          ! 1./dtime
  real(r8) :: qtend (pcols,pver)            ! moisture tendencies
  real(r8) :: ttend (pcols,pver)            ! temperature tendencies
  real(r8) :: ltend (pcols,pver)            ! cloud liquid water tendencies
  real(r8) :: evapheat(pcols,pver)          ! heating rate due to evaporation of precip
  real(r8) :: evapsnow(pcols,pver)          ! local evaporation of snow
  real(r8) :: prfzheat(pcols,pver)          ! heating rate due to freezing of precip (W/kg)
  real(r8) :: meltheat(pcols,pver)          ! heating rate due to phase change of precip
  real(r8) :: cmeheat (pcols,pver)          ! heating rate due to phase change of precip
  real(r8) :: prodsnow(pcols,pver)          ! local production of snow
  real(r8) :: totcw   (pcols,pver)          ! total cloud water mixing ratio
  real(r8) :: fice    (pcols,pver)          ! Fractional ice content within cloud
  real(r8) :: fsnow   (pcols,pver)          ! Fractional snow production
  real(r8) :: repartht(pcols,pver)          ! heating rate due to phase repartition of input precip
  real(r8) :: icimr(pcols,pver)             ! in cloud ice mixing ratio
  real(r8) :: icwmr(pcols,pver)             ! in cloud water mixing ratio
  real(r8) fwaut(pcols,pver)              
  real(r8) fsaut(pcols,pver)              
  real(r8) fracw(pcols,pver)              
  real(r8) fsacw(pcols,pver)              
  real(r8) fsaci(pcols,pver)              
  real(r8) cmeice(pcols,pver)   ! Rate of cond-evap of ice within the cloud
  real(r8) cmeliq(pcols,pver)   ! Rate of cond-evap of liq within the cloud
  real(r8) ice2pr(pcols,pver)   ! rate of conversion of ice to precip
  real(r8) liq2pr(pcols,pver)   ! rate of conversion of liquid to precip
  real(r8) liq2snow(pcols,pver)   ! rate of conversion of liquid to snow
  real(r8) temp(pcols)
  real(r8) res(pcols,pver)

  !======================================================================
  lchnk = state%lchnk
  ncol  = state%ncol


  call physics_state_copy(state,state1)   ! copy state to local state1.
  call physics_ptend_init(ptend_loc)  ! initialize local ptend type
  call physics_ptend_init(ptend_all)  ! initialize output ptend type
  call physics_tend_init(tend)   ! tend here is just a null place holder


  ! Associate pointers with physics buffer fields
  itim = pbuf_old_tim_idx()
  ifld = pbuf_get_fld_idx('QCWAT')
  qcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
  ifld = pbuf_get_fld_idx('TCWAT')
  tcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
  ifld = pbuf_get_fld_idx('LCWAT')
  lcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
  ifld = pbuf_get_fld_idx('CLD')
  cld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

   ifld = pbuf_get_fld_idx('QME')
   qme  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld = pbuf_get_fld_idx('PRAIN')
   prain  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld = pbuf_get_fld_idx('NEVAPR')
   nevapr  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

   ifld = pbuf_get_fld_idx('REL')
   rel  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld = pbuf_get_fld_idx('REI')
   rei  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)



  !+++sediment ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Allow the cloud liquid drops and ice particles to sediment
  ! Occurs before adding convectively detrained cloud water, because the phase of the
  ! of the detrained water is unknown.
  call t_startf('stratiform_sediment')


  ptend_loc%name         = 'pcwsediment'
  ptend_loc%ls           = .TRUE.
  ptend_loc%lq(1)        = .TRUE.
  ptend_loc%lq(ixcldice) = .TRUE.
  ptend_loc%lq(ixcldliq) = .TRUE.

  call cld_sediment_vel (ncol,                                    &
       icefrac, landfrac, ocnfrac, state1%pmid, state1%pdel, state1%t, &
       cld, state1%q(:,:,ixcldliq), state1%q(:,:,ixcldice), pvliq, pvice, landm, snowh)

  call cld_sediment_tend (ncol, dtime ,                                       &
       state1%pint, state1%pmid           , state1%pdel            , state1%t     , &
       cld  ,    state1%q(:,:,ixcldliq), state1%q(:,:,ixcldice) , pvliq, pvice, &
       ptend_loc%q(:,:,ixcldliq), ptend_loc%q(:,:,ixcldice), ptend_loc%q(:,:,1), ptend_loc%s  , &
       rain   , snow_sed   )

  ! convert rain and snow from kg/m2 to m/s
  snow_sed(:ncol) = snow_sed(:ncol)/1000.
  rain(:ncol) = rain(:ncol)/1000.
  ! compute total precip (m/s)
  prec_sed(:ncol) = rain(:ncol) + snow_sed(:ncol)

  ! record history variables
  lchnk = state1%lchnk
  call outfld('DQSED'   ,ptend_loc%q(:,:,1)       , pcols,lchnk)
  call outfld('DISED'   ,ptend_loc%q(:,:,ixcldice), pcols,lchnk)
  call outfld('DLSED'   ,ptend_loc%q(:,:,ixcldliq), pcols,lchnk)
  call outfld('HSED'    ,ptend_loc%s              , pcols,lchnk)
  call outfld('PRECSED' ,prec_sed                 , pcols,lchnk)
  call outfld('SNOWSED' ,snow_sed                 , pcols,lchnk)
  call outfld('RAINSED' ,rain                 , pcols,lchnk)

  ! add tendency from this process to tend from other processes here
  call physics_ptend_sum(ptend_loc,ptend_all, state)

  ! update physics state type state1 with ptend_loc 
  call physics_update(state1, tend, ptend_loc, dtime)
  call physics_ptend_init(ptend_loc)


  call t_stopf('stratiform_sediment')


  ! accumulate prec and snow
  prec_str(:ncol) = prec_sed(:ncol)
  snow_str(:ncol) = snow_sed(:ncol)


  !++detrain ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Put the detraining cloud water from convection into the cloud and environment. 

  ptend_loc%name         = 'pcwdetrain'
!!$    ptend_loc%ls           = .TRUE.
!!$    ptend_loc%lq(1)        = .TRUE.
!!$    ptend_loc%lq(ixcldice) = .TRUE.
  ptend_loc%lq(ixcldliq) = .TRUE.
  !
  ! Put all of the detraining cloud water from convection into the large scale cloud.
  ! It all goes in liquid for the moment.
  do k = 1,pver
     do i = 1,state1%ncol
!!$          ptend_loc%q(i,k,1)        = dlf(i,k) * (1.-cld(i,k))
!!$          ptend_loc%s(i,k)          =-dlf(i,k) * (1.-cld(i,k))*latvap
!!$          ptend_loc%q(i,k,ixcldice) = 0.
!!$          ptend_loc%q(i,k,ixcldliq) = dlf(i,k) * cld(i,k)
        ptend_loc%q(i,k,ixcldliq) = dlf(i,k)
     end do
  end do
  call outfld('ZMDLF' ,dlf  , pcols,state1%lchnk)

  ! add tendency from this process to tend from other processes here
  call physics_ptend_sum(ptend_loc,ptend_all, state)

  call physics_update(state1, tend, ptend_loc, dtime)
  call physics_ptend_init(ptend_loc)

  ! accumulate prec and snow, reserved liquid has now been used
  prec_str(:ncol) = prec_str(:ncol) - rliq(:ncol)  ! ( snow contribution is zero )

  !++++ cldfrc ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !
  ! cloud fraction after transport and convection,
  ! derive the relationship between rh and cld from 
  ! the employed c7loud scheme
  !
  pmid(:ncol,:pver) = state1%pmid(:ncol,:pver)
  t(:ncol,:pver) = state1%t(:ncol,:pver)
  q(:ncol,:pver) = state1%q(:ncol,:pver,1)
  omga(:ncol,:pver) = state1%omega(:ncol,:pver)
  pdel(:ncol,:pver) =  state1%pdel(:ncol,:pver)
  ps(:ncol) = state1%pint(:ncol,pverp)
  phis(:ncol) = state1%phis(:ncol)

  !
  call t_startf("cldfrc")
  call cldfrc(lchnk,   ncol,                                &
       pmid,      t,        q,     omga, phis, &
       cld,    clc,     pdel,   &
       cmfmc,   cmfmc2,  landfrac,snowh,   concld,  cldst,    &
       ts,      sst, ps,      zdu,     ocnfrac, rhu00, &
       relhum,  0  )    

  ! re-calculate cloud with perturbed rh             add call cldfrc  
  call cldfrc(lchnk,   ncol,                                &
       pmid,       t,      q,      omga, phis, &
       cld2,   clc,     pdel,   &
       cmfmc,   cmfmc2   ,landfrac,snowh,   concld2, cldst2,   &
       ts,      sst, ps,        zdu,   ocnfrac, rhu002, &
       relhum2, 1  )              

  call t_stopf("cldfrc")

  ! cldfrc does not define layer cloud for model layer at k=1
  ! so set rhu00(k=1)=2.0 to not calculate cme for this layer

  rhu00(:ncol,1) = 2.0 

  ! Add following to estimate rhdfda                       

  do k=1,pver
     do i=1,ncol
        if(relhum(i,k) < rhu00(i,k) ) then
           rhdfda(i,k)=0.0
        else if (relhum(i,k) >= 1.0 ) then
           rhdfda(i,k)=0.0
        else
           !under certain circumstances, rh+ cause cld not to changed
           !when at an upper limit, or w/ strong subsidence
           !need to further check whether this if-block is necessary

           if((cld2(i,k) - cld(i,k) ) < 1.e-4 ) then
              rhdfda(i,k) = 0.01*relhum(i,k)*1.e+4   !instead of 0.0
           else
              rhdfda(i,k)=0.01*relhum(i,k)/(cld2(i,k)-cld(i,k))
           endif
        endif
     enddo
  enddo


  call outfld('CONCLD  ',concld, pcols,lchnk)
  call outfld('CLDST   ',cldst,  pcols,lchnk)
  call outfld('CNVCLD  ',clc,    pcols,lchnk)

  !+ mp +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! cloud water and ice parameterizations
  call t_startf('stratiform_microphys')
  !drb   prec_pcw = 0.
  !drb   snow_pcw = 0.


  ! Initialize chunk id and size
  lchnk = state1%lchnk
  ncol  = state1%ncol

  ! associate local pointers with fields in the physics buffer
!!$    buffld1 => pbuf(ixbuffld1)%fld_ptr(1,1:pcols,1,     lchnk,1)
!!$    buffld2 => pbuf(ixbuffld2)%fld_ptr(1,1:pcols,1:pver,lchnk,1)

  ! Define fractional amount of cloud condensate in ice phase
  call cldwat_fice(ncol, state1%t, fice, fsnow)

  ! compute total cloud water
  totcw(:ncol,:pver) = state1%q(:ncol,:pver,ixcldice) + state1%q(:ncol,:pver,ixcldliq)

  ! save input cloud ice
  repartht(:ncol,:pver) = state1%q(:ncol,:pver,ixcldice)

  ! Repartition ice and liquid
!  state1%q(:ncol,:pver,ixcldice) = totcw(:ncol,:pver) * fice(:ncol,:pver)
!  state1%q(:ncol,:pver,ixcldliq) = totcw(:ncol,:pver) * (1.0_r8 - fice(:ncol,:pver))
  rdtime = 1./dtime
  ptend_loc%q(:ncol,:pver,ixcldice) = &
       ( totcw(:ncol,:pver) * fice(:ncol,:pver)            - state1%q(:ncol,:pver,ixcldice) ) * rdtime
  ptend_loc%q(:ncol,:pver,ixcldliq) = &
       ( totcw(:ncol,:pver) * (1.0_r8 - fice(:ncol,:pver)) - state1%q(:ncol,:pver,ixcldliq) ) * rdtime

  call outfld('REPARTICE'   ,ptend_loc%q(:,:,ixcldice), pcols,lchnk)
  call outfld('REPARTLIQ'   ,ptend_loc%q(:,:,ixcldliq), pcols,lchnk)

  ! Set output flags
  ptend_loc%name         = 'cldwat-repartition'
  ptend_loc%lq(ixcldice) = .true.
  ptend_loc%lq(ixcldliq) = .true.

  ! add tendency from this process to tend from other processes here
  call physics_ptend_sum(ptend_loc,ptend_all, state)

  ! update state for use below
  call physics_update (state1, tend, ptend_loc, dtime)
  call physics_ptend_init(ptend_loc)


  ! Set output flags for final cldwat update
  ptend_loc%name         = 'cldwat'
  ptend_loc%ls           = .true.
  ptend_loc%lq(1)        = .true.
  ptend_loc%lq(ixcldice) = .true.
  ptend_loc%lq(ixcldliq) = .true.


  ! Determine heating from change in cloud ice
  repartht(:ncol,:pver) = latice/dtime * (state1%q(:ncol,:pver,ixcldice) - repartht(:ncol,:pver))

  ! calculate the tendencies for moisture, temperature and cloud fraction
  qtend(:ncol,:pver) = (state1%q(:ncol,:pver,1) - qcwat(:ncol,:pver))*rdtime
  ttend(:ncol,:pver) = (state1%t(:ncol,:pver)   - tcwat(:ncol,:pver))*rdtime
  ltend(:ncol,:pver) = (totcw  (:ncol,:pver)   - lcwat(:ncol,:pver))*rdtime

  ! call microphysics package to calculate tendencies
  call t_startf('pcond')
  call pcond (lchnk,   ncol, &
       state1%t  , ttend      , state1%q(1,1,1), qtend   , state1%omega, &
       totcw    , state1%pmid , state1%pdel , cld       , fice       , fsnow, &
       qme      , prain   , prodsnow, nevapr   , evapsnow   , evapheat   , prfzheat, &
       meltheat , prec_pcw     , snow_pcw       , dtime         , fwaut      , &
       fsaut    , fracw      , fsacw      , fsaci      , ltend      , &
       rhdfda   , rhu00      , icefrac    , state1%zi   , ice2pr, liq2pr, liq2snow, snowh)
  call t_stopf('pcond')

  ! make it interactive
  do i = 1,ncol
     do k = 1,pver
        ptend_loc%s(i,k)          = qme(i,k) * (latvap + latice*fice(i,k)) &
             + evapheat(i,k) + prfzheat(i,k) + meltheat(i,k) + repartht(i,k)

        ptend_loc%q(i,k,1)        =-qme(i,k) + nevapr(i,k)

        ptend_loc%q(i,k,ixcldice) =qme(i,k)*fice(i,k)      - ice2pr(i,k)
        ptend_loc%q(i,k,ixcldliq) =qme(i,k)*(1.-fice(i,k)) - liq2pr(i,k)
     end do
  end do


#ifdef DEBUG
  if (lchnk.eq.248) then
     i = 12
     do k = 1,pver
        call debug_microphys_1(state1,ptend_loc,i,k, &
                dtime,qme,fice,snow_pcw,prec_pcw, &
                prain,nevapr,prodsnow, evapsnow, &
                ice2pr,liq2pr,liq2snow)
     end do
  endif
  call debug_microphys_2(state1,&
       snow_pcw,fsaut,fsacw ,fsaci, meltheat)
#endif

  ! Compute in cloud ice and liquid mixing ratios
  do k=1,pver
     do i = 1,ncol
        icimr(i,k) = (state1%q(i,k,ixcldice) + dtime*ptend_loc%q(i,k,ixcldice)) / max(0.01_r8,cld(i,k))
        icwmr(i,k) = (state1%q(i,k,ixcldliq) + dtime*ptend_loc%q(i,k,ixcldliq)) / max(0.01_r8,cld(i,k))
     end do
  end do


  ! convert precipitation from kg/m2 to m/s
  snow_pcw  (:ncol) = snow_pcw  (:ncol)/1000.
  prec_pcw(:ncol) = prec_pcw(:ncol)/1000.

  do i = 1,ncol
     do k = 1,pver
        cmeheat(i,k) = qme(i,k) * (latvap + latice*fice(i,k))
        cmeice (i,k) = qme(i,k) * fice(i,k)
        cmeliq (i,k) = qme(i,k) * (1. - fice(i,k))
     end do
  end do

  ! record history variables
  call outfld('FWAUT',fwaut,    pcols,lchnk)
  call outfld('FSAUT',fsaut,    pcols,lchnk)
  call outfld('FRACW',fracw,    pcols,lchnk)
  call outfld('FSACW',fsacw,    pcols,lchnk)
  call outfld('FSACI',fsaci,    pcols,lchnk)
  call outfld('ICIMR',icimr,    pcols,lchnk)
  call outfld('ICWMR',icwmr,    pcols,lchnk)

  call outfld('PCSNOW'  ,snow_pcw    , pcols,lchnk)
  call outfld('FICE'    ,fice    , pcols,lchnk)
  call outfld('CME'     ,qme     , pcols,lchnk)
  call outfld('CMEICE'  ,cmeice  , pcols,lchnk)
  call outfld('CMELIQ'  ,cmeliq  , pcols,lchnk)
  call outfld('ICE2PR'  ,ice2pr  , pcols,lchnk)
  call outfld('LIQ2PR'  ,liq2pr  , pcols,lchnk)
  call outfld('PRODPREC',prain, pcols,lchnk)
  call outfld('EVAPPREC',nevapr, pcols,lchnk)
  call outfld('EVAPSNOW',evapsnow, pcols,lchnk)
  call outfld('HPROGCLD',ptend_loc%s , pcols,lchnk)
  call outfld('HEVAP   ',evapheat, pcols,lchnk)
  call outfld('HMELT'   ,meltheat, pcols,lchnk)
  call outfld('HCME'    ,cmeheat , pcols,lchnk)
  call outfld('HFREEZ'  ,prfzheat, pcols,lchnk)
  call outfld('HREPART' ,repartht, pcols,lchnk)

  ! update boundary quantities
!!$    ptend_loc%hflx_srf = 0.
!!$    ptend_loc%hflx_top = 0.
!!$    ptend_loc%cflx_srf = 0.
!!$    ptend_loc%cflx_top = 0.


  ! add tendency from this process to tend from other processes here
  call physics_ptend_sum(ptend_loc,ptend_all, state)

  ! Set the name of the final package tendencies. Note that this
  ! is a special case in physics_update, so a change here must be 
  ! matched there.
  ptend_all%name = 'stratiform'


  ! used below
  call physics_update (state1, tend, ptend_loc, dtime)
  call physics_ptend_init(ptend_loc) 

  call t_stopf('stratiform_microphys')

  ! accumulate prec and snow
  prec_str(:ncol) = prec_str(:ncol) + prec_pcw(:ncol)
  snow_str(:ncol) = snow_str(:ncol) + snow_pcw(:ncol)

  ! Save off q and t after cloud water
  call cnst_get_ind('CLDLIQ', ixcldliq)
  call cnst_get_ind('CLDICE', ixcldice)

  do k=1,pver
     qcwat(:ncol,k) = state1%q(:ncol,k,1)
     tcwat(:ncol,k) = state1%t(:ncol,k)
     lcwat(:ncol,k) = state1%q(:ncol,k,ixcldice) + state1%q(:ncol,k,ixcldliq)

  end do

!
! Cloud water and ice particle sizes, saved in physics buffer for radiation
!
    call cldefr(lchnk, ncol, landfrac, state1%t, rel, rei, state1%ps, state1%pmid, landm, icefrac, snowh)


end subroutine stratiform_tend 


!===============================================================================

   subroutine debug_microphys_1(state1,ptend,i,k, &
        dtime,qme,fice,snow_pcw,prec_pcw, &
        prain,nevapr,prodsnow, evapsnow, &
        ice2pr,liq2pr,liq2snow)

     use physics_types, only: physics_state, physics_ptend
     use physconst,     only: tmelt

     implicit none
     
     integer, intent(in) :: i,k
     type(physics_state), intent(in) :: state1   ! local copy of the state variable
     type(physics_ptend), intent(in) :: ptend  ! local copy of the ptend variable
     real(r8), intent(in)  :: dtime                ! timestep
     real(r8), intent(in) :: qme(pcols,pver)          ! local condensation - evaporation of cloud water

     real(r8), intent(in) :: prain(pcols,pver)          ! local production of precipitation
     real(r8), intent(in) :: nevapr(pcols,pver)          ! local evaporation of precipitation
     real(r8), intent(in) :: prodsnow(pcols,pver)          ! local production of snow
     real(r8), intent(in) :: evapsnow(pcols,pver)          ! local evaporation of snow
     real(r8), intent(in) :: ice2pr(pcols,pver)   ! rate of conversion of ice to precip
     real(r8), intent(in) :: liq2pr(pcols,pver)   ! rate of conversion of liquid to precip
     real(r8), intent(in) :: liq2snow(pcols,pver)   ! rate of conversion of liquid to snow
     real(r8), intent(in) :: fice    (pcols,pver)          ! Fractional ice content within cloud
     real(r8), intent(in) :: snow_pcw(pcols)
     real(r8), intent(in) :: prec_pcw(pcols)

     real(r8) hs1, qv1, ql1, qi1, qs1, qr1, fice2, pr1, w1, w2, w3, fliq, res
     real(r8) w4, wl, wv, wi, wlf, wvf, wif, qif, qlf, qvf

     pr1 = 0
     hs1 = 0
     qv1 = 0
     ql1 = 0
     qi1 = 0
     qs1 = 0
     qr1 = 0
     w1 = 0
     wl = 0
     wv = 0
     wi = 0
     wlf = 0
     wvf = 0
     wif = 0


     write (6,*) 
     write (6,*) ' input state, t, q, l, i ', k, state1%t(i,k), state1%q(i,k,1), state1%q(i,k,ixcldliq),  state1%q(i,k,ixcldice)
     write (6,*) ' rain, snow, total from components before accumulation ', qr1, qs1, qr1+qs1
     write (6,*) ' total precip before accumulation                      ', k, pr1

     wv = wv + state1%q(i,k,1       )*state1%pdel(i,k)/gravit
     wl = wl + state1%q(i,k,ixcldliq)*state1%pdel(i,k)/gravit
     wi = wi + state1%q(i,k,ixcldice)*state1%pdel(i,k)/gravit

     qvf = state1%q(i,k,1) + ptend%q(i,k,1)*dtime
     qlf = state1%q(i,k,ixcldliq) + ptend%q(i,k,ixcldliq)*dtime
     qif = state1%q(i,k,ixcldice) + ptend%q(i,k,ixcldice)*dtime

     if (qvf.lt.0.) then
        write (6,*) ' qvf is negative *******', qvf
     endif
     if (qlf.lt.0.) then
        write (6,*) ' qlf is negative *******', qlf
     endif
     if (qif.lt.0.) then
        write (6,*) ' qif is negative *******', qif
     endif
     write (6,*) ' qvf, qlf, qif ', qvf, qlf, qif

     wvf = wvf + qvf*state1%pdel(i,k)/gravit
     wlf = wlf + qlf*state1%pdel(i,k)/gravit
     wif = wif + qif*state1%pdel(i,k)/gravit

     hs1 = hs1 + ptend%s(i,k)*state1%pdel(i,k)/gravit
     pr1 = pr1 + state1%pdel(i,k)/gravit*(prain(i,k)-nevapr(i,k))
     qv1 = qv1 - (qme(i,k)-nevapr(i,k))*state1%pdel(i,k)/gravit    ! vdot
     w1  = w1  + (qme(i,k)-prain(i,k))*state1%pdel(i,k)/gravit    ! cdot
     qi1 = qi1 + ((qme(i,k))*fice(i,k)        -ice2pr(i,k) )*state1%pdel(i,k)/gravit   ! idot
     ql1 = ql1 + ((qme(i,k))*(1._r8-fice(i,k))-liq2pr(i,k) )*state1%pdel(i,k)/gravit   ! ldot

     qr1 = qr1 &
          + ( liq2pr(i,k)-liq2snow(i,k)   &     ! production of rain
          -(nevapr(i,k)-evapsnow(i,k)) &     ! rain evaporation
          )*state1%pdel(i,k)/gravit
     qs1 = qs1 &
          + ( ice2pr(i,k) + liq2snow(i,k) &     ! production of snow.Note last term has phase change
          -evapsnow(i,k)               &     ! snow evaporation
          )*state1%pdel(i,k)/gravit

     if (state1%t(i,k).gt.tmelt) then
        qr1 = qr1 + qs1
        qs1 = 0.
     endif
     write (6,*) ' rain, snow, total after accumulation ', qr1, qs1, qr1+qs1
     write (6,*) ' total precip after accumulation      ', k, pr1
     write (6,*)
     write (6,*) ' layer prain, nevapr, pdel ', prain(i,k), nevapr(i,k), state1%pdel(i,k)
     write (6,*) ' layer prodsnow, ice2pr+liq2snow ', prodsnow(i,k), ice2pr(i,k)+liq2snow(i,k)
     write (6,*) ' layer prain-prodsnow, liq2pr-liq2snow ', prain(i,k)-prodsnow(i,k), liq2pr(i,k)-liq2snow(i,k)
     write (6,*) ' layer evapsnow, evaprain ', k, evapsnow(i,k), nevapr(i,k)-evapsnow(i,k)
     write (6,*) ' layer ice2pr, liq2pr, liq2snow ', ice2pr(i,k), liq2pr(i,k), liq2snow(i,k)
     write (6,*) ' layer ice2pr+liq2pr, prain ', ice2pr(i,k)+liq2pr(i,k), prain(i,k)
     write (6,*)
     write (6,*) ' qv1 vapor removed from col after accum  (vdot)   ', k, qv1
     write (6,*) ' - (precip produced - vapor removed) after accum  ', k, -pr1-qv1
     write (6,*) ' condensate produce after accum                   ', k, w1
     write (6,*) ' liq+ice tends accum                              ', k, ql1+qi1
     write (6,*) ' change in total water after accum                ', k, qv1+ql1+qi1
     write (6,*) ' imbalance in colum after accum                   ', k, qs1+qr1+qv1+ql1+qi1
     write (6,*) ' fice at this lev ', fice(i,k)
     write (6,*)

     res = abs((qs1+qr1+qv1+ql1+qi1)/max(abs(qv1),abs(ql1),abs(qi1),abs(qs1),abs(qr1),1.e-36))
     write (6,*) ' relative residual in column method 1             ', k, res

     write (6,*) ' relative residual in column method 2             ', k, abs((qs1+qr1+qv1+ql1+qi1)/max(abs(qv1+ql1+qi1),1.e-36))
     !            if (abs((qs1+qr1+qv1+ql1+qi1)/(qs1+qr1+1.e-36)).gt.1.e-14) then
     if (res.gt.1.e-14) then
        call endrun ('STRATIFORM_TEND')
     endif

     !             w3  = qme(i,k) * (latvap + latice*fice(i,k)) &
     !               + evapheat(i,k) + prfzheat(i,k) + meltheat(i,k)

     res = qs1+qr1-pr1
     w4 = max(abs(qs1),abs(qr1),abs(pr1)) 
     if (w4.gt.0.)  then
        if (res/w4.gt.1.e-14) then
           write (6,*) ' imbalance in precips calculated two ways '
           write (6,*) ' res/w4, pr1, qr1, qs1, qr1+qs1 ', &
                res/w4, pr1, qr1, qs1, qr1+qs1
           !                   call endrun()
        endif
     endif
     if (k.eq.pver) then
        write (6,*) ' pcond returned precip, rain and snow rates ', prec_pcw(i), prec_pcw(i)-snow_pcw(i), snow_pcw(i)
        write (6,*) ' I calculate ', pr1, qr1, qs1
        !               call endrun
        write (6,*) ' byrons water check ', wv+wl+wi-pr1*dtime, wvf+wlf+wif
     endif
     write (6,*)


   end subroutine debug_microphys_1

!===============================================================================


   subroutine debug_microphys_2(state1,&
        snow_pcw,fsaut,fsacw ,fsaci, meltheat)

     use ppgrid,        only:  pver
     use physics_types, only: physics_state
     
     implicit none

     type(physics_state), intent(in) :: state1   ! local copy of the state variable
     real(r8), intent(in) :: snow_pcw(pcols)
     real(r8), intent(in) :: fsaut(pcols,pver)              
     real(r8), intent(in) :: fsacw(pcols,pver)              
     real(r8), intent(in) :: fsaci(pcols,pver)              
     real(r8), intent(in) :: meltheat(pcols,pver)          ! heating rate due to phase change of precip


     integer  i,ncol,lchnk


     ncol = state1%ncol
     lchnk = state1%lchnk
     
     do i = 1,ncol
        if (snow_pcw(i).gt.0.01/8.64e4.and.state1%t(i,pver).gt.273.16) then
           write (6,*) ' stratiform: snow, temp, ', i, lchnk, &
                snow_pcw(i), state1%t(i,pver)
           write (6,*) ' t ', state1%t(i,:)
           write (6,*) ' fsaut ', fsaut(i,:)
           write (6,*) ' fsacw ', fsacw(i,:)
           write (6,*) ' fsaci ', fsaci(i,:)
           write (6,*) ' meltheat ', meltheat(i,:)
           call endrun ('STRATIFORM_TEND')
        endif
        
        if (snow_pcw(i)*8.64e4.lt.-1.e-5) then
           write (6,*) ' neg snow ', snow_pcw(i)*8.64e4
           write (6,*) ' stratiform: snow_pcw, temp, ', i, lchnk, &
                snow_pcw(i), state1%t(i,pver)
           write (6,*) ' t ', state1%t(i,:)
           write (6,*) ' fsaut ', fsaut(i,:)
           write (6,*) ' fsacw ', fsacw(i,:)
           write (6,*) ' fsaci ', fsaci(i,:)
           write (6,*) ' meltheat ', meltheat(i,:)
           call endrun ('STRATIFORM_TEND')
        endif
     end do
     

   end subroutine debug_microphys_2


!===============================================================================

end module stratiform

