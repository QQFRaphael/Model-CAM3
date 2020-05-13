#include <misc.h>
#include <params.h>

module vertical_diffusion

  !---------------------------------------------------------------------------------
  ! Module to compute vertical diffusion of momentum,  moisture, trace constituents
  ! and static energy. Separate modules compute 
  !   1. stresses associated with turbulent flow over orography
  !   2. eddy diffusivities, including nonlocal tranport terms
  !   3. molecular diffusivities
  !   4. coming soon... gravity wave drag
  ! Lastly, a implicit diffusion solver is called, and tendencies retrieved by 
  ! differencing the diffused and initial states.
  !
  ! calling sequence:
  !
  !  vertical_diffusion_init           initializes vertical diffustion constants and modules
  !        init_tms                      initializes turbulent mountain stress module
  !        init_hb_diff                  initializes eddy diffusivity module (includes PBL)
  !        init_molec_diff               initializes molecular diffusivity module 
  !        init_vdiff                    initializes diffusion solver module
  !  vertical_diffusion_ts_init        time step initialization (only used for upper boundary condition)
  !  vertical_diffusion_tend           computes vertical diffusion tendencies
  !        compute_tms                   computes turbulent mountain stresses
  !        compute_hb_diff               computes eddy diffusivities and countergradient terms
  !        compute_vdiff                 solves vertical diffusion equations (including molecular diffusivities)
  !
  !---------------------------Code history--------------------------------
  ! This is a reorganization of the original vertical diffusion module
  ! written by J. Rosinski in June 1992, and modified by many others.
  ! Initial coding for this version:  J. McCaa, September 2004.
  !---------------------------------------------------------------------------------

  use shr_kind_mod,     only: r8 => shr_kind_r8
  use ppgrid,           only: pcols, pver, pverp
  use constituents,     only: pcnst, pnats
  use diffusion_solver, only: vdiff_selector
  use abortutils,       only: endrun
  use physconst,        only: &
       cpair  , &     ! Specific heat of dry air
       gravit , &     ! Acceleration due to gravity
       rair   , &     ! Gas const for dry air
       zvir   , &     ! rh2o/rair - 1
       latvap , &     ! latent heat of vaporization
       karman , &     ! von Karman constant
       mwdry  , &     ! molecular weight of dry air
       avogad , &     ! Avogadro's number
       boltz          ! Boltzman's constant
       
  implicit none
  private          ! Make default type private to the module
  save
  
  !
  ! Public interfaces
  !
  public vd_register         ! Register multi-time-level variables with physics buffer
  public vertical_diffusion_init             ! Initialization
  public vertical_diffusion_ts_init          ! Time step initialization (only used for upper boundary condition)
  public vertical_diffusion_tend             ! Full routine
  ! 
  ! Private data
  !
  !   Local copies of physical constants
  !
  !  Other local data shared between routines
  logical :: do_molec_diff = .false.        ! switch for molecular diffusion
  logical :: do_tms        = .false.        ! switch for turbulent mountain stress
  type(vdiff_selector) :: fieldlist         ! Logical switches for moist vs dry mixing ratio diffusion
  integer :: ntop                           ! Top level to which vertical diffusion is applied (=1).
  integer :: nbot                           ! Bottom level to which vertical diffusion is applied (=pver).
  integer :: tke_idx                        ! indices for fields in the physics buffer
  character(len=8) :: vdiffnam(pcnst+pnats) ! names of v-diff tendencies
  integer :: ixcldice, ixcldliq             ! constituent indices for cloud liquid and ice water

contains

  !===============================================================================
  subroutine vd_register()
    !-----------------------------------------------------------------------
    ! Register physics buffer fields and constituents
    !-----------------------------------------------------------------------
    use phys_buffer, only: pbuf_times, pbuf_add
    
    ! Request physics buffer space for fields that persist across timesteps.
    call pbuf_add('tke', 'global', 1,pverp,pbuf_times, tke_idx) 

  end subroutine vd_register
  
  !===============================================================================
  subroutine vertical_diffusion_init()
    !-----------------------------------------------------------------------
    ! Initialization of time independent fields for vertical diffusion.
    ! Calls initialization routines for subsidiary modules
    !-----------------------------------------------------------------------
    use history,    only: addfld, add_default, phys_decomp
    use hb_diff, only: init_hb_diff
    use molec_diff, only: init_molec_diff
    use trb_mtn_stress, only: init_tms
    use diffusion_solver, only: init_vdiff, vdiff_select
    use constituents, only: cnst_get_ind, cnst_get_type_byind, cnst_name
    use pmgrid,       only: masterproc
    
    ! The next three lines should be replaced as soon as hypm is in a module
    use pmgrid,             only: plev, plevp
    ! include hypm (reference pressures at midpoints) 
#include <comhyb.h>

    !------------------------------Arguments--------------------------------
    ! none
    !---------------------------Local workspace-----------------------------
    integer :: ntop_eddy   ! Top level to which eddy vertical diffusion is applied.
    integer :: nbot_eddy   ! Bottom level to which eddy vertical diffusion is applied (=pver).
    integer :: ntop_molec  ! Top level to which molecular vertical diffusion is applied (=1).
    integer :: nbot_molec  ! Bottom level to which molecular vertical diffusion is applied.
    integer :: k           ! vertical loop index
    character(128) :: errstring                    ! error status for init_vdiff
    !-----------------------------------------------------------------------
    ! Get indices of cloud liquid and ice within the constituents array
    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)

    ! Initialize molecular diffusivity module
    ! Molecular diffusion turned on above ~60 km (50 Pa) if model top is above ~90 km (.1 Pa).
    ! Note that computing molecular diffusivities is a trivial expense, but constituent
    ! diffusivities depend on their molecular weights. Decomposing the diffusion matric
    ! for each constituent is a needless expense unless the diffusivity is significant.
    ntop_molec = 1       ! should always be 1
    nbot_molec = 0       ! should be set below about 70 km
    if (hypm(1) .lt. 0.1) then
       do_molec_diff = .true.
       do k = 1, pver
          if (hypm(k) .lt. 50.) nbot_molec  = k
       end do
       call init_molec_diff(r8, pcnst+pnats, rair, ntop_molec, nbot_molec, mwdry, &
            avogad, gravit, cpair, boltz)
       call addfld      ('TTPXMLC', 'K/S', 1, 'A','Top interf. temp. flux: molec. viscosity', phys_decomp)
       call add_default ('TTPXMLC', 1, ' ')
       if (masterproc) write (6,fmt='(a,i3,5x,a,i3)') 'NTOP_MOLEC =',ntop_molec, 'NBOT_MOLEC =',nbot_molec
    end if
    
    ! Initialize eddy diffusivity module
    ntop_eddy  = 1       ! no reason not to make this 1, if >1, must be <= nbot_molec
    nbot_eddy  = pver    ! should always be pver
    if(masterproc)write (6,fmt='(a,i3,5x,a,i3)') 'NTOP_EDDY  =',ntop_eddy,  'NBOT_EDDY  =',nbot_eddy
    call init_hb_diff(gravit, cpair, rair, zvir, ntop_eddy, nbot_eddy, hypm, karman)
    
    ! The vertical diffusion solver must operate over the full range of molecular and eddy diffusion
    ntop = min(ntop_molec,ntop_eddy)
    nbot = max(nbot_molec,nbot_eddy)
    
    ! Initialize turbulent mountain stress module
    do_tms = .false.         ! It is turned off for now
    if ( do_tms ) then
       call init_tms(r8, 4.0, karman, gravit, rair)
       call addfld ('TAUTMSX' ,'N/m2    ',1,    'A','Zonal turbulent mountain surface stress'  ,  phys_decomp)
       call addfld ('TAUTMSY' ,'N/m2    ',1,    'A','Meridional turbulent mountain surface stress', phys_decomp)
       call add_default ('TAUTMSX ', 1, ' ')
       call add_default ('TAUTMSY ', 1, ' ')
       if(masterproc)write(6,*)'Using turbulent mountain stress module'
    endif
    
    ! Initialize diffusion solver module
    call init_vdiff(r8, pcnst+pnats, rair, gravit, fieldlist, errstring)
    if (errstring.ne.'')call endrun(errstring)
    ! Select the fields which will be diffused using moist mixing ratios
    !       Default is everything -- use of dry mixing ratios for specific fields may
    !       be selected by changing them to false
    if(vdiff_select(fieldlist,'u').ne.'') call endrun(vdiff_select(fieldlist,'u'))
    if(vdiff_select(fieldlist,'v').ne.'') call endrun(vdiff_select(fieldlist,'v'))
    if(vdiff_select(fieldlist,'s').ne.'') call endrun(vdiff_select(fieldlist,'s'))
    do k = 1, pcnst+pnats
       if (cnst_get_type_byind(k).ne.'dry') then
          if (vdiff_select(fieldlist,'q',k).ne.'') call endrun(vdiff_select(fieldlist,'q',k))
       end if
    end do
    
    ! Diagnostic output fields
    do k = 1, pcnst+pnats
       vdiffnam(k) = 'VD'//cnst_name(k)
       if(k==1)vdiffnam(k) = 'VD01'    !**** compatibility with old code ****
       call addfld (vdiffnam(k),'kg/kg/s ',pver, 'A','Vertical diffusion of '//cnst_name(k),phys_decomp)
    end do
    call addfld('TKE     ','M2/S2   ',pverp, 'A','Turbulent Kinetic Energy',phys_decomp)
    call add_default (vdiffnam(1), 1, ' ')
    
  end subroutine vertical_diffusion_init
  
  !===============================================================================
  subroutine vertical_diffusion_ts_init
    !-----------------------------------------------------------------------
    ! timestep dependent setting, 
    ! at present only invokes upper bc code for molecular diffusion
    !-----------------------------------------------------------------------
    use molec_diff, only: init_timestep_molec_diff

    if (do_molec_diff) call init_timestep_molec_diff

  end subroutine vertical_diffusion_ts_init

  !===============================================================================
  subroutine vertical_diffusion_tend(                                     &
       ztodt    ,state    ,                               &
       taux     ,tauy     ,shflx    ,cflx     ,pblh     , &
       tpert    ,qpert    ,ustar    ,obklen   ,ptend    , &
       cldn     ,ocnfrac  ,landfrac ,sgh      ,pbuf     ,kvh ) 
    !-----------------------------------------------------------------------
    ! interface routine for vertical diffusion
    !-----------------------------------------------------------------------
    use physics_types,  only: physics_state, physics_ptend
    use history,        only: outfld
    use phys_buffer, only: pbuf_size_max, pbuf_fld,pbuf_old_tim_idx,pbuf_next_tim_idx,pbuf_times, pbuf_get_fld_idx
    use time_manager, only: is_first_step
    use geopotential, only: geopotential_dse
    ! The commented 'only' limiter from the following line acommodates broken pgf90 v.5.1.6
    use diffusion_solver ! , only: compute_vdiff, any, operator(.not.)
    use trb_mtn_stress, only: compute_tms
    use hb_diff, only: compute_hb_diff
    use wv_saturation, only: fqsatd
    use molec_diff, only: compute_molec_diff , vd_lu_qdecomp
    use constituents, only: qmincg
    use infnan

    !------------------------------Arguments--------------------------------
    real(r8), intent(in) :: taux(pcols)            ! x surface stress (N/m2)
    real(r8), intent(in) :: tauy(pcols)            ! y surface stress (N/m2)
    real(r8), intent(in) :: shflx(pcols)           ! surface sensible heat flux (w/m2)
    real(r8), intent(in) :: cflx(pcols,pcnst+pnats)! surface constituent flux (kg/m2/s)
    real(r8), intent(in) :: ztodt                  ! 2 delta-t
    real(r8), intent(in) :: cldn(pcols,pver)       ! new cloud fraction
    real(r8), intent(in) :: ocnfrac(pcols)         ! Ocean fraction
    real(r8), intent(in) :: landfrac(pcols)        ! Land fraction
    real(r8), intent(in) :: sgh(pcols)             ! standard deviation of orography
    type(physics_state), intent(in)  :: state      ! Physics state variables
    
    type(physics_ptend), intent(inout)                      :: ptend ! indivdual parameterization tendencies
    type(pbuf_fld), intent(inout), dimension(pbuf_size_max) :: pbuf  ! physics buffer
    
    real(r8), intent(out) :: pblh(pcols)              ! planetary boundary layer height
    real(r8), intent(out) :: tpert(pcols)             ! convective temperature excess
    real(r8), intent(out) :: qpert(pcols,pcnst+pnats) ! convective humidity and constituent excess
    real(r8), intent(out) :: ustar(pcols)             ! surface friction velocity
    real(r8), intent(out) :: obklen(pcols)            ! Obukhov length
    real(r8), intent(out) :: kvh(pcols,pverp)         ! eddy diffusivity for heat [m2/s]
    !
    !---------------------------Local storage-------------------------------
    integer :: lchnk                               ! chunk identifier
    integer :: ncol                                ! number of atmospheric columns
    integer :: i,k,m                               ! longitude,level,constituent indices

    real(r8) :: dtk(pcols,pver)                    ! T tendency from KE dissipation
    real(r8) :: tke(pcols,pverp)                   ! turbulent kinetic energy
    real(r8) :: cgs(pcols,pverp)                   ! counter-gradient star (cg/flux)
    real(r8) :: cgh(pcols,pverp)                   ! counter-gradient term for heat
    real(r8) :: rztodt                             ! 1./ztodt
    real(r8) :: ksrftms(pcols)                     ! turbulent mountain stress surface drag coefficient
    real(r8) :: tautmsx(pcols)                     ! u component of turbulent mountain stress
    real(r8) :: tautmsy(pcols)                     ! v component of turbulent mountain stress
    real(r8) :: tautotx(pcols)                     ! u component of total surface stress
    real(r8) :: tautoty(pcols)                     ! v component of total surface stress

    integer :: time_index                          ! time level index for physics buffer access
    real(r8) :: kvm(pcols,pverp)                   ! eddy diffusivity for momentum [m2/s]
    
    real(r8) :: kvq(pcols,pverp)                   ! diffusivity for constituents
    real(r8) :: th(pcols,pver)                     ! Potential temperature
    real(r8) :: topflx(pcols)                      ! molecular heat flux at top interface
    character(128) :: errstring                    ! error status for compute_vdiff

    !-----------------------------------------------------------------------
    rztodt = 1./ztodt
    lchnk = state%lchnk
    ncol  = state%ncol
    ! If first time step, initialize tke in pbuf at all time levels.
    ! This should be moved (to inti?) 
    if (is_first_step()) then
       do i = 1, pbuf_times
          pbuf(tke_idx)%fld_ptr(1,1:ncol,1:pverp,lchnk,i) = 0.01
       enddo
    endif
       
    ! All variables are modified by vertical diffusion

    ptend%name  = "vertical diffusion"
    ptend%lq(:) = .TRUE.
    ptend%ls    = .TRUE.
    ptend%lu    = .TRUE.
    ptend%lv    = .TRUE.

    !-----------------------------------------------------------------------
    !    Computation of turbulent mountain stress
    !-----------------------------------------------------------------------

    if ( do_tms ) then
       ! compute the turbulent mountain stress
       call compute_tms(  pcols   , pver       , ncol        , &
            state%u     , state%v  , state%t , state%pmid , state%exner , &
            state%zm    , sgh      , ksrftms, tautmsx, tautmsy, landfrac)
       tautotx(:ncol) = taux(:ncol) + tautmsx(:ncol)
       tautoty(:ncol) = tauy(:ncol) + tautmsy(:ncol)
    else
       ksrftms(:ncol) = 0.0
       tautotx(:ncol) = taux(:ncol)
       tautoty(:ncol) = tauy(:ncol)
    endif

    !-----------------------------------------------------------------------
    !    Computation of eddy diffusivities
    !-----------------------------------------------------------------------
    call t_startf('hb_diff')

    th(:ncol,:pver) = state%t(:ncol,:pver) * state%exner(:ncol,:pver)
    call compute_hb_diff(ncol      ,                   &
         th      ,state%t ,state%q ,state%zm,state%zi, &
         state%pmid,state%u,state%v,tautotx ,tautoty , &
         shflx   ,cflx    ,obklen  ,ustar   ,pblh    , &
         kvm     ,kvh     ,kvq     ,cgh     ,cgs     , &
         tpert   ,qpert   ,cldn    ,ocnfrac ,tke     )
    
    call t_stopf('hb_diff')
    !-----------------------------------------------------------------------
    !    Application of diffusivities
    !-----------------------------------------------------------------------

    ptend%q(:ncol,:,:) = state%q(:ncol,:,:)
    ptend%s(:ncol,:)   = state%s(:ncol,:)
    ptend%u(:ncol,:)   = state%u(:ncol,:)
    ptend%v(:ncol,:)   = state%v(:ncol,:)

    call t_startf('compute_vdiff')
    ! Call the diffusivity solver
    ! The final two arguments are optional function references to 
    ! constituent-independent and constituent-dependent moleculuar diffusivity routines
    ! There use allows the diffusion_solver module to be independent of CAM, and to be 
    ! used, for instance, by the Grenier-Bretherton PBL module.
    if (any(fieldlist)) call compute_vdiff(   state%lchnk ,                          &
         pcols         , pver               , pcnst+pnats , ncol      , state%pmid , &
         state%pint    , state%rpdel        , state%t     , ztodt     , taux       , &
         tauy          , shflx              , cflx        , ntop      , nbot       , &
         kvh           , kvm                , kvq         , cgs       , cgh        , &
         state%zi      , ksrftms            , qmincg      , fieldlist , &
         ptend%u       , ptend%v            , ptend%q     , ptend%s   , &
         tautmsx       , tautmsy            , dtk         , topflx    , errstring  , &
         do_molec_diff , compute_molec_diff , vd_lu_qdecomp )
    if (errstring.ne.'')call endrun(errstring)
 
    if (any(.not.fieldlist)) then
       if (do_molec_diff) then
          errstring = "Design flaw: dry vdiff not currently supported with molecular diffusion"
          call endrun(errstring)
       end if
       call compute_vdiff(  state%lchnk        , &
            pcols         , pver               , pcnst+pnats , ncol       , state%pmiddry , &
            state%pintdry , state%rpdeldry     , state%t     , ztodt      , taux          , &       
            tauy          , shflx              , cflx        , ntop       , nbot          , &       
            kvh           , kvm                , kvq         , cgs        , cgh           , &   
            state%zi      , ksrftms            , qmincg      , .not. fieldlist            , &
            ptend%u       , ptend%v            , ptend%q     , ptend%s    , &
            tautmsx       , tautmsy            , dtk         , topflx     , errstring     , &
            do_molec_diff , compute_molec_diff , vd_lu_qdecomp)
       if (errstring.ne.'')call endrun(errstring)
    end if
    call t_stopf('compute_vdiff')
    
    !-----------------------------------------------------------------------
    !    Diagnostics and output
    !-----------------------------------------------------------------------

    ! Convert the new profiles into vertical diffusion tendencies.
    ! Convert KE dissipative heat change into "temperature" tendency.
    ptend%s(:ncol,:)   = (ptend%s(:ncol,:)   - state%s(:ncol,:))   * rztodt
    ptend%u(:ncol,:)   = (ptend%u(:ncol,:)   - state%u(:ncol,:))   * rztodt
    ptend%v(:ncol,:)   = (ptend%v(:ncol,:)   - state%v(:ncol,:))   * rztodt
    ptend%q(:ncol,:,:) = (ptend%q(:ncol,:,:) - state%q(:ncol,:,:)) * rztodt

#ifdef PERGRO
    ! For pergro case, do not diffuse cloud water: replace with input values
!    ptend%lq(ixcldice) = .FALSE.
!    ptend%lq(ixcldliq) = .FALSE.
!    ptend%q(:ncol,:,ixcldice) = 0.
!    ptend%q(:ncol,:,ixcldliq) = 0.
#endif

    call outfld('TKE     ', tke  ,pcols,lchnk)
    call outfld ('PBLH    ',pblh ,pcols,lchnk)
    call outfld ('TPERT   ',tpert,pcols,lchnk)
    call outfld ('QPERT   ',qpert,pcols,lchnk)
    call outfld ('USTAR   ',ustar,pcols,lchnk)
    call outfld ('KVH     ',kvh,pcols,lchnk)
    call outfld ('KVM     ',kvm,pcols,lchnk)
    call outfld ('CGS     ',cgs,pcols,lchnk)
    dtk(:ncol,:) = dtk(:ncol,:)/cpair                ! normalize heating for history
    call outfld ('DTVKE   ',dtk,pcols,lchnk)
    dtk(:ncol,:) = ptend%s(:ncol,:)/cpair            ! normalize heating for history using dtk
    call outfld ('DTV     ',dtk  ,pcols,lchnk)
    call outfld ('DUV     ',ptend%u,pcols,lchnk)
    call outfld ('DVV     ',ptend%v,pcols,lchnk)
    do m = 1, pcnst+pnats
       call outfld(vdiffnam(m),ptend%q(1,1,m),pcols,lchnk)
    end do
    if ( do_tms ) then
       call outfld ('TAUTMSX', tautmsx, pcols, lchnk)
       call outfld ('TAUTMSY', tautmsy, pcols, lchnk)
    end if
    if (do_molec_diff) then
       call outfld ('TTPXMLC',topflx,pcols,lchnk)
    end if

    return
  end subroutine vertical_diffusion_tend

end module vertical_diffusion
