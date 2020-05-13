module molec_diff
!---------------------------------------------------------------------------------
! Module to compute molecular diffusivity for various constituents
!
! Public interfaces:
!    init_molec_diff           initializes time independent coefficients
!    init_timestep_molec_diff  time-step initialization for molecular diffusivity
!    compute_molec_diff   computes constituent-independent terms for moleculuar diffusivity
!    vd_lu_qdecomp        computes constituent-dependent terms for moleculuar diffusivity and 
!                         updates terms in the triadiagonal matrix used for the implicit 
!                         solution of the diffusion equation
!
!---------------------------Code history--------------------------------
! modularized:       J. McCaa, September 2004
!---------------------------------------------------------------------------------

  implicit none
  private          ! Make default type private to the module
  save
!
! Public interfaces
!
  public init_molec_diff     ! Initialization
  public init_timestep_molec_diff
  public compute_molec_diff  ! Full routine
  public vd_lu_qdecomp

  integer, parameter :: r8 = selected_real_kind(12) ! 8 byte real
  
  real(r8), parameter :: km_fac = 3.55E-7 ! molecular viscosity constant
  real(r8), parameter :: pr_num = 1.       ! Prandtl number
  real(r8), parameter :: pwr = 2._r8/3._r8    ! exponentiation factor
  real(r8), parameter :: d0     = 1.52E20  ! diffusion factor m-1 s-1 molec sqrt(kg/kmol/K)
                                           ! Aerononmy, Part B, Banks and Kockarts (1973), p39
                                           ! note text cites 1.52E18 cm-1 ...

  real(r8) :: rair       ! Gas const for dry air
  real(r8) :: mw_dry     ! molecular weight of dry air
  real(r8) :: n_avog     ! Avogadro's number (molec/kmol)
  real(r8) :: gravit     
  real(r8) :: cpair
  real(r8) :: kbtz       ! boltzman constant

  integer :: ntop_molec  ! Top level to which molecular vertical diffusion is applied (=1).
  integer :: nbot_molec  ! Bottom level to which molecular vertical diffusion is applied.
  real(r8), allocatable :: mw_fac(:)          ! sqrt(1/M_q + 1/M_d) in constituent diffusivity
  
contains

!===============================================================================

  subroutine init_molec_diff(kind, ncnst, rair_in, ntop_molec_in, nbot_molec_in, &
       mw_dry_in, n_avog_in, gravit_in, cpair_in, kbtz_in)
    
    use constituents, only: cnst_mw
    use upper_bc,     only: ubc_init
    
    integer, intent(in)   :: kind      ! kind of reals being passed in
    integer, intent(in)   :: ncnst     ! number of constituents
    real(r8), intent(in)  :: rair_in
    integer, intent(in)   :: ntop_molec_in  ! Top level to which molecular vertical diffusion is applied (=1).
    integer, intent(in)   :: nbot_molec_in  ! Bottom level to which molecular vertical diffusion is applied.
    real(r8), intent(in)  :: mw_dry_in      ! molecular weight of dry air
    real(r8), intent(in)  :: n_avog_in      ! Avogadro's number (molec/kmol)
    real(r8), intent(in)  :: gravit_in
    real(r8), intent(in)  :: cpair_in
    real(r8), intent(in) :: kbtz_in
    integer  :: m                                  ! constituent index
    
    if ( kind .ne. r8 ) then
       write(6,*) 'KIND of reals passed to init_molec_diff -- exiting.'
       stop 'init_molec_diff'
    endif
    
    rair   = rair_in
    mw_dry = mw_dry_in
    n_avog = n_avog_in
    gravit = gravit_in
    cpair  = cpair_in
    kbtz = kbtz_in
    ntop_molec = ntop_molec_in
    nbot_molec = nbot_molec_in
    
    ! Initialize upper boundary condition variables
    call ubc_init()

    ! Molecular weight factor in constitutent diffusivity
    ! ***** FAKE THIS FOR NOW USING MOLECULAR WEIGHT OF DRY AIR FOR ALL TRACERS ****
    allocate(mw_fac(ncnst))
    do m = 1, ncnst
       mw_fac(m) = d0 * mw_dry * sqrt(1._r8/mw_dry + 1._r8/cnst_mw(m)) / n_avog
    end do

  end subroutine init_molec_diff

!===============================================================================

  subroutine init_timestep_molec_diff()
    !-----------------------------------------------------------------------
    ! timestep dependent setting, 
    !-----------------------------------------------------------------------
    use upper_bc,     only: ubc_timestep_init

    call ubc_timestep_init
    
  end subroutine init_timestep_molec_diff

!===============================================================================

  integer function compute_molec_diff(  lchnk       , &
       pcols              , pver       , ncnst      , ncol    , t      , pmid   , pint        , &
       zi                 , ztodt      , kvh        , kvm     , tint   , rhoi   , tmpi2       , &
       kq_scal            , ubc_t      , ubc_mmr    , dse_top , cc_top , cd_top , cnst_mw_out , &
       cnst_fixed_ubc_out , mw_fac_out , ntop_molec_out, nbot_molec_out )
    
    use upper_bc,     only: ubc_get_vals
    use constituents, only: cnst_mw, cnst_fixed_ubc
    
    integer, intent(in) :: pcols
    integer, intent(in) :: pver
    integer, intent(in) :: ncnst
    integer, intent(in) :: ncol                      ! number of atmospheric columns
    integer, intent(in)  :: lchnk                    ! chunk identifier
    real(r8), intent(in) :: t(pcols,pver)            ! temperature input
    real(r8), intent(in) :: pmid(pcols,pver)         ! midpoint pressures
    real(r8), intent(in) :: pint(pcols,pver+1)       ! interface pressures
    real(r8), intent(in) :: zi(pcols,pver+1)         ! interface heights
    real(r8), intent(in) :: ztodt                    ! 2 delta-t
    
    real(r8), intent(inout) :: kvh(pcols,pver+1)     ! diffusivity for heat
    real(r8), intent(inout) :: kvm(pcols,pver+1)     ! viscosity (diffusivity for momentum)
    real(r8), intent(inout) :: tint(pcols,pver+1)    ! interface temperature
    real(r8), intent(inout) :: rhoi(pcols,pver+1)    ! density (rho) at interfaces
    real(r8), intent(inout) :: tmpi2(pcols,pver+1)   ! dt*(g*rho)**2/dp at interfaces
    real(r8), intent(out)   :: kq_scal(pcols,pver+1) ! kq_fac*sqrt(T)*m_d/rho for molecular diffusivity
    real(r8), intent(out)   :: ubc_mmr(pcols,ncnst)  ! upper bndy mixing ratios (kg/kg)
    real(r8), intent(out)   :: cnst_mw_out(ncnst)
    logical,  intent(out)   :: cnst_fixed_ubc_out(ncnst)
    real(r8), intent(out)   :: mw_fac_out(ncnst)
    real(r8), intent(out)   :: dse_top(pcols)       ! dse on top boundary
    real(r8), intent(out)   :: cc_top(pcols)        ! lower diagonal at top interface
    real(r8), intent(out)   :: cd_top(pcols)        ! cc_top * dse ubc value
    integer,  intent(out)   :: ntop_molec_out   
    integer,  intent(out)   :: nbot_molec_out   

    ! Local variables
    real(r8) :: mkvisc                             ! molecular kinematic viscosity c*tint**(2/3)/rho
    integer  :: m                                  ! constituent index
    integer  :: i                                  ! column index
    integer  :: k                                  ! level index
    real(r8) :: ubc_t(pcols)                       ! upper bndy temperature (K)

    ! get upper boundary values
    call ubc_get_vals (lchnk, ncol, ntop_molec, pint, zi, ubc_t, ubc_mmr)

    ! These are already computed, just need to be copied for output
    cnst_mw_out(:ncnst) = cnst_mw(:ncnst)
    cnst_fixed_ubc_out(:ncnst) = cnst_fixed_ubc(:ncnst)
    mw_fac_out(:ncnst) = mw_fac(:ncnst)
    ntop_molec_out = ntop_molec
    nbot_molec_out = nbot_molec
    
    ! Density and related factors for moecular diffusion and ubc
    ! Always have a fixed upper boundary T if molecular diffusion is active
    tint (:ncol,ntop_molec) = ubc_t(:ncol)
    rhoi (:ncol,ntop_molec) = pint(:ncol,ntop_molec) / (rair*tint (:ncol,ntop_molec))
    tmpi2(:ncol,ntop_molec) = ztodt * (gravit*rhoi(:ncol,ntop_molec))**2 &
         / (pmid(:ncol,ntop_molec) - pint(:ncol,ntop_molec))
    
    ! Compute molecular kinematic viscosity, heat diffusivity and factor for constituent diffusivity
    kq_scal(:ncol,1:ntop_molec-1) = 0.0
    do k = ntop_molec, nbot_molec
       do i = 1, ncol
          mkvisc   = km_fac * tint(i,k)**pwr / rhoi(i,k)
          kvm(i,k) = kvm(i,k) + mkvisc
          kvh(i,k) = kvh(i,k) + mkvisc*pr_num
          kq_scal(i,k) = sqrt(tint(i,k)) / rhoi(i,k)
       end do
    end do
    kq_scal(:ncol,nbot_molec+1:pver+1) = 0.0
    
    ! top boundary condition for dry static energy
    dse_top(:ncol) = cpair * tint(:ncol,ntop_molec) + gravit * zi(:ncol,ntop_molec)
    ! top value of cc for dry static energy
    do i=1, ncol
       cc_top(i) = ztodt * gravit**2 * rhoi(i,ntop_molec) * km_fac * ubc_t(i)**pwr   &
            / ( (pint(i,2)-pint(i,1)) * (pmid(i,1) - pint(i,1)) )
    enddo
    cd_top(:ncol) = cc_top(:ncol) *  dse_top(:ncol)
    
    compute_molec_diff = 1
    return
  end function compute_molec_diff

!==============================================================================
  integer function vd_lu_qdecomp(                                         &
       pcols      , pver , ncol       , fixed_ubc  ,mw         ,ubc_mmr    ,             &
       kv         ,kq_scal    ,mw_facm    ,tmpi       ,rpdel      , &
       ca         ,cc         ,dnom       ,ze         ,rhoi       , &
       tint       ,ztodt      ,ntop_molec, nbot_molec, cd_top     )
!-----------------------------------------------------------------------
! Add the molecular diffusivity to the turbulent diffusivity for a consitutent.

! Update the superdiagonal (ca(k)), diagonal (cb(k)) and subdiagonal (cc(k))
! coeffs of the tridiagonal diffusion matrix, also ze and denominator.
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    integer, intent(in) :: pcols
    integer, intent(in) :: pver
    integer, intent(in)     :: ncol                 ! number of atmospheric columns

    logical,  intent(in)    :: fixed_ubc            ! fixed upper bdny cond flag
    real(r8), intent(in)    :: kv(pcols,pver+1)      ! vertical diffusion coefficients
    real(r8), intent(in)    :: kq_scal(pcols,pver+1) ! kq_fac*sqrt(T)*m_d/rho for molecular diffusivity
    real(r8), intent(in)    :: mw                   ! molecular weight for this constituent
    real(r8), intent(in)    :: ubc_mmr(pcols)       ! upper bndy mixing ratios (kg/kg)
    real(r8), intent(in)    :: mw_facm              ! sqrt(1/M_q + 1/M_d) for this constituent
    real(r8), intent(in)    :: tmpi(pcols,pver+1)    ! dt*(g/R)**2/dp*pi(k+1)/(.5*(tm(k+1)+tm(k))**2
    real(r8), intent(in)    :: rpdel(pcols,pver)    ! 1./pdel  (thickness bet interfaces)
    real(r8), intent(in)    :: rhoi(pcols,pver+1)    ! density at interfaces
    real(r8), intent(in)    :: tint(pcols,pver+1)    ! interface temperature
    real(r8), intent(in)    :: ztodt                ! 2 delta-t
    integer, intent(in) :: ntop_molec
    integer, intent(in) :: nbot_molec

    real(r8), intent(inout) :: ca(pcols,pver)       ! -upper diagonal
    real(r8), intent(inout) :: cc(pcols,pver)       ! -lower diagonal
    real(r8), intent(inout) :: dnom(pcols,pver)     ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))
                                                    ! 1./(b(k) - c(k)*e(k-1))
    real(r8), intent(inout) :: ze(pcols,pver)       ! term in tri-diag. matrix system

    real(r8), intent(out)   :: cd_top(pcols)        ! term for updating top level with ubc

!---------------------------Local workspace-----------------------------
    integer :: i                                    ! longitude index
    integer :: k, kp1                               ! vertical indicies

    real(r8) :: rghd(pcols,pver+1)                   ! (1/H_i - 1/H) * (rho*g)^(-1)
    real(r8) :: kmq(ncol)                           ! molecular diffusivity for constituent
    real(r8) :: wrk0(ncol)                          ! work variable
    real(r8) :: wrk1(ncol)                          ! work variable

    real(r8) :: cb(pcols,pver)       ! -diagonal
    real(r8) :: kvq(pcols,pver+1)     ! output vert. diff. coeff
!-----------------------------------------------------------------------
! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the 
! tridiagonal diffusion matrix. The diagonal elements  (cb=1+ca+cc) are
! a combination of ca and cc; they are not required by the solver.
!-----------------------------------------------------------------------
    call t_startf('vd_lu_qdecomp')

    kvq(:,:)  = 0._r8
    cd_top(:) = 0._r8

! Compute difference between scale heights of consituent and dry air
    do k = ntop_molec, nbot_molec
       do i = 1, ncol
          rghd(i,k)  = gravit  / (kbtz * n_avog * tint(i,k)) * (mw - mw_dry)
          rghd(i,k)  = ztodt * gravit * rhoi(i,k) * rghd(i,k) 
       enddo
    enddo

!-----------------------------------------------------------------------
! Molecular diffusion
!-----------------------------------------------------------------------
    do k = nbot_molec-1, ntop_molec, -1
       kp1 = k + 1
       kmq(:ncol)  = kq_scal(:ncol,kp1) * mw_facm
       wrk0(:ncol) = (kv(:ncol,kp1) + kmq(:ncol)) * tmpi(:ncol,kp1)
       wrk1(:ncol) = kmq(:ncol) * .5_r8 * rghd(:ncol,kp1)
!add species separation term
       ca(:ncol,k  )  = (wrk0(:ncol) - wrk1(:ncol)) * rpdel(:ncol,k)
       cc(:ncol,kp1)  = (wrk0(:ncol) + wrk1(:ncol)) * rpdel(:ncol,kp1)
       kvq(:ncol,kp1) = kmq(:ncol)
    end do
    if( fixed_ubc ) then
       cc(:ncol,ntop_molec)  = kq_scal(:ncol,ntop_molec) * mw_facm &
            * (tmpi(:ncol,ntop_molec) + rghd(:ncol,ntop_molec)) &
            * rpdel(:ncol,ntop_molec)
    end if
! Calculate diagonal elements
    do k = nbot_molec-1, ntop_molec+1, -1
       kp1 = k + 1
       cb(:ncol,k) = 1._r8 + ca(:ncol,k) + cc(:ncol,k) &
            + rpdel(:ncol,k) * (kvq(:ncol,kp1)*rghd(:ncol,kp1) &
            - kvq(:ncol,k)*rghd(:ncol,k))
       kvq(:ncol,kp1) = kv(:ncol,kp1) + kvq(:ncol,kp1)
    end do

    k   = ntop_molec
    kp1 = k + 1
    if( fixed_ubc ) then
       cb(:ncol,k) = 1._r8 + ca(:ncol,k) &
            + rpdel(:ncol,k) * kvq(:ncol,kp1)*rghd(:ncol,kp1) &
            + kq_scal(:ncol,ntop_molec) * mw_facm &
            * (tmpi(:ncol,ntop_molec) - rghd(:ncol,ntop_molec)) &
            * rpdel(:ncol,ntop_molec)
    else
       cb(:ncol,k) = 1._r8 + ca(:ncol,k) &
            + rpdel(:ncol,k) * kvq(:ncol,kp1)*rghd(:ncol,kp1)
    end if

    k   = nbot_molec
    cb(:ncol,k) = 1._r8 + cc(:ncol,k) + ca(:ncol,k) &
         - rpdel(:ncol,k) * kvq(:ncol,k)*rghd(:ncol,k)
    do k = 1, nbot_molec+1, -1
       cb(:ncol,k) = 1._r8 + ca(:ncol,k) + cc(:ncol,k)
    end do


! Compute term for updating top level mixing ratio for ubc
    if( fixed_ubc ) then
       cd_top(:ncol) = cc(:ncol,ntop_molec)*ubc_mmr(:ncol)
    end if

!-----------------------------------------------------------------------
! Calculate e(k).  This term is 
! required in solution of tridiagonal matrix defined by implicit diffusion eqn.
!-----------------------------------------------------------------------
    do k = nbot_molec, ntop_molec+1, -1
       dnom(:ncol,k) = 1._r8/ (cb(:ncol,k) - ca(:ncol,k)*ze(:ncol,k+1))
       ze(:ncol,k)   = cc(:ncol,k)*dnom(:ncol,k)
    end do
    k = ntop_molec
    dnom(:ncol,k) = 1._r8/ (cb(:ncol,k) - ca(:ncol,k)*ze(:ncol,k+1))

    vd_lu_qdecomp = 1
    call t_stopf('vd_lu_qdecomp')
    return

  end function vd_lu_qdecomp

end module molec_diff
