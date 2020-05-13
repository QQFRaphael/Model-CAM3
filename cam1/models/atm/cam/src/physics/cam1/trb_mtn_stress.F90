module trb_mtn_stress

  implicit none

  private          ! Make default type private to the module
  save
!
! Public interfaces
!
  public init_tms     ! Initialization
  public compute_tms  ! Full routine
!
! Private data
!
  integer, parameter :: r8 = selected_real_kind(12) ! 8 byte real
  real(r8), parameter :: z0fac  = 0.025    ! factor determining z_0 from orographic standard deviation
  real(r8), parameter :: z0max  = 100.     ! max value of z_0 for orography
  real(r8), parameter :: horomin= 10.      ! min value of subgrid orographic height for mountain stress
  real(r8), parameter :: dv2min = 0.01     ! minimum shear squared
  real(r8) :: oroconst          ! converts from standard deviation to height
  real(r8) :: karman     ! von Karman constant
  real(r8) :: gravit     ! Acceleration due to gravity
  real(r8) :: rair       ! Gas const for dry air

contains

!===============================================================================
  subroutine init_tms(kind, oro_in, karman_in, gravit_in, rair_in)
    
    integer, intent(in)  :: kind      ! kind of reals being passed in
    real(r8), intent(in) :: oro_in, karman_in, gravit_in, rair_in
    
    if ( kind .ne. r8 ) then
       write(6,*) 'KIND of reals passed to init_tms -- exiting.'
       stop 'compute_tms'
    endif

    oroconst = oro_in
    karman = karman_in
    gravit = gravit_in
    rair = rair_in
    
    return
  end subroutine init_tms

  subroutine compute_tms( pcols    , pver    , ncol        , &
            u           , v        , t       , pmid       , exner       , &
            zm          , sgh      , ksrf    , taux       , tauy    , landfrac)
!-----------------------------------------------------------------------
! Turbulent mountain stress parameterization
! 
! Returns the surface stress associated with subgrid mountains
! For points where the orographic variance is small (including ocean),
! the returned stress is zero. 
!-----------------------------------------------------------------------


!------------------------------Arguments--------------------------------

    integer, intent(in) :: pcols                  ! number of columns dimensioned
    integer, intent(in) :: pver                   ! number of model layers
    integer, intent(in) :: ncol                   ! number of columns actually used

    real(r8), intent(in) :: u(pcols,pver)         ! midpoint zonal wind
    real(r8), intent(in) :: v(pcols,pver)         ! midpoint meridional wind
    real(r8), intent(in) :: t(pcols,pver)         ! midpoint temperatures
    real(r8), intent(in) :: pmid (pcols,pver)     ! midpoint pressures
    real(r8), intent(in) :: exner(pcols,pver)     ! exner function
    real(r8), intent(in) :: zm   (pcols,pver)     ! midpoint height
    real(r8), intent(in) :: sgh(pcols)            ! standard deviation of orography
    real(r8), intent(in)  :: landfrac(pcols)          ! surface "drag" coefficient
    
    real(r8), intent(out) :: ksrf(pcols)          ! surface "drag" coefficient
    real(r8), intent(out) :: taux(pcols)          ! surface "drag" coefficient
    real(r8), intent(out) :: tauy(pcols)          ! surface "drag" coefficient
    
    !---------------------------Local storage-------------------------------
    integer  :: i                                 ! loop indexes
    integer  :: kb,kt                             ! bottom and top of source region
    
    real(r8) :: horo                              ! orographic height
    real(r8) :: z0oro                             ! orographic z0 momentum
    real(r8) :: dv2                               ! (delta v)**2
    real(r8) :: ri                                ! richardson number
    real(r8) :: stabfri                           ! stability function of richardson number
    real(r8) :: rho                               ! density
    real(r8) :: cd                                ! drag coefficient
    real(r8) :: vmag                              ! velocity magnitude

!---------------------------------------------------------------------------
       
    do i = 1, ncol

! determine subgrid orgraphic height (mean to peak)
       horo  = oroconst *sgh(i)

! no mountain stress if horo is too small
       if (horo < horomin) then
          ksrf(i) = 0.
       else

! determine z0m for orography
          z0oro = min(z0fac * horo, z0max)

! calculate neutral drag coefficient
          cd = ( karman / log((zm(i,pver) + z0oro)/z0oro) )**2

! calculate the Richardson number over 1st 2 layers
          kt  = pver-1
          kb  = pver
          dv2 = max((u(i,kt) - u(i,kb))**2 + (v(i,kt) - v(i,kb))**2, dv2min)
          ri  = 2. * gravit * (t(i,kt)*exner(i,kt) - t(i,kb)*exner(i,kb)) * (zm(i,kt)-zm(i,kb)) &
               / ((t(i,kt) + t(i,kb)) * dv2)

! calculate the stability function and modify the neutral drag cofficient
! should probably follow Louis et al (1982) but for now just 1 for ri<0, 0 for ri>1, linear in 1-ri between 0 and 1
          stabfri = max(0._r8,min(1._r8, 1. - ri))
          cd  = cd * stabfri

! compute density, velocity magnitude and stress using bottom level properties
          rho     = pmid(i,pver) / (rair * t(i,pver)) 
          vmag    = sqrt(u(i,pver)**2 + v(i,pver)**2)
          ksrf(i) = rho * cd * vmag * landfrac(i)
          taux(i) =  - ksrf(i) * u(i,pver)
          tauy(i) =  - ksrf(i) * v(i,pver)
       end if
    end do
    
    return
  end subroutine compute_tms

end module trb_mtn_stress
