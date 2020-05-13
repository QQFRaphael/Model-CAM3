#include <misc.h>

!-------------------------------------------------------------------------------
!physics data types module
!-------------------------------------------------------------------------------
module physics_types

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, pver
  use constituents, only: ppcnst, qmin, cnst_name
  use geopotential, only: geopotential_dse
  use physconst,    only: zvir, gravit, cpair, rair
  use dycore,       only: dycore_is

  implicit none
  private          ! Make default type private to the module

  logical, parameter :: adjust_te = .FALSE.

! Public types:

  public physics_state
  public physics_tend
  public physics_ptend

! Public interfaces

  public physics_update
  public physics_ptend_reset
  public physics_ptend_init
  public physics_dme_adjust  ! adjust dry mass and energy for change in water
                             ! cannot be applied to eul or sld dycores
  public physics_state_copy  ! copy a state type
  public physics_ptend_sum   ! add 2 ptend types
  public physics_tend_init   ! initialize a tend type

  public set_state_pdry      ! calculate dry air masses in state variable
  public set_wet_to_dry
  public set_dry_to_wet

!-------------------------------------------------------------------------------
  type physics_state
     integer                                     :: &
          lchnk,   &! chunk index
          ncol      ! number of active columns
     real(r8), dimension(pcols)                  :: &
          lat,     &! latitude (radians)
          lon,     &! longitude (radians)
          ps,      &! surface pressure
          psdry,   &! dry surface pressure
          phis      ! surface geopotential
     real(r8), dimension(pcols,pver)             :: &
          t,       &! temperature (K)
          u,       &! zonal wind (m/s)
          v,       &! meridional wind (m/s)
          s,       &! dry static energy
          omega,   &! vertical pressure velocity (Pa/s) 
          pmid,    &! midpoint pressure (Pa) 
          pmiddry, &! midpoint pressure dry (Pa) 
          pdel,    &! layer thickness (Pa)
          pdeldry, &! layer thickness dry (Pa)
          rpdel,   &! reciprocal of layer thickness (Pa)
          rpdeldry,&! recipricol layer thickness dry (Pa)
          lnpmid,  &! ln(pmid)
          lnpmiddry,&! log midpoint pressure dry (Pa) 
          exner,   &! inverse exner function w.r.t. surface pressure (ps/p)^(R/cp)
          zm        ! geopotential height above surface at midpoints (m)

     real(r8), dimension(pcols,pver,ppcnst) :: &
          q         ! constituent mixing ratio (kg/kg moist or dry air depending on type)

     real(r8), dimension(pcols,pver+1)           :: &
          pint,    &! interface pressure (Pa)
          pintdry, &! interface pressure dry (Pa) 
          lnpint,  &! ln(pint)
          lnpintdry,&! log interface pressure dry (Pa) 
          zi        ! geopotential height above surface at interfaces (m)

     real(r8), dimension(pcols)                  :: &
          te_ini,  &! vertically integrated total (kinetic + static) energy of initial state
          te_cur,  &! vertically integrated total (kinetic + static) energy of current state
          tw_ini,  &! vertically integrated total water of initial state
          tw_cur    ! vertically integrated total water of new state
     integer :: count ! count of values with significant energy or water imbalances
     
  end type physics_state

!-------------------------------------------------------------------------------
  type physics_tend
     real(r8), dimension(pcols,pver)             :: dtdt, dudt, dvdt
     real(r8), dimension(pcols     )             :: flx_net
     real(r8), dimension(pcols)                  :: &
          te_tnd,  &! cumulative boundary flux of total energy
          tw_tnd    ! cumulative boundary flux of total water
  end type physics_tend

!-------------------------------------------------------------------------------
! This is for tendencies returned from individual parameterizations
  type physics_ptend
     character*24 :: name    ! name of parameterization which produced tendencies.

     logical ::             &
          ls,               &! true if dsdt is returned
          lu,               &! true if dudt is returned
          lv,               &! true if dvdt is returned
          lq(ppcnst)         ! true if dqdt() is returned

     integer ::             &
          top_level,        &! top level index for which nonzero tendencies have been set
          bot_level          ! bottom level index for which nonzero tendencies have been set

     real(r8), dimension(pcols,pver)             :: &
          s,                &! heating rate (J/kg/s)
          u,                &! u momentum tendency (m/s/s)
          v                  ! v momentum tendency (m/s/s)
     real(r8), dimension(pcols,pver,ppcnst) :: &
          q                  ! consituent tendencies (kg/kg/s)

! boundary fluxes
     real(r8), dimension(pcols) ::&
          hflux_srf,     &! net heat flux at surface (W/m2)
          hflux_top,     &! net heat flux at top of model (W/m2)
          taux_srf,      &! net zonal stress at surface (Pa)
          taux_top,      &! net zonal stress at top of model (Pa)
          tauy_srf,      &! net meridional stress at surface (Pa)
          tauy_top        ! net meridional stress at top of model (Pa)
     real(r8), dimension(pcols,ppcnst) ::&
          cflx_srf,      &! constituent flux at surface (kg/m2/s)
          cflx_top        ! constituent flux top of model (kg/m2/s)

  end type physics_ptend


!===============================================================================
contains
!===============================================================================

!===============================================================================
  subroutine physics_update(state, tend, ptend, dt)
!-----------------------------------------------------------------------
! Update the state and or tendency structure with the parameterization tendencies
!-----------------------------------------------------------------------
    use geopotential, only: geopotential_dse
    use physconst,    only: cpair, gravit, rair, zvir
    use constituents, only: cnst_get_ind
#if ( defined SCAM )
#include <max.h>
    use scamMod, only: switch
#endif
!------------------------------Arguments--------------------------------
    type(physics_ptend), intent(inout)  :: ptend   ! Parameterization tendencies

    type(physics_state), intent(inout)  :: state   ! Physics state variables
    type(physics_tend ), intent(inout)  :: tend    ! Physics tendencies

    real(r8), intent(in) :: dt                     ! time step
!
!---------------------------Local storage-------------------------------
    integer :: i,k,m                               ! column,level,constituent indices
    integer :: ixcldice, ixcldliq                  ! indices for CLDICE and CLDLIQ
    integer :: ncol                                ! number of columns
    character*40 :: name    ! param and tracer name for qneg3
!-----------------------------------------------------------------------
#if ( defined SCAM )
    ! The column radiation model does not update the state
    if(switch(CRM_SW+1)) return
#endif
    ncol = state%ncol

! Update u,v fields
    if(ptend%lu) then
       do k = ptend%top_level, ptend%bot_level
          do i = 1, ncol
             state%u  (i,k) = state%u  (i,k) + ptend%u(i,k) * dt
             tend%dudt(i,k) = tend%dudt(i,k) + ptend%u(i,k)
          end do
       end do
    end if

    if(ptend%lv) then
       do k = ptend%top_level, ptend%bot_level
          do i = 1, ncol
             state%v  (i,k) = state%v  (i,k) + ptend%v(i,k) * dt
             tend%dvdt(i,k) = tend%dvdt(i,k) + ptend%v(i,k)
          end do
       end do
    end if

! Update dry static energy
    if(ptend%ls) then
       do k = ptend%top_level, ptend%bot_level
          do i = 1, ncol
             state%s(i,k)   = state%s(i,k)   + ptend%s(i,k) * dt
             tend%dtdt(i,k) = tend%dtdt(i,k) + ptend%s(i,k)/cpair
          end do
       end do
    end if

! Update constituents, all schemes use time split q: no tendency kept
    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)
    do m = 1, ppcnst
       if(ptend%lq(m)) then
          do k = ptend%top_level, ptend%bot_level
             do i = 1,ncol
                state%q(i,k,m) = state%q(i,k,m) + ptend%q(i,k,m) * dt
             end do
          end do
! now test for mixing ratios which are too small
          name = trim(ptend%name) // '/' // trim(cnst_name(m))
          call qneg3(trim(name), state%lchnk, ncol, pcols, pver, m, m, qmin(m), state%q(1,1,m))
       end if
    end do

! special test for cloud water
    if(ptend%lq(ixcldliq)) then
       if (ptend%name == 'stratiform' .or. ptend%name == 'cldwat'  ) then
#ifdef PERGRO
          where (state%q(:ncol,:pver,ixcldliq) < 1.e-12)
             state%q(:ncol,:pver,ixcldliq) = 0.
          end where
#endif
       else if (ptend%name == 'convect_deep') then
          where (state%q(:ncol,:pver,ixcldliq) < 1.e-36)
             state%q(:ncol,:pver,ixcldliq) = 0.
          end where
       end if
    end if
    if(ptend%lq(ixcldice)) then
       if (ptend%name == 'stratiform' .or. ptend%name == 'cldwat'  ) then
#ifdef PERGRO
          where (state%q(:ncol,:pver,ixcldice) < 1.e-12)
             state%q(:ncol,:pver,ixcldice) = 0.
          end where
#endif
       else if (ptend%name == 'convect_deep') then
          where (state%q(:ncol,:pver,ixcldice) < 1.e-36)
             state%q(:ncol,:pver,ixcldice) = 0.
          end where
       end if
    end if

! Derive new temperature and geopotential fields if heating or water tendency not 0.
    if (ptend%ls .or. ptend%lq(1)) then
       call geopotential_dse(                                                                    &
            state%lnpint, state%lnpmid, state%pint  , state%pmid  , state%pdel  , state%rpdel  , &
            state%s     , state%q(1,1,1),state%phis , rair        , gravit      , cpair        , &
            zvir        , state%t     , state%zi    , state%zm    , ncol         )
    end if

! Reset all parameterization tendency flags to false
    call physics_ptend_reset(ptend)

  end subroutine physics_update


!===============================================================================
  subroutine physics_ptend_sum(ptend, ptend_sum, state)
!-----------------------------------------------------------------------
! Add ptend fields to ptend_sum for ptend logical flags = .true.
! Where ptend logical flags = .false, don't change ptend_sum
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    type(physics_ptend), intent(in)     :: ptend   ! New parameterization tendencies
    type(physics_ptend), intent(inout)  :: ptend_sum   ! Sum of incoming ptend_sum and ptend
    type(physics_state), intent(in)     :: state   ! New parameterization tendencies

!---------------------------Local storage-------------------------------
    integer :: i,k,m                               ! column,level,constituent indices
    integer :: ncol                                ! number of columns

!-----------------------------------------------------------------------
    ncol = state%ncol
    

! Update u,v fields
    if(ptend%lu) then
       ptend_sum%lu = .true.
       do i = 1, ncol
          do k = ptend%top_level, ptend%bot_level
             ptend_sum%u(i,k) = ptend_sum%u(i,k) + ptend%u(i,k)
          end do
          ptend_sum%taux_srf(i) = ptend_sum%taux_srf(i) + ptend%taux_srf(i)
          ptend_sum%taux_top(i) = ptend_sum%taux_top(i) + ptend%taux_top(i)
       end do
    end if

    if(ptend%lv) then
       ptend_sum%lv = .true.
       do i = 1, ncol
          do k = ptend%top_level, ptend%bot_level
             ptend_sum%v(i,k) = ptend_sum%v(i,k) + ptend%v(i,k)
          end do
          ptend_sum%tauy_srf(i) = ptend_sum%tauy_srf(i) + ptend%tauy_srf(i)
          ptend_sum%tauy_top(i) = ptend_sum%tauy_top(i) + ptend%tauy_top(i)
       end do
    end if


    if(ptend%ls) then
       ptend_sum%ls = .true.
       do i = 1, ncol
          do k = ptend%top_level, ptend%bot_level
             ptend_sum%s(i,k) = ptend_sum%s(i,k) + ptend%s(i,k)
          end do
          ptend_sum%hflux_srf(i) = ptend_sum%hflux_srf(i) + ptend%hflux_srf(i)
          ptend_sum%hflux_top(i) = ptend_sum%hflux_top(i) + ptend%hflux_top(i)
       end do
    end if

! Update constituents
    do m = 1, ppcnst
       if(ptend%lq(m)) then
          ptend_sum%lq(m) = .true.
          do i = 1,ncol
             do k = ptend%top_level, ptend%bot_level
                ptend_sum%q(i,k,m) = ptend_sum%q(i,k,m) + ptend%q(i,k,m)
             end do
             ptend_sum%cflx_srf(i,m) = ptend_sum%cflx_srf(i,m) + ptend%cflx_srf(i,m)
             ptend_sum%cflx_top(i,m) = ptend_sum%cflx_top(i,m) + ptend%cflx_top(i,m)
          end do
       end if
    end do


  end subroutine physics_ptend_sum

!===============================================================================
  subroutine physics_ptend_reset(ptend)
!-----------------------------------------------------------------------
! Reset the parameterization tendency structure to "empty"
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    type(physics_ptend), intent(inout)  :: ptend   ! Parameterization tendencies
!-----------------------------------------------------------------------
    integer :: m             ! Index for constiuent
!-----------------------------------------------------------------------

    if(ptend%ls) then
       ptend%s = 0.
       ptend%hflux_srf = 0.
       ptend%hflux_top = 0.
    endif
    if(ptend%lu) then
       ptend%u = 0.
       ptend%taux_srf = 0.
       ptend%taux_top = 0.
    endif
    if(ptend%lv) then
       ptend%v = 0.
       ptend%tauy_srf = 0.
       ptend%tauy_top = 0.
    endif
    do m = 1, ppcnst
       if(ptend%lq(m)) then
          ptend%q(:,:,m) = 0.
          ptend%cflx_srf(:,m) = 0.
          ptend%cflx_top(:,m) = 0.
       endif
    end do

    ptend%name  = "none"
    ptend%lq(:) = .FALSE.
    ptend%ls    = .FALSE.
    ptend%lu    = .FALSE.
    ptend%lv    = .FALSE.

    ptend%top_level = 1
    ptend%bot_level = pver

    return
  end subroutine physics_ptend_reset

!===============================================================================
  subroutine physics_ptend_init(ptend)
!-----------------------------------------------------------------------
! Initialize the parameterization tendency structure to "empty"
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    type(physics_ptend), intent(inout)  :: ptend   ! Parameterization tendencies
!-----------------------------------------------------------------------
    ptend%name  = "none"
    ptend%lq(:) = .true.
    ptend%ls    = .true.
    ptend%lu    = .true.
    ptend%lv    = .true.

    call physics_ptend_reset(ptend)

    return
  end subroutine physics_ptend_init

!===============================================================================
  subroutine physics_dme_adjust(state, tend, qini, dt)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Adjust the dry mass in each layer back to the value of physics input state
    ! 
    ! Method: Conserve the integrated mass, momentum and total energy in each layer
    !         by scaling the specific mass of consituents, specific momentum (velocity)
    !         and specific total energy by the relative change in layer mass. Solve for
    !         the new temperature by subtracting the new kinetic energy from total energy
    !         and inverting the hydrostatic equation
    !
    !         The mass in each layer is modified, changing the relationship of the layer 
    !         interfaces and midpoints to the surface pressure. The result is no longer in 
    !         the original hybrid coordinate. 
    !
    !         This procedure cannot be applied to the "eul" or "sld" dycores because they
    !         require the hybrid coordinate.
    ! 
    ! Author: Byron Boville

    ! !REVISION HISTORY:
    !   03.03.28  Boville    Created, partly from code by Lin in p_d_adjust
    ! 
    !-----------------------------------------------------------------------

    use constituents, only : cnst_get_type_byind

    implicit none
    !
    ! Arguments
    !
    type(physics_state), intent(inout) :: state
    type(physics_tend ), intent(inout) :: tend
    real(r8),            intent(in   ) :: qini(pcols,pver)    ! initial specific humidity
    real(r8),            intent(in   ) :: dt                  ! model physics timestep
    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: lchnk         ! chunk identifier
    integer  :: ncol          ! number of atmospheric columns
    integer  :: i,k,m         ! Longitude, level indices
    real(r8) :: fdq(pcols)    ! mass adjustment factor
    real(r8) :: te(pcols)     ! total energy in a layer
    real(r8) :: utmp(pcols)   ! temp variable for recalculating the initial u values
    real(r8) :: vtmp(pcols)   ! temp variable for recalculating the initial v values
    !
    !-----------------------------------------------------------------------
    ! verify that the dycore is not sld or eul
    if (dycore_is('SLD') .or. dycore_is('EUL')) return

    lchnk = state%lchnk
    ncol  = state%ncol

    ! adjust dry mass in each layer back to input value, while conserving
    ! constituents, momentum, and total energy
    do k = 1, pver

       ! adjusment factor is just change in water vapor
       fdq(:ncol) = 1. + state%q(:ncol,k,1) - qini(:ncol,k)

       ! adjust constituents to conserve mass in each layer
       do m = 1, ppcnst
          state%q(:ncol,k,m) = state%q(:ncol,k,m) / fdq(:ncol)
       end do

       if (adjust_te) then
          ! compute specific total energy of unadjusted state (J/kg)
          te(:ncol) = state%s(:ncol,k) + 0.5*(state%u(:ncol,k)**2 + state%v(:ncol,k)**2) 

          ! recompute initial u,v from the new values and the tendencies
          utmp(:ncol) = state%u(:ncol,k) - dt * tend%dudt(:ncol,k)
          vtmp(:ncol) = state%v(:ncol,k) - dt * tend%dvdt(:ncol,k)
          ! adjust specific total energy and specific momentum (velocity) to conserve each
          te     (:ncol)   = te     (:ncol)     / fdq(:ncol)
          state%u(:ncol,k) = state%u(:ncol,k  ) / fdq(:ncol)
          state%v(:ncol,k) = state%v(:ncol,k  ) / fdq(:ncol)
          ! compute adjusted u,v tendencies
          tend%dudt(:ncol,k) = (state%u(:ncol,k) - utmp(:ncol)) / dt
          tend%dvdt(:ncol,k) = (state%v(:ncol,k) - vtmp(:ncol)) / dt

          ! compute adjusted static energy
          state%s(:ncol,k) = te(:ncol) - 0.5*(state%u(:ncol,k)**2 + state%v(:ncol,k)**2)
       end if

! compute new total pressure variables
       state%pdel  (:ncol,k  ) = state%pdel(:ncol,k  ) * fdq(:ncol)
       state%pint  (:ncol,k+1) = state%pint(:ncol,k  ) + state%pdel(:ncol,k)
       state%lnpint(:ncol,k+1) = log(state%pint(:ncol,k+1))
       state%rpdel (:ncol,k  ) = 1./ state%pdel(:ncol,k  )
    end do

! compute new T,z from new s,q,dp
    if (adjust_te) then
       call geopotential_dse(state%lnpint, state%lnpmid  , state%pint ,  &
            state%pmid  , state%pdel    , state%rpdel,  &
            state%s     , state%q(1,1,1), state%phis , rair, gravit, cpair, zvir, &
            state%t     , state%zi      , state%zm   , ncol      )
    end if

  end subroutine physics_dme_adjust
!-----------------------------------------------------------------------

!===============================================================================
  subroutine physics_state_copy(state_in, state_out)
    
    use ppgrid,       only: pver, pverp
    use constituents, only: ppcnst, cnst_need_pdeldry

    implicit none

    !
    ! Arguments
    !
    type(physics_state), intent(in) :: state_in
    type(physics_state), intent(out) :: state_out

    !
    ! Local variables
    !
    integer i, k, m, ncol


    ncol = state_in%ncol
    
    state_out%lchnk = state_in%lchnk
    state_out%ncol  = state_in%ncol
    state_out%count = state_in%count 

    do i = 1, ncol
       state_out%lat(i)    = state_in%lat(i)
       state_out%lon(i)    = state_in%lon(i)
       state_out%ps(i)     = state_in%ps(i)
       state_out%phis(i)   = state_in%phis(i)
       state_out%te_ini(i) = state_in%te_ini(i) 
       state_out%te_cur(i) = state_in%te_cur(i) 
       state_out%tw_ini(i) = state_in%tw_ini(i) 
       state_out%tw_cur(i) = state_in%tw_cur(i) 
    end do

    do k = 1, pver
       do i = 1, ncol
          state_out%t(i,k)         = state_in%t(i,k) 
          state_out%u(i,k)         = state_in%u(i,k) 
          state_out%v(i,k)         = state_in%v(i,k) 
          state_out%s(i,k)         = state_in%s(i,k) 
          state_out%omega(i,k)     = state_in%omega(i,k) 
          state_out%pmid(i,k)      = state_in%pmid(i,k) 
          state_out%pdel(i,k)      = state_in%pdel(i,k) 
          state_out%rpdel(i,k)     = state_in%rpdel(i,k) 
          state_out%lnpmid(i,k)    = state_in%lnpmid(i,k) 
          state_out%exner(i,k)     = state_in%exner(i,k) 
          state_out%zm(i,k)        = state_in%zm(i,k)
       end do
    end do

    do k = 1, pverp
       do i = 1, ncol
          state_out%pint(i,k)      = state_in%pint(i,k) 
          state_out%lnpint(i,k)    = state_in%lnpint(i,k) 
          state_out%zi(i,k)        = state_in% zi(i,k) 
       end do
    end do


    if ( cnst_need_pdeldry ) then
       do i = 1, ncol
          state_out%psdry(i)  = state_in%psdry(i) 
       end do
       do k = 1, pver
          do i = 1, ncol
             state_out%lnpmiddry(i,k) = state_in%lnpmiddry(i,k) 
             state_out%pmiddry(i,k)   = state_in%pmiddry(i,k) 
             state_out%pdeldry(i,k)   = state_in%pdeldry(i,k) 
             state_out%rpdeldry(i,k)  = state_in%rpdeldry(i,k) 
          end do
       end do
       do k = 1, pverp
          do i = 1, ncol
             state_out%pintdry(i,k)   = state_in%pintdry(i,k)
             state_out%lnpintdry(i,k) = state_in%lnpintdry(i,k) 
          end do
       end do
    endif !cnst_need_pdeldry

    do m = 1, ppcnst
       do k = 1, pver
          do i = 1, ncol
             state_out%q(i,k,m) = state_in%q(i,k,m) 
          end do
       end do
    end do

  end  subroutine physics_state_copy
!===============================================================================

  subroutine physics_tend_init(tend)
    
    use ppgrid,       only:  pver
    use constituents, only: ppcnst
    
    implicit none
    
    !
    ! Arguments
    !
    type(physics_tend), intent(inout) :: tend

    !
    ! Local variables
    !

    tend%dtdt  = 0.
    tend%dudt = 0.
    tend%dvdt  = 0.
    tend%flx_net   = 0.
    tend%te_tnd  = 0.
    tend%tw_tnd   = 0.
    
end subroutine physics_tend_init

!===============================================================================

subroutine set_state_pdry (state,pdeld_calc)

  use ppgrid,  only: pver
  use pmgrid,  only: plev, plevp

  implicit none

#include <comhyb.h>

  type(physics_state), intent(inout) :: state
  logical, optional, intent(in) :: pdeld_calc    !  .true. do calculate pdeld [default]
                                                 !  .false. don't calculate pdeld 
  integer ncol
  integer i, k
  logical do_pdeld_calc

  if ( present(pdeld_calc) ) then
     do_pdeld_calc = pdeld_calc
  else
     do_pdeld_calc = .true.
  endif
  
  ncol = state%ncol

  state%psdry(:ncol) = ps0 * hyai(1) + state%ps(:ncol) * hybi(1)
  state%pintdry(:ncol,1) = ps0 * hyai(1) + state%ps(:ncol) * hybi(1)

  if (do_pdeld_calc)  then
     do k = 1, pver
        state%pdeldry(:ncol,k) = state%pdel(:ncol,k)*(1.-state%q(:ncol,k,1))
     end do
  endif
  do k = 1, pver
     state%pintdry(:ncol,k+1) = state%pintdry(:ncol,k)+state%pdeldry(:ncol,k)
     state%pmiddry(:ncol,k) = (state%pintdry(:ncol,k+1)+state%pintdry(:ncol,k))/2.
     state%psdry(:ncol) = state%psdry(:ncol) + state%pdeldry(:ncol,k)
  end do

  state%rpdeldry(:ncol,:) = 1./state%pdeldry(:ncol,:)
  state%lnpmiddry(:ncol,:) = log(state%pmiddry(:ncol,:))
  state%lnpintdry(:ncol,:) = log(state%pintdry(:ncol,:))

end subroutine set_state_pdry 

!===============================================================================

subroutine set_wet_to_dry (state)

  use constituents,  only: ppcnst, cnst_type

  type(physics_state), intent(inout) :: state

  integer m, ncol
  
  ncol = state%ncol

  do m = 1,ppcnst
     if (cnst_type(m).eq.'dry') then
        state%q(:ncol,:,m) = state%q(:ncol,:,m)*state%pdel(:ncol,:)/state%pdeldry(:ncol,:)
     endif
  end do

end subroutine set_wet_to_dry 

!===============================================================================

subroutine set_dry_to_wet (state)

  use constituents,  only: ppcnst, cnst_type

  type(physics_state), intent(inout) :: state

  integer m, ncol
  
  ncol = state%ncol

  do m = 1,ppcnst
     if (cnst_type(m).eq.'dry') then
        state%q(:ncol,:,m) = state%q(:ncol,:,m)*state%pdeldry(:ncol,:)/state%pdel(:ncol,:)
     endif
  end do

end subroutine set_dry_to_wet

!===============================================================================

end module physics_types
