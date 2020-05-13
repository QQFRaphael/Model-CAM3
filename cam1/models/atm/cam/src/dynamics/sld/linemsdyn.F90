#include <misc.h>
#include <params.h>

subroutine linemsdyn(lat     ,ps      ,u3      ,u3m1    ,v3      , &
                     v3m1    ,t3      ,t3m1    ,q3      ,etadot  , &
                     etamid  ,ztodt   ,                   &
                     vcour   ,vmax    ,vmaxt   ,detam   , &
                     ed1     ,fu      ,fv      ,lnpssld ,prhssld , &
                     tarrsld ,parrsld ,t2      ,div     ,tl      , &
                     tm      ,ql      ,qm      ,dpsl    ,dpsm    , &
                     phis    ,phisl   ,phism   ,omga    , &
                     urhs    ,vrhs    , &
                     trhs    ,prhs    ,nlon    ,cwava   ,flx_net)
!-----------------------------------------------------------------------
!
! Purpose:
! Driver for non-linear dynamics computations (in grid-point space)
!
! Original version:  CCM1
!
!-----------------------------------------------------------------------
!
! $Id: linemsdyn.F90,v 1.17.2.5 2005/01/27 18:02:31 olson Exp $
! $Author: olson $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: plon, plev, plevp
  use constituents, only: pcnst, pnats
  use pspect,       only: ptrn, ptrm
  use scanslt,      only: engy1lat
  use commap,       only: w, clat
  use grmult,       only: grmult_run
  use history, only: outfld
  use dynconst, only: omega
  use time_manager, only: get_step_size, is_first_step

  implicit none

#include <comctl.h>
#include <comhyb.h>
#include <comlun.h>

!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: lat                    ! latitude index for S->N storage
  real(r8), intent(in)   :: ps      (plon)         ! surface pressure (time n)
  real(r8), intent(in)   :: u3      (plon,plev)    ! u-wind (time n)
  real(r8), intent(inout):: u3m1    (plon,plev)    ! u-wind (time n-1)
  real(r8), intent(in)   :: v3      (plon,plev)    ! v-wind (time n)
  real(r8), intent(inout):: v3m1    (plon,plev)    ! v-wind (time n-1)
  real(r8), intent(in)   :: t3      (plon,plev)    ! temperature (time n)
  real(r8), intent(inout):: t3m1    (plon,plev)    ! temperature (time n-1)
  real(r8), intent(in)   :: q3      (plon,plev,pcnst+pnats)! constituents
  real(r8), intent(inout):: etadot  (plon,plevp)   ! vertical motion (3-d used by slt)
  real(r8), intent(in)   :: etamid  (plev)         ! midpoint values of eta (a+b)
  real(r8), intent(in)   :: ztodt                  ! 2*timestep unless nstep = 0
  real(r8), intent(out)  :: vcour   (plev)         ! maximum Courant number in vert.
  real(r8), intent(out)  :: vmax    (plev)         ! maximum wind speed squared (m^2/s^2)
  real(r8), intent(out)  :: vmaxt   (plev)         ! maximum truncated wind speed (m^2/s^2)
  real(r8), intent(in)   :: detam   (plev)         ! intervals between vert.levels
  real(r8), intent(inout):: ed1     (plon,plev)    ! etadot*dp/deta (for SLD)
  real(r8), intent(inout):: fu      (plon,plev)    ! nonlinear term - u momentum eqn.
  real(r8), intent(inout):: fv      (plon,plev)    ! nonlinear term - v momentum eqn.
  real(r8), intent(inout):: lnpssld (plon,plev)    ! RHS Ps term for SLD
  real(r8), intent(inout):: prhssld (plon,plev)    ! RHS Ps term for SLD
  real(r8), intent(inout):: tarrsld (plon,plev)    ! T  at arr. pt. (SLD)
  real(r8), intent(inout):: parrsld (plon,plev)    ! Ps at arr. pt. (SLD)
  real(r8), intent(inout):: t2      (plon,plev)    ! nonlinear term - temperature
  real(r8), intent(in)   :: div     (plon,plev)    ! divergence (time n)
  real(r8), intent(in)   :: tl      (plon,plev)    ! long derivative of T (n)
  real(r8), intent(in)   :: tm      (plon,plev)    ! lat  derivative of T (n)
  real(r8), intent(in)   :: ql      (plon,plev)    ! long derivative of q (n)
  real(r8), intent(in)   :: qm      (plon,plev)    ! lat  derivative of q (n-1)
  real(r8), intent(in)   :: dpsl    (plon)         ! long derivative of lnps (n)
  real(r8), intent(in)   :: dpsm    (plon)         ! lat  derivative of lnps (n)
  real(r8), intent(in)   :: phis    (plon)         ! surface geopotential
  real(r8), intent(in)   :: phisl   (plon)         ! long derivative of phis
  real(r8), intent(in)   :: phism   (plon)         ! lat  derivative of phis
  real(r8), intent(in)   :: omga    (plon,plev)    ! vertical velocity
  real(r8), intent(inout):: urhs    (plon,plev)    ! RHS of U  eqn valid for mid-point
  real(r8), intent(inout):: vrhs    (plon,plev)    ! RHS of V  eqn valid for mid-point
  real(r8), intent(inout):: trhs    (plon,plev)    ! RHS of T  eqn valid for mid-point
  real(r8), intent(inout):: prhs    (plon,plev)    ! RHS of Ps eqn valid for mid-point
  integer , intent(in)   :: nlon                   ! number of longitudes for this lat
  real(r8), intent(in)   :: cwava                  ! weight for global water vapor int.
  real(r8), intent(in)   :: flx_net(plon)          ! net flux from physics
!
!---------------------------Local workspace-----------------------------
!
  real(r8) :: dtime             ! timestep size
  real(r8) pmid  (plon,plev)    ! pressure at model levels (time n)
  real(r8) rpmid (plon,plev)    ! 1./pmid
  real(r8) pint  (plon,plevp)   ! pressure at model interfaces (n  )
  real(r8) pdel  (plon,plev)    ! pdel(k)   = pint  (k+1)-pint  (k)
  real(r8) ktoop (plon,plev)    ! (Kappa*T)*(omega/P)
  real(r8) logps (plon)         ! log(ps)
  real(r8) lvcour               ! local vertical courant number
  real(r8) dtdz                 ! dt/detam(k)

  real(r8) dpslat(plon,plev)    ! Pressure gradient term 
  real(r8) dpslon(plon,plev)    ! Pressure gradient term 
  real(r8) coslat               ! cosine(latitude)
  real(r8) rcoslat              ! 1./cosine(latitude)
  real(r8) coriol               ! Coriolis term
  real(r8) wind                 ! u**2 + v**2 (m/s)
  real(r8) utfac                ! asymmetric truncation factor for courant calculation
  real(r8) vtfac                ! asymmetric truncation factor for courant calculation
  real(r8) engy                 ! accumulator
  integer i,k                   ! indices
!
!----------------------------------------------------------------------
!
! Compute maximum wind speed this latitude (used in Courant number
! estimate)
!
  if (ptrm .lt. ptrn) then
     utfac = float(ptrm)/float(ptrn)
     vtfac = 1.
  else if (ptrn .lt. ptrm) then
     utfac = 1.
     vtfac = float(ptrn)/float(ptrm) 
  else if (ptrn .eq. ptrm) then
     utfac = 1.
     vtfac = 1.
  end if
  do k=1,plev
     vmax(k) = 0.
     vmaxt(k) = 0.
     do i=1,nlon
        wind = u3(i,k)**2 + v3(i,k)**2
        vmax(k) = max(wind,vmax(k))
!
! Change to Courant limiter for non-triangular truncations.
!
        wind = utfac*u3(i,k)**2 + vtfac*v3(i,k)**2
        vmaxt(k) = max(wind,vmaxt(k))
     end do
  end do
!
  coslat = cos(clat(lat))
  rcoslat = 1./coslat
  coriol  = 2.0*omega*sin(clat(lat))
!
! Set current time pressure arrays for model levels etc.
!
  call plevs0 (nlon    ,plon   ,plev    ,ps      ,pint    ,pmid    ,pdel)
  do k=1,plev
     do i=1,nlon
        rpmid(i,k) = 1./pmid(i,k)
     end do
  end do
!
! Accumulate statistics for diagnostic print
!
  call stats(lat     ,pint    ,pdel    ,ps      , &
             div     ,t3      ,q3(:,:,1),nlon   )
!     
! Compute log(surface pressure) for use by grmult and when adding
! tendency.
!     
  do i=1,nlon
     logps  (i) = log(ps  (i))
  end do

  if (adiabatic) t2     (:nlon,:)=0.
!     
! Compute (1/ps)etadot(dp/deta)
!     
  if (.not. is_first_step()) then
     call etadtn (lat   ,nlon    ,div     ,parrsld ,ed1     )
  end if
!     
! Compute integrals
!     
  call engy_te (cwava,w(lat),t3  ,u3  ,v3 ,phis    ,pdel, engy ,nlon)
  engy1lat(lat) = engy
!
! Include top/bottom flux integral to energy integral
!
  call flxint  (w(lat) ,flx_net ,engy ,nlon )
  engy1lat(lat) = engy1lat(lat) + engy*ztodt
!     
! Calculate non-linear terms in tendencies
!     
  call grmult_run(ztodt   ,rcoslat ,coslat  ,coriol  ,div     , &
              q3(1,1,1),t3     ,u3      ,v3               , &
                       ed1     ,ql      ,qm      ,tl      , &
              tm      ,phis    ,phisl   ,phism   ,dpsl    , &
              dpsm    ,omga    ,pmid    ,pdel    ,ps      , &
              logps   ,rpmid   ,etamid  ,fu      ,fv      , &
              t2      ,t3m1    ,u3m1    ,v3m1    ,etadot  , &
              dpslon  ,dpslat  ,urhs    ,vrhs    ,trhs    , &
              prhs    ,lnpssld ,prhssld ,tarrsld , &
              parrsld ,ktoop   ,nlon    )
  call outfld('ETADOT  ',etadot(1,1),plon,lat)
  call outfld('KTOOP   ',ktoop ,plon,lat)
!     
! Compute maximum vertical Courant number this latitude.
!     
  dtime = get_step_size()
  vcour(:) = 0.
  do k=2,plev
     dtdz = dtime/detam(k-1)
     do i=1,nlon
        lvcour = abs(etadot(i,k))*dtdz
        vcour(k) = max(lvcour,vcour(k))
     end do
  end do

  return
end subroutine linemsdyn

