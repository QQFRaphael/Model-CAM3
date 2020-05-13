#include <misc.h>
#include <params.h>

module iop
#if ( defined BFB_CAM_SCAM_IOP ) || ( defined SCAM )
!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: iop
! 
! !DESCRIPTION: 
! iop specific routines
!
! !USES:
!
  use shr_kind_mod, only: r8 => shr_kind_r8
  use constituents, only: pcnst, cnst_name,ppcnst
  use pmgrid,       only: plon,plev,plat,beglat,endlat
  use string_utils, only: to_lower
!
! !PUBLIC TYPES:
  implicit none

  private

   real(r8), allocatable :: betasav(:)
   real(r8), allocatable :: fixmassav(:)
   real(r8), allocatable :: alphasav(:,:)
   real(r8), allocatable :: clat_plon(:)          ! d(ps)/dt
   real(r8), allocatable :: alpha_plon(:,:)   
   real(r8), allocatable :: fixmas_plon(:)         
   real(r8), allocatable :: beta_plon(:)          
   real(r8), allocatable :: dqfx3sav(:,:,:)       
   real(r8), allocatable :: dqfx3savm1(:,:,:,:)       
   real(r8), allocatable :: divq3dsav(:,:,:,:)
   real(r8), allocatable :: divt3dsav(:,:,:)       
   real(r8), allocatable :: t3sav(:,:,:)       
   real(r8), allocatable :: u3sav(:,:,:)       
   real(r8), allocatable :: v3sav(:,:,:)       
   real(r8), allocatable :: t2sav(:,:,:)         ! temp tendency
   real(r8), allocatable :: q3sav(:,:,:,:)
   real(r8), allocatable :: pssav(:,:)
   real(r8), allocatable :: tssav(:,:)
   character(len=8) alphanam(pcnst)       ! alpha fixer
   character(len=8) dqfxnam(pcnst)       ! dq fixer

!
! !PUBLIC MEMBER FUNCTIONS:
  public :: init_iop_fields
!
! !PUBLIC DATA:
  public betasav, fixmassav, alphasav, clat_plon, alpha_plon, fixmas_plon, beta_plon,   &
         dqfx3sav, dqfx3savm1, divq3dsav, divt3dsav, t3sav, u3sav, v3sav, t2sav, q3sav, &
         pssav, tssav, alphanam, dqfxnam
!
! !REVISION HISTORY:
! Created by John Truesdale
!
!EOP
!
! !PRIVATE MEMBER FUNCTIONS:
!----------------------------------------------------------------------- 

contains
   subroutine init_iop_fields(ps, t3, u3, v3, q3, nocopy)
!------------------------------------------------------------------------------
! Coupler for converting dynamics output variables into physics input variables
! also writes dynamics variables (on physics grid) to history file
!------------------------------------------------------------------------------
   use pmgrid,       only: beglat,endlat,plon,numlats,plev
   use constituents, only: ppcnst

   implicit none

!------------------------------Arguments--------------------------------
    real(r8), intent(in) :: ps  (plon, beglat:endlat)               ! surface pressure
    real(r8), intent(in) :: t3  (plon, plev, beglat:endlat)         ! temperature
    real(r8), intent(in) :: u3  (plon, plev, beglat:endlat)         ! u-wind component
    real(r8), intent(in) :: v3  (plon, plev, beglat:endlat)         ! v-wind component
    real(r8), intent(in) :: q3  (plon, plev, ppcnst, beglat:endlat) ! constituents
    logical, intent(in), optional :: nocopy                         ! Flag to not copy variables

!-----------------------------------------------------------------------
        
   if(.not.allocated(betasav))     allocate (betasav(beglat:endlat))
   if(.not.allocated(fixmassav))   allocate (fixmassav(beglat:endlat))
   if(.not.allocated(alphasav))    allocate (alphasav(pcnst,beglat:endlat))
   if(.not.allocated(clat_plon))   allocate (clat_plon(plon))          ! d(ps)/dt
   if(.not.allocated(alpha_plon))  allocate (alpha_plon(plon,pcnst))
   if(.not.allocated(fixmas_plon)) allocate (fixmas_plon(plon))
   if(.not.allocated(beta_plon))   allocate (beta_plon(plon))
   if(.not.allocated(dqfx3sav))    allocate (dqfx3sav(plon,plev,pcnst))
   if(.not.allocated(dqfx3savm1))  allocate (dqfx3savm1(plon,plev,pcnst,beglat:endlat))
   if(.not.allocated(divq3dsav))   allocate (divq3dsav(plon,plev,ppcnst,beglat:endlat))
   if(.not.allocated(divt3dsav))   allocate (divt3dsav(plon,plev,beglat:endlat))
   if(.not.allocated(t3sav))       allocate (t3sav(plon,plev,beglat:endlat))
   if(.not.allocated(u3sav))       allocate (u3sav(plon,plev,beglat:endlat))
   if(.not.allocated(v3sav))       allocate (v3sav(plon,plev,beglat:endlat))
   if(.not.allocated(t2sav))       allocate (t2sav(plon,plev,beglat:endlat))  ! temp tendency
   if(.not.allocated(q3sav))       allocate (q3sav(plon,plev,ppcnst,beglat:endlat))
   if(.not.allocated(pssav))       allocate (pssav(plon,beglat:endlat))
   if(.not.allocated(tssav))       allocate (tssav(plon,beglat:endlat))
   if ( .not. present(nocopy) ) then
      t3sav(:,:,:)   = t3(:,:,:)
      u3sav(:,:,:)   = u3(:,:,:)
      v3sav(:,:,:)   = v3(:,:,:)
      q3sav(:,:,:,:) = q3(:,:,:,:)
      pssav(:,:)     = ps(:,:)
   end if
  end subroutine init_iop_fields

#endif

end module iop

