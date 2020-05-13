#include <misc.h>
#include <params.h>

module ghg_defaults

!------------------------------------------------------------------------------------------------
! Purpose:
! Provide default distributions of CH4, N2O, CFC11 and CFC12 to the radiation routines.
! **NOTE** CO2 is assumed by the radiation to a be constant value.  This value is
!          currently supplied directly by the chem_surfvals module.
!
! Revision history:
! 2004-08-29  B. Eaton        Create CAM interface to trcmix.
!------------------------------------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use ppgrid,         only: pcols, pver, begchunk, endchunk
use physics_types,  only: physics_state
use physconst,      only: mwdry, mwch4, mwn2o, mwf11, mwf12
use chem_surfvals,  only: chem_surfvals_get
use abortutils,     only: endrun
use error_messages, only: handle_err

implicit none
private
save

! Public interfaces
public ::&
   ghg_defaults_init,     &! initialization
   ghg_defaults_get_cnst   ! return pointer to constituent concentrations.

! Private variables

real(r8), parameter :: rmwn2o = mwn2o/mwdry      ! ratio of molecular weight n2o   to dry air
real(r8), parameter :: rmwch4 = mwch4/mwdry      ! ratio of molecular weight ch4   to dry air
real(r8), parameter :: rmwf11 = mwf11/mwdry      ! ratio of molecular weight cfc11 to dry air
real(r8), parameter :: rmwf12 = mwf12/mwdry      ! ratio of molecular weight cfc12 to dry air

integer, parameter :: ncnst = 4                        ! number of constituents
character(len=8), dimension(ncnst), parameter :: &
   cnst_names = (/'N2O  ', 'CH4  ', 'CFC11', 'CFC12'/) ! constituent names

! storage for constituents diagnosed by trcmix
real(r8), allocatable, target, dimension(:,:,:,:) :: diag_q  ! constituent mass mixing ratios

!================================================================================================
contains
!================================================================================================

subroutine ghg_defaults_init()
!-------------------------------------------------------------------------------
!
! Purpose:
! Allocate memory for constituent distributions.
!
!-------------------------------------------------------------------------------
   integer :: istat
!-------------------------------------------------------------------------------

   allocate(diag_q(pcols,pver,ncnst,begchunk:endchunk), stat=istat)
   call handle_err(istat, 'ghg_defaults_init: ERROR allocating diag_q')

end subroutine ghg_defaults_init

!================================================================================================

function ghg_defaults_get_idx(name)
!-------------------------------------------------------------------------------
! 
! Purpose:
! Return index of constituent diagnosed by trcmix.
! 
!-------------------------------------------------------------------------------
   character(len=*), intent(in) :: name                  ! constituent name
   integer                      :: ghg_defaults_get_idx  ! return value

! Local variable
   integer :: m
!-------------------------------------------------------------------------------

   ghg_defaults_get_idx = -1
   do m = 1, ncnst
      if (name == cnst_names(m)) then
         ghg_defaults_get_idx = m
         return
      end if
   end do
end function ghg_defaults_get_idx

!================================================================================================

subroutine ghg_defaults_get_cnst(name, state, q)
!-------------------------------------------------------------------------------
! 
! Purpose: 
! Return pointer to constituent concentrations.
!
!-------------------------------------------------------------------------------

   character(len=*),                  intent(in) :: name  ! constituent name
   type(physics_state),               intent(in) :: state
   real(r8), pointer, dimension(:,:)             :: q     ! constituent mass mixing ratio

   ! local variables
   integer :: idx
   integer :: lchnk            ! chunk identifier
!-------------------------------------------------------------------------------

   ! Make sure the requested constituent can be provided.
   idx = ghg_defaults_get_idx(name)
   if ( idx < 1 ) then
      call endrun ('ghg_defaults_get_cnst: '//name//' not implemented in ghg_defaults module')
   end if

   lchnk = state%lchnk
   call trcmix(name, state%ncol, state%lat, state%pmid, diag_q(:,:,idx,lchnk))
   q => diag_q(:,:,idx,lchnk)

end subroutine ghg_defaults_get_cnst

!================================================================================================

subroutine trcmix(name, ncol, clat, pmid, q)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Specify zonal mean mass mixing ratios of CH4, N2O, CFC11 and
! CFC12
! 
! Method: 
! Distributions assume constant mixing ratio in the troposphere
! and a decrease of mixing ratio in the stratosphere. Tropopause
! defined by ptrop. The scale height of the particular trace gas
! depends on latitude. This assumption produces a more realistic
! stratospheric distribution of the various trace gases.
! 
! Author: J. Kiehl
! 
!-----------------------------------------------------------------------

   ! Arguments
   character(len=*), intent(in)  :: name              ! constituent name
   integer,          intent(in)  :: ncol              ! number of columns
   real(r8),         intent(in)  :: clat(pcols)       ! latitude in radians for columns
   real(r8),         intent(in)  :: pmid(pcols,pver)  ! model pressures
   real(r8),         intent(out) :: q(pcols,pver)     ! constituent mass mixing ratio

   integer i                ! longitude loop index
   integer k                ! level index

   real(r8) coslat(pcols)   ! cosine of latitude
   real(r8) dlat            ! latitude in degrees
   real(r8) ptrop           ! pressure level of tropopause
   real(r8) pratio          ! pressure divided by ptrop
   real(r8) trop_mmr        ! tropospheric mass mixing ratio
   real(r8) scale           ! pressure scale height
!-----------------------------------------------------------------------

   do i = 1, ncol
      coslat(i) = cos(clat(i))
   end do

   if (name == 'CH4') then

      ! set tropospheric mass mixing ratios
      trop_mmr = rmwch4 * chem_surfvals_get('CH4VMR')

      do k = 1,pver
         do i = 1,ncol
            ! set stratospheric scale height factor for gases
            dlat = abs(57.2958 * clat(i))
            if(dlat.le.45.0) then
               scale = 0.2353
            else
               scale = 0.2353 + 0.0225489 * (dlat - 45)
            end if

            ! pressure of tropopause
            ptrop = 250.0e2 - 150.0e2*coslat(i)**2.0

            ! determine output mass mixing ratios
            if (pmid(i,k) >= ptrop) then
               q(i,k) = trop_mmr
            else
               pratio = pmid(i,k)/ptrop
               q(i,k) = trop_mmr * (pratio)**scale
            end if
         end do
      end do

   else if (name == 'N2O') then

      ! set tropospheric mass mixing ratios
      trop_mmr = rmwn2o * chem_surfvals_get('N2OVMR')

      do k = 1,pver
         do i = 1,ncol
            ! set stratospheric scale height factor for gases
            dlat = abs(57.2958 * clat(i))
            if(dlat.le.45.0) then
               scale = 0.3478 + 0.00116 * dlat
            else
               scale = 0.4000 + 0.013333 * (dlat - 45)
            end if

            ! pressure of tropopause
            ptrop = 250.0e2 - 150.0e2*coslat(i)**2.0

            ! determine output mass mixing ratios
            if (pmid(i,k) >= ptrop) then
               q(i,k) = trop_mmr
            else
               pratio = pmid(i,k)/ptrop
               q(i,k) = trop_mmr * (pratio)**scale
            end if
         end do
      end do

   else if (name == 'CFC11') then

      ! set tropospheric mass mixing ratios
      trop_mmr = rmwf11 * chem_surfvals_get('F11VMR')

      do k = 1,pver
         do i = 1,ncol
            ! set stratospheric scale height factor for gases
            dlat = abs(57.2958 * clat(i))
            if(dlat.le.45.0) then
               scale = 0.7273 + 0.00606 * dlat
            else
               scale = 1.00 + 0.013333 * (dlat - 45)
            end if

            ! pressure of tropopause
            ptrop = 250.0e2 - 150.0e2*coslat(i)**2.0

            ! determine output mass mixing ratios
            if (pmid(i,k) >= ptrop) then
               q(i,k) = trop_mmr
            else
               pratio = pmid(i,k)/ptrop
               q(i,k) = trop_mmr * (pratio)**scale
            end if
         end do
      end do

   else if (name == 'CFC12') then

      ! set tropospheric mass mixing ratios
      trop_mmr = rmwf12 * chem_surfvals_get('F12VMR')

      do k = 1,pver
         do i = 1,ncol
            ! set stratospheric scale height factor for gases
            dlat = abs(57.2958 * clat(i))
            if(dlat.le.45.0) then
               scale = 0.4000 + 0.00222 * dlat
            else
               scale = 0.50 + 0.024444 * (dlat - 45)
            end if

            ! pressure of tropopause
            ptrop = 250.0e2 - 150.0e2*coslat(i)**2.0

            ! determine output mass mixing ratios
            if (pmid(i,k) >= ptrop) then
               q(i,k) = trop_mmr
            else
               pratio = pmid(i,k)/ptrop
               q(i,k) = trop_mmr * (pratio)**scale
            end if
         end do
      end do

   end if

end subroutine trcmix

end module ghg_defaults
