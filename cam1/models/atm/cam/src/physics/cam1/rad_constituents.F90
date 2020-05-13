#include <misc.h>
#include <params.h>

module rad_constituents

!------------------------------------------------------------------------------------------------
! Purpose:
!
! Provide constituent distributions to the radiation routines.
! 
! By default the subroutine that returns constituent mixing ratios returns
! the distribution to be used for the interactive calculation.  That method
! also provides an optional argument for indicating that the requested
! distribution is passive and only used in the diagnostic radiative forcing
! calculation.
! 
! The logic to control which constituent distribution is interactive and
! which is passive is contained in this module.  By default, if a prognostic
! version of a constituent is found (by looking in the constituent array),
! then it is used for the interactive calculation, and a prescribed
! distribution will be used for the diagnostic calculation.
!
! Revision history:
! 2004-08-28  B. Eaton      Original version
!------------------------------------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use physics_types,  only: physics_state
use ozone_data,     only: ozone_data_get_cnst
use ghg_defaults,   only: ghg_defaults_get_cnst
use constituents,   only: cnst_get_ind
use abortutils,     only: endrun

implicit none
private
save

! Public interfaces

public ::&
   rad_constituents_defaultopts, &! set default namelist values
   rad_constituents_setopts,     &! set runtime namelist values
   rad_constituents_use_data_o3, &! returns value of use_data_o3 namelist variable
   rad_constituents_get           ! return pointer to constituent concentration in mmr.

! Private module data

! Namelist variables
logical :: use_data_o3 = .false.  ! true => use ozone dataset for interactive radiation calc even
                                  !         if ozone is available as a constituent

!================================================================================================
contains
!================================================================================================

subroutine rad_constituents_defaultopts( &
   use_data_o3_out)
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
!-----------------------------------------------------------------------

   logical,          intent(out), optional :: use_data_o3_out
!-----------------------------------------------------------------------

   if ( present(use_data_o3_out) ) then
      use_data_o3_out = use_data_o3
   endif
end subroutine rad_constituents_defaultopts

!================================================================================================

subroutine rad_constituents_setopts( &
   use_data_o3_in)
!----------------------------------------------------------------------- 
! Purpose: Set runtime options
!-----------------------------------------------------------------------

   logical,          intent(in), optional :: use_data_o3_in
!-----------------------------------------------------------------------

   if ( present(use_data_o3_in) ) then
      use_data_o3 = use_data_o3_in
   endif
end subroutine rad_constituents_setopts

!================================================================================================

logical function rad_constituents_use_data_o3()
   rad_constituents_use_data_o3 = use_data_o3
end function rad_constituents_use_data_o3

!================================================================================================

subroutine rad_constituents_get(name, state, q, passive)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Return pointer to constituent concentrations.
!
!-----------------------------------------------------------------------

   character(len=*),                  intent(in) :: name    ! constituent name
   type(physics_state), target,       intent(in) :: state   ! state contains prognostic constituents
   real(r8), pointer, dimension(:,:)             :: q       ! constituent mass mixing ratio
   logical, optional,                 intent(in) :: passive ! false for interactive radiation calc (default)

   ! local variables
   integer :: idx
   logical :: passive_rad_calc
   logical :: prognostic_constituent_present
   logical :: use_prognostic_constituent
!-----------------------------------------------------------------------

   if ( present(passive) ) then
      passive_rad_calc = passive
   else
      passive_rad_calc = .false.
   end if

   ! Check for requested constituent in constituent array.
   call cnst_get_ind(name, idx, abort=.false.)
   if ( idx == -1 ) then
      prognostic_constituent_present = .false.
   else
      prognostic_constituent_present = .true.
   end if

   ! Determine whether to use prognostic constituent or prescribed data
   !
   ! Default: If doing interactive radiative calculation, and a prognostic 
   ! constituent is present, then use it.
   if ( .not. passive_rad_calc  .and. prognostic_constituent_present ) then
      use_prognostic_constituent = .true.
   else   
      use_prognostic_constituent = .false.
   end if

   ! Override use of prognostic ozone for interactive calc
   if (name .eq. 'O3'  .and. use_data_o3) then
      use_prognostic_constituent = .false.
   end if

   ! Return prognostic constituent
   if (use_prognostic_constituent) then
      q => state%q(:,:,idx)
      return
   end if

   ! Return diagnosed or prescribed constituent
   select case (name)

   case ('O3')
      call ozone_data_get_cnst(state, q)

   case ('N2O', 'CH4', 'CFC11', 'CFC12')
      call ghg_defaults_get_cnst(name, state, q)
         
   case default
      call endrun('rad_constituents_get: unknown name:'//name)

   end select

end subroutine rad_constituents_get

!================================================================================================

end module rad_constituents
