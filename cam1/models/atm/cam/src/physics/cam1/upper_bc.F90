#include <misc.h>
#include <params.h>

module upper_bc

!---------------------------------------------------------------------------------
! Module to compute the upper boundary condition for temperature (dry static energy)
! and trace gases.   The standard CAM version does nothing.
!
! original code by Stacy Walters
! adapted by B. A. Boville
!---------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, pverp
  use constituents, only: ppcnst

  implicit none
  private
  save
!
! Public interfaces
!
  public :: ubc_defaultopts    ! set default values of namelist variables
  public :: ubc_setopts        ! get namelist input
  public :: ubc_init           ! global initialization
  public :: ubc_timestep_init  ! time step initialization
  public :: ubc_get_vals       ! get ubc values for this step

!================================================================================================
contains
!================================================================================================

subroutine ubc_defaultopts(tgcm_ubc_file_out, snoe_ubc_file_out)
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
!-----------------------------------------------------------------------

   character(len=*), intent(out), optional :: tgcm_ubc_file_out
   character(len=*), intent(out), optional :: snoe_ubc_file_out
!-----------------------------------------------------------------------

end subroutine ubc_defaultopts

!================================================================================================

subroutine ubc_setopts(tgcm_ubc_file_in, snoe_ubc_file_in)
!----------------------------------------------------------------------- 
! Purpose: Set runtime options
!-----------------------------------------------------------------------

   character(len=*), intent(in), optional :: tgcm_ubc_file_in
   character(len=*), intent(in), optional :: snoe_ubc_file_in
!-----------------------------------------------------------------------

end subroutine ubc_setopts

!===============================================================================

  subroutine ubc_init
!-----------------------------------------------------------------------
! Initialization of time independent fields for the upper boundary condition
! Calls initialization routine for MSIS, TGCM and SNOE
!-----------------------------------------------------------------------

  end subroutine ubc_init

!===============================================================================

  subroutine ubc_timestep_init
!-----------------------------------------------------------------------
! timestep dependent setting
!-----------------------------------------------------------------------

  end subroutine ubc_timestep_init

!===============================================================================

  subroutine ubc_get_vals (lchnk, ncol, ntop_molec, pint, zi, msis_temp, ubc_mmr)
!-----------------------------------------------------------------------
! interface routine for vertical diffusion and pbl scheme
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    integer,  intent(in)  :: lchnk                 ! chunk identifier
    integer,  intent(in)  :: ncol                  ! number of atmospheric columns
    integer,  intent(in)  :: ntop_molec            ! top of molecular diffusion region (=1)
    real(r8), intent(in)  :: pint(pcols,pverp)     ! interface pressures
    real(r8), intent(in)  :: zi(pcols,pverp)       ! interface geoptl height above sfc

    real(r8), intent(out) :: ubc_mmr(pcols,ppcnst) ! upper bndy mixing ratios (kg/kg)
    real(r8), intent(out) :: msis_temp(pcols)      ! upper bndy temperature (K)

  end subroutine ubc_get_vals

end module upper_bc
