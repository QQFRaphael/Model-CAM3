#include <misc.h>
#include <params.h>
subroutine esinti(epslon  ,latvap  ,latice  ,rh2o    ,cpair   ,tmelt   )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize es lookup tables
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Hack
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use wv_saturation, only: gestbl
   implicit none
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   real(r8), intent(in) :: epslon          ! Ratio of h2o to dry air molecular weights
   real(r8), intent(in) :: latvap          ! Latent heat of vaporization
   real(r8), intent(in) :: latice          ! Latent heat of fusion
   real(r8), intent(in) :: rh2o            ! Gas constant for water vapor
   real(r8), intent(in) :: cpair           ! Specific heat of dry air
   real(r8), intent(in) :: tmelt           ! Melting point of water (K)
!
!---------------------------Local workspace-----------------------------
!
   real(r8) tmn             ! Minimum temperature entry in table
   real(r8) tmx             ! Maximum temperature entry in table
   real(r8) trice           ! Trans range from es over h2o to es over ice
   logical ip           ! Ice phase (true or false)
!
!-----------------------------------------------------------------------
!
! Specify control parameters first
!
   tmn   = 173.16
   tmx   = 375.16
   trice =  20.00
   ip    = .true.
!
! Call gestbl to build saturation vapor pressure table.
!
   call gestbl(tmn     ,tmx     ,trice   ,ip      ,epslon  , &
               latvap  ,latice  ,rh2o    ,cpair   ,tmelt )
!
   return
end subroutine esinti

