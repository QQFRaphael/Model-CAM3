#include <params.h>
subroutine get_levels(nlev, press,ilev,ipress,number_only)

   use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
   use pmgrid
   use prognostics
!   use comsrf, only:
!-----------------------------------------------------------------------
#include <max.h>
#include <comhyb.h>

   real(r8) press(MAX_LEVELS)
   real(r8) ipress(MAX_LEVELS)
   integer nlev
   integer ilev
   integer number_only

!-----------------------------------------------------------------------      
!       
!  compute  model pressure levels 
!
   nlev = plev
   ilev=plevp
   if(number_only.eq.1)return
   do i = 1, nlev
      press(i) = 1000.0 * hyam(i) + ps(1,1,1) * hybm(i) / 100.0
   end do
   do i = 1, ilev
      ipress(i) = 1000.0 * hyai(i) + ps(1,1,1) * hybi(i) / 100.0
   end do

   return
end subroutine get_levels


