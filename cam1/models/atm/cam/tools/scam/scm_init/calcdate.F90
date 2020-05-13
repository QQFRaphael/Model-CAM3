!------------------------------------------------------------------------
! File: calcdate.F 
! Author: John Truesdale (jet@ucar.edu) 
! $Id: calcdate.F90,v 1.1.6.1 2004/05/18 20:51:26 jet Exp $
!
!------------------------------------------------------------------------
#include <params.h>
subroutine calcdate(inDate, inSecs,  outDate, outSecs)
!-----------------------------------------------------------------------
!  calcdate           Calculate Date from base date plus seconds
!
! INPUTS:
!
!	inDate	       Base date as YYMMDD.
!       inSecs         number of seconds the model has run
!
! OUTPUTS:
!       outDate        Current date as YYMMDD
!       outSecs        number of seconds into current date
!
!
!-----------------------------------------------------------------------
! Computational notes: 
!
! 86400 is the number of seconds in 1 day.
!
! Dividing an integer by 10**n has the effect of right-shifting the           
! decimal digits n positions (ex: 861231/100 = 008612).
!
! mod(integer,10**n) has the effect of extracting the rightmost     
! n decimal digits of the integer (ex: mod(861231,10000) = 1231).
!
use abortutils,   only: endrun
   implicit none
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer     inDate
   integer     inSecs       

!         
! Output arguments       
!                        
   integer     outSecs
   integer     outDate
!
!---------------------------Local workspace-----------------------------
!
   integer     YY
   integer     MM
   integer     DD
   integer     i
   integer     bmnth
   integer     bday
   integer     jday
   integer     jsec
   integer     jdcon(12)
   integer     ndm(12)
   integer     secs_per_year

   data ndm/31,28,31,30,31,30,31,31,30,31,30,31/
   data jdcon/0,31,59,90,120,151,181,212,243,273,304,334/
!
!-----------------------------------------------------------------------
!
! Check validity of input data
!
   bmnth = mod(inDate,10000)/100
   bday =  mod(inDate,100)
   if (bmnth.lt.1 .or. bmnth.gt.12) then
      write(6,*)' CALCDATE: Invalid base month input:',bmnth
      call endrun
   end if
   if (bday.lt.1 .or. bday.gt.ndm(bmnth)) then
      write(6,*)' CALCDATE: Invalid base day of base date input:',bday
      call endrun
   end if
!
!
!
   jday = jdcon(bmnth) + bday
   jsec  = (jday-1) * 86400 + insecs

   secs_per_year = 86400 * 365

   YY  = inDate/10000 + jsec/secs_per_year

   jsec = mod(jsec, secs_per_year)


   do i=1, 12
      if ((jsec - jdcon(i)*86400) .GE. 0) then
         MM = i
      end if
   end do

   jsec = jsec - jdcon(MM)*86400

   DD = jsec/86400 +1

   outSecs = mod(jsec,86400)

   outDate = YY*10000 + MM*100 + DD

!      write( *,* )'date =' , outDate
!
   return
end subroutine calcdate
