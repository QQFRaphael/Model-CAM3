!------------------------------------------------------------------------
! File: setiopupdate.F 
! Author: John Truesdale (jet@ucar.edu) 
! $Id: setiopupdate.F90,v 1.1.6.1 2004/05/18 20:51:32 jet Exp $
!
!------------------------------------------------------------------------
#include <params.h>
#include <max.h>

subroutine setiopupdate

!-----------------------------------------------------------------------
!   
! Open and read netCDF file to extract time information
!
!---------------------------Code history--------------------------------
!
! Written by John Truesdale    August, 1996
! 
!-----------------------------------------------------------------------
  use pmgrid
  use prognostics
  use time_manager, only: timemgr_init,nestep, get_curr_date, get_curr_calday,start_ymd,start_tod,get_nstep
  use scamMod, only :iopfile,use_userdata,doiopupdate,ioptimeidx
  implicit none
#if ( defined RS6000 )
  implicit automatic (a-z)
#endif

!------------------------------Includes-----------------------------------
#include <netcdf.inc>

!------------------------------Locals-----------------------------------

   integer NCID,i
   integer tsec_varID, time_dimID
   integer tsec(MAX_TIME_DIM), ntime 
   integer bdate, bdate_varID
   integer STATUS
   integer next_date, next_sec, last_date, last_sec 
   integer :: ncsec,ncdate                      ! current time of day,date
   integer :: yr, mon, day                      ! year, month, and day component
   save tsec, ntime, bdate
   save last_date, last_sec 
!------------------------------------------------------------------------------

   if ( get_nstep() .eq. 0 ) then
!     
!     Open  IOP dataset
!     
      STATUS = NF_OPEN( iopfile, NF_NOWRITE, NCID )
!     
!     Read time (tsec) variable 
!     
      STATUS = NF_INQ_VARID( NCID, 'tsec', tsec_varID )
      if ( STATUS .NE. NF_NOERR ) write(6,*)'ERROR - setiopupdate.F:', &
         'Cant get variable ID for tsec'

      STATUS = NF_INQ_VARID( NCID, 'bdate', bdate_varID )
      if ( STATUS .NE. NF_NOERR ) then
         STATUS = NF_INQ_VARID( NCID, 'basedate', bdate_varID )
         if ( STATUS .NE. NF_NOERR )         &
            write(6,*)'ERROR - setiopupdate.F:Cant get variable ID for bdate'
      endif

      STATUS = NF_INQ_DIMID( NCID, 'time', time_dimID )
      if ( STATUS .NE. NF_NOERR )  then
         STATUS = NF_INQ_DIMID( NCID, 'tsec', time_dimID )
         if ( STATUS .NE. NF_NOERR )  then
            write( 6,* )'ERROR - setiopupdate.F:Could not find variable dim ID for time'
            STATUS = NF_CLOSE ( NCID )
            return
         end if
      end if

      if ( STATUS .NE. NF_NOERR )  &
         write(6,*)'ERROR - setiopupdate.F:Cant get variable dim ID for time'

      STATUS = NF_INQ_DIMLEN( NCID, time_dimID, ntime )
      if ( STATUS .NE. NF_NOERR )then
         write(6,*)'ERROR - setiopupdate.F:Cant get time dimlen'
      endif

      STATUS = NF_GET_VAR_INT( NCID, tsec_varID, tsec )
      if ( STATUS .NE. NF_NOERR )then
         write(6,*)'ERROR - setiopupdate.F:Cant get variable tsec'
      endif
      STATUS = NF_GET_VAR_INT( NCID, bdate_varID, bdate )
      if ( STATUS .NE. NF_NOERR )then
         write(6,*)'ERROR - setiopupdate.F:Cant get variable bdate'
      endif
!     Close the netCDF file
      STATUS = NF_CLOSE( NCID )

!     
!     determine the last date in the dataset
!     
      call calcdate( bdate, tsec(ntime), last_date, last_sec )
!     
!     set the iop dataset index
!    
      do i=1,ntime           ! set the first ioptimeidx
         call calcdate( bdate, tsec(i), next_date, next_sec )
         if ( start_ymd .gt. next_date .or. (start_ymd .eq. next_date &
            .and. start_tod .ge. next_sec)) then
            iopTimeIdx = i
         endif
      enddo


      call get_curr_date(yr,mon,day,ncsec)
      ncdate=yr*10000 + mon*100 + day
      doiopupdate = .false.

!------------------------------------------------------------------------------
!     Check if iop data needs to be updated and set doiopupdate accordingly
!------------------------------------------------------------------------------
   else                      ! endstep > 1

      call calcdate( bdate, tsec(iopTimeIdx+1), next_date, next_sec )

      call get_curr_date(yr, mon, day, ncsec)
      ncdate = yr*10000 + mon*100 + day

      if ( ncdate .gt. next_date .or. (ncdate .eq. next_date &
         .and. ncsec .ge. next_sec)) then
         iopTimeIdx = iopTimeIdx + 1
         doiopupdate = .true.
#if DEBUG > 2
         print *, 'nstep = ',get_nstep()
         print *, 'ncdate=',ncdate,' ncsec=',ncsec
         print *, 'next_date=',next_date,' next_sec=',next_sec
         write(*,*)'******* do iop update'
#endif 
      else
         doiopupdate = .false.
      end if
   endif                     ! if (endstep .eq. 0 )
!
!     make sure we're
!     not going past end of iop data
!
   if ( ncdate .gt. last_date .or. (ncdate .eq. last_date &
      .and. ncsec .gt. last_sec))  then
      if ( .not. use_userdata ) then
         write(6,*)'ERROR - setiopupdate.c:Reached the end of the time varient dataset'
         stop
      else
         doiopupdate = .false.              
      end if
   endif

#if DEBUG > 1
   write(*,*)'iop time index = ' , ioptimeidx
#endif

   return

end subroutine setiopupdate

