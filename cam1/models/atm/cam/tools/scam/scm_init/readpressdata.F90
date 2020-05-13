
#include <params.h>
#include <max.h>
!------------------------------------------------------------------------
! File: readpressdata.F 
! Author: John Truesdale (jet@ucar.edu) 
! $Id: readpressdata.F90,v 1.1.6.2 2004/08/19 15:05:42 jmccaa Exp $
!
!------------------------------------------------------------------------
subroutine readpressdata( error_code )
!-----------------------------------------------------------------------
!   
! Open and read netCDF file containing hybrid pressure coordinates
!
!---------------------------Code history--------------------------------
!
! Written by John Truesdale, 11-25-97
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
   use pmgrid
   use getnetcdfdata
   use scamMod, only :pressfile
   implicit none
#if ( defined RS6000 )
   implicit automatic ( a-z )
#endif
#include <runtype.h>      
!-----------------------------------------------------------------------
#include <comhyb.h>
!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------

!------------------------------Inputs-----------------------------------

   integer error_code     ! out:  netcdf errors

!-----------------------------Externals---------------------------------


!------------------------------Locals-----------------------------------
!
!   real(r8) hyai( plev+1 )
!   real(r8) hybi( plev+1 )
!   real(r8) hyam( plev )
!   real(r8) hybm( plev )

   integer     NCID, STATUS, lev_dimID
   integer     nlev, i


!-----------------------------------------------------------------------

#if ( defined sun )
   external myhandler
   integer iexcept, ieee_handler, myhandler
#endif
!
!-----------------------------------------------------------------------
!
!     Trap ieee exceptions on SUN for debugging purposes
!
#if ( defined sun )
   iexcept = ieee_handler( 'set', 'common', myhandler )
   if ( iexcept .ne. 0 ) write( 6,* )'ieee trapping not supported here'
#endif
   error_code = PRES

!
!     Open pressure dataset
!

   STATUS = NF_OPEN( pressfile, NF_NOWRITE, NCID )
   if( STATUS .NE. NF_NOERR )  then
      write( 6,* )'ERROR - readpressdata.F:Cant open pressure data netcdf file', &
         pressfile
      STATUS = NF_CLOSE( NCID )
      return
   endif

! 
! ---------------------------------------------------------------------
! get pressure level variables
! 

   STATUS =  NF_INQ_DIMID( NCID, 'lev', lev_dimID )
   if ( STATUS .NE.  NF_NOERR )  then
      write( 6,* )'ERROR - readpressdata.F:Cant get variable dimension for lev'
      STATUS = NF_CLOSE( NCID )
      return
   endif
!
!     check that number of levels is correct
!
   STATUS = NF_INQ_DIMLEN( NCID, lev_dimID, nlev )
   if ( nlev .NE. plev ) then
      write(*,*) 'ERROR: readpressdata.F:'
      write(*,*) pressfile, 'has wrong number of levels.'
      write(*,*) 'Model is compiled for ', plev, ' levels, ' 
      write(*,*) 'but you  have ',nlev,' levels!'
      STATUS = NF_CLOSE( NCID )
      call exit(-1)
      return
   endif
!
! ---------------------------------------------------------------
!
   call getncdata( NCID, 1, 1, 1, 'hyam', hyam, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readpressdata.F:Cant get variable hyam'
      STATUS = NF_CLOSE( NCID )
      return
   endif

   call getncdata( NCID, 1, 1, 1, 'hyai', hyai, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readpressdata.F:Cant get variable hyai'
      STATUS = NF_CLOSE( NCID )
      return
   endif

   call getncdata ( NCID, 1, 1, 1, 'hybm', hybm, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readpressdata.F:Cant get variable hybm'
      STATUS = NF_CLOSE( NCID )
      return
   endif

   call getncdata( NCID, 1, 1, 1, 'hybi', hybi, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readpressdata.F:Cant get variable hybi'
      STATUS = NF_CLOSE( NCID )
      return
   endif

!     Close the netCDF file
   STATUS = NF_CLOSE( NCID )
   error_code = 0

   return
end subroutine readpressdata
