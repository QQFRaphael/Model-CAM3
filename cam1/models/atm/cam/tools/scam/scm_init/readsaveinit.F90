!------------------------------------------------------------------------
! File: readsaveinit.F 
! Author: John Truesdale (jet@ucar.edu) 
! $Id: readsaveinit.F90,v 1.1.6.1 2004/05/18 20:51:29 jet Exp $
!------------------------------------------------------------------------
#include <params.h>
#include <max.h>
#include <runtype.h>      

subroutine readsaveinit(error_code)

!-----------------------------------------------------------------------
!   
! Open and read netCDF file containing saved initial conditions
!
!
!---------------------------Code history--------------------------------
!
! Written by John Truesdale
! 10-96
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
   use pmgrid
   use prognostics
   use buffer
   use prognostics, only: phis
   use getnetcdfdata
   use scamMod, only :userfile,sicfile,uobs,vobs,tobs,qobs
   implicit none
#if ( defined RS6000 )
   implicit automatic (a-z)
#endif


!-----------------------------------------------------------------------
#include <comfrc.h>
!-----------------------------------------------------------------------
#include <netcdf.inc>
!------------------------------Inputs-----------------------------------

   integer error_code     ! returns netcdf errors


!-----------------------------Externals---------------------------------

!------------------------------Locals-----------------------------------
!

   integer NCID, lev_dimid, nlev
   integer STATUS
   integer k

   real(r8) dplevs( MAX_DATASET_LEVS )
   real(r8) surfdat
   real(r8) tmp( MAX_DATASET_LEVS )
   real(r8) tmpdata

   logical have_surfdat


!-----------------------------------------------------------------------

#if (defined sun)
   external myhandler
   integer iexcept, ieee_handler, myhandler
#endif
!
!-----------------------------------------------------------------------
!
!     Trap ieee exceptions on SUN for debugging purposes
!
#if (defined sun)
   iexcept = ieee_handler( 'set', 'common', myhandler )
   if ( iexcept .ne. 0 ) write(6,*)'ieee trapping not supported here'
#endif

   error_code = SIC
!
!     Open initial dataset
!
   STATUS = NF_OPEN( sicfile, NF_NOWRITE, NCID )
   if( STATUS .NE. NF_NOERR ) then
      write(6,*)'ERROR - readsaveinit:', &
         'Cant open saved initial conditions file', &
         userfile
      return
   end if

! ====================================================================
!
!     get size of level dimension
! 
   STATUS =  NF_INQ_DIMID ( NCID, 'lev', lev_dimID )
   if ( STATUS .NE.  NF_NOERR )  then
      write( 6,* )'ERROR - interplevs.F:', &
         'Cant get variable dim for lev'
      return
   endif

   STATUS = NF_INQ_DIMLEN( NCID, lev_dimID, nlev )
!
!     extract the pressure level data - in the saved-initial-conditions
!       datasets, lev is expressed in millibars, not pascals
!
!
   call getncdata ( NCID, 1,1,1, &
      'lev', dplevs, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - :', &
         'Cant get variable lev'
      return
   endif

   call getncdata( NCID, 1, 1, 1, 'PS', tmpdata, STATUS)
   if ( STATUS .NE. NF_NOERR ) then
      write(6,*)'ERROR - readsaveinit:', &
         'Cant get variable Ps'
      STATUS = NF_CLOSE (NCID)
      return
   endif
   ps(1,1,n3)=tmpdata
!=====================================================================

   call getncdata( NCID, 1, 1, 1, 'PHIS', phis, STATUS)
   if ( STATUS .NE. NF_NOERR ) then
      write(6,*)'ERROR - readsaveinit:', &
         'Cant get variable phis'
      STATUS = NF_CLOSE (NCID)
      return
   endif

   have_surfdat = .false.

   call getinterpncdata( NCID, 1, 1, 1, 'U',    &
      have_surfdat, surfdat, .false., &
      dplevs, nlev, uobs, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readinitdata:', &
         'Cant get variable U'
      STATUS = NF_CLOSE( NCID )
      return
   endif
   do k=1, PLEV
      u3(1,k,1,n3) = uobs(k)        !     set u to uobs 
   end do



   call getinterpncdata( NCID, 1, 1, 1, 'V',    &
      have_surfdat, surfdat, .false., &
      dplevs, nlev, vobs, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readinitdata:', &
         'Cant get variable V'
      STATUS = NF_CLOSE( NCID )
      return
   endif
   do k=1, PLEV
      v3(1,k,1,n3) = vobs(k)        !     set v to vobs 
   end do



   call getinterpncdata( NCID, 1, 1, 1, 'T',    &
      have_surfdat, surfdat, .false., &
      dplevs, nlev, tobs, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readinitdata:', &
         'Cant get variable T'
      STATUS = NF_CLOSE( NCID )
      return
   endif
   do k=1, plev
      t3(1,k,1,n3) = tobs(k)
   end do

   call getinterpncdata( NCID, 1, 1, 1, 'Q',    &
      have_surfdat, surfdat, .false., &
      dplevs, nlev, qobs, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readinitdata:', &
         'Cant get variable Q '
      STATUS = NF_CLOSE( NCID )
      return
   endif
   do k=1, plev
      q3(1,k,1,1,n3) = qobs(k)
   end do

   call getinterpncdata( NCID, 1, 1, 1, 'OMEGA',    &
      have_surfdat, surfdat, .false., &
      dplevs, nlev, wfld, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readinitdata:', &
         'Cant get variable OMEGA'
      STATUS = NF_CLOSE( NCID )
      return
   endif

   call getinterpncdata( NCID, 1, 1, 1, 'DIVQ',    &
      have_surfdat, surfdat, .false., &
      dplevs, nlev, divq, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readinitdata:', &
         'Cant get variable DIVQ'
      STATUS = NF_CLOSE( NCID )
      return
   endif

   call getinterpncdata( NCID, 1, 1, 1, 'DIVT',    &
      have_surfdat, surfdat, .false., &
      dplevs, nlev, divt, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readinitdata:', &
         'Cant get variable DIVT'
      STATUS = NF_CLOSE( NCID )
      return
   endif

   call getinterpncdata( NCID, 1, 1, 1, 'DIVU',    &
      have_surfdat, surfdat, .false., &
      dplevs, nlev, divu, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readinitdata:', &
         'Cant get variable DIVU'
      STATUS = NF_CLOSE( NCID )
      return
   endif

   call getinterpncdata( NCID, 1, 1, 1, 'DIVV',    &
      have_surfdat, surfdat, .false., &
      dplevs, nlev, divv, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readinitdata:', &
         'Cant get variable DIVV'
      STATUS = NF_CLOSE( NCID )
      return
   endif


   STATUS = NF_CLOSE( NCID )
   error_code = 0

   return
end subroutine readsaveinit
