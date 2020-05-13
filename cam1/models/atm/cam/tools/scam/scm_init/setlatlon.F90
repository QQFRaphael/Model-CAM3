#include <params.h>
#include <max.h>
!------------------------------------------------------------------------
! File: setlatlon.F 
! Author: John Truesdale (jet@ucar.edu)
! $Id: setlatlon.F90,v 1.1.6.1 2004/05/18 20:51:32 jet Exp $
!
!------------------------------------------------------------------------
subroutine setlatlonidx( )
   use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
   use pmgrid
   use prognostics
   use buffer
   use comsrf
   use scamMod, only :modelfile,columnlat,columnlon,initlonidx,initlatidx
   implicit none
!-----------------------------------------------------------------------
#include <comfrc.h>
!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------


!------------------------------locals-----------------------------------

   real(r8) datalats( MAX_LAT_DIM )
   real(r8) datalons( MAX_LON_DIM )
   real(r8) prev, next             !
   integer  i                ! 
   integer  status
   integer  nlat_dimID, nlat, lat_varID 
   integer  nlon_dimID, nlon, lon_varID
   integer  NCID
   logical  use_nf_real

!
! Check mode: double or single precision.
! 

#if USE_4BYTE_REAL
   use_nf_real = .true.
#else
   use_nf_real = .false.
#endif


   STATUS = NF_OPEN( modelfile, NF_NOWRITE, NCID )
   if( STATUS .NE. NF_NOERR )  then
      write( 6,* )'ERROR - setlatlon.F:', &
         'Cant open model data netcdf file', &
         modelfile
      STATUS = NF_CLOSE( NCID )
      return
   endif

   STATUS =  NF_INQ_DIMID ( NCID, 'lat', nlat_dimID )
   if ( STATUS .NE.  NF_NOERR )  then
      write( 6,* )'ERROR - setlatlon.F:', &
         'Cant get variable dim for lat'
      STATUS = NF_CLOSE( NCID )
      return
   endif
   STATUS =  NF_INQ_DIMID ( NCID, 'lon', nlon_dimID )
   if ( STATUS .NE.  NF_NOERR )  then
      write( 6,* )'ERROR - setlatlon.F:', &
         'Cant get variable dim for lon'
      STATUS = NF_CLOSE( NCID )
      return
   endif

   STATUS = NF_INQ_DIMLEN( NCID, nlat_dimID, nlat )
   STATUS = NF_INQ_DIMLEN( NCID, nlon_dimID, nlon )

   STATUS = NF_INQ_VARID( NCID, 'lat', lat_varID )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - setlatlon.F:', &
         'Cant get variable ID for lat'
      STATUS = NF_CLOSE( NCID )
      return
   endif

   if (use_nf_real) then
      STATUS = NF_GET_VAR_REAL( NCID, lat_varID, datalats )      
   else
      STATUS = NF_GET_VAR_DOUBLE( NCID, lat_varID, datalats )
   endif

   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - setlatlon:', &
         'Cant get variable lat'
      STATUS = NF_CLOSE ( NCID )
      return
   endif

   STATUS = NF_INQ_VARID( NCID, 'lon', lon_varID )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - setlatlon.F:', &
         'Cant get variable ID for lon'
      STATUS = NF_CLOSE( NCID )
      return
   endif

   if (use_nf_real) then
      STATUS =  NF_GET_VAR_REAL( NCID, lon_varID, datalons )         
   else
      STATUS = NF_GET_VAR_DOUBLE( NCID, lon_varID, datalons )
   endif

   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - setlatlon.F:', &
         'Cant get variable lon'
      STATUS = NF_CLOSE ( NCID )
      return
   endif

   initLatIdx = 1
   do i = 1, nlat
      next = datalats(i)
      prev = datalats(initLatIdx)
      if ( abs(columnLat - next) .LT.  &
         abs(columnLat - prev) ) then
         initLatIdx = i
      endif
   enddo

   initLonIdx = 1
   do i = 1, nlon
      if ( datalons(i) .le. 180.0 ) then 
         next = datalons(i)       
      else
         next = datalons(i) - 360.0
      endif
      if ( datalons(initLonIdx) .le. 180.0 ) then 
         prev = datalons(initLonIdx)
      else
         prev = datalons(initLonIdx) - 360.0
      endif


      if ( abs(columnLon - next) .LT. &
         abs(columnLon - prev) ) then
         initLonIdx = i
      endif
   enddo

   return
end subroutine setlatlonidx




