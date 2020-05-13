!------------------------------------------------------------------------
! File: getinterpncdata.F 
! Author: John Truesdale (jet@ucar.edu) 
! $Id: getinterpncdata.F90,v 1.1.6.2 2004/08/19 15:05:40 jmccaa Exp $
!
!------------------------------------------------------------------------
#include <params.h>
#include <max.h>
subroutine getinterpncdata( NCID, latIdx, lonIdx, TimeIdx, &
   varName, have_surfdat, surfdat, fill_ends, &
   press, npress,  outData, STATUS )



!     getinterpncdata: extracts the entire level dimension for a 
!     particular lat,lon,time from a netCDF file
!     and interpolates it onto the input pressure levels, placing
!     result in outData, and the error status inx STATUS

   use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
   use pmgrid
   use scamMod, only: switch
   implicit none
!-----------------------------------------------------------------------
#include <netcdf.inc>


!     ---------- inputs ------------

   integer     nlev          ! number of levels in dataset
   integer     NCID          ! NetCDF ID
   integer     latIdx        ! latitude index
   integer     lonIdx        ! longitude index
   integer     TimeIdx       ! time index
   real(r8)        surfdat       ! surface value to be added before interpolation
   logical     have_surfdat  ! is surfdat provided
   logical     fill_ends ! extrapolate the end values
   integer     npress        ! number of dataset pressure levels
   real(r8)        press(npress) ! dataset pressure levels

!     ---------- outputs ----------

   real(r8)        outData(*)    ! interpolated output data
   integer     STATUS        ! return status of netcdf calls

!     -------  locals ---------

   real(r8)        tmp( MAX_DATASET_LEVS )
   real(r8)        dx, dy, m             ! slope for interpolation of surfdat
   integer     varID
   integer     var_ndims
   integer     dims_set
   integer     i
   integer     var_dimIDs( NF_MAX_VAR_DIMS )
   integer     start( NF_MAX_VAR_DIMS ) 
   integer     count( NF_MAX_VAR_DIMS )

   character   varName*(*)
   character   dim_name*( NF_MAX_NAME )
   real(r8)        missing_val
   logical     usable_var
   logical     use_nf_real
!
! Check mode: double or single precision
!

#if USE_4BYTE_REAL
   use_nf_real = .true.
#else
   use_nf_real = .false.
#endif

!
! Get var ID.  Check to make sure it's there.
!
   STATUS = NF_INQ_VARID( NCID, varName, varID )

   if ( STATUS .NE. NF_NOERR )  return

!
! Check the var variable's information with what we are expecting
! it to be.
!

   STATUS = NF_INQ_VARNDIMS( NCID, varID, var_ndims )
   if ( var_ndims .GT. 4 ) then
      write( 6,* ) 'ERROR - extractdata.F: The input var',varName, &
         'has', var_ndims, 'dimensions'
      STATUS = -1
   endif

!
!     surface variables
!
   if ( var_ndims .EQ. 0 ) then
      if (use_nf_real) then
         STATUS = NF_GET_VAR_REAL( NCID, varID, outData )
      else
         STATUS = NF_GET_VAR_DOUBLE( NCID, varID, outData )
      endif

      return
   endif

   STATUS = NF_INQ_VARDIMID( NCID, varID, var_dimIDs )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* ) 'ERROR - extractdata.F:Cant get dimension IDs for', varName
      return
   endif
!     
!     Initialize the start and count arrays 
!     
   dims_set = 0
   nlev = 1
   do i =  var_ndims, 1, -1

      usable_var = .false.
      STATUS = NF_INQ_DIMNAME( NCID, var_dimIDs( i ), dim_name )

      if ( dim_name .EQ. 'lat' ) then
         start( i ) =  latIdx
         count( i ) = 1           ! Extract a single value 
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'lon' ) then
         start( i ) = lonIdx
         count( i ) = 1           ! Extract a single value
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'lev' ) then
         STATUS = NF_INQ_DIMLEN( NCID, var_dimIDs( i ), nlev )
         start( i ) = 1
         count( i ) = nlev       ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'ilev' ) then
         STATUS = NF_INQ_DIMLEN( NCID, var_dimIDs( i ), nlev )
         start( i ) = 1
         count( i ) = nlev        ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'time' .OR. dim_name .EQ. 'tsec' ) then 
         start( i ) = TimeIdx
         count( i ) = 1           ! Extract a single value 
         dims_set = dims_set + 1   
         usable_var = .true.
      endif

      if ( usable_var .EQV. .false. ) then
         write( 6,* )'ERROR - extractdata.F: The input var ',varName, &
            ' has an unusable dimension ', dim_name
         STATUS = 1
      endif
   end do

   if ( dims_set .NE. var_ndims ) then
      write( 6,* )'ERROR - extractdata.F: Could not find all the', &
         ' dimensions for input var ', varName
      write( 6,* )'Found ',dims_set, ' of ',var_ndims
      STATUS = 1
   endif

   if (use_nf_real) then
      STATUS = NF_GET_VARA_REAL( NCID, varID, start, count, tmp )
   else
      STATUS = NF_GET_VARA_DOUBLE( NCID, varID, start, count, tmp )
   endif

   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - extractdata.F: Could not get data for input var ', varName
      return
   endif

   if ( nlev .eq. 1 ) then
      outdata(1) = tmp(1)
      return                 ! no need to do interpolation 
   endif
!   if ( use_ccmiop .and. nlev.eq.plev) then
   if ( nlev.eq.plev) then
      outData(:nlev)= tmp(:nlev)! no need to do interpolation 
   else
!
!     add the surface data if available, else
!     fill in missing surface data by extrapolation
!
      if(.not.switch(CRM_SW+1)) then
      if ( have_surfdat ) then
         tmp(npress) = surfdat
      else
         dy = press(npress-1) - press(npress-2)
         dx = tmp(npress-1) - tmp(npress-2)
         if ( dx .ne. 0.0 ) then
            m = dy/dx
            tmp(npress) = ((press(npress) - press(npress-1)) / m ) + tmp(npress-1)
         else
            tmp(npress) = tmp(npress-1)
         endif
         surfdat = tmp(npress)
         endif
      endif

#if DEBUG > 1
!
!     check data for missing values
!

      STATUS = NF_GET_ATT_DOUBLE( NCID, varID, 'missing_value', missing_val )
      if ( STATUS .NE. NF_NOERR ) then
         missing_val = -9999999.0
      endif
!
! reset status to zero
!     
      STATUS = 0
!
      do i=1, npress
         if ( tmp(i) .eq. missing_val ) then
            print *, 'ERROR - missing value found in ', varname
            print *, 'time,lat,lon,lev = ' ,timeidx, latidx, lonidx, i
            stop
         endif
      enddo
#endif
!
      call interplevs( tmp, press, npress, fill_ends,outdata )
   endif
   return
 end subroutine getinterpncdata

