#include <params.h>
#include <max.h>
!------------------------------------------------------------------------
! File: interplevs.F 
! Author: John Truesdale (jet@ucar.edu) 
! $Id: interplevs.F90,v 1.1.6.1 2004/05/18 20:51:28 jet Exp $
!
!------------------------------------------------------------------------
subroutine interplevs( inputdata,   dplevs,   nlev, &
                       fill_ends,   outdata)

   use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
   use pmgrid
   use buffer
   use prognostics
   implicit none

!
!     WARNING: ps, siga and sigb must be initialized before calling this routine
!

!------------------------------Commons----------------------------------
#include <comhyb.h>
!-----------------------------------------------------------------------


!     ------- inputs -----------
   integer nlev                 ! num press levels in dataset

   real(r8) inputdata(nlev)     ! data from netcdf dataset
   real(r8) dplevs(nlev)         ! input data pressure levels 

   logical fill_ends            ! fill in missing end values(used for
! global model datasets)


! ------- outputs ----------
   real(r8) outdata( PLEV )      ! interpolated column data

! ------- locals -----------

   integer wksplen 
   parameter ( wksplen = 3 * MAX_DATASET_LEVS) !workspace length

   real(r8) mplevs( PLEV )
   real(r8) workspace( wksplen )
   real(r8) interpdata( PLEV )


   integer dstart_lev, dend_lev 
   integer mstart_lev, mend_lev
   integer data_nlevs, model_nlevs, i
   integer wkspminlen , STATUS

!
!     Initialize  model_pressure_levels.  ps should be set in the calling
!     routine to the value in the dataset
!
   do i = 1, plev
      mplevs( i ) = 1000.0 * hyam( i ) + ps(1,1,n3) * hybm( i ) / 100.0
   end do
!     
!     the following algorithm assumes that pressures are increasing in the
!     arrays
!     
!     
!     Find the data pressure levels that are just outside the range
!     of the model pressure levels, and that contain valid values
!     
   dstart_lev = 1
   do i= 1, nlev
      if ( dplevs(i) .LE. mplevs(1) ) dstart_lev  = i
   end do

   dend_lev = nlev
   do i= nlev, 1, -1
      if ( dplevs(i) .GE. mplevs(plev) ) then
         dend_lev  = i
      endif
   end do
!         
!     Find the model pressure levels that are just inside the range
!     of the data pressure levels
!
   mstart_lev = 1
   do i=plev, 1, -1
      if ( mplevs( i ) .GE. dplevs( dstart_lev ) )  mstart_lev = i
   end do

   mend_lev = plev
   do i=1,plev
      if ( mplevs( i ) .LE. dplevs( dend_lev ) ) mend_lev = i
   end do

   data_nlevs = dend_lev - dstart_lev +1
   model_nlevs = mend_lev - mstart_lev +1
!
!     interpolate data onto the model pressure levels
!
   call lininterp (inputdata,dplevs(dstart_lev),1,data_nlevs, &
		   interpdata,mplevs(mstart_lev),model_nlevs)
   do i=1 , model_nlevs
      outdata( i+mstart_lev-1 ) = interpdata( i )
   end do
!
!     fill in the missing end values 
!           (usually  done if this is global model dataset)
!
   if ( fill_ends ) then 
      do i=1, mstart_lev
         outdata(i) = inputdata(1)
      end do
      do i= mend_lev, plev
         outdata(i) = inputdata(nlev)
      end do
   end if

   return
end subroutine interplevs
