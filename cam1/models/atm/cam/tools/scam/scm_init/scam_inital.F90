#include <misc.h>
#include <params.h>
#include <max.h>

subroutine scam_inital(error_code)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Define initial conditions for first run of case
! 
! Method: 
! 
! Author: 
! Original version:  CCM1
! Standardized:      L. Bath, June 1992
!                    T. Acker, March 1996
! Reviewed:          B. Boville, April 1996
!
! $Id: scam_inital.F90,v 1.1.6.3 2004/12/10 16:25:44 jeff Exp $
! $Author: jeff $
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,       only: pcols,begchunk,endchunk,pver
   use pmgrid,       only: plon, masterproc, plev
   use inidat,       only: read_inidat
   use prognostics, only: ps, u3, v3, t3, q3, vort, div, ptimelevels, initialize_prognostics,n3,n3m1,n3m2
   use buffer,       only: initialize_buffer
   use scanslt,      only: scanslt_alloc
   use comsrf,       only: initialize_comsrf,frzmlt,Focn,tsocn
   use chem_surfvals, only: chem_surfvals_init
   use ioFileMod,    only: getfil
   use radae,        only: initialize_radbuffer
   use phys_buffer,  only: pbuf_allocate,pbuf, pbuf_times, pbuf_get_fld_idx
   use phys_grid,    only: phys_grid_init
   use time_manager, only: timemgr_init,nestep, get_curr_date, get_curr_calday
   use filenames,    only: ncdata
   use constituents, only: cnst_get_ind
   use history,      only: initialize_iop_history
   use scamMod, only :isrestart,use_iop,use_analysis,use_saveinit, &
	use_pert_init,have_u,have_v

#if (defined COUP_CSM)
   use ccsm_msg, only: initialize_ccsm_msg
#endif
use infnan
!-----------------------------------------------------------------------
   implicit none
!------------------------------Parameters-------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comlun.h>
!-----------------------------------------------------------------------
#include <runtype.h>
!-----------------------------------------------------------------------
   include 'netcdf.inc'
!---------------------------Local variables-----------------------------
!
   integer n,m,k,i,lat      ! index
   integer :: yr, mon, day   ! components of a date
   integer :: ncdate         ! current date in integer format [yyyymmdd]
   integer :: ncsec          ! current time of day [seconds]

   real(r8) :: calday            ! current calendar day
   real(r8) :: tprt,qprt,pertval
   integer :: ixcldice, ixcldliq
   real(r8), pointer, dimension(:,:,:,:) :: qcwat, lcwat, tcwat, cld
!
!-----------------------------------------------------------------------
   integer error_code     ! returns netcdf errors
   character(len=256) locfn ! local filename
!
!-----------------------------------------------------------------------
!
!
!     Trap ieee exceptions on SUN for debugging purposes
!
#if ( defined sun )
   iexcept = ieee_handler( 'set', 'common', myhandler )
   if ( iexcept .ne. 0 ) write(6,*)'ieee trapping not supported here'
#endif

! Initialize time manager.

   nestep=20
   call timemgr_init()
   calday = get_curr_calday()
   call get_curr_date(yr, mon, day, ncsec)
   ncdate = yr*10000 + mon*100 + day

! Initialize chem surface values before default initial distributions 
! are set in inidat.
   if (.not. isrestart)   call chem_surfvals_init()
!
! Initialize prognostics variables 
!
   if (.not. isrestart) call initialize_prognostics
   if (.not. isrestart) call scanslt_alloc()

   call readpressdata(error_code)
   if (error_code .ne. 0) then
      print *, 'Got error code in readpressdata returning'
      return
   end if

!
! Set commons
!
   call initcom
!
! Define physics data structures
!
!
! Define physics data structures
!
   if (.not.isrestart)    call phys_grid_init

#if (defined COUP_CSM)
!
! Initialize ccsm arrays (must be done after phys_grid_init where
! begchunk and endchunk are defined
!
   call initialize_ccsm_msg
#endif
!
! Initialize buffer, comsrf, and radbuffer variables 
! (which must occur after the call to phys_grid_init)
!
   if (.not. isrestart) call pbuf_allocate('global')
   if (.not. isrestart) call initialize_buffer
   if (.not. isrestart) call initialize_comsrf
   if (.not. isrestart) call initialize_radbuffer

   if (.not. isrestart) call initialize_iop_history
!
! Read in initial data
!
      if (use_iop) then
         call setiopupdate
         call readiopdata(error_code)
      endif

      call read_inidat(MODEL) 

      if (use_analysis) then
	call read_inidat(ANAL) 
     end if

      if (use_iop) then
         call setiopupdate
         call readiopdata(error_code)
!
! read_inidat reads all values in at (:,:,1) and these values
! are copied into all three time levels below. Readiopdata is called
! in the time loop and reads values into the n3 position so we must
! update the (:,:,1) position to let initial readiopcall work with
! model code.
!
         ps(:,:,1)     = ps(:,:,n3)
         if(have_u)u3(:,:,:,1)   = u3(:,:,:,n3)
         if(have_v)v3(:,:,:,1)   = v3(:,:,:,n3)
         t3(:,:,:,1)   = t3(:,:,:,n3)
         q3(:,:,:,:,1) = q3(:,:,:,:,n3)
         vort(:,:,:,1) = vort(:,:,:,n3)
         div(:,:,:,1)  = div(:,:,:,n3)


      endif

   if (error_code .ne. 0) then
      print *, 'Got error code in readiopdate and returning'
      return
   end if

   if (use_saveinit) call readsaveinit(error_code) 
   if (error_code .ne. 0) then
      print *, 'Got error code in readsaveinit and returning'
      return
   end if
!
! add random perturbation to the initial condition for T and q
!

   tprt = 0.9                ! +/- .9 degrees
   qprt = .06                ! +/- 6 %
   if ( use_pert_init ) then
!   if (seedval==0.) then
      call random_seed()
      write(6,*)'Initializing f90 pseudo random generator with random value'
!   else
!      call random_seed(put=seedval)
!      write(6,*)'Initializing f90 pseudo random generator with value ',seedval
!   end if
      write(6,*)'INIDAT: Adding random perturbation bounded by +/- ', &
         tprt,'degrees to initial temperature field'
      write(6,*)'INIDAT: Adding random perturbation bounded by +/- ', &
         qprt*100  ,'% to initial Q field'
      do k=1,plev
	 call random_number(pertval)
         pertval = tprt * (2.*(0.5 - pertval))
!         t3(1,k,1,n3) = t3(1,k,1,n3) + pertval
         t3(1,k,1,1) = t3(1,k,1,1) + pertval
      end do
      do k=1,plev
	 call random_number(pertval)
         pertval = qprt * (2.*(0.5 - pertval))
!         q3(1,k,1,1,n3) = q3(1,k,1,1,n3) * (1. + pertval)
         q3(1,k,1,1,1) = q3(1,k,1,1,1) * (1. + pertval)
      end do
   end if
!
! Make all time levels of prognostics contain identical data.
! Fields to be convectively adjusted only *require* n3 time
! level since copy gets done in linems.
!

!       
! Make all time levels of prognostics contain identical data.
! Fields to be convectively adjusted only *require* n3 time
! level since copy gets done in linems.
!
   do n=2,ptimelevels
      ps(:,:,n)     = ps(:,:,1)
      u3(:,:,:,n)   = u3(:,:,:,1)
      v3(:,:,:,n)   = v3(:,:,:,1)
      t3(:,:,:,n)   = t3(:,:,:,1)
      q3(:,:,:,:,n) = q3(:,:,:,:,1)
      vort(:,:,:,n) = vort(:,:,:,1)
      div(:,:,:,n)  = div(:,:,:,1)
   end do
   
#if ( defined COUP_SOM )
! define an initial ocean fraction
      ocnfrac(:pcols,:) = 1. - landfrac(:pcols,:) - icefrac(:pcols,:)
      write(6,*)'INIDAT: ocnfrac=',ocnfrac(1,begchunk)
!
!JR Could read in Focn from initial dataset if available
      Focn(:,:) = 0.
#else
      Focn(:,:) = inf
      frzmlt(:,:) = 0.  ! needs to be 0, otherwise test in tstm always true
      tsocn(:,:) = inf
#endif

   return
end subroutine scam_inital
