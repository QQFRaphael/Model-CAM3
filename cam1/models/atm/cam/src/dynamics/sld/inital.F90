#include <misc.h>
#include <params.h>

subroutine inital

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
! $Id: inital.F90,v 1.13.2.14 2005/03/08 17:06:54 eaton Exp $
! $Author: eaton $
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, masterproc
   use inidat,       only: read_inidat
   use prognostics,  only: ps, u3, v3, t3, q3, qm1, div, ptimelevels, initialize_prognostics
   use buffer,       only: initialize_buffer
   use scanslt,      only: slt_alloc
   use comsrf,       only: initialize_comsrf
   use chem_surfvals, only: chem_surfvals_init
   use ioFileMod,    only: getfil
   use radae,        only: initialize_radbuffer
   use phys_buffer,  only: pbuf_allocate
   use phys_grid,    only: phys_grid_init
#if (defined SPMD)
   use spmd_dyn,     only: spmdbuf
#endif
   use time_manager, only: timemgr_init
   use filenames,    only: ncdata, bnd_topo
#if (defined COUP_CSM)
   use ccsm_msg,     only: initialize_ccsm_msg
#endif
!-----------------------------------------------------------------------
   implicit none
!------------------------------Parameters-------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comlun.h>
!-----------------------------------------------------------------------
   include 'netcdf.inc'
!---------------------------Local variables-----------------------------
!
   integer k,n              ! indices
   character(len=256) locfn ! local filename
!
!-----------------------------------------------------------------------
!
!
! Obtain initial dataset and topography dataset
!
   if (masterproc) then
      call getfil (ncdata, locfn)
      call wrap_open (locfn, NF_NOWRITE, ncid_ini)

! Backward compatibility: look for boundary data on initial file if boundary file name not provided.
      if (trim(bnd_topo) /= 'bnd_topo') then
         call getfil(bnd_topo, locfn)
         call wrap_open(locfn, NF_NOWRITE, ncid_topo)
      else
         ncid_topo = ncid_ini
      end if
   end if
!
! Check for consistent settings on initial dataset
!
   call readinitial (ncid_ini)

! Initialize time manager.

   call timemgr_init()

! Initialize ghg surface values before default initial distributions 
! are set in inidat.
   call chem_surfvals_init()
!
! Initialize comslt and prognostics variables 
!
   call initialize_prognostics
   call slt_alloc()
!
! Set commons
!
   call initcom
!
! Define physics data structures
!
   call phys_grid_init

#if (defined SPMD)
! Allocate communication buffers for
! collective communications in realloc
! routines and in dp_coupling
!
   call spmdbuf ()
#endif

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
   call pbuf_allocate('global')
   call initialize_buffer
   call initialize_comsrf
   call initialize_radbuffer
!
! Read in initial data
!
   call read_inidat

! Close the topographic boundary dataset
! Backward compatibility: don't close if ncid_topo = ncid_ini
   if (masterproc .and. (ncid_topo /= ncid_ini)) call wrap_close(ncid_topo)
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
      div(:,:,:,n)  = div(:,:,:,1)
   end do
   qm1(:,:,:,:) = q3(:,:,:,:,1)
   call print_memusage ('post-inidat')

   return
end subroutine inital
