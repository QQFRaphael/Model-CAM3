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
! $Id: inital.F90,v 1.15.2.17 2005/03/08 17:06:46 eaton Exp $
! $Author: eaton $
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plev, plevp, beglat, endlat, masterproc
   use inidat,       only: read_inidat
   use prognostics,  only: ps, u3, v3, t3, q3, vort, div, ptimelevels, initialize_prognostics, pdeld
   use buffer,       only: initialize_buffer
   use comspe,       only: alp, dalp
   use comsrf,       only: initialize_comsrf
   use chem_surfvals, only: chem_surfvals_init
   use ioFileMod,    only: getfil
   use radae,        only: initialize_radbuffer
   use phys_buffer,  only: pbuf_allocate
   use phys_grid,    only: phys_grid_init
   use constituents, only: ppcnst, cnst_need_pdeldry
   use rgrid,   only: nlon

#if (defined SPMD)
   use spmd_dyn,     only: spmdbuf
#endif
   use time_manager, only: timemgr_init
   use filenames,    only: ncdata, bnd_topo
   use scanslt,      only: scanslt_alloc
#if (defined BFB_CAM_SCAM_IOP )
   use history,         only: initialize_iop_history
#endif
#if (defined COUP_CSM)
   use ccsm_msg, only: initialize_ccsm_msg
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
   integer n,i,k,lat              ! index
   real(r8) pdel(plon,plev)     ! pressure arrays needed to calculate
   real(r8) pint(plon,plevp)    !     pdeld
   real(r8) pmid(plon,plev)

   character(len=256) locfn ! local filename
!
!-----------------------------------------------------------------------
!
!
! Obtain initial and topography datasets
!
   if (masterproc) then
      call getfil (ncdata, locfn)
      call wrap_open (locfn, NF_NOWRITE, ncid_ini)

! Backward compatibility: look for topography data on initial file if topo file name not provided.
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
! Initialize prognostics variables 
!
   call initialize_prognostics
   call scanslt_alloc()
!
! Set commons
!
   call initcom
!
! Define physics data structures
!
   call phys_grid_init
!
#if (defined SPMD)
! Allocate communication buffers for
! collective communications in realloc
! routines and in dp_coupling
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

#if (defined BFB_CAM_SCAM_IOP )
   call initialize_iop_history
#endif
!
! Read in initial data
!
   call read_inidat

! Close the topography dataset
! Backward compatibility: don't close if ncid_topo = ncid_ini
   if (masterproc .and. (ncid_topo /= ncid_ini)) call wrap_close(ncid_topo)

! Recover space used for ALP and DALP arrays
! (no longer needed after spectral truncations
! inside of read_inidat)
   deallocate ( alp )
   deallocate ( dalp )

!
! If dry-type tracers are present, initialize pdeld
! First, set current time pressure arrays for model levels etc. to get pdel
!
      do lat=beglat,endlat
         call plevs0(nlon(lat), plon, plev, ps(1,lat,1), pint, pmid, pdel)
         if (  cnst_need_pdeldry ) then
            do k=1,plev
               do i=1,nlon(lat)
                  pdeld(i,k,lat,1) = pdel(i,k)*(1.-q3(i,k,1,lat,1))
               end do !i
            end do !k
         endif !cnst_need_pdeldry
      end do !lat
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
      q3(1:plon,:,:,:,n) = q3(1:plon,:,:,:,1)
      vort(:,:,:,n) = vort(:,:,:,1)
      div(:,:,:,n)  = div(:,:,:,1)
      if (  cnst_need_pdeldry )  pdeld(1:plon,:,:,n) = pdeld(1:plon,:,:,1)
   end do
   call print_memusage ('post-inidat')

   return
end subroutine inital
