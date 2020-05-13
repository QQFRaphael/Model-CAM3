#include <misc.h>
#include <params.h>

module spmd_dyn

!----------------------------------------------------------------------- 
! 
! Purpose: SPMD implementation of CAM spectral Eulerian dynamics.
! 
! Author: CCM Core Group
! Modified: P. Worley, September 2002, November 2003, December 2003,
!                      November 2004, January 2005
! 
!-----------------------------------------------------------------------

#if (defined SPMD)

   use shr_kind_mod, only: r8 => shr_kind_r8
   use rgrid,        only: nlon
   use pmgrid,       only: plat, masterproc, iam, numlats, &
                           beglat, endlat, begirow, endirow, plev
   use scanslt,      only: beglatex, endlatex, numbnd, numlatsex
   use mpishorthand, only: mpir8, mpicom
   use infnan,       only: inf
   use abortutils,   only: endrun

   implicit none

   private
   save

   public spmdinit_dyn, compute_gsfactors, spmdbuf
   public spmd_dyn_defaultopts, spmd_dyn_setopts

   logical, public :: local_dp_map=.false. ! flag indicates that mapping between dynamics 
                                           !  and physics decompositions does not require 
                                           !  interprocess communication
   integer, public :: block_buf_nrecs      ! number of local grid points (lon,lat,lev)
                                           !  in dynamics decomposition (including level 0)
   integer, public :: chunk_buf_nrecs      ! number of local grid points (lon,lat,lev)
                                           !  in physics decomposition (including level 0)

   integer, public :: npes                 ! Total number of MPI tasks
   integer, public :: nsmps                ! Total number of SMP nodes
   integer, public, allocatable ::        &
    cut(:,:),                             &! partition for MPI tasks
    cutex(:,:)                             ! extended partition 
   integer, public :: proc(plat)           ! MPI task id associated with a given lat.
   integer, public :: neighs               ! number of south neighbors to comm guardcells
   integer, public, allocatable :: neighs_proc(:)    ! sorted south process neighbors
   integer, public :: neighn               ! number of north neighbors to comm guardcells
   integer, public, allocatable :: neighn_proc(:)    ! sorted north process neighbors
   integer, public :: npessp               ! number of MPI tasks in spectral space
   integer, public :: maxlats              ! max number of lats on any MPI task
   integer, public :: maxcols              ! max number of columns on any MPI task
   integer, public, allocatable :: nlat_p(:)    ! number of latitudes per MPI task
   integer, public, allocatable :: ncol_p(:)    ! number of columns per MPI task
   integer, public, allocatable :: proc_smp_map(:) ! map of process/SMP node assignments
   integer, public :: realloc4_steps       ! number of swaps in realloc4 algorithms
   integer, public, allocatable :: realloc4_proc(:)
                                           ! swap partner in each step of 
                                           ! realloc4 algorithms
   integer, public, allocatable :: realloc4_step(:)
                                           ! step in realloc4 algorithms
                                           ! in which communicate with a given
                                           ! process
   integer, public :: allgather_steps      ! number of swaps in allgather algorithm
   integer, public, allocatable :: allgather_proc(:)
                                           ! swap partner in each step of 
                                           ! allgather (realloc5/7) algorithm
   integer, public, allocatable :: allgather_step(:)
                                           ! step in allgather (realloc5/7) algorithm
                                           ! in which communicate with a given
                                           ! process
!
   logical, private, parameter :: def_equi_by_col = .true.          ! default
   logical, private :: dyn_equi_by_col = def_equi_by_col 
                                           ! flag indicating whether to assign
                                           ! latitudes to equidistribute columns or
                                           ! latitudes. This only matters when using a
                                           ! reduced grid.
!
   logical, private, parameter :: def_mirror = .false.          ! default
   logical, private :: mirror = def_mirror ! flag indicating whether latitudes and their
                                           ! reflections across the equator should assigned 
                                           ! to consecutive processes
!
! Dynamics communication transpose algorithm option:
!  0: use mpi_alltoallv
!  1: use point-to-point MPI-1 two-sided implementation
!  2: use point-to-point MPI-2 one-sided implementation if supported, 
!       otherwise use MPI-1 implementation
!  3: use Co-Array Fortran implementation if supported, 
!       otherwise use MPI-1 implementation
   integer, private, parameter :: min_alltoall = 0
   integer, private, parameter :: max_alltoall = 3
   integer, private, parameter :: def_alltoall = 0         ! default
   integer, public :: dyn_alltoall  = def_alltoall
!
! Dynamics communication allgather (realloc5/7) algorithm option:
!  0: use mpi_allgatherv
!  1: use point-to-point MPI-1 two-sided implementation
!  2: use point-to-point MPI-2 one-sided implementation if supported, 
!       otherwise use MPI-1 implementation
!  3: use Co-Array Fortran implementation if supported, 
!       otherwise use MPI-1 implementation
   integer, private, parameter :: min_allgather = 0
   integer, private, parameter :: max_allgather = 3
   integer, private, parameter :: def_allgather = 0         ! default
   integer, public :: dyn_allgather = def_allgather
!
! Collective communication send/receive buffers
#if (defined CAF)
   real(r8), public, allocatable :: buf1(:)[:],buf2(:)[:] ! buffers for packing MPI msgs
#else
   real(r8), public, allocatable :: buf1(:),buf2(:) ! buffers for packing MPI msgs
#endif
   integer, public :: spmdbuf_siz = 0        ! buffer size (in r8s)
   integer, public :: buf1win                ! buf1 Window id
   integer, public :: buf2win                ! buf2 Window id

CONTAINS

!========================================================================

  subroutine spmd_dyn_defaultopts(npr_yz_out, geopktrans_out,    &
               tracertrans_out, ompnest_out, force_2d_out,       &
               modcomm_transpose_out, modcomm_geopk_out,         &
               dyn_alltoall_out, dyn_allgather_out,              &
               dyn_equi_by_col_out )
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
! Author: Art Mirin / Pat Worley
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
! FV-only arguments
     integer, intent(out), optional :: npr_yz_out(4)
     integer, intent(out), optional :: geopktrans_out
     integer, intent(out), optional :: tracertrans_out
     integer, intent(out), optional :: ompnest_out
     integer, intent(out), optional :: force_2d_out
     integer, intent(out), optional :: modcomm_transpose_out
     integer, intent(out), optional :: modcomm_geopk_out
! EUL/SLD arguments
     integer, intent(out), optional :: dyn_alltoall_out
     integer, intent(out), optional :: dyn_allgather_out
     logical, intent(out), optional :: dyn_equi_by_col_out
!-----------------------------------------------------------------------
     if ( present(dyn_alltoall_out) ) then
       dyn_alltoall_out = def_alltoall
     endif
     if ( present(dyn_allgather_out) ) then
       dyn_allgather_out = def_allgather
     endif
     if ( present(dyn_equi_by_col_out) ) then
       dyn_equi_by_col_out = def_equi_by_col
     endif
  end subroutine spmd_dyn_defaultopts

!========================================================================

  subroutine spmd_dyn_setopts(npr_yz_in, geopktrans_in,       &
               tracertrans_in, ompnest_in, force_2d_in,       &
               modcomm_transpose_in, modcomm_geopk_in,        &
               dyn_alltoall_in, dyn_allgather_in,             &
               dyn_equi_by_col_in )
!----------------------------------------------------------------------- 
! Purpose: Set runtime options
! Author: Art Mirin / Pat Worley
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
! FV-only arguments
     integer, intent(in), optional :: npr_yz_in(4)
     integer, intent(in), optional :: geopktrans_in
     integer, intent(in), optional :: tracertrans_in
     integer, intent(in), optional :: ompnest_in
     integer, intent(in), optional :: force_2d_in
     integer, intent(in), optional :: modcomm_transpose_in
     integer, intent(in), optional :: modcomm_geopk_in
! EUL/SLD arguments
     integer, intent(in), optional :: dyn_alltoall_in
     integer, intent(in), optional :: dyn_allgather_in
     logical, intent(in), optional :: dyn_equi_by_col_in
!-----------------------------------------------------------------------
     if ( present(dyn_alltoall_in) ) then
       dyn_alltoall = dyn_alltoall_in
       if ((dyn_alltoall.lt.min_alltoall).or. &
           (dyn_alltoall.gt.max_alltoall)) then
         write(6,*)                                          &
           'SPMD_DYN_SETOPTS:  ERROR:  dyn_alltoall=', &
           dyn_alltoall_in,                              &
           '  is out of range.  It must be between ',        &
           min_alltoall,' and ',max_alltoall
         call endrun
       endif
     endif
     if ( present(dyn_allgather_in) ) then
       dyn_allgather = dyn_allgather_in
       if ((dyn_allgather.lt.min_allgather).or. &
           (dyn_allgather.gt.max_allgather)) then
         write(6,*)                                          &
           'SPMD_DYN_SETOPTS:  ERROR:  dyn_allgather=', &
           dyn_allgather_in,                              &
           '  is out of range.  It must be between ',        &
           min_allgather,' and ',max_allgather
         call endrun
       endif
     endif
     if ( present(dyn_equi_by_col_in) ) then
       dyn_equi_by_col = dyn_equi_by_col_in
     endif
  end subroutine spmd_dyn_setopts

!========================================================================

   subroutine spmdinit_dyn ()
!----------------------------------------------------------------------- 
! 
! Purpose: Distribute latitudes among available processes
! 
! Method: Distribution is S->N for processes 0->npes
! 
! Author: CCM Core Group
! Modified: P. Worley, November 2003 to improve SMP load balance, and to
!           change distribution to 
!             S->E for processes 0,2,..,npes-2
!           and 
!             N->E for processes 1,3,..,npes-1
!           when mirror flag is set (at request of physics)
! Modified: P. Worley, November 2004 to improve load balance for 
!           reduced grid by equidistributing columns (not latitudes)
!           in latitude decomposition. Used when equi_by_col flag is set.
!           On by default, and gives identical decomposition as
!           equidistributing by latitude when using a full grid.
! 
!-----------------------------------------------------------------------
      use comspe, only: numm
      use spmd_utils
#if (defined MODCM_DP_TRANSPOSE)
      use parutilitiesmodule, only : parinit
#endif
!-----------------------------------------------------------------------
!
! Local workspace
!
      integer i         ! loop index
      integer tot_cols  ! total number of columns in computational grid
      integer m2,m3,m5  ! 2, 3, 5 prime factors for problem decomposition
      integer tot_nx    ! total number of latitudes/columns in 
                        ! computational grid
      integer nx_base   ! approx. number of latitudes/columns per proc
      integer nx_p(0:npes-1)      ! number of latitudes/columns per process
      integer nx_smp(0:npes-1)    ! number of latitudes/columns per SMP
      integer nproc_smp(0:npes-1) ! number of MPI processes per SMP
      integer workleft  ! amount of work still to be parcelled out

      integer smpid     ! SMP id
      integer smpids    ! SMP id for SH process
      integer smpidn    ! SMP id for NH process
      integer procj     ! process offset loop index
      integer procid    ! process id
      integer procids   ! process id SH
      integer procidn   ! process id NH

      integer max_ncols ! maximum number of columns assigned to a process
      integer min_max_ncols    ! minmax number of columns assigned 
                               ! to a process over all latitude assignments
      integer ncol      ! number of columns assigned to current process
      integer ncol_curtot  ! current total number of columns assigned
      integer ncol_curgoal ! target number of columns to be assigned to process
      integer lat       ! latitude index
      integer iend      ! ending latitude band of work for a given proc
      integer neighn_minlat(plat)    ! minimum latitude in north neighbor
      integer neighs_maxlat(plat)    ! maximum latitude in south neighbor

      real(r8) avgnx_proc(0:npes-1) ! average number of latitudes/columns per 
                                    ! MPI process in a given SMP node
      real(r8) minavgnx_proc        ! minimum average number of 
                                    ! latitudes/columns per 
                                    ! MPI process over SMP nodes
      real(r8) alpha    ! slop factor in assigning latitudes to processes
      real(r8) opt_alpha! best slop factor in assigning latitudes to processes

      logical done      ! exit flag for latitude assignment loop
!
!-----------------------------------------------------------------------
!
! Initialize Pilgrim library
!
#if (defined MODCM_DP_TRANSPOSE)
      call parinit(mpicom)
#endif
!
! Initialize mirror flag
!
      mirror = phys_mirror_decomp_req
!
! Allocate memory for MPI task partition array
! and extended partition
!
      allocate (cut  (2,0:npes-1))
      allocate (cutex(2,0:npes-1))
!
! Allocate memory for number of lats per proc
!
      allocate (nlat_p (0:npes-1))
      nlat_p(0:npes-1) = 0
!
! Allocate memory for number of columns per proc
!
      allocate (ncol_p (0:npes-1))
      ncol_p(0:npes-1) = 0
!
! determine total number of columns
!
      tot_cols = 0
      do lat=1,plat
         tot_cols = tot_cols + nlon(lat)
      enddo
!
! Make sure number of PEs, latitudes, and columns are kosher
!
      call factor (plat, m2, m3, m5)

      if (m2 < 1) then
         call endrun ('SPMDINIT_DYN: Problem size is not divisible by 2')
      end if

      if (masterproc) then
         write (6,*) 'Problem factors: 2**',m2,' * 3**',m3,' * 5**',m5
      end if
      call factor (npes, m2, m3, m5)
      
      if (mod(npes,2) /= 0) then
         write(6,*)'SPMDINIT_DYN: nprocs(',npes,') must be a multiple of 2'
         call endrun
      end if

      if ((dyn_equi_by_col) .and. (mod(tot_cols,2) /= 0)) then
         write(6,*)'SPMDINIT_DYN: Total number of columns(', &
                   tot_cols,') must be a multiple of 2'
         call endrun
      end if
!
! Determine approximate number of columns or latitudes per process
!
      if (dyn_equi_by_col) then
         tot_nx = tot_cols
      else
         tot_nx = plat
      endif
      nx_base = tot_nx/npes
      do procid=0,npes-1
         nx_p(procid) = nx_base
      enddo
!
! Calculate initial distribution of columns or latitudes and 
! distribution of processes by SMP
!
      nx_smp(0:npes-1) = 0
      nproc_smp(0:npes-1) = 0
      do procid=0,npes-1
         smpid = proc_smp_map(procid)
         nproc_smp(smpid) = nproc_smp(smpid) + 1
      enddo
!
      do smpid=0,nsmps-1
         nx_smp(smpid)     = nx_base*nproc_smp(smpid)
         avgnx_proc(smpid) = float(nx_base)
      enddo
!
! Equi-distribute remaining columns or latitudes across SMPs
! without increasing per process imbalance beyond minimum
!
      workleft = tot_nx - npes*nx_base
      do while (workleft > 0)
!
! (a) Find minimun number of columns or latitudes assigned to an SMP
!
         minavgnx_proc = avgnx_proc(0)
         do smpid=1,nsmps-1
            if (minavgnx_proc > avgnx_proc(smpid)) then
               minavgnx_proc = avgnx_proc(smpid)
            endif
         enddo
!
! (b) Assign an additional column or latitude to processes with 
!     nx_base latitudes/columns in SMPs with the minimum 
!     average number of latitudes/columns
!
         do procid=npes/2-1,0,-1
            if (mirror) then
               procids = 2*procid
               procidn = procids + 1
            else
               procids = procid
               procidn = npes - procids - 1
            endif
!
            smpids = proc_smp_map(procids)
            smpidn = proc_smp_map(procidn)
            if ((nx_p(procids) .eq. nx_base)  .and. &
                ((avgnx_proc(smpids) .eq. minavgnx_proc) .or. &
                 (avgnx_proc(smpidn) .eq. minavgnx_proc)) .and. &
                (workleft > 0)) then
!
               nx_p(procids) = nx_p(procids) + 1
               nx_smp(smpids) = nx_smp(smpids) + 1
               avgnx_proc(smpids) = &
                  float(nx_smp(smpids))/float(nproc_smp(smpids))
!
               nx_p(procidn) = nx_p(procids)
               nx_smp(smpidn) = nx_smp(smpidn) + 1
               avgnx_proc(smpidn) = &
                  float(nx_smp(smpidn))/float(nproc_smp(smpidn))
!
               workleft = workleft - 2
            endif
         enddo
      end do
!
! Partition latitudes over processes, equidistributing either 
! a) columns, or
! b) latitudes
!
      if (dyn_equi_by_col) then
!
! Evaluate different latitude assignments
!
         min_max_ncols = tot_cols
         do i=0,10
            alpha = .05*i
            max_ncols = 0
!
            iend = 0
            ncol_curtot  = 0
            ncol_curgoal = 0
            do procid=0,npes/2-1
               if (mirror) then
                  procids = 2*procid
               else
                  procids = procid
               endif
               ncol_curgoal = ncol_curgoal + nx_p(procids)
               ncol = 0
!
               done = .false.
!
! Add latitudes until near column per process goal for current process
!
               do while ((.not. done) .and. &
                  (ncol_curtot < ncol_curgoal))
                  if (iend .ge. plat/2) then
                     write(6,*)'SPMDINIT_DYN: error in assigning latitudes to processes'
                     call endrun
                  endif
                  if (ncol_curtot + nlon(iend+1) .le. &
                      ncol_curgoal + alpha*nlon(iend+1)) then
                     iend = iend + 1
                     ncol = ncol + nlon(iend)
                     ncol_curtot = ncol_curtot + nlon(iend)
                  else
                     done = .true.
                  endif
               enddo
               if (ncol > max_ncols) max_ncols = ncol
!
            enddo
            if (max_ncols < min_max_ncols) then
               min_max_ncols = max_ncols
               opt_alpha = alpha
            endif
         enddo
!
! Determine latitude assignments when equidistributing columns
!
         iend = 0
         ncol_curtot = 0
         ncol_curgoal = 0
         do procid=0,npes/2-1
            if (mirror) then
               procids = 2*procid
               procidn = procids + 1
            else
               procids = procid
               procidn = npes - procids - 1
            endif
            ncol_curgoal = ncol_curgoal + nx_p(procids)
            ncol_p(procids) = 0
!
            cut(1,procids) = iend + 1
            cut(2,procids) = iend
            done = .false.
!
! Add latitudes until near column per process goal for current process
!
            do while ((.not. done) .and. &
               (ncol_curtot < ncol_curgoal))
               if (ncol_curtot + nlon(iend+1) .le. &
                  ncol_curgoal + opt_alpha*nlon(iend+1)) then
                  iend = iend + 1
                  cut(2,procids) = iend
                  ncol_p(procids) = ncol_p(procids) + nlon(iend)
                  ncol_curtot = ncol_curtot + nlon(iend)
                  nlat_p(procids) = nlat_p(procids) + 1
               else
                  done = .true.
               endif
            enddo
!
! Assign mirror latitudes
!
            cut(1,procidn) = plat - cut(2,procids) + 1
            cut(2,procidn) = plat - cut(1,procids) + 1
            ncol_p(procidn) = ncol_p(procids)
            nlat_p(procidn) = nlat_p(procids)
!
! Save local information
!
            if (iam == procids .or. iam == procidn) then
               beglat = cut(1,iam)
               endlat = cut(2,iam)
               numlats = nlat_p(iam)
               begirow = cut(1,procids)
               endirow = cut(2,procids)
            end if
!
         enddo
!
      else
!
! Determine latitude assignments when
! equidistributing latitudes
!
         iend = 0
         do procid=0,npes/2-1
            if (mirror) then
               procids = 2*procid
               procidn = procids + 1
            else
               procids = procid
               procidn = npes - procids - 1
            endif
!
            nlat_p(procids) = nx_p(procids)
            cut(1,procids) = iend + 1
            cut(2,procids) = iend + nlat_p(procids)
            iend = iend + nlat_p(procids)
!
            ncol_p(procids) = 0
            do lat=cut(1,procids),cut(2,procids)
               ncol_p(procids) = ncol_p(procids) + nlon(lat)
            enddo
!
! Assign mirror latitudes
!
            nlat_p(procidn) = nx_p(procidn)
            cut(1,procidn) = plat - cut(2,procids) + 1
            cut(2,procidn) = plat - cut(1,procids) + 1
!
            ncol_p(procidn) = 0
            do lat=cut(1,procidn),cut(2,procidn)
               ncol_p(procidn) = ncol_p(procidn) + nlon(lat)
            enddo
!
! Save local information
!
            if (iam == procids .or. iam == procidn) then
               beglat = cut(1,iam)
               endlat = cut(2,iam)
               numlats = nlat_p(iam)
               begirow = cut(1,procids)
               endirow = cut(2,procids)
            end if
!
         enddo
      endif
!
! Calculate maximum number of latitudes and columns assigned to a process
!
      maxlats = maxval(nlat_p)
      maxcols = maxval(ncol_p)
!
      do procid=0,npes-1
         if (masterproc) then
            write(6,*)'procid ',procid,' assigned ', &
                      cut(2,procid)-cut(1,procid)+1,' latitude values from', &
                      cut(1,procid),' through ',cut(2,procid),' containing', &
                      ncol_p(procid),' vertical columns'
         end if
!
! Determine which process is responsible for the defined latitudes
!
         do lat=cut(1,procid),cut(2,procid)
            proc(lat) = procid
         end do
!
! The extended regions are simply "numbnd" wider at each
! side. The extended region do not go beyond 1 and plat, though
!
         cutex(1,procid) = cut(1,procid) - numbnd
         cutex(2,procid) = cut(2,procid) + numbnd
         if (iam == procid) then
            beglatex = cutex(1,procid) + numbnd
            endlatex = cutex(2,procid) + numbnd
            numlatsex = endlatex - beglatex + 1
         end if
      end do
!
! Determine neighbor processes needed for boundary communication.  
! North first.
!
      neighn = 0
      neighn_minlat(:) = -1
      do procid=0,npes-1
         if (procid /= iam) then
            if ((cut(1,procid) > cut(2,iam)) .and. &
                (cut(1,procid) <= cut(2,iam)+numbnd)) then
               neighn_minlat(cut(1,procid)) = procid
               neighn = neighn + 1
            endif
         endif
      enddo
!
! Sort north processes by increasing latitude
!
      allocate (neighn_proc (neighn))
      neighn = 0
      do lat=1,plat
         if (neighn_minlat(lat) /= -1) then
            neighn = neighn + 1
            neighn_proc(neighn) = neighn_minlat(lat)
         endif
      enddo
!
! South next.
!
      neighs = 0
      neighs_maxlat(:) = -1
      do procid=0,npes-1
         if (procid /= iam) then
            if ((cut(2,procid) < cut(1,iam)) .and. &
                (cut(2,procid) >= cut(1,iam)-numbnd)) then
               neighs_maxlat(cut(2,procid)) = procid
               neighs = neighs + 1
            endif
         endif
      enddo
!
! Sort south processes by decreasing latitude
!
      allocate (neighs_proc (neighs))
      neighs = 0
      do lat=plat,1,-1
         if (neighs_maxlat(lat) /= -1) then
            neighs = neighs + 1
            neighs_proc(neighs) = neighs_maxlat(lat)
         endif
      enddo
!
      if (masterproc) then
         write(6,*)'-----------------------------------------'
         write(6,*)'Number of lats passed north & south = ',numbnd
         write(6,*)'Node  Partition  Extended Partition'
         write(6,*)'-----------------------------------------'
         do procid=0,npes-1
            write(6,200) procid,cut(1,procid),cut(2,procid) ,cutex(1,procid), cutex(2,procid)
200         format(i3,4x,i3,'-',i3,7x,i3,'-',i3)
         end do
      end if
!      write(6,*)'iam=',iam,'Number of south neighbors needed for bndry exchange = ',neighs
!      write(6,*)'iam=',iam,'Number of north neighbors needed for bndry exchange = ',neighn

      call decomp_wavenumbers ()
!
! Precompute swap partners and number of steps in realloc4 alltoall algorithm.
! First, determine number of swaps.
!
      realloc4_steps = 0
      do procj=1,ceil2(npes)-1
         procid = pair(npes,procj,iam)
         if (procid >= 0) then
           if ((numm(iam) > 0 .or. numm(procid) > 0)) then
             realloc4_steps = realloc4_steps + 1
           end if
         end if
      end do
!
! Second, determine swap partners.
!
      allocate( realloc4_proc(realloc4_steps) )
      allocate( realloc4_step(0:npes-1) )
      realloc4_step(:) = -1
      realloc4_steps = 0
      do procj=1,ceil2(npes)-1
         procid = pair(npes,procj,iam)
         if (procid >= 0) then
            if ((numm(iam) > 0 .or. numm(procid) > 0)) then
              realloc4_steps = realloc4_steps + 1
              realloc4_proc(realloc4_steps) = procid
              realloc4_step(procid) = realloc4_steps
            end if
         end if
      end do
!
! Precompute swap partners in realloc5/7 allgather algorithm.
      allocate( allgather_proc(npes-1) )
      allocate( allgather_step(0:npes-1) )
      allgather_step(:) = -1
      allgather_steps = 0
      do procj=1,ceil2(npes)-1
         procid = pair(npes,procj,iam)
         if (procid >= 0) then
            allgather_steps = allgather_steps + 1
            allgather_proc(allgather_steps) = procid
            allgather_step(procid) = allgather_steps
         end if
      end do
!
      return
   end subroutine spmdinit_dyn

!========================================================================

   subroutine factor (nitems, m2, m3, m5)
!----------------------------------------------------------------------- 
! 
! Purpose: Factor a given number into powers of 2,3,5
! 
! Method: Brute force application of "mod" function
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: nitems      ! Number to be factored into powers of 2,3,5
      integer, intent(out) :: m2,m3,m5   ! Powers of 2, 3, and 5 respectively
!
! Local workspace
!
      integer num                        ! current number to be factored
!
!-----------------------------------------------------------------------
!
      num = nitems
      m2 = 0
      m3 = 0
      m5 = 0
      
2     if (mod(num,2) == 0) then
         m2 = m2 + 1
         num = num/2
         goto 2
      end if
      
3     if (mod(num,3) == 0) then
         m3 = m3 + 1
         num = num/3
         goto 3
      end if
      
5     if (mod(num,5) == 0) then
         m5 = m5 + 1
         num = num/5
         goto 5
      end if
      
      if (num /= 1) then
         write(6,*) 'FACTOR: ',nitems,' has a prime factor other than 2, 3, or 5.  Aborting...'
         call endrun
      end if
      
      return
   end subroutine factor

!========================================================================

   subroutine decomp_wavenumbers
!----------------------------------------------------------------------- 
! 
! Purpose: partition the spectral work among the given number of processes
! 
! Method: Approximately equidistribute both the number of spectral 
!         coefficients and the number of wavenumbers assigned to each 
!         MPI task using a modified version of the mapping due to
!         Barros and Kauranne. 
! 
! Author: P. Worley, September 2002
! 
!-----------------------------------------------------------------------
      use pspect, only: pmmax
      use comspe, only: numm, maxm, locm, locrm, nlen, lpspt, lnstart
      use infnan, only: bigint
!
! Local workspace
!
      integer procid      ! process id
      integer m, lm       ! global and local fourier wavenumber indices
      integer mstride     ! Stride over wavenumbers used in decomposition
      integer begm1       ! Starting Fourier wavenumbers owned by an MPI task
      integer begm2       !  when using Barros & Kauranne decomposition
      integer speccount(0:npes-1)
                          ! number of spectral coefficients assigned to
                          ! each MPI task
!-----------------------------------------------------------------------
!
! determine upper bound on number of wavenumbers to be assigned to each 
! process
      if (mod(pmmax,npes) .eq. 0) then
         maxm = pmmax/npes
      else
         maxm = (pmmax/npes) + 1
      endif
      allocate ( locm(1:maxm, 0:npes-1) )
      allocate ( locrm(1:2*maxm, 0:npes-1) )
!
! assign wavenumbers to approximately equidistribute the number 
! of spectral coefficients assigned to each process
      mstride = 2*npes
      npessp = 0
      do procid = 0,npes-1
         numm(procid) = 0
         speccount(procid) = 0
         begm1 = procid + 1
         begm2 = mstride - procid
         do m=begm1,pmmax,mstride
            numm(procid) = numm(procid) + 1
            locm(numm(procid),procid) = m
            speccount(procid) = speccount(procid) + nlen(m)
         enddo
         do m=begm2,pmmax,mstride
            numm(procid) = numm(procid) + 1
            locm(numm(procid),procid) = m
            speccount(procid) = speccount(procid) + nlen(m)
         enddo
!
         if (numm(procid) .gt. 0) then
            npessp = npessp + 1
         endif
!
      enddo
!
      do procid = 0,npes-1
         if (masterproc) then
            write(6,*)'procid ',procid,' assigned ', speccount(procid), &
                      ' spectral coefficients and ', numm(procid), &
                      ' m values: ', (locm(lm,procid),lm=1,numm(procid))
         end if
         do lm=1,numm(procid)
            locrm(2*lm-1,procid) = 2*locm(lm,procid)-1
            locrm(2*lm  ,procid) = 2*locm(lm,procid)
         enddo
         do lm=numm(procid)+1,maxm
            locm(lm,procid) = bigint
            locrm(2*lm-1,procid) = bigint
            locrm(2*lm  ,procid) = bigint
         enddo
      enddo
!
! Calculate number of local spectral coefficients
      lpspt = 0
      do lm=1,numm(iam)
         lpspt = lpspt + nlen(locm(lm,iam))
      enddo
!
! Evaluate displacement info based on truncation params and
! wavenumber assignment
      allocate ( lnstart(1:maxm) )
      lnstart(1) = 0
      do lm=2,numm(iam)
         lnstart(lm) = lnstart(lm-1) + nlen(locm(lm-1,iam))
      enddo
!   
      return
   end subroutine decomp_wavenumbers

!========================================================================

  subroutine spmdbuf 
!----------------------------------------------------------------------- 
! 
! Purpose: allocate spmd pack buffers used in collective communications
! 
! Author: CCM Core Group
!
! Note: Call after phys_grid_init
! 
!-----------------------------------------------------------------------
     use error_messages, only: alloc_err
     use comspe,         only: nlen, maxm
     use constituents,   only: pcnst, ppcnst
!-----------------------------------------------------------------------
!
! Local workspace
!
     integer :: maxcount(5),m
     integer :: length,i,lm,istat1,istat2
     integer :: bsiz, glb_bsiz       ! buffer size (in bytes)
!
! realloc4a max: 8  2 plev*numm*numlats (e.g. tdyn)
!                1  2     *numm*numlats (bpstr)
!
     maxcount(1) = (npes-1)*maxlats*(2*maxm*(plev*8 + 1))
!
! realloc4b max: 8  2 plev*numm*numlats (e.g. vort)
!                4  2     *numm*numlats (e.g. dps)
!
     maxcount(2) = (npes-1)*maxlats*(2*maxm*(plev*8 + 4))
!
! realloc5 max: 6 numlats         (e.g. tmass)
!               5 numlats  *pcnst (e.g. hw1lat)
!               2 4*numlats*pcnst (e.g. hw2al)
!
     maxcount(3) = npes*maxlats*(6 + (5 + 2*4)*pcnst)
!
! realloc7 max: 3 plev *numlats    (e.g. vmax2d)
!               5      *numlats    (e.g. psurf)
!
     maxcount(4) = npes*maxlats*(3*plev + 5)
!
! dp_coupling max:
!
     if (.not. local_dp_map) then
        maxcount(5) = (4 + ppcnst)*max(block_buf_nrecs,chunk_buf_nrecs)
     else
        maxcount(5) = 0
     endif
!
     m = maxval(maxcount)
     call mpipack_size (m, mpir8, mpicom, bsiz)
     call mpiallmaxint(bsiz, glb_bsiz, 1, mpicom)
     write(6,*) 'SPMDBUF: Allocating SPMD buffers of size ',glb_bsiz
     spmdbuf_siz = glb_bsiz/8 + 1
#if (defined CAF)
     allocate(buf1(spmdbuf_siz)[*], stat=istat1)
     allocate(buf2(spmdbuf_siz)[*], stat=istat2)
#else
     allocate(buf1(spmdbuf_siz), stat=istat1)
     allocate(buf2(spmdbuf_siz), stat=istat2)
#endif
     call alloc_err( istat1, 'spmdbuf', 'buf1', spmdbuf_siz )
     call alloc_err( istat2, 'spmdbuf', 'buf2', spmdbuf_siz )
     call mpiwincreate(buf1,spmdbuf_siz*8,mpicom,buf1win)
     call mpiwincreate(buf2,spmdbuf_siz*8,mpicom,buf2win)
     return
  end subroutine spmdbuf

!========================================================================

  subroutine compute_gsfactors (numperlat, numtot, numperproc, displs)
!----------------------------------------------------------------------- 
! 
! Purpose: Compute arguments for gatherv, scatterv
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!
! Input arguments
!
     integer, intent(in) :: numperlat    ! number of elements per latitude
!
! Output arguments
!
     integer, intent(out) :: numtot               ! total number of elements (to send or recv)
     integer, intent(out) :: numperproc(0:npes-1) ! per-PE number of items to receive
     integer, intent(out) :: displs(0:npes-1)     ! per-PE displacements
!
! Local variables
!
     integer :: p                    ! index
   
     numtot = numperlat*numlats
   
     do p=0,npes-1
        numperproc(p) = numperlat*nlat_p(p)
     end do
     
     displs(0) = 0
     do p=1,npes-1
        displs(p) = numperlat*(cut(1,p)-1)
     end do
     
  end subroutine compute_gsfactors

#endif

end module spmd_dyn
