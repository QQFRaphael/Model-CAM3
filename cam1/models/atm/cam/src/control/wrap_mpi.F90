#include <misc.h>
!----------------------------------------------------------------------- 
! 
! Purpose:
!
! 	Wrapper routines for the MPI (Message Passing) library for the
!	distributed memory (SPMD) version of the code. Also data with
!	"shorthand" names for the MPI data types.
!
! Entry points:
!      mpibarrier             Calls mpi_barrier
!      mpifinalize            Calls mpi_finalize
!      mpipack_size           Calls mpi_pack
!      mpipack                Calls mpi_pack
!      mpiunpack              Calls mpi_unpack
!      mpisendrecv            Calls mpi_sendrecv
!      mpiisend               Calls mpi_isend
!      mpiirsend              Calls mpi_irsend
!      mpiissend              Calls mpi_issend
!      mpiirecv               Calls mpi_irecv
!      mpiwait                Calls mpi_wait
!      mpiwaitall             Calls mpi_waitall
!      mpisend                Calls mpi_send
!      mpirsend               Calls mpi_rsend
!      mpissend               Calls mpi_ssend
!      mpirecv                Calls mpi_recv
!      mpigather              Calls mpi_gather
!      mpigatherv             Calls mpi_gatherv
!      mpisum                 Calls mpi_sum
!      mpiscatter             Calls mpi_scatter
!      mpiscatterv            Calls mpi_scatterv
!      mpibcast               Calls mpi_bcast
!      mpiallmaxint           Calls mpi_allreduce on integer vector with mpi_max operator
!      mpialltoallv           Calls mpi_alltoallv
!      mpialltoallint         Calls mpi_alltoall for integer data
!      mpiallgatherv          Calls mpi_allgatherv
!      mpiallgatherint        Calls mpi_allgatherv for integer data
!      altalltoallv           Calls one of a number of alternative
!                              implementations of alltoallv; currently 
!                              includes MPI-1, MPI-2 and
!                              Co-Array Fortran-based implementations
!      mpiwincreate           Calls mpi_win_create and mpi_win_fence
!
! Author: Many
! 
!-----------------------------------------------------------------------
!
! Compile these routines only when SPMD is defined
!
#if (defined SPMD)

!****************************************************************

   subroutine mpibarrier (comm)
!
! MPI barrier, have threads wait until all threads have reached this point
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   integer, intent(in):: comm
 
   integer ier   !MP error code
 
   call mpi_barrier (comm, ier)
   if (ier.ne.mpi_success) then
      write(6,*)'mpi_barrier failed ier=',ier
      call endrun
   end if
 
   return
   end subroutine mpibarrier
 
!****************************************************************
 
   subroutine mpifinalize
!
! End of all MPI communication
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   integer ier   !MP error code
 
   call mpi_finalize (ier)
   if (ier.ne.mpi_success) then
      write(6,*)'mpi_finalize failed ier=',ier
      call endrun
   end if
 
   return
   end subroutine mpifinalize
 
!****************************************************************
 
   subroutine mpipack_size (incount, datatype, comm, size)
!
! Returns the size of the packed data
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   integer, intent(in):: incount
   integer, intent(in):: datatype
   integer, intent(in):: comm
   integer, intent(out):: size
 
   integer ier   !MP error code
 
   call mpi_pack_size (incount, datatype, comm, size, ier)
   if (ier.ne.mpi_success) then
      write(6,*)'mpi_pack_size failed ier=',ier
      call endrun
   end if
 
   return
   end subroutine mpipack_size
 
!****************************************************************
 
   subroutine mpipack (inbuf, incount, datatype, outbuf, outsize,    &
                       position, comm)
!
! Pack the data and send it.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   real(r8), intent(in):: inbuf(*)
   real(r8), intent(out):: outbuf(*)
   integer, intent(in):: incount
   integer, intent(in):: datatype
   integer, intent(out):: outsize
   integer, intent(inout):: position
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
   call mpi_pack (inbuf, incount, datatype, outbuf, outsize,         &
                  position, comm, ier)
   if (ier.ne.mpi_success) then
      write(6,*)'mpi_pack failed ier=',ier
      call endrun
   end if
 
   return
   end subroutine mpipack
 
!****************************************************************
 
   subroutine mpiunpack (inbuf, insize, position, outbuf, outcount,  &
                         datatype, comm)
!
! Un-packs the data from the packed receive buffer
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   real(r8), intent(in):: inbuf(*)
   real(r8), intent(out):: outbuf(*)
   integer, intent(in):: insize
   integer, intent(inout):: position
   integer, intent(in):: outcount
   integer, intent(in):: datatype
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
   call mpi_unpack (inbuf, insize, position, outbuf, outcount,       &
                    datatype, comm, ier)
   if (ier.ne.mpi_success) then
      write(6,*)'mpi_unpack failed ier=',ier
      call endrun
   end if
 
   return
   end subroutine mpiunpack
 
!****************************************************************
 
   subroutine mpisendrecv (sendbuf, sendcount, sendtype, dest, sendtag,  &
                           recvbuf, recvcount, recvtype, source,recvtag, &
                           comm)
!
! Blocking send and receive.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   real(r8), intent(in):: sendbuf(*)
   real(r8), intent(out):: recvbuf(*)
   integer, intent(in):: sendcount
   integer, intent(in):: sendtype
   integer, intent(in):: dest
   integer, intent(in):: sendtag
   integer, intent(in):: recvcount
   integer, intent(in):: recvtype
   integer, intent(in):: source
   integer, intent(in):: recvtag
   integer, intent(in):: comm
 
   integer :: status(MPI_STATUS_SIZE)
   integer ier   !MP error code
 
   call t_startf ('mpi_sendrecv')
   call mpi_sendrecv (sendbuf, sendcount, sendtype, dest, sendtag,   &
                      recvbuf, recvcount, recvtype, source, recvtag, &
                      comm, status, ier)
   if (ier.ne.mpi_success) then
      write(6,*)'mpi_sendrecv failed ier=',ier
      call endrun
   end if
!
! ASSUME nrecv = nsend for stats gathering purposes.  This is not actually
! correct, but its the best we can do since recvcount is a Max number
!
   nsend = nsend + 1
   nrecv = nrecv + 1
   nwsend = nwsend + sendcount
   nwrecv = nwrecv + sendcount

   call t_stopf ('mpi_sendrecv')
 
   return
   end subroutine mpisendrecv
 
!****************************************************************
 
   subroutine mpiisend (buf, count, datatype, dest, tag, comm, request)
!
! Does a non-blocking send.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   real (r8), intent(in):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: dest
   integer, intent(in):: tag
   integer, intent(in):: comm
   integer, intent(out):: request
 
   integer ier   !MP error code
 
   call t_startf ('mpi_isend')
   call mpi_isend (buf, count, datatype, dest, tag, comm, request, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_isend failed ier=',ier
      call endrun
   end if
   nsend = nsend + 1
   nwsend = nwsend + count
   call t_stopf ('mpi_isend')
 
   return
   end subroutine mpiisend
 
!****************************************************************
 
   subroutine mpiirsend (buf, count, datatype, dest, tag, comm, request)
!
! Does a non-blocking ready send.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   real (r8), intent(in):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: dest
   integer, intent(in):: tag
   integer, intent(in):: comm
   integer, intent(out):: request
 
   integer ier   !MP error code
 
   call t_startf ('mpi_irsend')
   call mpi_irsend (buf, count, datatype, dest, tag, comm, request, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_irsend failed ier=',ier
      call endrun
   end if
   nsend = nsend + 1
   nwsend = nwsend + count
   call t_stopf ('mpi_irsend')
 
   return
   end subroutine mpiirsend
 
!****************************************************************
 
   subroutine mpiissend (buf, count, datatype, dest, tag, comm, request)
!
! Does a non-blocking synchronous send.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   real (r8), intent(in):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: dest
   integer, intent(in):: tag
   integer, intent(in):: comm
   integer, intent(out):: request
 
   integer ier   !MP error code
 
   call t_startf ('mpi_issend')
   call mpi_issend (buf, count, datatype, dest, tag, comm, request, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_issend failed ier=',ier
      call endrun
   end if
   nsend = nsend + 1
   nwsend = nwsend + count
   call t_stopf ('mpi_issend')
 
   return
   end subroutine mpiissend
 
!****************************************************************
 
   subroutine mpiirecv (buf, count, datatype, source, tag, comm, request)
!
! Does a non-blocking receive.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   real (r8), intent(out):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: source
   integer, intent(in):: tag
   integer, intent(in):: comm
   integer, intent(out):: request
 
   integer ier   !MP error code
 
   call t_startf ('mpi_irecv')
   call mpi_irecv (buf, count, datatype, source, tag, comm, request, ier )
   if (ier/=mpi_success) then
      write(6,*)'mpi_irecv failed ier=',ier
      call endrun
   end if
   nrecv = nrecv + 1
   nwrecv = nwrecv + count
   call t_stopf ('mpi_irecv')
 
   return
   end subroutine mpiirecv
 
!****************************************************************
 
   subroutine mpiwait (request, status)
!
! Waits for a nonblocking operation to complete.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   integer, intent(inout):: request
   integer, intent(out):: status
 
   integer ier   !MP error code
 
   call t_startf ('mpi_wait')
   call mpi_wait (request, status, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_wait failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_wait')
 
   return
   end subroutine mpiwait
 
!****************************************************************
 
   subroutine mpiwaitall (count, array_of_requests, array_of_statuses)
!
! Waits for a collection of nonblocking operations to complete.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   integer, intent(in):: count
   integer, intent(inout):: array_of_requests(*)
   integer, intent(out):: array_of_statuses(*)
 
   integer ier   !MP error code
 
   call t_startf ('mpi_waitall')
   call mpi_waitall (count, array_of_requests, array_of_statuses, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_waitall failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_waitall')
 
   return
   end subroutine mpiwaitall
 
!****************************************************************
 
   subroutine mpisend (buf, count, datatype, dest, tag, comm)
!
! Does a blocking send
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   real (r8), intent(in):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: dest
   integer, intent(in):: tag
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
   call t_startf ('mpi_send')
   call mpi_send (buf, count, datatype, dest, tag, comm, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_send failed ier=',ier
      call endrun
   end if
   nsend = nsend + 1
   nwsend = nwsend + count
   call t_stopf ('mpi_send')
 
   return
   end subroutine mpisend
 
!****************************************************************
 
   subroutine mpirsend (buf, count, datatype, dest, tag, comm)
!
! Does a blocking ready send
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   real (r8), intent(in):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: dest
   integer, intent(in):: tag
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
   call t_startf ('mpi_rsend')
   call mpi_rsend (buf, count, datatype, dest, tag, comm, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_rsend failed ier=',ier
      call endrun
   end if
   nsend = nsend + 1
   nwsend = nwsend + count
   call t_stopf ('mpi_rsend')
 
   return
   end subroutine mpirsend
 
!****************************************************************
 
   subroutine mpissend (buf, count, datatype, dest, tag, comm)
!
! Does a blocking synchronous send
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   real (r8), intent(in):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: dest
   integer, intent(in):: tag
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
   call t_startf ('mpi_ssend')
   call mpi_ssend (buf, count, datatype, dest, tag, comm, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_ssend failed ier=',ier
      call endrun
   end if
   nsend = nsend + 1
   nwsend = nwsend + count
   call t_stopf ('mpi_ssend')
 
   return
   end subroutine mpissend
 
!****************************************************************
 
   subroutine mpirecv (buf, count, datatype, source, tag, comm)
!
! Does a blocking receive
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   real (r8), intent(out):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: source
   integer, intent(in):: tag
   integer, intent(in):: comm
 
   integer status (MPI_STATUS_SIZE) ! Status of message
   integer ier   !MP error code
 
   call t_startf ('mpi_recv')
   call mpi_recv (buf, count, datatype, source, tag, comm, status, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_recv failed ier=',ier
      call endrun
   end if
   nrecv = nrecv + 1
   nwrecv = nwrecv + count
   call t_stopf ('mpi_recv')
 
   return
   end subroutine mpirecv
 
!****************************************************************
 
   subroutine mpigather (sendbuf, sendcnt, sendtype, recvbuf, recvcnt, &
                         recvtype, root, comm)
!
! Collects different messages from each thread on masterproc
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   real (r8), intent(in):: sendbuf(*)
   real (r8), intent(out):: recvbuf(*)
   integer, intent(in):: sendcnt
   integer, intent(in):: sendtype
   integer, intent(in):: recvcnt
   integer, intent(in):: recvtype
   integer, intent(in):: root
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
   call t_startf ('mpi_gather')
   call mpi_gather (sendbuf, sendcnt, sendtype,                      &
                    recvbuf, recvcnt, recvtype, root, comm, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_gather failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_gather')
 
   return
   end subroutine mpigather
 
!****************************************************************
 
   subroutine mpigatherv (sendbuf, sendcnt, sendtype, recvbuf, recvcnts, &
                          displs, recvtype, root, comm)
!
! Collects different messages from each thread on masterproc
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   real (r8), intent(in)  :: sendbuf(*)
   real (r8), intent(out) :: recvbuf(*)
   integer, intent(in) :: displs(*)
   integer, intent(in) :: sendcnt
   integer, intent(in) :: sendtype
   integer, intent(in) :: recvcnts(*)
   integer, intent(in) :: recvtype
   integer, intent(in) :: root
   integer, intent(in) :: comm
 
   integer ier   ! MPI error code
 
   call t_startf ('mpi_gather')
   call mpi_gatherv (sendbuf, sendcnt, sendtype, recvbuf, recvcnts, displs, recvtype, &
                     root, comm, ier)
   if (ier /= mpi_success) then
      write(6,*)'mpi_gather failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_gather')
 
   return
   end subroutine mpigatherv
 
!****************************************************************
 
   subroutine mpisum (sendbuf, recvbuf, cnt, datatype, root, comm)
!
! Sums sendbuf across all processors on communicator, returning 
! result to root.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   real (r8), intent(in):: sendbuf(*)
   real (r8), intent(out):: recvbuf(*)
   integer, intent(in):: cnt
   integer, intent(in):: datatype
   integer, intent(in):: root
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
   call t_startf ('mpi_reduce')
   call mpi_reduce (sendbuf, recvbuf, cnt, datatype, mpi_sum, &
                    root, comm, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_reduce failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_reduce')
 
   return
   end subroutine mpisum
 
!****************************************************************
 
   subroutine mpiscatter (sendbuf, sendcnt, sendtype, recvbuf, recvcnt, &
                          recvtype, root, comm)
!
! Sends different messages from masterproc to each thread
! 
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   real (r8),intent(in):: sendbuf(*)
   real (r8), intent(out):: recvbuf(*)
   integer,intent(in):: sendcnt
   integer,intent(in):: sendtype
   integer,intent(in):: recvcnt
   integer,intent(in):: recvtype
   integer,intent(in):: root
   integer,intent(in):: comm
 
   integer ier   !MP error code
 
   call t_startf ('mpi_scatter')
   call mpi_scatter (sendbuf, sendcnt, sendtype, recvbuf, recvcnt, &
                     recvtype, root, comm, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_scatter failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_scatter')
 
   return
   end subroutine mpiscatter
 
!****************************************************************
 
   subroutine mpiscatterv (sendbuf, sendcnts, displs, sendtype, recvbuf, &
                           recvcnt, recvtype, root, comm)
!
! Sends different messages from masterproc to each thread
! 
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   real (r8), intent(in)  :: sendbuf(*)
   real (r8), intent(out) :: recvbuf(*)
   integer, intent(in) :: displs(*)
   integer, intent(in) :: sendcnts(*)
   integer, intent(in) :: sendtype
   integer, intent(in) :: recvcnt
   integer, intent(in) :: recvtype
   integer, intent(in) :: root
   integer, intent(in) :: comm
 
   integer ier   !MP error code
 
   call t_startf ('mpi_scatter')
   call mpi_scatterv (sendbuf, sendcnts, displs, sendtype, recvbuf, recvcnt, &
                      recvtype, root, comm, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_scatter failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_scatter')
 
   return
   end subroutine mpiscatterv
 
!****************************************************************
 
   subroutine mpibcast (buffer, count, datatype, root, comm )
!
! Broadcasts a message from masterproc to all threads
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   real (r8), intent(inout):: buffer(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: root
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
   call t_startf ('mpi_bcast')
   call mpi_bcast (buffer, count, datatype, root, comm, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_bcast failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_bcast')
 
   return
   end subroutine mpibcast
!****************************************************************
 
   subroutine mpiallmaxint (sendbuf, recvbuf, count, comm)
!
! Allreduce integer vector maximum
! 
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   integer, intent(in)  :: sendbuf(*)
   integer, intent(out) :: recvbuf(*)
   integer, intent(in)  :: count
   integer, intent(in)  :: comm
 
   integer :: ier              ! MPI error code

   call t_startf ('mpi_allreduce')
   call mpi_allreduce (sendbuf, recvbuf, count, mpiint, &
                       mpimax, comm, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_allreduce failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_allreduce')

   return
   end subroutine mpiallmaxint

!****************************************************************
 
   subroutine mpialltoallv (sendbuf, sendcnts, sdispls, sendtype, &
                            recvbuf, recvcnts, rdispls, recvtype, &
                            comm)
!
! All-to-all scatter/gather
! 
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   real (r8), intent(in)  :: sendbuf(*)
   real (r8), intent(out) :: recvbuf(*)
   integer, intent(in) :: sdispls(*)
   integer, intent(in) :: sendcnts(*)
   integer, intent(in) :: sendtype
   integer, intent(in) :: recvcnts(*)
   integer, intent(in) :: rdispls(*)
   integer, intent(in) :: recvtype
   integer, intent(in) :: comm
 
   integer :: ier              ! MPI error code

   call t_startf ('mpi_alltoallv')
   call mpi_alltoallv (sendbuf, sendcnts, sdispls, sendtype, &
                       recvbuf, recvcnts, rdispls, recvtype, &
                       comm, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_alltoallv failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_alltoallv')

   return
   end subroutine mpialltoallv
!****************************************************************
 
   subroutine mpialltoallint (sendbuf, sendcnt, recvbuf, recvcnt, &
                              comm)
!
! All-to-all scatter/gather
! 
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   integer, intent(in)  :: sendbuf(*)
   integer, intent(in)  :: sendcnt
   integer, intent(out) :: recvbuf(*)
   integer, intent(in)  :: recvcnt
   integer, intent(in)  :: comm
 
   integer :: ier              ! MPI error code

   call t_startf ('mpi_alltoallint')
   call mpi_alltoall (sendbuf, sendcnt, mpiint, &
                      recvbuf, recvcnt, mpiint, &
                      comm, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_alltoallint failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_alltoallint')

   return
   end subroutine mpialltoallint

!****************************************************************
 
   subroutine mpiallgatherv (sendbuf, sendcnt, sendtype, &
                             recvbuf, recvcnts, rdispls, recvtype, &
                             comm)
!
! Collect data from each task and broadcast resulting
! vector to all tasks
! 
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   real (r8), intent(in)  :: sendbuf(*)
   real (r8), intent(out) :: recvbuf(*)
   integer, intent(in) :: sendcnt
   integer, intent(in) :: sendtype
   integer, intent(in) :: recvcnts(*)
   integer, intent(in) :: rdispls(*)
   integer, intent(in) :: recvtype
   integer, intent(in) :: comm
 
   integer ier   !MP error code
 
   call t_startf ('mpi_allgatherv')
   call mpi_allgatherv (sendbuf, sendcnt, sendtype, &
                        recvbuf, recvcnts, rdispls, recvtype, &
                        comm, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_allgatherv failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_allgatherv')
 
   return
   end subroutine mpiallgatherv
!****************************************************************
 
   subroutine mpiallgatherint (sendbuf, scount, recvbuf, rcount, comm)
!
! Collects integer data from each task and broadcasts resulting
! vector to all tasks
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   integer, intent(in)  :: sendbuf(*)
   integer, intent(out) :: recvbuf(*)
   integer, intent(in)  :: scount
   integer, intent(in)  :: rcount
   integer, intent(in)  :: comm
 
   integer ier   !MP error code

   call t_startf ('mpi_allgather')
   call mpi_allgather (sendbuf, scount, mpiint, recvbuf, rcount, &
                       mpiint, comm, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_allgather failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_allgather')
 
   return
   end subroutine mpiallgatherint

!****************************************************************
   subroutine altalltoallv (option, mytid, nprocs, steps, dests, &
                 sendbuf, sbuf_siz, sendcnts, sdispls, sendtype, &
                 recvbuf, rbuf_siz, recvcnts, rdispls, recvtype, &
                 msgtag, pdispls, desttype, recvwin, comm)
!
! All-to-all scatter/gather implemented using Co-Array
! Fortran one-sided commands, MPI-2 one sided commands,
! SWAP module MPI-1 commands, or MPI_SENDRECV.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils,   only: endrun
   use pmgrid
   use swap_comm,    only: swapm

   implicit none

   integer, intent(in) :: option               ! 0: sendrecv
                                               ! 1: swap package
                                               ! 2: mpi2 
                                               ! 3: co-array fortran
   integer, intent(in) :: mytid
   integer, intent(in) :: nprocs
   integer, intent(in) :: steps
   integer, intent(in) :: dests(steps)
   integer, intent(in) :: sbuf_siz
   integer, intent(in) :: sendcnts(0:nprocs-1)
   integer, intent(in) :: sdispls(0:nprocs-1)
   integer, intent(in) :: sendtype
   integer, intent(in) :: rbuf_siz
   integer, intent(in) :: recvcnts(0:nprocs-1)
   integer, intent(in) :: rdispls(0:nprocs-1)
   integer, intent(in) :: recvtype
   integer, intent(in) :: msgtag
   integer, intent(in) :: pdispls(0:nprocs-1)   ! displacement at 
                                                !  destination
   integer, intent(in) :: desttype
   integer, intent(in) :: recvwin
   integer, intent(in) :: comm

#if (defined CAF)
   real (r8), intent(in)  :: sendbuf(sbuf_siz)[*]
   real (r8), intent(out) :: recvbuf(rbuf_siz)[*]

   integer :: istart, iend, jstart, jend
#else
   real (r8), intent(in)  :: sendbuf(sbuf_siz)
   real (r8), intent(out) :: recvbuf(rbuf_siz)
#endif

   integer :: loption          ! local copy of option
   integer :: dest             ! MPI remote process id
   integer :: ier              ! MPI error code
   integer :: i                ! loop index
   integer :: sndids(steps)    ! nonblocking MPI send request ids
   integer :: rcvids(steps)    ! nonblocking MPI recv request ids
   integer :: status(MPI_STATUS_SIZE)
#if ( defined MPI2)
   integer(kind=MPI_ADDRESS_KIND) :: ddispls
#endif

!-----------------------------------------------------------------------
   loption = option
!
!  Co-Array Fortran implementation of alltoallv
!
   if (loption .eq. 3) then
#if ( defined CAF )
      call t_startf ('caf_alltoallv')
      if (this_image() .ne. (mytid+1)) then
         call endrun('altalltoallv (caf_alltoallv) failed: MPI id .ne. CAF id')
      endif

      call sync_images()

!DIR$ CONCURRENT
      do i = 1, steps
         dest = dests(i)
         if (sendcnts(dest) > 0) then
            istart = sdispls(dest)+1
            iend   = istart+sendcnts(dest)-1
            jstart = pdispls(dest)+1
            jend   = jstart+sendcnts(dest)-1
            recvbuf(jstart:jend)[dest+1] = sendbuf(istart:iend)
         end if
      end do

      call sync_images()
      call t_stopf ('caf_alltoallv')
#else
      loption = 0
#endif
!
!  MPI-2 one-sided implementation of alltoallv
!
   elseif (loption .eq. 2) then
#ifdef MPI2
      call t_startf ('mpi2_alltoallv')
      call mpi_win_fence(0,recvwin,ier)
      do i=1, steps
         dest = dests(i)
         if (sendcnts(dest) > 0) then
            ddispls = pdispls(dest)
            call mpi_put(sendbuf(sdispls(dest)+1), sendcnts(dest), sendtype, &
                         dest, ddispls, sendcnts(dest), desttype, &
                         recvwin, ier)
         endif
      end do
!
! wait for completion
      call mpi_win_fence(0,recvwin,ier)
      if (ier/=mpi_success) then
         write(6,*)'altalltoallv (mpi2_alltoallv) failed ier=',ier
         call endrun
      end if
      call t_stopf ('mpi2_alltoallv')
#else
      loption = 0
#endif
!
!  MPI-1 two-sided implementation of alltoallv
!  using SWAP routines
!
   elseif (loption .eq. 1) then
      call t_startf ('swap_alltoallv')
!
      call swapm(steps, msgtag, dests, sendcnts, sdispls, &
                 sendbuf, sbuf_siz, recvcnts, rdispls, &
                 recvbuf, rbuf_siz)
!
      call t_stopf ('swap_alltoallv')
!
!  Anything else defined to be MPI_SENDRECV implementation
!
   else
!
      loption = 0
!
   endif
!
!  MPI_SENDRECV implementation of alltoallv
!
   if (loption .eq. 0) then
      call t_startf ('mpi1_alltoallv')
      do i=1, steps
         dest = dests(i)
         call mpi_sendrecv (sendbuf(sdispls(dest)+1), sendcnts(dest), &
                            sendtype, dest, msgtag,                   &
                            recvbuf(rdispls(dest)+1), recvcnts(dest), &
                            recvtype, dest, msgtag,                   &
                            comm, status, ier)
      end do
!
! test for error
      if (ier/=mpi_success) then
         write(6,*)'altalltoallv (mpi1_alltoallv) failed ier=',ier
         call endrun
      end if
      call t_stopf ('mpi1_alltoallv')
   endif
!
!  Local copy (if necessary)
   if (sendcnts(mytid) > 0) then
      do i=1,sendcnts(iam)
         recvbuf(rdispls(mytid)+i) = sendbuf(sdispls(mytid)+i)
      enddo
   endif
!
   nsend = nsend + steps
   nrecv = nrecv + steps
   do i=1,steps
      dest = dests(i)
      nwsend = nwsend + sendcnts(dest)
      nwrecv = nwrecv + recvcnts(dest)
   enddo
!
   return
   end subroutine altalltoallv

!****************************************************************

   subroutine mpiwincreate(base,size,comm,win)
!
! Creates window for MPI2 one-sided commands
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use abortutils, only: endrun

   implicit none

   real(r8), intent(in)  :: base(*)
   integer,  intent(in)  :: size
   integer,  intent(in)  :: comm
   integer,  intent(out) :: win
!
#ifdef MPI2
   integer(kind=MPI_ADDRESS_KIND) :: size8
   integer :: ier, info
!
   call t_startf ('mpi_win_create')
   info = MPI_INFO_NULL
   size8 = size
   call mpi_win_create(base,size8,8,info,comm,win,ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_win_create failed ier=',ier
      call endrun
   end if
   call mpi_win_fence(0,win,ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_win_fence failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_win_create')
#endif

   return
   end subroutine mpiwincreate
!****************************************************************
!
! If SPMD is not turned on
!
#else
   subroutine wrap_mpi
   use abortutils, only: endrun
   implicit none
!
! A unused stub routine to make the compiler happy when SPMD is
! turned off (which means you don't need anything in this file).
!
   call endrun ('(WRAP_MPI): This should not be called at all')
   end subroutine wrap_mpi
#endif

