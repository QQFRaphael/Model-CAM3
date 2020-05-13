#include <misc.h>
#include <params.h>

module swap_comm

!----------------------------------------------------------------------- 
! 
! Purpose: swap communication routines used in performance portable 
!          distributed transpose algorithms.
! 
! Entry points:
!      swap_comm_init          Initialize swap module.
!
!      swap_comm_defaultopts   Get default runtime options.
!      swap_comm_setopts       Set runtime options.
!      
!      swapm                   Implementation of swap for multiple
!                              messages using MPI point-to-point routines.
!
! Author: P. Worley
!-----------------------------------------------------------------------

#if (defined SPMD)
!-----------------------------------------------------------------------
!- use statements ------------------------------------------------------
!-----------------------------------------------------------------------

   use abortutils, only: endrun
   use mpishorthand

!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
   private                   ! Make the default access private
   save

!-----------------------------------------------------------------------
! Public interfaces ----------------------------------------------------
!-----------------------------------------------------------------------
   public swap_comm_init
   public swap_comm_defaultopts 
   public swap_comm_setopts 
   public swapm

!-----------------------------------------------------------------------
! Private data ---------------------------------------------------------
!-----------------------------------------------------------------------
! Swap communication order option:
!  0: simple swap: send/recv
!  1: ordered swap: [send/recv]|[recv/send]
!  2: delayed-recv swap: send ... recv
   integer, private, parameter :: min_comm_order = 0
   integer, private, parameter :: max_comm_order = 2
   integer, private, parameter :: def_comm_order = 0              ! default
   integer, private :: swap_comm_order = def_comm_order

! Swap communication protocol option:
!  1, 3, 5, 7, 9:                  nonblocking send
!  2, 3, 4, 5, 8, 9:               nonblocking receive
!  4, 5:                           ready send
!  6 .and. swap_comm_order .eq. 0: sendrecv
!  6 .and. swap_comm_order .eq. 1: explicitly synchronous  
!  7, 8, 9, .or. 10:               synchronous send          
   integer, private, parameter :: min0_comm_protocol =  1
   integer, private, parameter :: max0_comm_protocol =  9
   integer, private, parameter :: min1_comm_protocol =  0
   integer, private, parameter :: max1_comm_protocol = 10
   integer, private, parameter :: def_comm_protocol  =  6        ! default
   integer, private :: swap_comm_protocol = def_comm_protocol

! Swap communication maximum request count:
!  <=0: do not limit number of outstanding send/receive requests
!   >0: do not allow more than swap_comm_maxreq outstanding
!       nonblocking send requests or nonblocking receive requests
   integer, private, parameter :: def_comm_maxreq = -1           ! default
   integer, private :: swap_comm_maxreq = def_comm_maxreq

! Swap communicators
   integer, private :: swap_com = mpi_comm_null
                                      ! primary MPI communicator
   integer, private :: handshake_com  = mpi_comm_null
                                      ! MPI communicator for 
                                      !  handshaking messages

!-----------------------------------------------------------------------
! Subroutines and functions --------------------------------------------
!-----------------------------------------------------------------------
contains

!
!========================================================================
!
   subroutine swap_comm_init()

!----------------------------------------------------------------------- 
! 
! Purpose: Create communicators to be used in swap communication.
! 
! Method: 
!
! Author:  P. Worley
!-----------------------------------------------------------------------
   use mpishorthand, only: mpicom
!-----------------------------------------------------------------------
   implicit none
!---------------------------Input arguments-----------------------------
!
!---------------------------Local workspace-----------------------------
!
   integer ier               ! return error status    
!
!-----------------------------------------------------------------------
!
   call mpi_comm_dup(mpicom, swap_com, ier)
   if (ier /= mpi_success) then
      write(6,*)                                         &
         'SWAP_COMM_INIT:  ERROR:  mpi_comm_dup failed with IER=', ier
      call endrun
   endif
   call mpi_comm_dup(mpicom, handshake_com, ier)
   if (ier /= mpi_success) then
      write(6,*)                                         &
         'SWAP_COMM_INIT:  ERROR:  mpi_comm_dup failed with IER=', ier
      call endrun
   endif
!
   return
   end subroutine swap_comm_init
!
!========================================================================
!
   subroutine swap_comm_defaultopts(swap_comm_order_out, &
                                    swap_comm_protocol_out, &
                                    swap_comm_maxreq_out  )
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
! Author: P. Worley (modelled after Tom Henderson's code)
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
   ! swap module communication order option
   integer, intent(out), optional :: swap_comm_order_out
   ! swap module communication protocol option
   integer, intent(out), optional :: swap_comm_protocol_out
   ! swap module communication nonblocking request maximum
   integer, intent(out), optional :: swap_comm_maxreq_out
!-----------------------------------------------------------------------
   if ( present(swap_comm_order_out) ) then
      swap_comm_order_out = def_comm_order
   endif
   if ( present(swap_comm_protocol_out) ) then
      swap_comm_protocol_out = def_comm_protocol
   endif
   if ( present(swap_comm_maxreq_out) ) then
      swap_comm_maxreq_out = def_comm_maxreq
   endif
!
   return
   end subroutine swap_comm_defaultopts
!
!========================================================================
!
   subroutine swap_comm_setopts(swap_comm_order_in, &
                                swap_comm_protocol_in, &
                                swap_comm_maxreq_in )
!----------------------------------------------------------------------- 
! Purpose: Set runtime options
! Author: P. Worley (modelled after Tom Henderson's code)
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
     ! swap module communication order option
     integer, intent(in), optional :: swap_comm_order_in
     ! swap module communication protocol option
     integer, intent(in), optional :: swap_comm_protocol_in
     ! swap module communication nonblocking request maximum
     integer, intent(in), optional :: swap_comm_maxreq_in
!-----------------------------------------------------------------------
     if ( present(swap_comm_order_in) ) then
        swap_comm_order = swap_comm_order_in
        if ((swap_comm_order < min_comm_order) .or. &
            (swap_comm_order > max_comm_order)) then
           write(6,*)                                         &
              'SWAP_COMM_SETOPTS:  ERROR:  swap_comm_order=', &
              swap_comm_order,                                &
              ' is out of range.  It must be between ',       &
              min_comm_order,' and ',max_comm_order
           call endrun
        endif
     endif
!
     if ( present(swap_comm_protocol_in) ) then
        swap_comm_protocol = swap_comm_protocol_in
        if ((swap_comm_order .eq. 0) .or. &
            (swap_comm_order .eq. 2)) then
           if ((swap_comm_protocol < min0_comm_protocol) .or. &
               (swap_comm_protocol > max0_comm_protocol)) then
              write(6,*)                                            &
                 'SWAP_COMM_SETOPTS:  ERROR:  swap_comm_protocol=', &
                 swap_comm_protocol,                                &
                 ' is out of range.  It must be between ',          &
                 min0_comm_protocol,' and ',max0_comm_protocol,     &
                 ' when swap_comm_order= ', swap_comm_order
              call endrun
           endif
        else
           if ((swap_comm_protocol < min1_comm_protocol) .or. &
               (swap_comm_protocol > max1_comm_protocol)) then
              write(6,*)                                            &
                 'SWAP_COMM_SETOPTS:  ERROR:  swap_comm_protocol=', &
                 swap_comm_protocol,                                &
                 ' is out of range.  It must be between ',          &
                 min1_comm_protocol,' and ',max1_comm_protocol,     &
                 ' when swap_comm_order= ', swap_comm_order
              call endrun
           endif
        endif
     endif
!
     if ( present(swap_comm_maxreq_in) ) then
        swap_comm_maxreq = swap_comm_maxreq_in
     endif
!
     return
   end subroutine swap_comm_setopts
!
!========================================================================
!
   subroutine swapm(cnt, mtag, swapnodes, sndlths, sdispls, &
                    sndbuf, sbuf_siz, rcvlths, rdispls, &
                    rcvbuf, rbuf_siz)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Implementation of swap for multiple message using MPI point-to-point
! routines. 
! 
! Method: 
! if (swap_comm_order .eq. 0) simple swap: send/recv
! if (swap_comm_order .eq. 1) ordered swap: [send/recv]|[recv/send]
! if (swap_comm_order .eq. 2) delayed-recv swap: send ... recv
! if (swap_comm_protocol .eq. 1, 3, 5, 7, .or. 9) nonblocking send
! if (swap_comm_protocol .eq. 2, 3, 4, 5, 8, .or. 9)  nonblocking receive
! if (swap_comm_protocol .eq. 4 .or. 5)  ready send
! if (swap_comm_protocol .eq. 6 .and. swap_comm_order .eq. 0) sendrecv
! if (swap_comm_protocol .eq. 6 .and. swap_comm_order .eq. 1) explicitly synchronous  
! if (swap_comm_protocol .eq. 7, 8, 9, .or. 10)  synchronous send          
!
! Author of original version:  P. Worley
! Ported to CAM: P. Worley, December 2003
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: iam
   use spmd_dyn, only: npes
!-----------------------------------------------------------------------
   implicit none
!---------------------------Input arguments--------------------------
!
   integer, intent(in)   :: cnt              ! number of swaps to initiate
   integer, intent(in)   :: mtag             ! MPI message tag
   integer, intent(in)   :: swapnodes(cnt)   ! MPI process id of swap partners
   integer, intent(in)   :: sndlths(0:npes-1)! length of outgoing message
   integer, intent(in)   :: sdispls(0:npes-1)! offset from beginning of send
                                             !  buffer where outgoing messages
                                             !  should be sent from
   integer, intent(in)   :: sbuf_siz         ! size of send buffer
   integer, intent(in)   :: rcvlths(0:npes-1)! length of incoming messages
   integer, intent(in)   :: rdispls(0:npes-1)! offset from beginning of receive 
                                             !  buffer where incoming messages
                                             !  should be placed
   integer, intent(in)   :: rbuf_siz         ! size of receive buffer
   real(r8), intent(in)  :: sndbuf(sbuf_siz) ! outgoing message buffer
   real(r8), intent(out) :: rcvbuf(rbuf_siz) ! incoming message buffer
!
!---------------------------Local workspace-----------------------------
!
   integer  i,jb,je                        ! loop indices
   integer  max_cnt                        ! maximum number of outstanding
                                           ! send/receive requests
   integer  lcnt                           ! current loop count
   integer  p                              ! process index
   integer  offset_s                       ! index of message beginning in 
                                           !  send buffer
   integer  offset_r                       ! index of message beginning in 
                                           !  receive buffer
   integer  sndids(cnt)                    ! send request ids
   integer  rcvids(cnt)                    ! receive request ids
   real(r8) signal                         ! ready send signal
   integer  ier                            ! return error status    
   integer  status(MPI_STATUS_SIZE,cnt)    ! MPI status integer
!
!-------------------------------------------------------------------------------------
!
! Determine how many outstanding send/receive requests to allow.
   if (swap_comm_maxreq < 1) then
      max_cnt = cnt
   else
      max_cnt = swap_comm_maxreq
   endif
!
   signal = 1.0
   ier = mpi_success
!
   do jb=1,cnt,max_cnt
      je = min(jb+max_cnt-1,cnt)
      lcnt = (je-jb+1)
!
      if (swap_comm_order .eq. 0) then
!
         if (swap_comm_protocol <= 5) then
!
            if (swap_comm_protocol <= 1) then
!
! do not block for the send
               do i=jb,je
                  p = swapnodes(i)
                  offset_s = sdispls(p)+1
                  call mpi_isend( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                  swap_com, sndids(i), ier )
               enddo
               do i=jb,je
                  p = swapnodes(i)
                  offset_r = rdispls(p)+1
                  call mpi_recv ( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag, &
                                  swap_com, status(1,i), ier )
               enddo
               call mpi_waitall ( lcnt, sndids(jb), status(1,jb), ier )
!
            elseif (swap_comm_protocol <= 3) then
!
! post the receive before the send, increasing odds that the
! receive will be posted before the message arrives.
               do i=jb,je
                  p = swapnodes(i)
                  offset_r = rdispls(p)+1
                  call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag,&
                                 swap_com, rcvids(i), ier )
               enddo
!
               if (swap_comm_protocol .eq. 2) then
! complete outstanding receives
                  do i=jb,je
                     p = swapnodes(i)
                     offset_s = sdispls(p)+1
                     call mpi_send( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                    swap_com, ier )
                  enddo
                  call mpi_waitall( lcnt, rcvids(jb), status(1,jb), ier )
               else
! also do not block for the send
                  do i=jb,je
                     p = swapnodes(i)
                     offset_s = sdispls(p)+1
                     call mpi_isend( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                     swap_com, sndids(i), ier )
                  enddo
                  call mpi_waitall ( lcnt, rcvids(jb), status(1,jb), ier )
                  call mpi_waitall ( lcnt, sndids(jb), status(1,jb), ier )
               endif
!
            else
!
! post the receive before send to allow use of ready send.
               do i=jb,je
                  p = swapnodes(i)
                  offset_r = rdispls(p)+1
                  call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag,&
                                  swap_com, rcvids(i), ier )
                  call mpi_send( signal, 1, mpir8, p, mtag, handshake_com, &
                                 ier )
               enddo
!    
               if (swap_comm_protocol .eq. 4) then
! complete receive of ready send
                  do i=jb,je
                     p = swapnodes(i)
                     offset_s = sdispls(p)+1
                     call mpi_recv ( signal, 1, mpir8, p, mtag, handshake_com, &
                                     status(1,i), ier )
                     call mpi_rsend( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                     swap_com, ier )
                  enddo
                  call mpi_waitall ( lcnt, rcvids(jb), status(1,jb), ier )
               else
! also do not block for send, enabling overlap of communication with computation.
                  do i=jb,je
                     p = swapnodes(i)
                     offset_s = sdispls(p)+1
                     call mpi_recv  ( signal, 1, mpir8, p, mtag, handshake_com, &
                                      status(1,i), ier )
                     call mpi_irsend( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                      swap_com, sndids(i), ier )
                  enddo
                  call mpi_waitall ( lcnt, rcvids(jb), status(1,jb), ier )
                  call mpi_waitall ( lcnt, sndids(jb), status(1,jb), ier )
               endif
!
            endif
!
         elseif (swap_comm_protocol <= 9) then
!
            if (swap_comm_protocol <= 7) then
!
               if (swap_comm_protocol .eq. 6) then
! native sendrecv
                  do i=jb,je
                     p = swapnodes(i)
                     offset_s = sdispls(p)+1
                     offset_r = rdispls(p)+1
                     call mpi_sendrecv( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                        rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag, &
                                        swap_com, status(1,i), ier )
                  enddo
!
               else
! do not block for the synchronous send, enabling overlap of 
! communication with computation.
                  do i=jb,je
                     p = swapnodes(i)
                     offset_s = sdispls(p)+1
                     call mpi_issend( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                      swap_com, sndids(i), ier )
                  enddo
                  do i=jb,je
                     p = swapnodes(i)
                     offset_r = rdispls(p)+1
                     call mpi_recv ( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag, &
                                     swap_com, status(1,i), ier )
                  enddo
                  call mpi_waitall ( lcnt, sndids(jb), status(1,jb), ier )
               endif
!
            else
!
! post the receive before the synchronous send, increasing odds that the
! receive will be posted before the message arrives.
               do i=jb,je
                  p = swapnodes(i)
                  offset_r = rdispls(p)+1
                  call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag,&
                                  swap_com, rcvids(i), ier )
               enddo
!
               if (swap_comm_protocol .eq. 8) then
! complete outstanding receive,
                  do i=jb,je
                     p = swapnodes(i)
                     offset_s = sdispls(p)+1
                     call mpi_ssend( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                     swap_com, ier )
                  enddo
                  call mpi_waitall ( lcnt, rcvids(jb), status(1,jb), ier )
               else
! also do not block for the synchronous send, enabling overlap of
! communication with computation
                  do i=jb,je
                     p = swapnodes(i)
                     offset_s = sdispls(p)+1
                     call mpi_issend( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                      swap_com, sndids(i), ier )
                  enddo
                  call mpi_waitall ( lcnt, rcvids(jb), status(1,jb), ier )
                  call mpi_waitall ( lcnt, sndids(jb), status(1,jb), ier )
               endif
!
            endif
!
         else
!
            write (0,901) swap_comm_order, swap_comm_protocol
  901      format(/,' fatal error in subroutine swapm:', &
                  /,' unknown communication protocol specified',/,  &
                    ' swap_comm_order = ',i6, ' swap_comm_protocol = ',i6)
           call endrun
!
         endif
!
      elseif (swap_comm_order .eq. 1) then
! ordered swap:
! if (iam <= swapnode) send/recv
! if (iam >= swapnode) recv/send
!
         if (swap_comm_protocol <= 5) then
!
            if (swap_comm_protocol <= 1) then
!
               if (swap_comm_protocol .eq. 0) then
!
                  do i=jb,je
                     p = swapnodes(i)
                     offset_s = sdispls(p)+1
                     offset_r = rdispls(p)+1
                     if (iam <= p) then
                        call mpi_send ( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                        swap_com, ier )
                        call mpi_recv ( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag, &
                                        swap_com, status(1,i), ier )
                     else
                        call mpi_recv ( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag, &
                                        swap_com, status(1,i), ier )
                        call mpi_send ( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                        swap_com, ier )
                     endif
                  enddo
!
               else
!
! do not block for the send, enabling overlap of communication with computation.
                  do i=jb,je
                     p = swapnodes(i)
                     offset_s = sdispls(p)+1
                     offset_r = rdispls(p)+1
                     if (iam <= p) then
                        call mpi_isend( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                        swap_com, sndids(i), ier )
                        call mpi_recv ( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag, &
                                        swap_com, status(1,i), ier )
                     else
                        call mpi_recv ( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag, &
                                        swap_com, status(1,i), ier )
                        call mpi_isend( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                        swap_com, sndids(i), ier )
                     endif
                  enddo
                  call mpi_waitall ( lcnt, sndids(jb), status(1,jb), ier )
!
               endif
!
            elseif (swap_comm_protocol <= 3) then
!
! post the receive before the initial send, increasing odds 
! that the receive will be posted before the message arrives.
               do i=jb,je
                  p = swapnodes(i)
                  offset_r = rdispls(p)+1
                  call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag,&
                                  swap_com, rcvids(i), ier )
               enddo
!
               if (swap_comm_protocol .eq. 2) then
!
! complete outstanding receive.
                  do i=jb,je
                     p = swapnodes(i)
                     offset_s = sdispls(p)+1
                     if (iam <= p) then
                        call mpi_send ( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                        swap_com, ier )
                        call mpi_wait( rcvids(i), status(1,i), ier )
                     else
                        call mpi_wait( rcvids(i), status(1,i), ier )
                        call mpi_send ( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                        swap_com, ier )
                     endif
                  enddo
!
               else
!
! also do not block for the send, enabling overlap of communication with computation.
                  do i=jb,je
                     p = swapnodes(i)
                     offset_s = sdispls(p)+1
                     if (iam <= p) then
                        call mpi_isend( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                        swap_com, sndids(i), ier )
                        call mpi_wait( rcvids(i), status(1,i), ier )
                     else
                        call mpi_wait( rcvids(i), status(1,i), ier )
                        call mpi_isend( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                        swap_com, sndids(i), ier )
                     endif
                  enddo
                  call mpi_waitall ( lcnt, sndids(jb), status(1,jb), ier )
    
               endif
!
            else
!
! post the receive before the send to allow use of forcetypes. 
               do i=jb,je
                  p = swapnodes(i)
                  offset_r = rdispls(p)+1
                  if (iam <= p) then
                     call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag,&
                                     swap_com, rcvids(i), ier )
                  else
                     call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag,&
                                     swap_com, rcvids(i), ier )
                     call mpi_send( signal, 1, mpir8, p, mtag, handshake_com, ier )
                  endif
               enddo
!
               if (swap_comm_protocol .eq. 4) then
!
! complete forcetype receive.
                  do i=jb,je
                     p = swapnodes(i)
                     offset_s = sdispls(p)+1
                     if (iam <= p) then
                        call mpi_recv ( signal, 1, mpir8, p, mtag, handshake_com, &
                                        status(1,i), ier )
                        call mpi_rsend( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                        swap_com, ier )
                        call mpi_wait ( rcvids(i), status(1,i), ier )
                     else
                        call mpi_wait ( rcvids(i), status(1,i), ier )
                        call mpi_rsend( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                        swap_com, ier )
                     endif
                  enddo
!
               else
!
! also do not block for the send, enabling overlap of communication with computation.
                  do i=jb,je
                     p = swapnodes(i)
                     offset_s = sdispls(p)+1
                     if (iam <= p) then
                        call mpi_recv ( signal, 1, mpir8, p, mtag, handshake_com, &
                                        status(1,i), ier )
                        call mpi_irsend( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                         swap_com, sndids(i), ier )
                        call mpi_wait ( rcvids(i), status(1,i), ier )
                     else
                        call mpi_wait ( rcvids(i), status(1,i), ier )
                        call mpi_irsend( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                         swap_com, sndids(i), ier )
                     endif
                  enddo
                  call mpi_waitall ( lcnt, sndids(jb), status(1,jb), ier )
!
               endif
!
            endif
!
         elseif (swap_comm_protocol <= 10) then
!
            if (swap_comm_protocol <= 7) then
!
               if (swap_comm_protocol .eq. 6) then
!
! synchronous ordered swap 
                  do i=jb,je
                     p = swapnodes(i)
                     offset_s = sdispls(p)+1
                     offset_r = rdispls(p)+1
                     if (iam <= p) then
                        call mpi_recv ( signal, 1, mpir8, p, mtag, handshake_com, &
                                        status(1,i), ier )
                        call mpi_send ( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                        swap_com, ier )
                        call mpi_recv ( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag, &
                                        swap_com, status(1,i), ier )
                     else
                        call mpi_send ( signal, 1, mpir8, p, mtag, handshake_com, ier )
                        call mpi_recv ( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag, &
                                        swap_com, status(1,i), ier )
                        call mpi_send ( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                        swap_com, ier )
                     endif
                  enddo
!
               else
!
! do not block for the synchronous send, enabling overlap of communication
! with computation.
                  do i=jb,je
                     p = swapnodes(i)
                     offset_s = sdispls(p)+1
                     offset_r = rdispls(p)+1
                     if (iam <= p) then
                        call mpi_issend( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                         swap_com, sndids(i), ier )
                        call mpi_recv ( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag, &
                                        swap_com, status(1,i), ier )
                     else
                        call mpi_recv ( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag, &
                                        swap_com, status(1,i), ier )
                        call mpi_issend( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                         swap_com, sndids(i), ier )
                     endif
                  enddo
                  call mpi_waitall ( lcnt, sndids(jb), status(1,jb), ier )
!
               endif
!
            elseif (swap_comm_protocol <= 9) then
!
! post the receive before the initial synchronous send, increasing odds 
! that the receive will be posted before the message arrives.
               do i=jb,je
                  p = swapnodes(i)
                  offset_r = rdispls(p)+1
                  call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag,&
                                 swap_com, rcvids(i), ier )
               enddo
!
               if (swap_comm_protocol .eq. 8) then
!
! complete outstanding receive.
                  do i=jb,je
                     p = swapnodes(i)
                     offset_s = sdispls(p)+1
                     if (iam <= p) then
                        call mpi_ssend( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                        swap_com, ier )
                        call mpi_wait( rcvids(i), status(1,i), ier )
                     else
                        call mpi_wait( rcvids(i), status(1,i), ier )
                        call mpi_ssend( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                        swap_com, ier )
                     endif
                  enddo
!
               else
!
! also do not block for the synchronous send, enabling overlap of communication
! with computation.
                  do i=jb,je
                     p = swapnodes(i)
                     offset_s = sdispls(p)+1
                     if (iam <= p) then
                        call mpi_issend( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                         swap_com, sndids(i), ier )
                        call mpi_wait( rcvids(i), status(1,i), ier )
                     else
                        call mpi_wait( rcvids(i), status(1,i), ier )
                        call mpi_issend( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                         swap_com, sndids(i), ier )
                     endif
                  enddo
                  call mpi_waitall ( lcnt, sndids(jb), status(1,jb), ier )
!
               endif
!
            else
! ordered swap using synchronous sends
               do i=jb,je
                  p = swapnodes(i)
                  offset_s = sdispls(p)+1
                  offset_r = rdispls(p)+1
                  if (iam <= p) then
                     call mpi_ssend ( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                      swap_com, ier )
                     call mpi_recv ( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag, &
                                     swap_com, status(1,i), ier )
                  else
                     call mpi_recv ( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag, &
                                     swap_com, status(1,i), ier )
                     call mpi_ssend ( sndbuf(offset_s), sndlths(p), mpir8, p, mtag, &
                                      swap_com, ier )
                  endif
               enddo
!
            endif
!
         else
!
! protocol error
            write (0,901) swap_comm_order, swap_comm_protocol
            call endrun
!
         endif
!
      else
! undefined swap option
!
          write (0,900) swap_comm_order
  900     format(/,' fatal error in subroutine swapm:', &
                 /,' unknown communication option specified',/, &
                   ' swap_comm_order = ',i6)                                 
          call endrun                                            
!
      endif
!
      if (ier /= mpi_success) then
         write(6,*)'MPI command in swapm failed, ier=',ier
         call endrun
      end if
!
   enddo
!
   return
   end subroutine swapm
!
#endif

end module swap_comm
