program fmain
!
! Purpose: Main program for reduced<->full grid interpolation.
!
! Method:  Parse cmdline, open files, set dimension sizes, and call main driver
!
   use mkrgrid, only: driver
   use globals, only: maxnlat, ncidi, ncido, londimid, latdimid, levdimid, &
                      ilevdimid, unlimdimid
   implicit none
!
! Local workspace
!
   include 'netcdf.inc'

   integer, parameter :: maxmonosiz = 1000  ! max size of monolist (cmdline args)
   character(len=16) monolist(maxmonosiz)   ! list of vars for monotonic cubic interp
   integer :: monosiz = 0                   ! initial size of monolist

   integer :: i, n               ! loop indices
   integer :: nargs              ! number of cmd line args
   integer :: plon = -1          ! longitude dimension length (init to bad value)
   integer :: nlat = -1          ! latitude dimension length (init to bad value)
   integer :: nlev = -1          ! level dimension length (init to bad value)
   integer :: ntime              ! length of unlimited dimension (if present)
   integer :: size = -1          ! dimension length
   integer :: cmdlen             ! length of cmd line
   integer :: hislen             ! length of input history attribute
   integer :: totlen             ! length of output history attribute
   integer :: ngatts             ! number of global attributes
   integer :: ndims              ! number of dimensions
   integer :: dimid              ! dimension id
   integer :: nvars              ! number of variables
   integer :: ret                ! return code
   integer :: old_mode           ! returned from nf_set_fill

   logical :: itexists           ! return fron "inquire"
   logical :: reverse = .false.  ! true => convert reduced->full
   logical :: verbose = .false.  ! verbose printing
   logical :: silent = .false.   ! unverbose printing

   character(len=8) :: interp = 'linear'    ! interpolation type (default linear)
   character(len=80) :: arg                 ! cmdline argument
   character(len=80) :: file1 = ' '         ! input netcdf file to be converted
   character(len=80) :: file2 = ' '         ! output converted netcdf file
   character(len=256) :: cmdline            ! input command line
   character(len=nf_max_name) :: attname    ! attribute name
   character(len=nf_max_name) :: name       ! dimension name
   character, allocatable :: history(:)     ! history attribute
!
! Externals
!
   integer, external :: iargc               ! number of cmdline arguments
!
! Parse argument list
!
   nargs = iargc()
   n = 1
   cmdline = char(10) // 'mkrgridnew '
   do while (n <= nargs)
      arg = ' '
      call getarg (n, arg)
      n = n + 1
      select case (arg)
      case ('-c')
         interp = 'cubic'
         cmdline = trim(cmdline) // ' -c'
      case ('-f')
         interp = 'fourier'
         cmdline = trim(cmdline) // ' -f'
      case ('-l')
         interp = 'linear'
         cmdline = trim(cmdline) // ' -l'
      case ('-m')
         call getarg (n, arg)
         n = n + 1
         monosiz = monosiz + 1
         if (monosiz > maxmonosiz) then
            write(6,*)'parameter too small. Recompile with maxmonosiz > ',maxmonosiz
            stop 999
         end if
         monolist(monosiz) = arg(1:len_trim(arg))
         cmdline = trim(cmdline) // ' -m ' // trim(monolist(monosiz))
      case ('-r')
         reverse = .true.
         cmdline = trim(cmdline) // ' -r'
      case ('-s')
         silent = .true.
         cmdline = trim(cmdline) // ' -s'
      case ('-v')
         verbose = .true.
         cmdline = trim(cmdline) // ' -v'
      case default
         if (file1 == ' ') then
            file1 = arg
         else if (file2 == ' ') then
            file2 = arg
         else
            write (6,*) 'Argument ', arg,' is not known'
            call usage_exit (' ')
         end if
         cmdline = trim(cmdline) // ' ' // trim(arg)
      end select
   end do
   
   if (file1 == ' ' .or. file2 == ' ') then
      call usage_exit ('Must enter an input file and an output file')
   else if (silent .and. verbose) then
      call usage_exit ('-s cannot be specified with -v')
   end if
   
   inquire (file=file1, exist=itexists)
   if (.not.itexists) then
      write(6,*)'Unable to find input file1 =', file1
      stop 999
   end if
!
! Open input and output netcdf files
!
   call wrap_open (file1, NF_NOWRITE, ncidi)
   call wrap_create (file2, NF_CLOBBER, ncido)

   ret = nf_set_fill (ncido, NF_NOFILL, old_mode)
   if (ret /= NF_NOERR) then
      write(6,*)'Error calling nf_set_fill:'
      write(6,*)nf_strerror (ret)
      stop 999
   end if

!
! Copy dimension and attribute information from input file to output file
!
   call wrap_inq (ncidi, ndims, nvars, ngatts, unlimdimid)
   if (unlimdimid < 0) then
      write(6,*)'INFO: Input file has no unlimited dimension'
   end if

   if (nvars < 1) then
      write(6,*)'nvars =', nvars, ' must be > 0'
      stop 999
   end if

   ret = nf_inq_dimlen (ncidi, unlimdimid, ntime)
   if (ret /= NF_NOERR) then
      ntime = -1
      write(6,*)'INFO: Input file has no unlimited dimension'
   end if
   
   do n=1,ndims
      call wrap_inq_dim (ncidi, n, name, size)
      if (n == unlimdimid) then
         call wrap_def_dim (ncido, name, NF_UNLIMITED, dimid)
      else
         call wrap_def_dim (ncido, name, size, dimid)
      end if
      
      if (name == 'lon') then
         londimid = dimid
         call wrap_inq_dimlen (ncidi, londimid, plon)
      else if (name == 'lat') then
         latdimid = dimid
         call wrap_inq_dimlen (ncidi, latdimid, nlat)
      else if (name == 'lev') then
         levdimid = dimid
         call wrap_inq_dimlen (ncidi, levdimid, nlev)
      else if (name == 'ilev') then
         ilevdimid = dimid
      end if
!
! This code relies on input and output file dimension ids matching, so
! check to ensure it is true 
!
      if (dimid /= n) then
         write(6,*)'Input dimid not equal to output dimid'
         stop 999
      end if
   end do
!
! Ensure both lon and lat dimsizes are ok (implicitly ensures that lon and lat
! dims are there) 
!
   if (plon < 0) then
      write(6,*)'dimension var lon not found'
      stop 999
   end if

   if (nlat < 0) then
      write(6,*)'dimension var lat not found'
      stop 999
   end if

   if (nlat > maxnlat) then
      write(6,*)'maxnlat too small: recompile with this parameter > ',nlat
      stop 999
   end if
!
! Copy global attributes
!
   do n=1,ngatts
      call wrap_inq_attname (ncidi, NF_GLOBAL, n, attname)
      call wrap_copy_att (ncidi, NF_GLOBAL, attname, ncido, NF_GLOBAL)
   end do
!
! Add to or define history attribute.
!
   cmdlen = len_trim (cmdline)
   
   if (nf_inq_attlen (ncido, nf_global, 'history', hislen) == nf_noerr) then
      totlen = cmdlen + hislen
      allocate (history(totlen))
      if (nf_get_att_text (ncido, nf_global, 'history', history) /= nf_noerr) then
         write(6,*)'mkrgridnew: bad attempt to get history attribute'
         stop 999
      end if
   else
      hislen = 0
      totlen = cmdlen
      allocate (history(totlen))
   end if
   
   do i=1,cmdlen
      history(hislen+i) = cmdline(i:i)
   end do
  
   call wrap_put_att_text (ncido, nf_global, 'history', totlen, history)
   deallocate (history)
!
! Call driving routine to convert the data
!
   call driver (plon, nlat, nlev, ntime, nvars, &
                interp, monolist, reverse, silent, verbose, &
                monosiz)
!
! Close files and exit
!
   call wrap_close (ncidi)
   call wrap_close (ncido)

   write(6,*)'Successful completion of mkrgridnew.  Output file is ', trim(file2)

   stop 0
end program fmain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine usage_exit (string)

   implicit none

   character(len=*), intent(in) :: string  ! string to print
   
   if (string /= ' ') write (6,*) string
   write (6,*) 'Usage: mkrgridnew [-c] [-f] [-l] [-m field]... ', &
               '[-r] [-s] [-v] infile outfile [ < namelist]'
   write (6,*)
   write (6,*) '-c: Use cubic interpolation'
   write (6,*) '-f: Use fourier interpolation'
   write (6,*) '-l: Use linear interpolation (default)'
   write (6,*) '-m field ...: Use monotonic cubic interpolation for "field"'
   write (6,*) '-r: Convert reduced->full (default is full->reduced)'
   write (6,*) '-s: Silent output (opposite of verbose)'
   write (6,*) '-v: Verbose output'

   stop 999
end subroutine usage_exit

#ifdef UNICOSMP
subroutine getarg (n, arg)
  implicit none
  integer, intent(in) :: n
  character(len=*), intent(out) :: arg

  integer :: ilen
  integer :: ierr

  call pxfgetarg (n, arg, ilen, ierr)
  if (ierr /= 0) then
    write(6,*)'getarg: ierr not 0:', ierr
    stop 999
  end if
  return
end subroutine getarg

integer function iargc ()
  integer, external :: ipxfargc

  iargc = ipxfargc ()
  return
end function iargc
#endif
