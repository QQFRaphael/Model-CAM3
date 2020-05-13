program fmain
!----------------------------------------------------------------------------------
!
! Purpose: Bin, neighbor-fill, cap, and smooth 1x1 yearly averaged mixed layer depths 
! to grid of output file specified by -o argument.  The input mixed 
! layer depths were presumeably created by running definemld1x1, the code for which is 
! also in this directory.  It is assumed that the binning is to an output grid coarser 
! than the 1x1 degree resolution of the MLD data.  If data are desired on a finer grid,
! interpolation rather than binning would probably be more appropriate.
!
! Step sequence:
!
!   1. Bin from 1x1 degree grid to desired output resolution, accounting for missing 
!      (i.e. land) values.
!   2. For points on the coarse grid which the (coarse) dataset says are non-land but for 
!      which no points on the 1x1 grid fell in the bin, perform a nearest neighbor fill to 
!      define values at those points.  Currently implemented as avg of all non-land N/S and 
!      E/W neighbors.  So does not work on reduced grid.
!   3. Cap the output at 200 meters.  Mixed layer depths deeper than this can cause
!      excessively long adjustment times for the Slab Ocean Model (SOM). See Kiehl, et. al,
!      Description of the NCAR Community Climate Model (CCM3).
!   4. Apply a 1-2-1 smoothing simultaneously in the x and y directions to prevent
!      computational noise.  This smoothing was applied 10 times in the original
!      implementation of SOM.  The desired number of applications can be specified by the 
!      -s argument to this program.  In addition to the neighbor fill algorithm, this also
!      does not currently work on a reduced grid.
!
! Usage is: definemldbdy -i input_1x1_MLD_file -o input/output_file_on_output_grid [-v]
!                        [-s number_of_smoothing_iterations] 
!
!----------------------------------------------------------------------------------
   use precision

   implicit none

   include 'netcdf.inc'
!
! Local workspace
!
   real(r8), parameter :: fillvaluein = -99.9_r8     ! land flag for 1x1 grid
   real(r8), parameter :: fillvalueout = 1.e36_r8    ! land flag for output grid

   character*80 :: filein =  ' '          ! input filename for 1x1 MLD data
   character*80 :: fileout = ' '          ! input (grid) /output (MLD) filename for output data
   character*80 :: arg
   character*256 :: cmdline               ! input command line
   character, allocatable :: history(:)   ! history attribute

   logical verbose                        ! Add print statements

   integer cmdlen, hislen, totlen         ! character array lengths
!
! Netcdf file, dimension, and variable id's for input file
!
   integer :: ncidin = -1                 ! input file
   integer :: ncidout = -1                ! output file
   integer :: londimidin = -1             ! longitude dimension
   integer :: latdimidin = -1             ! latitude dimension
   integer :: timedimidin = -1            ! time dimension
   integer :: lonidin = -1                ! longitude variable
   integer :: rlonid = -1                 ! longitude variable (reduced grid)
   integer :: latidin = -1                ! latitude variable
   integer :: nlonid = -1                 ! number of longitudes (reduced grid)
   integer :: mldidin = -1                ! MLD variable
!
! Netcdf file, dimension, and variable id's for output file
!
   integer :: londimidout = -1            ! longitude dimension
   integer :: latdimidout = -1            ! latitude dimension
   integer :: timedimidout = -1           ! time dimension
   integer :: lonidout = -1               ! longitude variable
   integer :: latidout = -1               ! latitude variable
   integer :: timeidout = -1              ! time variable
   integer :: mldidout = -1               ! MLD variable
   integer :: landfracid = -1             ! land fraction

   integer :: dimid(3)                    ! dimension id's for variable to be defined
   integer :: start(3) = (/1, 1, 1/)      ! starting position for variable written
   integer :: count(3) = (/-1, -1, 1/)    ! number of values for variable written
   integer :: plonin = -1                 ! longitude dimension (input grid)
   integer :: nlatin = -1                 ! latitude dimension (input grid)
   integer :: plonout = -1                ! longitude dimension (output grid)
   integer :: nlatout = -1                ! latitude dimension (output grid)
   integer :: ntimeout = -1               ! number of time slices (output grid)
   integer :: sm121count = 10             ! number of times to apply 1-2-1 smoothing (default 10)
   integer i, j                           ! longitude, latitude indices
   integer ret                            ! return code
   integer nargs                          ! input arg
   integer n                              ! index loops thru input args

   real(r8), allocatable :: lonin(:)      ! longitude on input grid (degrees)
   real(r8), allocatable :: latin(:)      ! latitude on input grid (degrees)
   real(r8), allocatable :: mldin(:,:)    ! mixed layer depths on input grid

   integer , allocatable :: nlon(:)       ! number of lons per lat (reduced output grid)
   real(r8), allocatable :: latout(:)     ! latitude on output grid
   real(r8), allocatable :: rlon(:,:)     ! longitude on (reduced) output grid
   real(r8), allocatable :: mldout(:,:)   ! mixed layer depths on output grid
   real(r8), allocatable :: landfrac(:,:) ! land/ocean/sea ice flag on output grid
   real(r8), allocatable :: fillpts(:,:)  ! diagnostic

   integer iargc
   external iargc
!
! Default settings before parsing argument list
!
   verbose = .false.

! parse command line arguments, saving them to be written to history attribute

   nargs = iargc()
   n = 1
   cmdline = char(10) // 'definemldbdy '
   do while (n <= nargs)
      arg = ' '
      call getarg (n, arg)
      n = n + 1
      select case (arg)
      case ('-i')
         call getarg (n, arg)
         n = n + 1
         filein = arg
         cmdline = trim(cmdline) // ' -i ' // trim(filein)
      case ('-o')
         call getarg (n, arg)
         n = n + 1
         fileout = arg
         cmdline = trim(cmdline) // ' -o ' // trim(fileout)
      case ('-s')
         call getarg (n, arg)
         n = n + 1
         read(arg,'(i)') sm121count
         write(6,*)'Data will be run through 1-2-1 filter ',sm121count,' times'
         cmdline = trim(cmdline) // ' -s ' // trim(arg)
      case ('-v')
         verbose = .true.
         cmdline = trim(cmdline) // ' -v'
      case default
         write (6,*) 'Argument ', arg,' is not known'
         call usage_exit (' ')
      end select
   end do

   cmdline = trim(cmdline) // char(10)

   if (filein == ' ') then
      call usage_exit ('No input file specified')
   else if (fileout == ' ') then
      call usage_exit ('No output file specified')
   end if
!
! Open 1x1 file for reading and get dimension info
!
   call wrap_open (trim(filein), nf_nowrite, ncidin)

   call wrap_inq_dimid (ncidin, 'lon', londimidin)
   call wrap_inq_dimlen (ncidin, londimidin, plonin)

   call wrap_inq_dimid (ncidin, 'lat', latdimidin)
   call wrap_inq_dimlen (ncidin, latdimidin, nlatin)

   call wrap_inq_dimid (ncidin, 'time', timedimidin)
!
! Allocate spece for needed data from 1x1 file
!
   allocate (lonin(plonin))
   allocate (latin(nlatin))
   allocate (mldin(plonin,nlatin))
!
! Read needed grid info and mixed layer depths from 1x1 file
!
   call wrap_inq_varid (ncidin, 'lon', lonidin)
   call wrap_inq_varid (ncidin, 'lat', latidin)
   call wrap_inq_varid (ncidin, 'MLDANN', mldidin)

   call wrap_get_var_double (ncidin, lonidin, lonin)
   call wrap_get_var_double (ncidin, latidin, latin)

   count(1) = plonin
   count(2) = nlatin
   call wrap_get_vara_double (ncidin, mldidin, start, count, mldin)
!
! Open output file for read/write and get dimension info
!
   call wrap_open (trim(fileout), nf_write, ncidout)

   call wrap_inq_dimid (ncidout, 'lon', londimidout)
   call wrap_inq_dimlen (ncidout, londimidout, plonout)

   call wrap_inq_dimid (ncidout, 'lat', latdimidout)
   call wrap_inq_dimlen (ncidout, latdimidout, nlatout)

   call wrap_inq_dimid (ncidout, 'time', timedimidout)
   call wrap_inq_dimlen (ncidout, timedimidout, ntimeout)
!
! Allocate space for needed input and output data on output file
!
   allocate (rlon(plonout,nlatout))
   allocate (nlon(nlatout))
   allocate (latout(nlatout))
   allocate (mldout(plonout,nlatout))
   allocate (landfrac(plonout,nlatout))
   allocate (fillpts(plonout,nlatout))
!
! Read needed grid info from output file
!
   call wrap_inq_varid (ncidout, 'lon', lonidout)
   call wrap_inq_varid (ncidout, 'lat', latidout)
   call wrap_inq_varid (ncidout, 'time', timeidout)
!
! Enter define mode for history attribute and possible creation of mixed layer depth variable
!
   if (nf_redef (ncidout) /= nf_noerr) then
      write(6,*)'Entering define mode failed on ', fileout
      stop 999
   end if
!
! Add to or define history attribute.
!
   cmdlen = len_trim (cmdline)

   if (nf_inq_attlen (ncidout, nf_global, 'history', hislen)  == nf_noerr) then
      totlen = cmdlen + hislen
      allocate (history(totlen))
      if (nf_get_att_text (ncidout, nf_global, 'history', history) /= nf_noerr) then
         write(6,*)'nf_get_att_text() failed for history attribute'
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
   
   call wrap_put_att_text (ncidout, nf_global, 'history', totlen, history)
   call wrap_put_att_int (ncidout, nf_global, 'sm121count', NF_INT, 1, sm121count)
   deallocate (history)
!
! If MLD already exists on the output file, overwrite it.  Otherwise define it for output
!
   ret = nf_inq_varid (ncidout, 'MLDANN', mldidout)
   if (ret /= nf_noerr) then
      dimid(1) = londimidout
      dimid(2) = latdimidout

      ret = nf_def_var (ncidout, 'MLDANN', nf_double, 2, dimid, mldidout) 
      if (ret /= nf_noerr) then
         write(6,*)'nf_def_var failed for MLDANN'
         write(6,*)nf_strerror (ret)
         stop 999
      end if

      call wrap_put_att_double (ncidout, mldidout, "_FillValue", NF_DOUBLE, 1, fillvalueout);

   end if
!
! Done with data definition
!
   if (nf_enddef (ncidout) /= nf_noerr) then
      write(6,*)'Ending define mode failed on ', fileout
      stop 999
   end if
!
! Define nlon and rlon if output file is on reduced grid
!
   ret = nf_inq_varid (ncidout, 'nlon', nlonid)
   if (ret == nf_noerr) then                            ! probably reduced grid
      ret = nf_inq_varid (ncidout, 'rlon', rlonid)
      if (ret == nf_noerr) then
         call wrap_get_var_int (ncidout, nlonid, nlon)
         call wrap_get_var_double (ncidout, rlonid, rlon)
      else
         write(6,*)filein, ' has nlon but not rlon => is screwed up'
         stop 999
      end if
   else                                                 ! full grid: define nlon and rlon
      do j=1,nlatout
         nlon(j) = plonout
         do i=1,plonout
            rlon(i,j) = (i-1)*360./plonout
         end do
      end do
   end if
!
! Stop if reduced grid since neighbor fill and 1-2-1 filter algorithms only work on
! full grid.
!
   do j=1,nlatout
      if (nlon(j) /= plonout) then
         write(6,*) 'Output file is on a reduced grid: 1-2-1 smoothing and neighbor fill'
         write(6,*) 'algorithms cannot handle that.  Best to define on full grid and'
         write(6,*) 'then run mkrgrid'
         stop 999
      end if
   end do

   call wrap_get_var_double (ncidout, latidout, latout)
!
! Bin 1x1 mixed layer depths to output grid
!
   call binf2c (plonin,  nlatin,  lonin, latin, mldin, fillvaluein, &
                plonout, nlatout, nlon,  rlon,  latout, mldout, fillvalueout, &
                verbose)
!
! Nearest neighbor fill: use landfrac to define points that need to be filled.
!
   ret = nf_inq_varid (ncidout, 'LANDFRAC', landfracid)
   if (ret == nf_noerr) then
      write(6,*)'Found LANDFRAC: doing neighborfill'
      call wrap_get_var_double (ncidout, landfracid, landfrac)
   else
      write(6,*)'NOTE: LANDFRAC not found on ', trim(fileout), '. Filling MLD over entire globe'
   end if

   landfrac(:,:) = 0.

   call neighborfill (plonout, nlatout, nlon, landfrac, mldout, &
                      fillpts, fillvalueout)
!
! Cap the mixed layer depths at 200 meters
!
   do j=1,nlatout
      do i=1,nlon(j)
         if (mldout(i,j) /= fillvalueout) then
            mldout(i,j) = min (200._r8, mldout(i,j))
         end if
      end do
   end do
!
! Smooth the data if sm121count > 0
!
   do n=1,sm121count
      call sm121 (plonout, nlatout, nlon, mldout, fillvalueout)
   end do
!
! Write out mixed layer depths.
!
   count(1) = plonout
   count(2) = nlatout

   start(3) = 1
   call wrap_put_vara_double (ncidout, mldidout, start, count, mldout)

   if (nf_close (ncidout) /= nf_noerr) then
      write(6,*)'WARNING: BAD CLOSE of ', trim(fileout),'. Output is probably screwed up'
      stop 999
   end if

   ret = nf_close (ncidin)
   ret = nf_close (ncidout)

   stop
end program fmain

subroutine usage_exit (arg)
   implicit none
   character*(*) arg
   
   if (arg /= ' ') write (6,*) arg
   write (6,*) 'Usage: definemldbdy -i inputfile -o outputfile [-v]'
   write (6,*) '       -i inputfile is 1x1 netcdf MLD file'
   write (6,*) '       -o outputfile is input/output netcdf file on which to bin input MLD'
   write (6,*) '       -s number is number of smoothing iterations (default 0)'
   write (6,*) '       -v verbose mode'
   stop 999
end subroutine usage_exit
