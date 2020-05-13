module mkrgrid
!
! Routines to drive the interpolation between full and reduced grids
!
   use prec,          only: r8
   use globals,       only: def_fillvalue, ncidi, ncido, londimid, latdimid, levdimid, &
                            ilevdimid, unlimdimid
   use interpolation, only: lininterp, cubinterp, fourier

PRIVATE

   public :: driver

   include 'netcdf.inc'

CONTAINS

   subroutine driver (plon, nlat, nlev, ntime, nvars, &
                      interp, monolist, reverse, silent, verbose, &
                      monosiz)
!
! Purpose: Driver for the interpolation or copying of variables from input
!          file to output file.
!
! Method:  Loop over input variables, and copy to output file or call the
!          appropriate interpolation routine.
!
      implicit none
!
! Arguments
!
      integer, intent(in) :: plon                 ! longitude dimension
      integer, intent(in) :: nlat                 ! latitude dimension
      integer, intent(in) :: nlev                 ! level dimension
      integer, intent(in) :: ntime                ! size of unlimited dimension (if present)
      integer, intent(in) :: nvars                ! number of variables
      character(len=*), intent(in) :: interp      ! interpolation type
      character(len=*), intent(in) :: monolist(:) ! vars to be interpolated monotonic cubic
      logical, intent(in) :: reverse              ! true => converting reduced->full
      logical, intent(in) :: silent               ! unverbose printing
      logical, intent(in) :: verbose              ! verbose printing
      integer, intent(in) :: monosiz              ! size of monolist
!
! Local workspace
!
      character(len=nf_max_name) :: name          ! variable name
      character(len=nf_max_name) :: dimname       ! dimension name
      character(len=nf_max_name) :: attname       ! attribute name

      integer :: natts                            ! number of attributes for a given var
      integer :: nvdims                           ! number of dimensions for this variable
      integer :: vardids(nf_max_var_dims)         ! variable dimension ids
      integer :: numlev                           ! number of levels
      integer :: i, n                             ! indices
      integer :: dimlen                           ! dimension length
      integer :: t                                ! index over unlimited dimension
      integer :: nt                               ! number of time samples for this variable
      integer :: v                                ! loop index over variable id input file
      integer :: vo                               ! variable id on output file
      integer :: xtype                            ! variable type (netcdf)
      integer :: start(nf_max_var_dims) = -1      ! start array for netcdf calls
      integer :: kount(nf_max_var_dims) = -1      ! kount array for netcdf calls
      integer :: size = -1                        ! size of variable
      integer :: attnum                           ! attribute number

      logical :: mono                             ! monotonic cubic interpolation
!
! Dynamic
!
      integer :: varmap(nvars)                    ! map from input varid to output varid
      real(r8) :: fillvaluein(nvars)              ! fillvalue for input vars
      real(r8) :: fillvalueot(nvars)              ! fillvalue for output vars
!
! Allocatables
!
      character, allocatable :: cbuf(:)      ! character data to be copied

      integer, allocatable :: ibuf(:)        ! integer data to be copied

      real(r8), allocatable :: buf(:)        ! floating point data to be copied
      real(r8), allocatable :: arrin(:,:,:)  ! input array for interpolation
      real(r8), allocatable :: arrot(:,:,:)  ! output array for interpolation

      if (reverse) then
         write(6,*)'Reduced -> full'
      else
         write(6,*)'Full -> reduced'
      end if
!
! Determine interpolation type
!
      if (interp == 'fourier' .or. interp == 'linear' .or. interp == 'cubic') then
         write(6,*) 'Using ', interp, ' for grid interpolation'
         do i=1,monosiz
            write(6,*)'field ',monolist(i),' will be interpolated cubic monotonic'
         end do
      else
         write(6,*)interp, ' is not a valid interpolation type'
         stop 99
      end if
!
! Handle special case variables, then loop over variables to interp or copy
!
      call special_cases (plon, nlat, interp, reverse, silent)

      varmap(:) = -1         ! init to bad value
      do v=1,nvars
         call wrap_inq_var (ncidi, v, name, xtype, nvdims, vardids, natts)
!
! Skip any special case variables: they have already been dealt with
!
         if (is_special_case (name, reverse)) then
            cycle
         end if
!
! Normal case variables
!
         call wrap_def_var (ncido, name, xtype, nvdims, vardids, vo)
         varmap(v) = vo
!
! Copy attributes
!
         do n=1,natts
            call wrap_inq_attname (ncidi, v, n, attname)
            call wrap_copy_att (ncidi, v, attname, ncido, vo)
         end do
!
! If converting full->reduced, add _FillValue attribute to fields which will
! be interpolated if the attribute is not already present.  Use NF_DOUBLE attribute
! type for both float and double vars. per netcdf manual pg. 85.
!
         if (nf_inq_attid (ncidi, v, '_FillValue', attnum) == NF_NOERR) then
            call wrap_get_att_double (ncidi, v, '_FillValue', fillvaluein(v))
            fillvalueot(v) = fillvaluein(v)
         else
            fillvaluein(v) = 0.
            fillvalueot(v) = def_fillvalue
            if (.not.reverse .and. dointerp (nvdims, vardids, xtype)) then
               call wrap_put_att_double (ncido, vo, '_Fillvalue', NF_FLOAT, 1, fillvalueot(v))
            end if
         end if
      end do

      if (nf_enddef (ncido) /= NF_NOERR) stop 999
!
! 2nd loop over variables not folded into first solely for efficiency:
! nf_enddef/nf_redef sequence for each variable can be expensive
!
      do v=1,nvars
         vardids(:) = -999       ! init to bad value other than -1
         call wrap_inq_var (ncidi, v, name, xtype, nvdims, vardids, natts)
!
! Skip any special case variables: they have already been dealt with
!
         if (is_special_case (name, reverse)) then
            cycle
         end if
!
! More efficient to use varmap than inq_varid
!
         vo = varmap(v)
         if (vo < 1) then
            write(6,*)'driver: bad vo=',vo
            stop 999
         end if
!
! Determine length of 3rd dimension
!
         if (dointerp (nvdims, vardids, xtype)) then
            numlev = 1           ! if no level dimension present
            if (nvdims > 2) then
               if (levdimid > 0 .and. vardids(3) == levdimid) then
                  numlev = nlev
               else if (ilevdimid > 0 .and. vardids(3) == ilevdimid) then
                  numlev = nlev+1
               else if (unlimdimid > 0 .and. vardids(3) == unlimdimid) then
                  numlev = 1
               else
                  call wrap_inq_dimlen (ncidi, vardids(3), numlev)
                  if (nf_inq_dimname (ncidi, vardids(3), dimname) == NF_NOERR) then
                     write(6,*)trim(name),' has unknown dimension name ',trim(dimname), &
                          ' of length ',numlev
                     write(6,*)'Using that as number of "levels"'
                  else
                     write(6,*)'driver: nf_inq_dimname failure for dimid ',vardids(3)
                     stop 999
                  end if
               end if
            end if

            allocate (arrin(plon,nlat,numlev))
            allocate (arrot(plon,nlat,numlev))

            if (unlimdimid /= -1 .and. vardids(nvdims) == unlimdimid) then
               nt = ntime
            else
               nt = 1
            end if
!
! Loop over the length of the time coordinate for this variable.
!
            do t=1,nt
               if (vardids(3) == unlimdimid) then
                  start(1:3) = (/1,   1,   t/)
                  kount(1:3) = (/plon,nlat,1/)
               else if (vardids(4) == unlimdimid) then
                  start(1:4) = (/1,   1,   1,     t/)
                  kount(1:4) = (/plon,nlat,numlev,1/)
               else if (nvdims /= 2 .or. nvdims /= 3) then
                  write(6,*)'driver: choking on variable with ', nvdims, ' dims'
                  stop 999
               else
!
! 2-d or 3-d variable w/o unlimited dim.
!
                  start(1:4) = (/1,   1,   1,     -1/)
                  kount(1:4) = (/plon,nlat,numlev,-1/)
               end if
!
! Retrieve the input array
!
               call wrap_get_vara_double (ncidi, v, start, kount, arrin)
!
! Determine whether cubic monotonic flag has been set for this field
!
               mono = .false.
               do i=1,monosiz
                  if (monolist(i) == name) mono = .true.
               end do

               if (.not. silent) then
                  if (mono) then
                     write(6,*)'Interpolating time sample ',t,' of ',trim(name),' cubic monotonic'
                  else 
                     write(6,*)'Interpolating time sample ',t,' of ',trim(name), ' ', interp
                  end if
               end if
!
! Call the appropriate interpolation routine
!         
               if (interp == 'mono' .or. interp == 'cubic') then         
                  call cubinterp (plon, nlat, numlev, reverse, mono, &
                                  arrin, arrot, fillvaluein(v), fillvalueot(v))
               else if (interp == 'fourier') then
                  call fourier (plon, nlat, numlev, reverse, arrin, &
                                arrot, name, fillvaluein(v), fillvalueot(v))
               else if (interp == 'linear') then
                  call lininterp (plon, nlat, numlev, reverse, arrin, &
                                  arrot, fillvaluein(v), fillvalueot(v))
               else
                  write(6,*)interp, ' is not a valid interpolation type'
                  stop 999
               end if
!
! Write to output file
!
               call wrap_put_vara_double (ncido, vo, start, kount, arrot)
            end do

            deallocate (arrin)
            deallocate (arrot)
         else
!
! Copy the variable since we are not interpolating
!
            size = 1
            start(:) = -1
            kount(:) = -1
            do n=1,nvdims
               call wrap_inq_dimlen (ncidi, vardids(n), dimlen)
               start(n) = 1
               kount(n) = dimlen
               size = size*dimlen
            end do

            if (xtype == NF_FLOAT .or. xtype == NF_DOUBLE) then
               if (.not.silent) then
                  write(6,*)'copying floating point var ',trim(name), ' of size ',size
               end if
               allocate (buf(size))
               call wrap_get_vara_double (ncidi, v, start, kount, buf)
               call wrap_put_vara_double (ncido, vo, start, kount, buf)
               deallocate (buf)
            else if (xtype == NF_INT) then
               if (.not.silent) then
                  write(6,*)'copying integer var ',trim(name), ' of size ',size
               end if
               allocate (ibuf(size))
               call wrap_get_vara_int (ncidi, v, start, kount, ibuf)
               call wrap_put_vara_int (ncido, vo, start, kount, ibuf)
               deallocate (ibuf)
            else if (xtype == NF_CHAR) then
               if (.not.silent) then
                  write(6,*)'copying character var ',trim(name), ' of size ',size
               end if
               allocate (cbuf(size))
               call wrap_get_vara_text (ncidi, v, start, kount, cbuf)
               call wrap_put_vara_text (ncido, vo, start, kount, cbuf)
               deallocate (cbuf)
            else
               write(6,*)'Unknown type for variable ',name
               stop 999
            end if
         end if
      end do
      
      return
   end subroutine driver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   logical function is_special_case (name, reverse)
!
! Purpose: Determine whether "name" is a special case variable.  Variables
!          such as "lon" "nlon" and "rlon" are hendled differently.
!
! Method:  Use a hardwired set of names, which differs depending on the
!          direction of the interpolation.
!
      implicit none
!
! Arguments
!
      character(len=*), intent(in) :: name  ! variable name
      logical, intent(in) :: reverse        ! true => converting reduced->full
!
! Local workspace
!
      logical bad      ! true => this variable should not exist on this file

      if (reverse) then

         is_special_case = name == 'lon' .or. name == 'rlon' .or. name == 'nlon' .or. &
                           name == 'wnummax'
         bad = name == 'old_rlon' .or. name == 'old_nlon' .or. &
               name == 'old_wnummax'

      else

         bad = name == 'rlon'
         is_special_case = name == 'nlon' .or. name == 'wnummax' .or. &
                           name == 'lon' .or.name == 'old_rlon' .or. &
                           name == 'old_nlon' .or. name == 'old_wnummax'
      end if

      if (bad) then
         write(6,*)'Variable ',trim(name),' is inconsistent with conversion'
         stop 999
      end if

      return
   end function is_special_case

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   logical function dointerp (nvdims, vardids, xtype)
!
! Purpose: Determine whether to interpolate.
!
! Method:  Return true for anything which is of a floating point type, and
!          dimensioned lon x lat (or lon x lat x time or lon x lat x
!          something x time).
!
      implicit none
!
! Arguments
!
      integer, intent(in) :: nvdims
      integer, intent(in) :: vardids(*)
      integer, intent(in) :: xtype

      dointerp = .false.
      if (nvdims > 1 .and. nvdims < 5) then
!
! Nesting is so strict error checking under lf95 passes
!
         if (vardids(1) == londimid .and. vardids(2) == latdimid .and. &
             (xtype == NF_FLOAT .or. xtype == NF_DOUBLE)) then
            dointerp = .true.
         end if

         if (nvdims == 4 .and. vardids(4) /= unlimdimid) then
            dointerp = .false.
         end if
      end if

      return
   end function dointerp
end module mkrgrid
