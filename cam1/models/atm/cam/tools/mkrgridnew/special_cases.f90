subroutine special_cases (plon, nlat, interp, reverse, silent)
!
! Purpose: Read namelist (if applicable) and handle special case variables
!
! Method:  Set or copy special case variables on the output file
!
   use prec,    only: r8
   use globals, only: def_fillvalue, nlon, wnummax, ncidi, ncido, londimid, latdimid

   implicit none

   include 'netcdf.inc'
!
! Arguments
!
   integer, intent(in) :: plon            ! longitude dimension
   integer, intent(in) :: nlat            ! latitude dimension
   character(len=*), intent(in) :: interp ! interpolation type
   logical, intent(in) :: reverse         ! true => convert reduced->full
   logical, intent(in) :: silent          ! unverbose printing
!
! Local workspace
!
   character(len=nf_max_name) :: name     ! variable name
   character(len=80) :: text              ! attribute text

   integer :: ret                         ! return code
   integer :: nvdims                      ! number of variable dimensions
   integer :: i, irow, j                  ! spatial indices
   integer :: length                      ! string length
   integer :: xtype                       ! variable type (netcdf)
   integer :: xtype_rlon                  ! rlon variable type (netcdf)
   integer :: ndims                       ! number of dimensions
   integer :: natts                       ! number of attributes
   integer :: lonid                       ! id for lon
   integer :: rlonid                      ! id for rlon
   integer :: nlonid                      ! id for nlon
   integer :: vardids(nf_max_var_dims)    ! variable dimension ids
   integer :: wnummaxid                   ! id for wnummax

   real(r8) :: alon(plon)                 ! lon dimension var (full grid)
   real(r8) :: rlon(plon,nlat)            ! lon dimension var (reduced grid)

   logical :: haswnummax                  ! whether wnummax found on input file
!
! Namelist values cannot be dynamically allocated
!
   namelist /reduced/ nlon, wnummax
!
! Ensure that special case variables are there, and define appropriate output
! variables
!
   if (reverse) then
      call wrap_inq_varid (ncidi, 'nlon', nlonid)
      call wrap_inq_varid (ncidi, 'rlon', rlonid)

      call wrap_inq_var (ncidi, nlonid, name, xtype, ndims, vardids, natts)

      if (xtype /= NF_INT .or. ndims /= 1 .or. vardids(1) /= latdimid) then
         write(6,*)'Variable nlon on input file is not as expected'
         stop 999
      end if

      ret = nf_inq_varid (ncidi, 'wnummax', wnummaxid)
      if (ret == NF_NOERR) then
         haswnummax = .true.
         call wrap_inq_var (ncidi, wnummaxid, name, xtype, ndims, vardids, natts)
         if (xtype /= NF_INT .or. ndims /= 1 .or. vardids(1) /= latdimid) then
            write(6,*)'Variable wnummax on input file is not as expected'
            stop 999
         end if
         call wrap_get_var_int (ncidi, wnummaxid, wnummax)
         nvdims = 1
         vardids(1) = latdimid
         call wrap_def_var (ncido, 'old_wnummax', nf_int, nvdims, vardids, wnummaxid)

         text = 'old cutoff Fourier wavenumber'
         length = len_trim (text)
         call wrap_put_att_text (ncido, wnummaxid, 'long_name', length, text)
      else
         haswnummax = .false.
         if (interp == 'fourier') then
            write(6,*)'Cannot do fourier interpolation when wnummax not on input dataset'
            stop 999
         end if
      end if

      call wrap_inq_var (ncidi, rlonid, name, xtype_rlon, ndims, vardids, natts)
      if ((xtype_rlon /= nf_double .and. xtype_rlon /= nf_float) .or. &
           ndims /= 2 .or. vardids(1) /= londimid .or. vardids(2) /= latdimid) then
         write(6,*)'Variable rlon on input file is not as expected'
         stop 999
      end if
!
! Get the variables from the input file
!
      call wrap_get_var_int (ncidi, nlonid, nlon)
      call wrap_get_var_double (ncidi, rlonid, rlon)
!
! Define special case variables on the output file.
!
! nlon
!
      nvdims = 1
      vardids(1) = latdimid
      call wrap_def_var (ncido, 'old_nlon', nf_int, nvdims, vardids, nlonid)

      text = 'old number of longitudes'
      length = len_trim (text)
      call wrap_put_att_text (ncido, nlonid, 'long_name', length, text)
!
! rlon
!
      nvdims = 2
      vardids(1) = londimid
      vardids(2) = latdimid
      call wrap_def_var (ncido, 'old_rlon', xtype_rlon, nvdims, vardids, rlonid)

      text = 'old reduced longitude'
      length = len_trim (text)
      call wrap_put_att_text (ncido, rlonid, 'long_name', length, text)

      text = 'degrees_east'
      length = len_trim (text)
      call wrap_put_att_text (ncido, rlonid, 'units', length, text)
!
! lon
!
      nvdims = 1
      vardids(1) = londimid
      call wrap_def_var (ncido, 'lon', xtype_rlon, nvdims, vardids, lonid)

      text = 'longitude'
      length = len_trim (text)
      call wrap_put_att_text (ncido, lonid, 'long_name', length, text)

      text = 'degrees_east'
      length = len_trim (text)
      call wrap_put_att_text (ncido, lonid, 'units', length, text)
      
      do i=1,plon
         alon(i) = (i-1) * 360.0 / plon
      end do
!
! End define mode, write the variables, and go back to define mode
!
      if (nf_enddef (ncido) /= nf_noerr) then
         write(6,*)'bad call to nf_enddef'
         stop 999
      end if

      call wrap_put_var_int (ncido, nlonid, nlon)
      if (haswnummax) then
         call wrap_put_var_int (ncido, wnummaxid, wnummax)
      end if
      call wrap_put_var_double (ncido, rlonid, rlon)
      call wrap_put_var_double (ncido, lonid, alon)

      if (nf_redef (ncido) /= nf_noerr) then
         write(6,*)'special_cases: bad call to nf_redef'
         stop 999
      end if

   else
!
! Full -> reduced.
!
      nlon(:) = plon
      wnummax(:) = -1
      
      write(6,*)'If the code appears to be hung, it might be trying to read a namelist'
      write(6,*)'from stdin, required to define the reduced grid when converting'
      write(6,*)'full->reduced. If so, kill the job and rerun with "< nl.r1up"'
      write(6,*)'(or whatever) appended to the cmd line'
      read (5,reduced)

      do irow=1,nlat/2
         if (wnummax(irow) == -1) then
            wnummax(irow) = nlon(irow)/2
         end if
         wnummax(nlat-irow+1) = wnummax(irow)
      end do
!
! Ensure valid input
!
      do j=1,nlat
         if (nlon(j) < 1 .or. nlon(j) > plon) then
            write(6,*)'special_cases: bad nlon(',j,')=',nlon(j)
         end if
         if (wnummax(j) > nlon(j)/2) then
            write(6,*)'special_cases: bad wnummax(',j,')=',wnummax(j)
         end if
      end do

      do j=1,nlat
         rlon(:,j) = def_fillvalue
         do i=1,nlon(j)
            rlon(i,j) = (i-1) * 360.0 / nlon(j)
         end do
      end do
!
! Make output "rlon" the same type as input "lon"
!
      nvdims = 2
      vardids(1) = londimid
      vardids(2) = latdimid
      call wrap_inq_varid (ncidi, 'lon', lonid)
      if (nf_inq_vartype (ncidi, lonid, xtype) /= nf_noerr) stop 999
      call wrap_def_var (ncido, 'rlon', xtype, nvdims, vardids, rlonid)
      call wrap_put_att_double (ncido, rlonid, '_Fillvalue', NF_DOUBLE, 1, def_fillvalue)

      nvdims = 1
      vardids(1) = latdimid
      call wrap_def_var (ncido, 'nlon', nf_int, nvdims, vardids, nlonid)

      nvdims = 1
      vardids(1) = latdimid
      call wrap_def_var (ncido, 'wnummax', nf_int, nvdims, vardids, wnummaxid)
!
! End define mode, write the variables, and go back to define mode
!
      if (nf_enddef (ncido) /= nf_noerr) then
         write(6,*)'bad call to nf_enddef'
         stop 999
      end if

      call wrap_put_var_double (ncido, rlonid, rlon)
      call wrap_put_var_int (ncido, nlonid, nlon)
      call wrap_put_var_int (ncido, wnummaxid, wnummax)

      if (nf_redef (ncido) /= nf_noerr) then
         write(6,*)'bad call to nf_redef'
         stop 999
      end if

   end if

   if (.not. silent) then
      if (reverse) then
         write(6,*)'Input reduced grid:'
      else
         write(6,*)'Output reduced grid:'
      end if
      write(6,*)'nlon=', nlon(:nlat)

      if (reverse) then
         write(6,*)'Input fourier wavenumber truncation'
      else
         write(6,*)'Output fourier wavenumber truncation'
      end if
      write(6,*)'wnummax=', wnummax(:nlat/2)
   end if

   return
end subroutine special_cases
