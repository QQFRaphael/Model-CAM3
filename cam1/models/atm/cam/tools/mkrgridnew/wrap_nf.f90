!-------------------------------------------------------------------------------
!
! Purpose:
!
! Wrapper routines for the netCDF library for input and output data.
! NOTE: This code assumes that if r4 = r8, the proper floating point netcdf
! put/get routine is nf...real instead of nf...double.  Not perfect but
! avoids machine-specific preprocessor ifdefs.
!
! Author: Jim Rosinski
!
!-------------------------------------------------------------------------------
subroutine wrap_create (path, cmode, ncid)
  implicit none
  include 'netcdf.inc'
      
  character(len=*) path
  integer cmode, ncid

  integer ret

  ret =      nf_create (path, cmode, ncid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_inq (nfid, ndims, nvars, ngatts, unlimdimid)
  implicit none
  include 'netcdf.inc'
  
  integer nfid, ndims, nvars, ngatts, unlimdimid
  
  integer ret

  ret =      nf_inq (nfid, ndims, nvars, ngatts, unlimdimid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_inq_dim (nfid, dimid, name, length)
  implicit none
  include 'netcdf.inc'
  
  integer nfid, dimid, length
  character(len=*) name
  
  integer ret

  ret =      nf_inq_dim (nfid, dimid, name, length)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_inq_dimid (nfid, dimname, dimid)
  implicit none
  include 'netcdf.inc'
  
  integer nfid, dimid
  character(len=*) dimname
  
  integer ret

  ret =      nf_inq_dimid (nfid, dimname, dimid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_inq_unlimdim (nfid, unlimdimid)
  implicit none
  include 'netcdf.inc'
  
  integer nfid, unlimdimid
  
  integer ret

  ret =      nf_inq_unlimdim (nfid, unlimdimid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_inq_dimlen (nfid, dimid, dimlen)
  implicit none
  include 'netcdf.inc'
  
  integer nfid, dimid, dimlen
  
  integer ret
  
  ret =      nf_inq_dimlen (nfid, dimid, dimlen)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

!------------------------------------------------------------------------------
subroutine wrap_inq_vardimid (nfid, varid, dimids)
  implicit none
  include 'netcdf.inc'
  
  integer nfid, varid, dimids(*)
  
  integer ret
  
  ret =      nf_inq_vardimid (nfid, varid, dimids)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_inq_varid (nfid, varname, varid)
  implicit none
  include 'netcdf.inc'
  
  integer nfid, varid
  character(len=*) varname
  
  integer ret
  
  ret =      nf_inq_varid (nfid, varname, varid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine 

!------------------------------------------------------------------------------

subroutine wrap_inq_var (nfid, varid, varname, xtype, ndims, dimids, natts)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, xtype, ndims, dimids(*), natts
  character(len=*) varname

  integer ret

  ret =      nf_inq_var (nfid, varid, varname, xtype, ndims, dimids, natts)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_inq_varname (nfid, varid, varname)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  character(len=*) varname

  integer ret

  ret =      nf_inq_varname (nfid, varid, varname)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_inq_attname (nfid, varid, num, attname)
  implicit none
      
  include 'netcdf.inc'

  integer nfid, varid, num
  character(len=*) attname

  integer ret

  ret =      nf_inq_attname (nfid, varid, num, attname)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_get_att_text (nfid, varid, attname, atttext)
  implicit none
      
  include 'netcdf.inc'

  integer nfid, varid
  character(len=*) attname, atttext

  integer ret

  ret =      nf_get_att_text (nfid, varid, attname, atttext)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_put_att_text (nfid, varid, attname, length, atttext)
  implicit none
      
  include 'netcdf.inc'

  integer nfid, varid, length
  character(len=*) attname, atttext

  integer ret

  ret =      nf_put_att_text (nfid, varid, attname, length, atttext)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_get_att_double (nfid, varid, attname, dvals)
  use prec
  implicit none
      
  include 'netcdf.inc'

  integer nfid, varid
  character(len=*) attname
  real(r8) dvals

  integer ret

  if (r4 == r8) then
     ret = nf_get_att_real (nfid, varid, attname, dvals)
  else
     ret = nf_get_att_double (nfid, varid, attname, dvals)
  end if
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_put_att_double (nfid, varid, attname, xtype, len, attval)
  use prec
  implicit none
      
  include 'netcdf.inc'

  integer nfid, varid, xtype, len
  character(len=*) attname
  real(r8) attval

  integer ret

  if (r4 == r8) then
     ret = nf_put_att_real (nfid, varid, attname, xtype, len, attval)
  else
     ret = nf_put_att_double (nfid, varid, attname, xtype, len, attval)
  end if
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_copy_att (nfid, varid, attname, nfido, varido)
  implicit none
      
  include 'netcdf.inc'

  integer nfid, varid, nfido, varido
  character(len=*) attname

  integer ret

  ret =      nf_copy_att (nfid, varid, attname, nfido, varido)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_def_dim (nfid, dimname, len, dimid)
  implicit none
  include 'netcdf.inc'

  integer nfid, len, dimid
  character(len=*) dimname
      
  integer ret

  ret =      nf_def_dim (nfid, dimname, len, dimid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_def_var (nfid, name, xtype, nvdims, vdims, varid)
  implicit none
  include 'netcdf.inc'

  integer nfid, xtype, nvdims, varid
  integer vdims(nvdims)
  character(len=*) name
  
  integer ret

  ret =      nf_def_var (nfid, name, xtype, nvdims, vdims, varid)
!      write(6,*)'WRAP_DEF_VAR: ',name,' has varid ',varid,nvdims,
!     $'dimensions of ids ',(vdims(i),i=1,nvdims)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_get_var_double (nfid, varid, arr)
  use prec
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  real(r8) arr(*)

  integer ret

  if (r4 == r8) then
     ret = nf_get_var_real (nfid, varid, arr)
  else
     ret = nf_get_var_double (nfid, varid, arr)
  end if
  if (ret.ne.NF_NOERR) then
    write(6,*)'WRAP_GET_VAR_DOUBLE: error reading varid =', varid
    call handle_error (ret)
  end if
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_get_var_int (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, arr(*)

  integer ret

  ret =      nf_get_var_int (nfid, varid, arr)
  if (ret.ne.NF_NOERR) then
    write(6,*)'WRAP_GET_VAR_INT: error reading varid =', varid
    call handle_error (ret)
  end if
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_get_var_text (nfid, varid, text)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  character(len=*) text

  integer ret

  ret =      nf_get_var_text (nfid, varid, text)
  if (ret.ne.NF_NOERR) then
    write(6,*)'WRAP_GET_VAR_TEXT: error reading varid =', varid
    call handle_error (ret)
  end if
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_put_var_text (nfid, varid, text)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  character(len=*) text

  integer ret

  ret =      nf_put_var_text (nfid, varid, text)
  if (ret.ne.NF_NOERR) then
    write(6,*)'WRAP_PUT_VAR_TEXT: error writing varid =', varid
    call handle_error (ret)
  end if
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_put_vara_text (nfid, varid, start, kount, text)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  integer start(*), kount(*)
  character(len=*) text

  integer ret

  ret =      nf_put_vara_text (nfid, varid, start, kount, text)
  if (ret.ne.NF_NOERR) then
    write(6,*)'WRAP_PUT_VARA_TEXT: error writing varid =', varid
    call handle_error (ret)
  end if
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_get_vara_text (nfid, varid, start, kount, text)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  integer start(*), kount(*)
  character(len=*) text

  integer ret

  ret =      nf_get_vara_text (nfid, varid, start, kount, text)
  if (ret.ne.NF_NOERR) then
    write(6,*)'WRAP_GET_VARA_TEXT: error writing varid =', varid
    call handle_error (ret)
  end if
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_get_vara_double (nfid, varid, start, kount, arr)
  use prec
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, start(*), kount(*)
  real(r8) arr(*)

  integer ret
  
  if (r4 == r8) then
     ret = nf_get_vara_real (nfid, varid, start, kount, arr)
  else
     ret = nf_get_vara_double (nfid, varid, start, kount, arr)
  end if
  if (ret.ne.NF_NOERR) then
    write(6,*)'WRAP_GET_VARA8: error reading varid =', varid
    call handle_error (ret)
  end if
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_get_vara_int (nfid, varid, start, kount, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, start(*), kount(*)
  integer arr(*)

  integer ret

  ret =      nf_get_vara_int (nfid, varid, start, kount, arr)
  if (ret.ne.NF_NOERR) then
    write(6,*)'WRAP_GET_VARA8: error reading varid =', varid
    call handle_error (ret)
  end if
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_open (path, omode, ncid)
  implicit none
  include 'netcdf.inc'

  character(len=*) path
  integer omode
  integer ncid
  integer ret
  
  ret =      nf_open (path, omode, ncid)
  if (ret.ne.NF_NOERR) then
    write(6,*)'WRAP_OPEN: nf_open failed for file ',path
    call handle_error (ret)
  end if
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_close (ncid)
  implicit none
  include 'netcdf.inc'

  integer ncid
  integer ret

  ret =      nf_close (ncid)
  if (ret.ne.NF_NOERR) then
    write(6,*)'WRAP_CLOSE: nf_close failed for id ',ncid
    call handle_error (ret)
  end if
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_put_var_int (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  integer arr(*)

  integer ret
  ret =      nf_put_var_int (nfid, varid, arr)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_put_vara_int (nfid, varid, start, kount, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, start(*), kount(*)
  integer arr(*)

  integer ret
  ret =      nf_put_vara_int (nfid, varid, start, kount, arr)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_put_var_double (nfid, varid, arr)
  use prec
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  real(r8) arr(*)

  integer ret

  if (r4 == r8) then
     ret = nf_put_var_real (nfid, varid, arr)
  else
     ret = nf_put_var_double (nfid, varid, arr)
  end if
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_put_vara_double (nfid, varid, start, kount, arr)
  use prec
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  integer start(*), kount(*)
  real(r8) arr(*)
  
  integer ret

  if (r4 == r8) then
     ret = nf_put_vara_real (nfid, varid, start, kount, arr)
  else
     ret = nf_put_vara_double (nfid, varid, start, kount, arr)
  end if
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

!------------------------------------------------------------------------------

subroutine handle_error(ret)
  implicit none
  include 'netcdf.inc'

  integer ret
  
  write(6,*)'handle_error:', nf_strerror(ret)
  call abort
end subroutine
