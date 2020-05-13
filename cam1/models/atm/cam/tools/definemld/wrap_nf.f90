subroutine wrap_open (file, mode, nfid)
   implicit none
   include 'netcdf.inc'

   character(len=*) file
   integer mode, nfid

   integer ret

   ret = nf_open (file, mode, nfid)
   if (ret /= nf_noerr) then
      write(6,*)nf_strerror(ret)
      write(6,*)'Unable to open file ', file
      stop 999
   end if
end subroutine wrap_open

subroutine wrap_inq_varid (nfid, varname, varid)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  character*(*) varname

  integer ret

  ret = nf_inq_varid (nfid, varname, varid)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     stop 999
  end if
end subroutine wrap_inq_varid

subroutine wrap_inq_dimlen (nfid, dimid, dimlen)
  implicit none
  include 'netcdf.inc'

  integer nfid, dimid, dimlen

  integer ret

  ret = nf_inq_dimlen (nfid, dimid, dimlen)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     stop 999
  end if
end subroutine wrap_inq_dimlen

subroutine wrap_inq_dimid (nfid, dimname, dimid)
  implicit none
  include 'netcdf.inc'

  integer nfid, dimid
  character*(*) dimname

  integer ret

  ret = nf_inq_dimid (nfid, dimname, dimid)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     stop 999
  end if
end subroutine wrap_inq_dimid

subroutine wrap_inq_var (nfid, varid, varname, xtype, ndims, dimids, natts)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, xtype, ndims, dimids(nf_max_dims), natts
  character*(*) varname

  integer ret

  ret = nf_inq_var (nfid, varid, varname, xtype, ndims, dimids, natts)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     stop 999
  end if
end subroutine wrap_inq_var

subroutine wrap_def_dim (nfid, dimname, len, dimid)
  implicit none
  include 'netcdf.inc'

  integer nfid, len, dimid
  character*(*) dimname

  integer ret

  ret = nf_def_dim (nfid, dimname, len, dimid)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     stop 999
  end if
end subroutine wrap_def_dim

subroutine wrap_get_var_double (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  double precision arr(*)

  integer ret

  ret = nf_get_var_double (nfid, varid, arr)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     stop 999
  end if
end subroutine wrap_get_var_double

subroutine wrap_get_var_int (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  integer arr(*)

  integer ret

  ret = nf_get_var_int (nfid, varid, arr)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     stop 999
  end if
end subroutine wrap_get_var_int

subroutine wrap_put_var_double (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  double precision arr(*)

  integer ret
  ret = nf_put_var_double (nfid, varid, arr)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     stop 999
  end if
end subroutine wrap_put_var_double

subroutine wrap_get_vara_double (nfid, varid, start, count, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, start(*), count(*)
  double precision arr(*)

  integer ret

  ret = nf_get_vara_double (nfid, varid, start, count, arr)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     stop 999
  end if
end subroutine wrap_get_vara_double

subroutine wrap_put_vara_double (nfid, varid, start, count, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  integer start(*), count(*)
  double precision arr(*)

  integer ret
  ret = nf_put_vara_double (nfid, varid, start, count, arr)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     stop 999
  end if
end subroutine wrap_put_vara_double

subroutine wrap_put_att_text (nfid, varid, name, len, text)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, len
  character*(*) name, text

  integer ret

  ret = nf_put_att_text (nfid, varid, name, len, text)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     stop 999
  end if
end subroutine wrap_put_att_text

subroutine wrap_put_att_double (nfid, varid, name, xtype, len, dvals)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, xtype, len
  character*(*) name
  double precision dvals

  integer ret

  ret = nf_put_att_double (nfid, varid, name, xtype, len, dvals)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     stop 999
  end if
end subroutine wrap_put_att_double

subroutine wrap_put_att_int (nfid, varid, name, xtype, len, ivals)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, xtype, len
  character*(*) name
  integer ivals

  integer ret

  ret = nf_put_att_int (nfid, varid, name, xtype, len, ivals)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     stop 999
  end if
end subroutine wrap_put_att_int

