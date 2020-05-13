#include <misc.h>
#include <params.h>
#if ( defined SCAM )
#include <max.h>
#endif

module ncdio_atm

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: ncdio_atm
!
! !DESCRIPTION: 
! Generic interfaces to write fields to netcdf files
!
! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  use shr_sys_mod  , only: shr_sys_flush         ! Standardized system subroutines
  use phys_grid    , only: scatter_field_to_chunk
  use pmgrid       , only: masterproc
  use abortutils   , only: endrun
#if ( defined SPMD )
   use mpishorthand
#endif
#if ( defined SCAM )
   use getnetcdfdata
   use scamMod, only :initlonidx, initlatidx, initTimeIdx
#endif

!
! !PUBLIC TYPES:
  implicit none

PRIVATE

  include 'netcdf.inc'
  save
  public :: check_ret   ! checks return status of netcdf calls
  public :: check_var   ! determine if variable is on netcdf file
  public :: check_dim   ! validity check on dimension
!
!EOP
!
  interface infld
     module procedure infld_real_2d
     module procedure infld_real_3d
  end interface

  public :: infld

#if ( defined SCAM )
  real(r8) dplevs( MAX_DATASET_LEVS +1 )
  integer dplev
  integer STATUS
  real(r8) surfdat
#endif
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_dim
!
! !INTERFACE:
  subroutine check_dim(ncid, dimname, value)
!
! !DESCRIPTION: 
! Validity check on dimension
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: dimname
    integer, intent(in) :: value
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: dimid, dimlen    ! temporaries
!-----------------------------------------------------------------------
    
    call check_ret(nf_inq_dimid (ncid, trim(dimname), dimid), 'check_dim')
    call check_ret(nf_inq_dimlen (ncid, dimid, dimlen), 'check_dim')
    if (dimlen /= value) then
       write (6,*) 'CHECK_DIM error: mismatch of input dimension ',dimlen, &
            ' with expected value ',value,' for variable ',trim(dimname)
       call endrun()
    end if

  end subroutine check_dim

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_var
!
! !INTERFACE:
  subroutine check_var(ncid, varname, varid, readvar)
!
! !DESCRIPTION: 
! Check if variable is on netcdf file
!
! !ARGUMENTS:
    implicit none
    integer, intent(in)          :: ncid
    character(len=*), intent(in) :: varname
    integer, intent(out)         :: varid
    logical, intent(out)         :: readvar 
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ret     ! return value
!-----------------------------------------------------------------------

    readvar = .true.
    if (masterproc) then
       ret = nf_inq_varid (ncid, varname, varid)
       if (ret/=NF_NOERR) then
          write(6,*)'CHECK_VAR Warning:  variable ',trim(varname),' is not on initial dataset'
          call shr_sys_flush(6)
          readvar = .false.
       end if
    end if
  end subroutine check_var

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_ret
!
! !INTERFACE:
  subroutine check_ret(ret, calling)
!
! !DESCRIPTION: 
! Check return status from netcdf call
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ret
    character(len=*) :: calling
!
!EOP
!-----------------------------------------------------------------------

    if (ret /= NF_NOERR) then
       write(6,*)'netcdf error from ',trim(calling)
       call endrun()
    end if

  end subroutine check_ret

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: infld_real_2d
!
! !INTERFACE:
  subroutine infld_real_2d(varname ,ncid    ,lonnam  ,latnam  , dim1b   , &
                           dim1e   ,dim2b   ,dim2e   ,field   , readvar , &
                           grid_map)
!
! !DESCRIPTION: 
! Netcdf i/o of 2d initial real field from netCDF file
!
! !USES
!
   use string_utils, only: to_upper
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: varname  ! variable name
    integer         , intent(in)  :: ncid     ! input unit
    character(len=*), intent(in)  :: lonnam   ! name of longitude dimension of field on file
    character(len=*), intent(in)  :: latnam   ! name of latitude  dimension of field on file
    integer         , intent(in)  :: dim1b    ! start of first  dimension of array to be returned
    integer         , intent(in)  :: dim1e    ! end   of first  dimension of array to be returned
    integer         , intent(in)  :: dim2b    ! start of second dimension of array to be returned
    integer         , intent(in)  :: dim2e    ! end   of second dimension of array to be returned
    real(r8)        , intent(out) :: field(dim1b:dim1e,dim2b:dim2e) ! array to be returned (decomposed or global)
    logical         , intent(out) :: readvar  ! true => variable is on initial dataset
    character(len=*), intent(in)  :: grid_map ! flag indicating which grid to map data to
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i                        ! index
    integer :: ier                      ! error status
    integer :: varid                    ! variable id
    integer :: dimlon, dimlat           ! lon, lat, lev dimension lengths
    integer tmptype
    integer ndims                       ! number of dimensions
    integer dims(NF_MAX_VAR_DIMS)       ! variable shape
    integer londimid, latdimid          ! Dimension ID's
    integer strt(3)                     ! start lon, lat, time indices for netcdf 2-d
    integer cnt (3)                     ! lon, lat, time counts for netcdf 2-d
    data strt/3*1/                      ! 
    data cnt /1,1,1/                    ! 2-d arrs
    real(r8), pointer :: tmp(:,:)       ! input data
    logical :: readvar_tmp              ! if true, variable is on tape
    character*(NF_MAX_NAME) tmpname
    character(len= 8) :: grid_map_tmp   ! Local character string
    character(len=32) :: subname='INFLD_REAL_2D' ! subroutine name
!
!-----------------------------------------------------------------------
!
    grid_map_tmp = trim( to_upper(grid_map) )
!
! Error conditions
!
    if ( grid_map_tmp /= 'GLOBAL' .and. grid_map_tmp /= 'PHYS' ) then
       if(masterproc) then
          write(6,*) trim(subname), ' Error:  invalid grid-map flag specified for field ',trim(varname)
          write(6,*) '                            grid_map = ', grid_map_tmp
       end if
       call endrun()
    end if
    if(.not. masterproc .and. grid_map_tmp == 'GLOBAL') then
       write(6,*) trim(subname),' Error:  while reading field ', trim(varname)
       write(6,*) 'If mapping data to global grid, then call this routine from masterproc only'
       call endrun()
    end if
!
! Read netCDF file
!
    if (masterproc) then
!
! Check if field is on file; get netCDF variable id
!
       call check_var(ncid, varname, varid, readvar_tmp)
!
! If field is on file:
!
       if (readvar_tmp) then

#if ( defined SCAM )
          allocate ( tmp(1,1) )
          call getncdata( ncid, initLatIdx, initLonIdx,  initTimeIdx, &
                          varname, tmp, STATUS )
#else
!
! Get dimension id's and sizes
!
          call wrap_inq_dimid  (ncid, lonnam  , londimid)
          call wrap_inq_dimid  (ncid, latnam  , latdimid)
          call wrap_inq_dimlen (ncid, londimid, dimlon)
          call wrap_inq_dimlen (ncid, latdimid, dimlat)
!
! Check order of dimensions in variable
!
          call wrap_inq_var (ncid, varid, tmpname, tmptype, ndims, dims , ier)
          if (dims(1) /= londimid .or. dims(2) /= latdimid .or. ndims > 3) then
             write(6,*) trim(subname), ' Error: Bad number of dims or ordering while reading field ', trim(varname)
             call endrun()
          end if
!
! Allocate memory and read variable
!
          cnt(1) = dimlon
          cnt(2) = dimlat
          allocate ( tmp(dimlon,dimlat) )
          call wrap_get_vara_realx (ncid, varid, strt, cnt, tmp)
#endif

       end if  ! end of readvar_tmp

    end if  ! end masterproc
!
! Passing data back to calling routine:
!
    if(grid_map_tmp == 'GLOBAL') then
!
! If passing global array back to calling routine, do so only on masterproc
!
       if(masterproc) then
          if(readvar_tmp) then
                
#if ( defined SCAM )
             field(1,1) = tmp(1,1)
#else
             field(dim1b:dim1b+dimlon-1,dim2b:dim2b+dimlat-1) = tmp(:dimlon,:dimlat)
#endif
             deallocate (tmp)
          end if
       end if

    else
!
! Mapping across processors
!
#ifdef SPMD
       call mpi_bcast(readvar_tmp, 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= 0) then
          write(6,*) trim(subname),' Error:  broadcast error while reading ', trim(varname)
          call endrun()
       end if
#endif

       if(readvar_tmp) then

#ifdef SPMD
          call mpi_bcast(dimlon, 1, MPI_INTEGER, 0, mpicom, ier)
          if (ier /= 0) then
             write(6,*) trim(subname),' Error:  broadcast error while reading ', trim(varname)
             call endrun()
          end if
          call mpi_bcast(dimlat, 1, MPI_INTEGER, 0, mpicom, ier)
          if (ier /= 0) then
             write(6,*) trim(subname),' Error:  broadcast error while reading ', trim(varname)
             call endrun()
          end if
#endif

#if ( defined SCAM )
          if(.not. masterproc) allocate ( tmp(1,1) )
          if(grid_map_tmp == 'PHYS') call scatter_field_to_chunk(1,1,1,     1,tmp,field(dim1b,dim2b))
#else
          if(.not. masterproc) allocate ( tmp(dimlon,dimlat) )
          if(grid_map_tmp == 'PHYS') call scatter_field_to_chunk(1,1,1,dimlon,tmp,field(dim1b,dim2b))
#endif
          deallocate (tmp)

       end if

    end if

    readvar = readvar_tmp

    return

  end subroutine infld_real_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: infld_real_3d
!
! !INTERFACE:
  subroutine infld_real_3d(varname ,ncid    ,lonnam  ,levnam  ,latnam  , &
                           dim1b   ,dim1e   ,dim2b   ,dim2e   ,dim3b   , &
                           dim3e   ,field   ,readvar ,grid_map)
!
! !DESCRIPTION: 
! Netcdf i/o of 3d initial real field from netCDF file
!
! !USES
!
   use string_utils, only: to_upper
#if ( defined SCAM )
   use constituents, only: ppcnst, cnst_name
#endif
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: varname  ! variable name
    integer         , intent(in)  :: ncid     ! input unit
    character(len=*), intent(in)  :: lonnam   ! name of longitude dimension of field on file
    character(len=*), intent(in)  :: levnam   ! name of level     dimension of field on file
    character(len=*), intent(in)  :: latnam   ! name of latitude  dimension of field on file
    integer         , intent(in)  :: dim1b    ! start of first  dimension of array to be returned
    integer         , intent(in)  :: dim1e    ! end   of first  dimension of array to be returned
    integer         , intent(in)  :: dim2b    ! start of second dimension of array to be returned
    integer         , intent(in)  :: dim2e    ! end   of second dimension of array to be returned
    integer         , intent(in)  :: dim3b    ! start of third  dimension of array to be returned
    integer         , intent(in)  :: dim3e    ! end   of third  dimension of array to be returned
    real(r8)        , intent(out) :: field(dim1b:dim1e,dim2b:dim2e,dim3b:dim3e) ! array to be returned (decomposed or global)
    logical         , intent(out) :: readvar  ! true => variable is on initial dataset
    character(len=*), intent(in)  :: grid_map ! flag indicating which grid to map data to
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i, k, lat                ! indices
    integer :: ier                      ! error status
    integer :: varid                    ! variable id
    integer :: dimlon, dimlat, dimlev   ! lon, lat, lev dimension lengths
    integer tmptype
    integer ndims                       ! number of dimensions
    integer dims(NF_MAX_VAR_DIMS)       ! variable shape
    integer londimid, latdimid, levdimid ! Dimension ID's
    integer strt(4)                     ! start lon, lat, time indices for netcdf 2-d
    integer cnt (4)                     ! lon, lat, time counts for netcdf 2-d
    data strt/4*1/                      ! 
    data cnt /1,1,1,1/                  ! 3-d arrs
    real(r8), pointer :: tmp    (:,:,:) ! input data
    real(r8), pointer :: tmp_rev(:,:,:) ! input data with lev/lat dimensions reversed
    logical :: readvar_tmp              ! if true, variable is on tape
    character*(NF_MAX_NAME) tmpname
    character(len= 8) :: grid_map_tmp   ! Local character string
    character(len=32) :: subname='INFLD_REAL_3D' ! subroutine name
#if ( defined SCAM )
    integer m
    logical l_interp
    real(r8), pointer :: coldat(:)
#endif
!
!-----------------------------------------------------------------------
!
    grid_map_tmp = trim( to_upper(grid_map) )
!
! Error conditions
!
    if ( grid_map_tmp /= 'GLOBAL' .and. grid_map_tmp /= 'PHYS' ) then
       if(masterproc) then
          write(6,*) trim(subname), ' Error:  invalid grid-map flag specified for field ',trim(varname)
          write(6,*) '                            grid_map = ', grid_map_tmp
       end if
       call endrun()
    end if
    if(.not. masterproc .and. grid_map_tmp == 'GLOBAL') then
       write(6,*) trim(subname),' Error:  while reading field ', trim(varname)
       write(6,*) 'If mapping data to global grid, then call this routine from masterproc only'
       call endrun()
    end if
!
! Read netCDF file
!
    if (masterproc) then
!
! Check if field is on file; get netCDF variable id
!
       call check_var(ncid, varname, varid, readvar_tmp)
!
! If field is on file:
!
       if (readvar_tmp) then
!
! Get dimension id's and sizes
!
          call wrap_inq_dimid  (ncid, lonnam  , londimid)
          call wrap_inq_dimid  (ncid, levnam  , levdimid)
          call wrap_inq_dimid  (ncid, latnam  , latdimid)
          call wrap_inq_dimlen (ncid, londimid, dimlon)
          call wrap_inq_dimlen (ncid, latdimid, dimlat)
          call wrap_inq_dimlen (ncid, levdimid, dimlev)

#if ( defined SCAM )
          allocate ( tmp(1,dim2e+1-dim2b,1) )
          call getncdata( ncid, initLatIdx, initLonIdx,  initTimeIdx, &
               'lev', dplevs, STATUS )
          call getinterpncdata( ncid, initLatIdx, initLonIdx, initTimeIdx, &
               varname, .FALSE.,surfdat, .TRUE., &
               dplevs, dimlev, tmp(1,:,1), STATUS )
#else
!
! Check order of dimensions in variable, allocate memory, and read
! (reverse lat/lev dimensions if field is dimensioned backwards on file)
!
          call wrap_inq_var (ncid, varid, tmpname, tmptype, ndims, dims , ier)

          allocate ( tmp(dimlon,dimlev,dimlat) )
          if     (dims(1) == londimid .and. dims(2) == latdimid .and. &
                  dims(3) == levdimid .and. ndims   <= 4) then
                     cnt(1) = dimlon
                     cnt(2) = dimlat
                     cnt(3) = dimlev
                     allocate ( tmp_rev(dimlon,dimlat,dimlev) )
                     call wrap_get_vara_realx (ncid, varid, strt, cnt, tmp_rev)
                     do lat=1,dimlat
                        do k=1,dimlev
                           do i=1,dimlon
                              tmp(i,k,lat) = tmp_rev(i,lat,k)
                           end do
                        end do
                     end do
                     deallocate (tmp_rev)
          elseif (dims(1) == londimid .and. dims(2) == levdimid .and. &
                  dims(3) == latdimid .and. ndims   <= 4) then
                     cnt(1) = dimlon
                     cnt(2) = dimlev
                     cnt(3) = dimlat
                     call wrap_get_vara_realx (ncid, varid, strt, cnt, tmp)
          else
             write(6,*) trim(subname), ' Error: Bad number of dims or ordering while reading field ', trim(varname)
             call endrun()
          end if
#endif
       end if  ! end of readvar_tmp

    end if  ! end masterproc
!
! Passing data back to calling routine:
!
    if(grid_map_tmp == 'GLOBAL') then
!
! If passing global array back to calling routine, do so only on masterproc
!
       if(masterproc) then
          if(readvar_tmp) then

#if ( defined SCAM )
             field(1,dim2b:dim2e,1) = tmp(1,dim2b:dim2e,1)
#else
             field(dim1b:dim1b+dimlon-1,dim2b:dim2b+dimlev-1,dim3b:dim3b+dimlat-1) = tmp(:dimlon,:dimlev,:dimlat)
#endif
             deallocate (tmp)
          end if
       end if

    else
!
! Mapping across processors
!
#ifdef SPMD
       call mpi_bcast(readvar_tmp, 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= 0) then
          write(6,*) trim(subname),' Error:  broadcast error while reading ', trim(varname)
          call endrun()
       end if
#endif

       if(readvar_tmp) then

#ifdef SPMD
          call mpi_bcast(dimlon, 1, MPI_INTEGER, 0, mpicom, ier)
          if (ier /= 0) then
             write(6,*) trim(subname),' Error:  broadcast error while reading ', trim(varname)
             call endrun()
          end if
          call mpi_bcast(dimlev, 1, MPI_INTEGER, 0, mpicom, ier)
          if (ier /= 0) then
             write(6,*) trim(subname),' Error:  broadcast error while reading ', trim(varname)
             call endrun()
          end if
          call mpi_bcast(dimlat, 1, MPI_INTEGER, 0, mpicom, ier)
          if (ier /= 0) then
             write(6,*) trim(subname),' Error:  broadcast error while reading ', trim(varname)
             call endrun()
          end if
#endif

#if ( defined SCAM )
          if(.not. masterproc) allocate ( tmp(1,dimlev,1) )
          if(grid_map_tmp == 'PHYS') call scatter_field_to_chunk(1,dim2e+1-dim2b,1,     1,tmp,field(dim1b,dim2b,dim3b))
#else
          if(.not. masterproc) allocate ( tmp(dimlon,dimlev,dimlat) )
          if(grid_map_tmp == 'PHYS') call scatter_field_to_chunk(1,dimlev,1,dimlon,tmp,field(dim1b,dim2b,dim3b))
#endif
          deallocate (tmp)

       end if

    end if

    readvar = readvar_tmp

    return

  end subroutine infld_real_3d

end module ncdio_atm






