#include <misc.h>
#include <params.h>

module inidat

!----------------------------------------------------------------------- 
! 
! Purpose: Read initial dataset and process fields as appropriate
!
! Method: Read and process one field at a time
! 
! Author: J. Olson  May 2004
! 
!-----------------------------------------------------------------------

   use pmgrid
   use rgrid
   use prognostics
   use ncdio_atm
   use shr_kind_mod, only: r8 => shr_kind_r8
   use abortutils  , only: endrun
#if ( defined SPMD )
   use mpishorthand
#endif

   implicit none

PRIVATE
   include 'netcdf.inc'

   integer ixcldice,ixcldliq  ! indices into q3 array for cloud liq and cloud ice
   real(r8), allocatable :: ps_tmp  (:,:  )
   real(r8), allocatable :: phis_tmp(:,:  )
   real(r8), allocatable :: q3_tmp  (:,:,:)
   real(r8), allocatable :: cld_liq (:,:,:)
   real(r8), allocatable :: cld_ice (:,:,:)
   real(r8), allocatable :: t3_tmp  (:,:,:)
   real(r8), allocatable :: arr3d_a (:,:,:)
   real(r8), allocatable :: arr3d_b (:,:,:)

   logical readvar            ! inquiry flag:  true => variable exists on netCDF file

   public :: read_inidat

contains

   subroutine read_inidat
!
!-----------------------------------------------------------------------
!
! Purpose:
! Read initial dataset and spectrally truncate as appropriate.
!
!-----------------------------------------------------------------------
!
! $Id: inidat.F90,v 1.30.4.32 2005/03/08 17:06:48 eaton Exp $
! $Author: eaton $
!
!-----------------------------------------------------------------------
!
    use buffer
    use comsrf
    use phys_grid,    only: scatter_field_to_chunk, gather_chunk_to_field
    use phys_buffer,  only: pbuf, pbuf_times, pbuf_get_fld_idx
    use constituents, only: ppcnst, cnst_name, cnst_read_iv, cnst_get_ind

#include <comctl.h>
#include <comlun.h>
!
!---------------------------Local workspace-----------------------------
!
    integer i,m,n                           ! indices
    integer ncol

    real(r8), pointer, dimension(:,:,:,:) :: cldptr
    real(r8), pointer, dimension(:,:    ) :: arr2d_tmp
    real(r8), pointer, dimension(:,:    ) :: arr2d
    character*16 fieldname                  ! field name

    character*16 :: subname='READ_INIDAT'   ! subroutine name
!
!-----------------------------------------------------------------------
!     May 2004 revision described below (Olson)
!-----------------------------------------------------------------------
!
! This routine reads and processes fields one at a time to minimize 
! memory usage.
!
!   State fields (including PHIS) are read into a global array on 
!     masterproc, processed, and scattered to all processors on the
!     appropriate grid 
!
!   Physics fields are read in and scattered immediately to all
!     processors on the physics grid.
!
!-----------------------------------------------------------------------
!

!
!-------------------------------------
! Allocate memory for temporary arrays
!-------------------------------------
!
! Note if not masterproc still might need to allocate array for spmd case
! since each processor calls MPI_scatter 
!
    allocate ( ps_tmp  (plon,plat     ) )
    allocate ( phis_tmp(plon,plat     ) )
    allocate ( q3_tmp  (plon,plev,plat) )
    allocate ( cld_liq (plon,plev,plat) )
    allocate ( cld_ice (plon,plev,plat) )
    allocate ( t3_tmp  (plon,plev,plat) )
!
!---------------------
! Read required fields
!---------------------

!
!-----------
! 3-D fields
!-----------
!
    allocate ( arr3d_a (plon,plev,plat) )
    allocate ( arr3d_b (plon,plat,plev) )

#if ( ! defined OFFLINE_DYN )

    if(masterproc) then

       fieldname = 'US'
       call infld(fieldname, ncid_ini, 'lon', 'lev', 'slat', 1, plon, 1, plev, 1, plat, &
                                                      arr3d_a, readvar, grid_map='global')
       if(.not. readvar) call endrun()
    end if
    call process_inidat('U')

    if(masterproc) then

       fieldname = 'VS'
       call infld(fieldname, ncid_ini, 'slon', 'lev', 'lat', 1, plon, 1, plev, 1, plat, &
                                                      arr3d_a, readvar, grid_map='global')
       if(.not. readvar) call endrun()
    end if
    call process_inidat('V')

#endif

    if(masterproc) then

       fieldname = 'T'
       call infld(fieldname, ncid_ini, 'lon', 'lev', 'lat', 1, plon, 1, plev, 1, plat, &
                                                      t3_tmp, readvar, grid_map='global')
       if(.not. readvar) call endrun()

    end if
    call process_inidat('T')
!
! Get indices for cloud liquid and cloud ice so that these
! fields can be saved off for later use
!
    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)
!
! Constituents (read and process one at a time)
!
    do m = 1,ppcnst

       if(masterproc) then
          readvar   = .false.
          fieldname = cnst_name(m)
          if(cnst_read_iv(m)) &
             call infld(fieldname, ncid_ini, 'lon', 'lev', 'lat', 1, plon, 1, plev, 1, plat, &
                                                           arr3d_a, readvar, grid_map='global')
       end if
       call process_inidat('CONSTS', m_cnst=m)

    end do

    deallocate ( arr3d_a  )
    deallocate ( arr3d_b  )
!
!-----------
! 2-D fields
!-----------
!
    if(masterproc) then

       fieldname = 'PS'
       call infld(fieldname, ncid_ini, 'lon', 'lat', 1, plon, 1, plat, &
                                    ps_tmp  , readvar, grid_map='global')
       if(.not. readvar) call endrun()

    end if
    call process_inidat('PS')

    if(masterproc) then

       fieldname = 'PHIS'
       readvar   = .false.
       if (ideal_phys .or. aqua_planet) then
          phis_tmp(:,:) = 0.
       else
          call infld(fieldname, ncid_topo, 'lon', 'lat', 1, plon, 1, plat, &
                                       phis_tmp, readvar, grid_map='global')
          if (.not. readvar) &
             call endrun(trim(subname)//' ERROR: PHIS not found on topo dataset.')
       end if

    end if
    call process_inidat('PHIS')

    fieldname = 'SGH'
    if (ideal_phys .or. aqua_planet) then
       sgh(:,:) = 0.
    else
       call infld(fieldname, ncid_topo, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                                      sgh, readvar, grid_map='phys')
       if (.not. readvar) &
          call endrun(trim(subname)//' ERROR: SGH not found on topo dataset.')
    end if

    fieldname = 'SGH30'
    if (ideal_phys .or. aqua_planet) then
       sgh30(:,:) = 0.
    else
       call infld(fieldname, ncid_topo, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                                      sgh30, readvar, grid_map='phys')
       if(.not. readvar) then
          if (masterproc) write (6,*) trim(subname), ' Warning: SGH30 not found on topo dataset.'
          if (masterproc) write (6,*) '        The field SGH30 will be filled using data from SGH.'
          do i = begchunk,endchunk
             ncol = get_ncols_p(i)
             sgh30(:ncol,i) = sgh(:ncol,i)
          end do
       end if
    end if

    fieldname = 'LANDM_COSLAT'
    if (aqua_planet) then
       landm(:,:) = 0.
    else
       call infld(fieldname, ncid_topo, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                                    landm, readvar, grid_map='phys')
       if (.not. readvar) &
          call endrun(trim(subname)//' ERROR: LANDM_COSLAT not found on topo dataset.')
    end if

    fieldname = 'LANDFRAC'
    if (aqua_planet) then
       landfrac(:,:) = 0.
    else
       call infld(fieldname, ncid_topo, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                                 landfrac, readvar, grid_map='phys')
       if (.not. readvar) &
          call endrun(trim(subname)//' ERROR: LANDFRAC not found on topo dataset.')
    end if

#if ( ! defined COUP_CSM )

    fieldname = 'TSICE'
    call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                              tsice   , readvar, grid_map='phys')
    if(.not. readvar) call endrun()

    fieldname = 'SNOWHICE'
    call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                              snowhice, readvar, grid_map='phys')
    if(.not. readvar) call endrun()
!
! read TS1, TS2, TS3, TS4 in "plevmx" loop
!
    allocate ( arr2d(1:pcols,begchunk:endchunk) )
    do m = 1,plevmx
       fieldname = tsnam(m)
       call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                                    arr2d, readvar, grid_map='phys')
       if(.not. readvar) call endrun()
       do i = begchunk,endchunk
          surface_state2d(i)%tssub(:,m) = arr2d(:,i)
       end do
    end do
    deallocate ( arr2d )

#if ( defined COUP_SOM )

    fieldname = 'SICTHK'
    call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                              sicthk  , readvar, grid_map='phys')
    if(.not. readvar) call endrun()

    fieldname = 'ICEFRAC'
    call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                              icefrac , readvar, grid_map='phys')
    if(.not. readvar) call endrun()

    fieldname = 'TSOCN'
    call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                              tsocn   , readvar, grid_map='phys')
    if(.not. readvar) call endrun()

! define an initial ocean fraction and non-land ice fraction
! The 1st "where" stmt used to be done in update_srf_fractions (dev45)

    do i = begchunk,endchunk
       ncol = get_ncols_p(i)
       where (icefrac(:ncol,i) + landfrac(:ncol,i) > 1.0)
          icefrac(:ncol,i) = 1. - landfrac(:ncol,i)
       end where

       where (landfrac(:ncol,i) < 1.)
          aice(:ncol,i) = icefrac(:ncol,i)/(1. - landfrac(:ncol,i))
       elsewhere
          aice(:ncol,i) = 0.
       end where
       ocnfrac(:ncol,i) = 1. - landfrac(:ncol,i) - icefrac(:ncol,i)
    end do

    write(6,*)'INIDAT: ocnfrac=',ocnfrac(1,begchunk)
!
! Master needs global landfrac
!
    call gather_chunk_to_field(1,1,1,plon,landfrac,landfrac_field)
!   write(6,*)'INIDAT iam=',iam,' landfrac=',landfrac
!   write(6,*)'INIDAT iam=',iam,' landfrac_field=',landfrac_field
!
!JR Could read in Focn from initial dataset if available
!
    Focn(:,:) = 0.
#else
    Focn(:,:) = inf
    frzmlt(:,:) = 0.  ! needs to be 0, otherwise test in tstm always true
    tsocn(:,:) = inf
#endif
#endif
!
! Global integral of geopotential height (diagnostic for now)
!
    call global_int

    deallocate ( ps_tmp   )
    deallocate ( phis_tmp )
!
! --------------------------------------------------------------------
! Read optional fields (if not found, will be arbitrarily initialized)
! --------------------------------------------------------------------
!

!
! 2-D fields
!
    fieldname = 'PBLH'
    call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                                 pblht, readvar, grid_map='phys')
    if(.not. readvar) then
       pblht(:,:) = 0.
       if (masterproc) write(6,*) trim(fieldname), ' initialized to 0.'
    end if

    fieldname = 'TPERT'
    call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                                 tpert, readvar, grid_map='phys')
    if(.not. readvar) then
       tpert(:,:) = 0.
       if (masterproc) write(6,*) trim(fieldname), ' initialized to 0.'
    end if

    fieldname = 'QPERT'
    qpert(:,:,:) = 0.
    allocate ( arr2d(1:pcols,begchunk:endchunk) )
    call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                                 arr2d, readvar, grid_map='phys')
    if(.not. readvar) then
       if (masterproc) write(6,*) trim(fieldname), ' initialized to 0.'
       arr2d(:,:) = 0.
    end if
    do i = begchunk,endchunk
       qpert(:,1,:) = arr2d(:,:)
    end do
    deallocate ( arr2d )

    fieldname = 'TBOT'
    allocate ( arr2d(1:pcols,begchunk:endchunk) )
    call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                                 arr2d, readvar, grid_map='phys')
    if(.not. readvar) then
       allocate ( arr2d_tmp(1:plon,1:plat) )
       if (masterproc) write(6,*) trim(fieldname), ' initialized with lowest level of T'
       if (masterproc) arr2d_tmp(:plon,:) = t3_tmp(:plon,plev,:)
       call scatter_field_to_chunk(1,1,1,plon,arr2d_tmp,arr2d)
       deallocate ( arr2d_tmp )
    end if
    do i = begchunk,endchunk
       surface_state2d(i)%tbot(:) = arr2d(:,i)
    end do
    deallocate ( arr2d )

    fieldname = 'TSICERAD'
    call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                                             tsice_rad, readvar, grid_map='phys')
    if(.not. readvar) then
       tsice_rad(:,:) = tsice(:,:)
       if (masterproc) write(6,*) trim(fieldname), ' initialized with TSICE.'
    end if
!
! 3-D fields
!
    fieldname = 'CLOUD'
    m = pbuf_get_fld_idx('CLD')
    cldptr => pbuf(m)%fld_ptr(1,1:pcols,1:pver,begchunk:endchunk,1:pbuf_times)
    call infld(fieldname, ncid_ini, 'lon', 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
                                                       cldptr(:,:,:,1), readvar, grid_map='phys')
    if(.not. readvar) then
       cldptr(:,:,:,1) = 0.
       if (masterproc) write(6,*) trim(fieldname), ' initialized to 0.'
    end if
    do n = 2, pbuf_times
       cldptr(:,:,:,n) = cldptr(:,:,:,1)
    end do

    fieldname = 'QCWAT'
    m = pbuf_get_fld_idx(fieldname)
    cldptr => pbuf(m)%fld_ptr(1,1:pcols,1:pver,begchunk:endchunk,1:pbuf_times)
    call infld(fieldname, ncid_ini, 'lon', 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
                                                       cldptr(:,:,:,1), readvar, grid_map='phys')
    if(.not. readvar) then
       allocate ( arr3d_a(1:plon,1:plev,1:plat) )
       if (masterproc) write(6,*) trim(fieldname), ' initialized with Q'
       if (masterproc) arr3d_a(:plon,:,:) = q3_tmp(:plon,:,:)
       call scatter_field_to_chunk(1,plev,1,plon,arr3d_a,cldptr(:,:,:,1))
       deallocate ( arr3d_a )
    end if
    do n = 2, pbuf_times
       cldptr(:,:,:,n) = cldptr(:,:,:,1)
    end do

    fieldname = 'LCWAT'
    m = pbuf_get_fld_idx(fieldname)
    cldptr => pbuf(m)%fld_ptr(1,1:pcols,1:pver,begchunk:endchunk,1:pbuf_times)
    call infld(fieldname, ncid_ini, 'lon', 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
                                                       cldptr(:,:,:,1), readvar, grid_map='phys')
    if(.not. readvar) then
       allocate ( arr3d_a(1:plon,1:plev,1:plat) )
       if (masterproc) write(6,*) trim(fieldname), ' initialized with CLDICE + CLDLIQ'
       if (masterproc) arr3d_a(:plon,:,:) = cld_ice(:plon,:,:) + cld_liq(:plon,:,:)
       call scatter_field_to_chunk(1,plev,1,plon,arr3d_a,cldptr(:,:,:,1))
       deallocate ( arr3d_a )
    end if
    do n = 2, pbuf_times
       cldptr(:,:,:,n) = cldptr(:,:,:,1)
    end do

    fieldname = 'TCWAT'
    m = pbuf_get_fld_idx(fieldname)
    cldptr => pbuf(m)%fld_ptr(1,1:pcols,1:pver,begchunk:endchunk,1:pbuf_times)
    call infld(fieldname, ncid_ini, 'lon', 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
                                                       cldptr(:,:,:,1), readvar, grid_map='phys')
    if(.not. readvar) then
       allocate ( arr3d_a(1:plon,1:plev,1:plat) )
       if (masterproc) write(6,*) trim(fieldname), ' initialized with T'
       if (masterproc) arr3d_a(:plon,:,:) = t3_tmp(:plon,:,:)
       call scatter_field_to_chunk(1,plev,1,plon,arr3d_a,cldptr(:,:,:,1))
       deallocate ( arr3d_a )
    end if
    do n = 2, pbuf_times
       cldptr(:,:,:,n) = cldptr(:,:,:,1)
    end do

    deallocate ( q3_tmp  )
    deallocate ( cld_liq )
    deallocate ( cld_ice )
    deallocate ( t3_tmp  )

    return

  end subroutine read_inidat

!*********************************************************************

  subroutine process_inidat(fieldname, m_cnst)
!
!-----------------------------------------------------------------------
!
! Purpose:
! Post-process input fields
!
!-----------------------------------------------------------------------
!
! $Id: inidat.F90,v 1.30.4.32 2005/03/08 17:06:48 eaton Exp $
! $Author: eaton $
!
!-----------------------------------------------------------------------
!
    use history,      only: fillvalue
    use constituents, only: cnst_name, qmin
    use chemistry   , only: chem_implements_cnst, chem_init_cnst
    use aerosol_intr, only: aerosol_implements_cnst, aerosol_init_cnst
    use tracers     , only: tracers_implements_cnst, tracers_init_cnst
    use stratiform,      only: stratiform_implements_cnst, stratiform_init_cnst
#if ( defined SPMD )
    use spmd_dyn          , only: comm_y, comm_z
    use parutilitiesmodule, only: parcollective2d, BCSTOP
#endif

#include <comlun.h>
#include <perturb.h>

!
! Input arguments
!
    character(len=*), intent(in)           :: fieldname ! fields to be processed
    integer         , intent(in), optional :: m_cnst    ! constituent index
!
!---------------------------Local workspace-----------------------------
!
    integer  i,j,k,lat                     ! grid and constituent indices
    real(r8) pertval                       ! perturbation value
    integer  varid                         ! netCDF variable id
    integer  ret, attlen                   ! netcdf return values
    logical  phis_hires                    ! true => PHIS came from hi res topo
    real(r8), allocatable :: uv_local (:,:,:)
    character*256 text
    character*80 trunits                   ! tracer untis

    real(r8), pointer, dimension(:,:,:) :: q_tmp

    character*16 :: subname='PROCESS_INIDAT' ! subroutine name

    select case (fieldname)

!----------
! Process U
!----------

    case ('U')

!
! SJL: initialize j=1 because later on arr3d_b will be copied to u3s using f90 array syntax
!
       if(masterproc) then
          arr3d_b(:plon,1,:plev) = fillvalue
!
! Transpose array
!

!$omp parallel do private(i, j, k)
          do k = 1, plev
             do j = 1, plat-1
                do i = 1, plon
                   arr3d_b(i,j+1,k) = arr3d_a(i,k,j)
                enddo
             enddo
          enddo
       end if
!
! Scatter u3
!

#if ( defined SPMD )
       allocate( uv_local(plon,beglat:endlat,beglev:endlev) )
       call scatter( arr3d_b, strip3dxyz, uv_local, mpicom )

!$omp parallel do private(i,j,k)
       do k = beglev,endlev
          do j = beglat,endlat
             do i = 1,plon
                u3s(i,j,k) = uv_local(i,j,k)
             enddo
          enddo
       enddo
       deallocate( uv_local )
#else

!$omp parallel do private(i, j, k)
       do k = beglev,endlev
          do j = beglat,endlat
             do i = 1,plon
                u3s(i,j,k) = arr3d_b(i,j,k)
             enddo
          enddo
       enddo
#endif

!----------
! Process V
!----------

    case ('V')
       if(masterproc) then
!
! Transpose array
!

!$omp parallel do private(i, j, k)
          do k = 1, plev
             do j = 1, plat
                do i = 1, plon
                   arr3d_b(i,j,k) = arr3d_a(i,k,j)
                enddo
             enddo
          enddo
       end if
!
! Scatter v3
!

#if ( defined SPMD )
       allocate( uv_local(plon,beglat:endlat,beglev:endlev) )
       call scatter( arr3d_b, strip3dxyz, uv_local, mpicom )
!$omp parallel do private(i,j,k)
      do k = beglev,endlev
         do j = beglat,endlat
            do i = 1,plon
               v3s(i,j,k) = uv_local(i,j,k)
            enddo
         enddo
      enddo
#else
!$omp parallel do private(i, j, k)
      do k = beglev,endlev
         do j = beglat,endlat
            do i = 1,plon
               v3s(i,j,k) = arr3d_b(i,j,k)
            enddo
         enddo
      enddo
#endif

!----------
! Process T
!----------

    case ('T')

!
! Add random perturbation to temperature if required
!
       if(masterproc) then

          if (pertlim.ne.0.0) then
             write(6,*)trim(subname), ':  Adding random perturbation bounded by +/-', &
                       pertlim,' to initial temperature field'
             do lat = 1,plat
                do k = 1,plev
                   do i = 1,nlon(lat)
                      call random_number (pertval)
                      pertval = 2.*pertlim*(0.5 - pertval)
                      t3_tmp(i,k,lat) = t3_tmp(i,k,lat)*(1. + pertval)
                   end do
                end do
             end do
          end if
!
! Average T at the poles.
!

!$omp parallel do private(k)
          do k = 1, plev
             call xpavg(t3_tmp(1,k,   1), plon)
             call xpavg(t3_tmp(1,k,plat), plon)
          enddo

       end if   ! end of if-masterproc

! Scatter t3

#if ( defined SPMD )
      call scatter( t3_tmp, strip3dxzy, t3, mpicom )
#else
!$omp parallel do private(i, j, k)
      do j = beglat,endlat
         do k = beglev,endlev
            do i = 1,plon
               t3(i,k,j) = t3_tmp(i,k,j)
            enddo
         enddo
      enddo
#endif

!---------------------
! Process Constituents
!---------------------

    case ('CONSTS')

       if (.not. present(m_cnst)) then
          call endrun('  '//trim(subname)//' Error:  m_cnst needs to be present in the'// &
                      ' argument list')
       end if

       if(masterproc) then
!
! Check that all tracer units are in mass mixing ratios
!
          if(readvar) then
             call wrap_inq_varid    (NCID_INI,cnst_name(m_cnst), varid)
             call wrap_get_att_text (NCID_INI,varid,'units',trunits)
             if (trunits(1:5) .ne. 'KG/KG' .and. trunits(1:5) .ne. 'kg/kg') then
                call endrun('  '//trim(subname)//' Error:  Units for tracer ' &
                            //trim(cnst_name(m_cnst))//' must be in KG/KG')
             end if
!
! Constituents not read from initial file are initialized by the package that implements them.
!
          else
             if(m_cnst == 1) then
                call endrun('  '//trim(subname)//' Error:  Q must be on Initial File')
             end if

             allocate ( q_tmp(plon,plev,plat) )
             q_tmp(:,:,:) = 0.
             write(6,*) 'Warning:  Not reading ',cnst_name(m_cnst), ' from IC file.'
             if (stratiform_implements_cnst(cnst_name(m_cnst))) then
                call stratiform_init_cnst(cnst_name(m_cnst), q_tmp)
                write(6,*) '          ', cnst_name(m_cnst), ' initialized by "stratiform_init_cnst"'
             else if (chem_implements_cnst(cnst_name(m_cnst))) then
                call chem_init_cnst(cnst_name(m_cnst), q_tmp)
                write(6,*) '          ', cnst_name(m_cnst), ' initialized by "chem_init_cnst"'
             else if (tracers_implements_cnst(cnst_name(m_cnst))) then
                call tracers_init_cnst(cnst_name(m_cnst), q_tmp)
                write(6,*) '          ', cnst_name(m_cnst), ' initialized by "tracers_init_cnst"'
             else if (aerosol_implements_cnst(cnst_name(m_cnst))) then
                call aerosol_init_cnst(cnst_name(m_cnst), q_tmp)
                write(6,*) '          ', cnst_name(m_cnst), ' initialized by "aerosol_init_cnst"'
             else
                write(6,*) '          ', cnst_name(m_cnst), ' set to 0.'
             end if

             arr3d_a(:plon,:,:) = q_tmp(:,:,:)

             deallocate ( q_tmp )
          end if

!$omp parallel do private(lat)
          do lat = 1,plat
             call qneg3(trim(subname), lat   ,nlon(lat),plon   ,plev    , &
                  m_cnst, m_cnst, qmin(m_cnst) ,arr3d_a(1,1,lat))
          end do

!$omp parallel do private(k)
          do k = 1, plev
             call xpavg(arr3d_a(:,k,   1), plon)
             call xpavg(arr3d_a(:,k,plat), plon)
          enddo
!
! if "Q", "CLDLIQ", or "CLDICE", save off for later use
!
          if(m_cnst == 1       ) q3_tmp (:plon,:,:) = arr3d_a(:plon,:,:)
          if(m_cnst == ixcldliq) cld_liq(:plon,:,:) = arr3d_a(:plon,:,:)
          if(m_cnst == ixcldice) cld_ice(:plon,:,:) = arr3d_a(:plon,:,:)

!$omp parallel do private(j,k)
          do k = 1, plev
             do j = 1, plat
                arr3d_b(:,j,k) = arr3d_a(:,k,j)
             enddo
          enddo

       end if   ! end of if-masterproc

       call scatter_q_field_to_block(arr3d_b, m_cnst)

!-----------
! Process PS
!-----------

    case ('PS')

!
! Average PS at the poles.
!
       if(masterproc) then
          call xpavg(ps_tmp(1,   1), plon)
          call xpavg(ps_tmp(1,plat), plon)
       end if

#if ( defined SPMD )
       if (myid_z .eq. 0) then
          call scatter( ps_tmp, strip2d, ps, comm_y )
       endif
       if (twod_decomp .eq. 1) then
          call parcollective2d( comm_z, BCSTOP, plon, endlat-beglat+1, ps ) 
       endif
#else
       ps(:,:) = ps_tmp(:,:)
#endif

!-------------
! Process PHIS
!-------------

    case ('PHIS')

       if(masterproc) then
!
! Check for presence of 'from_hires' attribute to decide whether to filter
!
          if(readvar) then
             call wrap_inq_varid (ncid_topo, 'PHIS', varid)
             ret = nf_inq_attlen (ncid_topo, varid, 'from_hires', attlen)
             if (ret.eq.NF_NOERR .and. attlen.gt.256) then
                call endrun('  '//trim(subname)//' Error:  from_hires attribute length is too long')
             end if
             ret = nf_get_att_text (ncid_topo, varid, 'from_hires', text)
             if (ret.eq.NF_NOERR .and. text(1:4).eq.'true') then
                phis_hires = .true.
!                write(6,*) trim(subname), ': Will filter input PHIS: attribute from_hires is true'
             else
                phis_hires = .false.
!                write(6,*)trim(subname), ': Will not filter input PHIS: attribute ', &
!                          'from_hires is either false or not present'
             end if
!
! Average PHIS at the poles.
!
             call xpavg(phis_tmp(1,   1), plon)
             call xpavg(phis_tmp(1,plat), plon)

          end if

       end if    ! end of if-masterproc

#if ( defined SPMD )
      if (myid_z .eq. 0) then
         call scatter( phis_tmp, strip2d, phis, comm_y )
      endif
      if (twod_decomp .eq. 1) then
         call parcollective2d( comm_z, BCSTOP, plon, endlat-beglat+1, phis ) 
      endif
#else
      phis(:,:) = phis_tmp(:,:)
#endif

    end select

    return

  end subroutine process_inidat

!*********************************************************************

  subroutine global_int
!
!-----------------------------------------------------------------------
!
! Purpose:
! Compute global integral of geopotential height
!
!-----------------------------------------------------------------------
!
! $Id: inidat.F90,v 1.30.4.32 2005/03/08 17:06:48 eaton Exp $
! $Author: eaton $
!
!-----------------------------------------------------------------------
!
    use commap
    use physconst,    only: gravit

#include <comctl.h>
#include <comhyb.h>
#include <comqfl.h>

!
!---------------------------Local workspace-----------------------------
!
    integer i,lat              ! grid indices
    real(r8) zgssum            ! partial sums of phis
    real(r8) zgsint_tmp        ! Geopotential integral
!
!-----------------------------------------------------------------------
!
    if(masterproc) then

       zgsint_tmp = 0.
!              
! Accumulate average geopotential
!
       do lat = 1, plat
          zgssum = 0.
          do i = 1,nlon(lat)
             zgssum = zgssum + phis_tmp(i,lat)
          end do
          zgsint_tmp = zgsint_tmp + w(lat)*zgssum/nlon(lat)
       end do
!
! Normalize average height
!
       zgsint_tmp = zgsint_tmp*.5/gravit
       zgsint     = zgsint_tmp
!
! Globally avgd sfc. partial pressure of dry air (i.e. global dry mass):
! (appears to have no relevance to FV dycore but left in for now)
!
       tmass0 = 98222./gravit
       if (ideal_phys) tmass0 = 100000./gravit
       write(6,800) zgsint
800    format(/72('*')//'INIDAT:  Globally averaged geopotential height = ',f16.10,' meters'//72('*')/)

    end if    ! end of if-masterproc

#if ( defined SPMD )
    call mpibcast (zgsint,1,mpir8,0,mpicom)
#endif

    return

  end subroutine global_int

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: scatter_q_field_to_block --- scatter a 3D constituent array to prognostic array q3
!
! !INTERFACE: 
   subroutine scatter_q_field_to_block(xyz, cnst_idx)

! !USES:
      use shr_kind_mod, only: r8 => shr_kind_r8
#if ( defined SPMD )
      use mpishorthand, only: mpicom
#endif

! !INPUT PARAMETERS:
      real(r8), dimension(plon,plat,plev), intent(in) :: &
         xyz        ! 3D constituent field
      integer, intent(in) :: &
         cnst_idx   ! constituent index in prognostic array q3

! !DESCRIPTION:
!
! Scatter a 3D constituent array from the master processor to the 
! q3 array in the prognostics module.
!
! !REVISION HISTORY:
!
!   02.07.31   Eaton      Initial version
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer j, k
      real(r8), allocatable :: xyz_local(:,:,:)

!-----------------------------------------------------------------------

#if ( defined SPMD )

      allocate( xyz_local(plon,beglat:endlat,beglev:endlev) )
      call scatter( xyz, strip3dxyz, xyz_local, mpicom )
      do k = beglev,endlev
         do j = beglat,endlat
            q3(:,j,k,cnst_idx) = xyz_local(:,j,k)
         enddo
      enddo
      deallocate( xyz_local )

#else

      do j = beglat,endlat
         do k = beglev,endlev
            q3(:,j,k,cnst_idx) = xyz(:,j,k)
         enddo
      enddo

#endif

      return

!EOC
   end subroutine scatter_q_field_to_block

!-----------------------------------------------------------------------
end module inidat

