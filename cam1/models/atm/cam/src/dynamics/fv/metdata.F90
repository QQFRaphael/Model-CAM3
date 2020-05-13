#include <misc.h>
#include <params.h>

module metdata
!----------------------------------------------------------------------- 
!
! BOP
!
! !MODULE: metdata
!
! !DESCRIPTION
! Handles reading and interpolating offline meteorological data which
! is used to drive the dynamics.
!
! !USES
  use shr_kind_mod, only: r8 => shr_kind_r8,r4 => shr_kind_r4
  use time_manager, only: get_curr_date, get_step_size, calendar
  use pmgrid,       only: plon, plat, plev, masterproc, iam
  use pmgrid,       only: beglat, endlat, beglev, endlev
  use pmgrid,       only: myid_z, strip3dxyz, twod_decomp, strip2d
  use dynamics_vars,only: ng_d, ng_s
  use ppgrid,       only: pcols, pver, begchunk, endchunk
  use time_manager, only: get_curr_calday, get_curr_date, get_step_size
  use abortutils,   only: endrun

  use phys_grid,    only: scatter_field_to_chunk
#if ( defined SPMD )
  use mpishorthand, only: mpicom, mpir8, mpiint,mpichar
  use spmd_dyn          , only: comm_y, comm_z
  use parutilitiesmodule, only: parcollective2d, BCSTOP, parscatter
#endif

  implicit none

  private  ! all unless made public
  save 

! !PUBLIC MEMBERS

  public init_met        ! subroutine to open files, allocate arrays, etc
  public advance_met     ! subroutine to read more data and interpolate
  public get_met_fields  ! interface to set meteorology fields
  public get_us_vs       
  public offline_met_defaultopts
  public offline_met_setopts
  public met_winds_on_walls
  public write_met_restart
  public read_met_restart
  public met_ps_next  ! PS interpoloted to the next timestep

  !------------------------------------------------------------------
  ! Interface to access the meteorology fields.  Possible invocations
  ! are as follows:
  !   call get_met_fields( physics_state, us, vs , tend, dt )
  !   call get_met_fields( u, v )
  !   call get_met_fields( srfflx_state )
  !------------------------------------------------------------------
  Interface get_met_fields                       ! overload accessors
     Module Procedure get_dyn_flds
     Module Procedure get_uv_centered
     Module Procedure get_srf_flds
  End Interface
  
  real(r8), allocatable :: met_ps_next(:,:)   ! PS interpolated to next timestep

  logical :: cell_wall_winds = .false.  ! true => met data winds are cell centered
  logical :: remove_met_file = .false.  ! delete metdata file when finished with it

! !REVISION HISTORY:
!   31 Oct 2003  Francis Vitt     Creation
!   05 Feb 2004  F Vitt  Removed reading/inperpolating PS for current timestep
!                        -- only met_ps_next is needed
!   10 Nov 2004  F Vitt  Implemented ability to read from series of files
!   16 Dec 2004  F Vitt  Added offline_met_defaultopts and offline_met_setopts
!
! EOP
!----------------------------------------------------------------------- 
! $Id: metdata.F90,v 1.1.4.1 2004/12/21 17:11:57 fvitt Exp $
! $Author: fvitt $
!----------------------------------------------------------------------- 

  type input2d
     real(r8), dimension(:,:), pointer :: data
  endtype input2d

  type input3d
     real(r8), dimension(:,:,:), pointer :: data
  endtype input3d

  real(r8), allocatable :: met_t(:,:,:)  ! interpolated temperature 
  real(r8), allocatable :: met_u(:,:,:)  ! interpolated zonal wind
  real(r8), allocatable :: met_v(:,:,:)  ! interpolated meridional wind
  real(r8), allocatable :: met_us(:,:,:) ! interpolated zonal wind -staggered
  real(r8), allocatable :: met_vs(:,:,:) ! interpolated meridional wind -staggered
  real(r8), allocatable :: met_q(:,:,:)  ! interpolated water vapor

  real(r8), allocatable :: met_shflx(:,:)! interpolated surface pressure
  real(r8), allocatable :: met_qflx(:,:) ! interpolated water vapor flux
  real(r8), allocatable :: met_taux(:,:) ! interpolat
  real(r8), allocatable :: met_tauy(:,:) ! interpolated

  type(input3d) :: met_ti(2)
  type(input3d) :: met_ui(2)
  type(input3d) :: met_vi(2)
  type(input3d) :: met_usi(2)
  type(input3d) :: met_vsi(2)
  type(input3d) :: met_qi(2)

  type(input2d) :: met_psi_next(2)
  type(input2d) :: met_shflxi(2)
  type(input2d) :: met_qflxi(2)
  type(input2d) :: met_tauxi(2)
  type(input2d) :: met_tauyi(2)

  integer :: t_id, u_id, v_id, q_id, ps_id ! var ids of the data in the netCDF
  integer :: us_id, vs_id ! ids for the staggered winds
  integer :: shflx_id, qflx_id, taux_id, tauy_id

  integer :: dateid           ! var id of the date in the netCDF
  integer :: secid            ! var id of the sec data 
  real(r8) :: datatimem = -1.e36     ! time of prv. values read in
  real(r8) :: datatimep = -1.e36     ! time of nxt. values read in
  real(r8) :: datatimemn = -1.e36    ! time of prv. values read in for next timestep
  real(r8) :: datatimepn  = -1.e36   ! time of nxt. values read in for next timestep

  integer, parameter :: nm=1    ! array index for previous (minus) data
  integer, parameter :: np=2    ! array indes for next (plus) data

  real(r8) :: curr_mod_time ! model time - calendar day
  real(r8) :: next_mod_time ! model time - calendar day - next time step

  character(len=256) :: curr_filename, next_filename, metdata_file
  integer :: curr_fileid, next_fileid     ! the id of the NetCDF file
  real(r8), pointer, dimension(:)  :: curr_data_times, next_data_times

  real(r8) :: alpha = 1.0 ! don't read in water vapor  
  !   real(r8), private :: alpha = 0.0  ! read in water vaper each time step

  logical :: online_test = .false.
  logical :: debug = .false.

contains

!--------------------------------------------------------------------------
! Get the default runtime options
!--------------------------------------------------------------------------
  subroutine offline_met_defaultopts( met_data_file_out, &
                                      met_remove_file_out, & 
                                      met_cell_wall_winds_out )

    implicit none

    character(len=256), intent(out), optional :: met_data_file_out
    logical, intent(out), optional :: met_remove_file_out
    logical, intent(out), optional :: met_cell_wall_winds_out

    if ( present( met_data_file_out ) ) then
       met_data_file_out = metdata_file
    endif

    if ( present( met_remove_file_out ) ) then
       met_remove_file_out = remove_met_file
    endif

    if ( present( met_cell_wall_winds_out ) ) then
       met_cell_wall_winds_out =  cell_wall_winds
    endif

  end subroutine offline_met_defaultopts

!--------------------------------------------------------------------------
! Set runtime options
!--------------------------------------------------------------------------
  subroutine offline_met_setopts( met_data_file_in, &
                                  met_remove_file_in, &
                                  met_cell_wall_winds_in )

    implicit none

    character(len=256), intent(in), optional :: met_data_file_in
    logical, intent(in), optional :: met_remove_file_in
    logical, intent(in), optional :: met_cell_wall_winds_in

    if ( present( met_data_file_in ) ) then
       metdata_file = met_data_file_in 
    endif

    if ( present( met_remove_file_in ) ) then
       remove_met_file = met_remove_file_in
    endif

    if ( present( met_cell_wall_winds_in ) ) then
       cell_wall_winds = met_cell_wall_winds_in  
    endif

    if (masterproc) then
       write(6,*)'Time-variant meteorological dataset (metdata_file) is: ', trim(metdata_file)
       write(6,*)'Meteorological data file will be removed (remove_met_file): ', remove_met_file
       write(6,*)'Meteorological winds are on cell walls (cell_wall_winds): ', cell_wall_winds
    endif

  end subroutine offline_met_setopts

!--------------------------------------------------------------------------
! Opens file, allocates arrays
!--------------------------------------------------------------------------
  subroutine init_met()
    
    implicit none

#include <comctl.h>
    character(len=256) :: filen

    curr_fileid = -1
    next_fileid = -1

    if (.not. nlres ) then ! initial run (see comctl.h)
       curr_filename = metdata_file
       next_filename = ''
    else
       ! restart run
       ! curr_filename & next_filename already set by restart_dynamics
    endif

    call open_met_datafile( curr_filename, curr_fileid, curr_data_times, check_dims=.true. )

    if ( len_trim(next_filename) > 0 ) &
         call open_met_datafile( next_filename, next_fileid, next_data_times )

    if (masterproc) then

       if (online_test) then
          call wrap_inq_varid( curr_fileid, 'arch_T', t_id )
          call wrap_inq_varid( curr_fileid, 'arch_US', us_id )
          call wrap_inq_varid( curr_fileid, 'arch_VS', vs_id )
          call wrap_inq_varid( curr_fileid, 'arch_Q', q_id )
          call wrap_inq_varid( curr_fileid, 'arch_PS', ps_id )
       else
          call wrap_inq_varid( curr_fileid, 'T', t_id )
          if ( .not. cell_wall_winds ) then
             call wrap_inq_varid( curr_fileid, 'U', u_id )
             call wrap_inq_varid( curr_fileid, 'V', v_id )
          else
             call wrap_inq_varid( curr_fileid, 'US', us_id )
             call wrap_inq_varid( curr_fileid, 'VS', vs_id )
          endif
          call wrap_inq_varid( curr_fileid, 'Q', q_id )
          call wrap_inq_varid( curr_fileid, 'PS', ps_id )
       endif

       call wrap_inq_varid( curr_fileid, 'SHFLX', shflx_id )
       call wrap_inq_varid( curr_fileid, 'QFLX', qflx_id )
       call wrap_inq_varid( curr_fileid, 'TAUX', taux_id )
       call wrap_inq_varid( curr_fileid, 'TAUY', tauy_id )

    endif ! masterproc

!
! allocate space for data arrays ...
! 

    ! physics grid

    allocate( met_ti(nm)%data(pcols,begchunk:endchunk,plev) )
    allocate( met_ti(np)%data(pcols,begchunk:endchunk,plev) )
    allocate( met_t(pcols,begchunk:endchunk,plev) )

    allocate( met_qi(nm)%data(pcols,begchunk:endchunk,plev) )
    allocate( met_qi(np)%data(pcols,begchunk:endchunk,plev) )
    allocate( met_q(pcols,begchunk:endchunk,plev) )

    allocate( met_shflxi(nm)%data(pcols,begchunk:endchunk) )
    allocate( met_shflxi(np)%data(pcols,begchunk:endchunk) )
    allocate( met_shflx(pcols,begchunk:endchunk) )

    allocate( met_qflxi(nm)%data(pcols,begchunk:endchunk) )
    allocate( met_qflxi(np)%data(pcols,begchunk:endchunk) )
    allocate( met_qflx(pcols,begchunk:endchunk) )

    allocate( met_tauxi(nm)%data(pcols,begchunk:endchunk) )
    allocate( met_tauxi(np)%data(pcols,begchunk:endchunk) )
    allocate( met_taux(pcols,begchunk:endchunk) )

    allocate( met_tauyi(nm)%data(pcols,begchunk:endchunk) )
    allocate( met_tauyi(np)%data(pcols,begchunk:endchunk) )
    allocate( met_tauy(pcols,begchunk:endchunk) )

    ! dynamics grid

    allocate( met_psi_next(nm)%data(plon, beglat:endlat) )
    allocate( met_psi_next(np)%data(plon, beglat:endlat) )
    allocate( met_ps_next(plon, beglat:endlat) )

    allocate( met_us(plon, beglat-ng_d:endlat+ng_s, beglev:endlev) )
    allocate( met_vs(plon, beglat-ng_s:endlat+ng_d, beglev:endlev) )

    if (cell_wall_winds) then
       allocate( met_usi(nm)%data(plon, beglat:endlat, beglev:endlev) )
       allocate( met_usi(np)%data(plon, beglat:endlat, beglev:endlev) )
       allocate( met_vsi(nm)%data(plon, beglat:endlat, beglev:endlev) )
       allocate( met_vsi(np)%data(plon, beglat:endlat, beglev:endlev) )
    endif

    if (.not. cell_wall_winds) then

       allocate( met_u(plon, beglat-ng_d:endlat+ng_d, beglev:endlev) )
       allocate( met_ui(nm)%data(plon, beglat:endlat, beglev:endlev) )
       allocate( met_ui(np)%data(plon, beglat:endlat, beglev:endlev) )

       allocate( met_v(plon, beglat-ng_s:endlat+ng_d, beglev:endlev) )
       allocate( met_vi(nm)%data(plon, beglat:endlat, beglev:endlev) )
       allocate( met_vi(np)%data(plon, beglat:endlat, beglev:endlev) )

    endif

  end subroutine init_met


!-----------------------------------------------------------------------
! Reads more data if needed and interpolates data to current model time 
!-----------------------------------------------------------------------
  subroutine advance_met()

    implicit none

    call t_startf('MET__advance')

!
! master processor reads the data and scatters
! the data to each mpi task
! all mpi tasks need to call read routines since they all need
! to receive the scatter data
!
    call get_model_time()

    if ( ( curr_mod_time > datatimep ) .or. &
         ( next_mod_time > datatimepn ) ) then
       call check_files()
    endif

    if ( curr_mod_time > datatimep ) then
       call read_next_metdata()
    end if

    if ( next_mod_time > datatimepn ) then
       call read_next_ps()
    end if

! need to inperpolate the data, regardless !
! each mpi tasks needs to interpolate
    call interpolate_metdata()

    call t_stopf('MET__advance')

  end subroutine advance_met

!-------------------------------------------------------------------
! Method to get some the meteorology data. 
! Sets the following srfflx_state member fields to the 
! meteorology data :
!   qflx
!   shflx
!   taux
!   tauy
!-------------------------------------------------------------------
  subroutine get_srf_flds( srfflx_state2d )
    use comsrf,    only: srfflx_state
    use phys_grid,      only: get_ncols_p
!    use history,        only: outfld

    implicit none

    type(srfflx_state), intent(inout), dimension(begchunk:endchunk) :: srfflx_state2d

    integer :: c,ncol

    call t_startf('MET__GET_SRF_FLDS')

    do c=begchunk,endchunk
       ncol = get_ncols_p(c)
       srfflx_state2d(c)%wsx(:ncol)    = met_taux(:ncol,c)
       srfflx_state2d(c)%wsy(:ncol)    = met_tauy(:ncol,c)
       srfflx_state2d(c)%shf(:ncol)    = met_shflx(:ncol,c)
       srfflx_state2d(c)%cflx(:ncol,1) = met_qflx(:ncol,c)
    end do                    ! Chunk loop

    if (debug) then
      print*,'METDATA maxval(met_taux),minval(met_taux): ',maxval(met_taux),minval(met_taux)
      print*,'METDATA maxval(met_tauy),minval(met_tauy): ',maxval(met_tauy),minval(met_tauy)
      print*,'METDATA maxval(met_shflx),minval(met_shflx): ',maxval(met_shflx),minval(met_shflx)
      print*,'METDATA maxval(met_qflx),minval(met_qflx): ',maxval(met_qflx),minval(met_qflx)
    endif

!!$    if ( debug ) then
!!$       do c = begchunk, endchunk
!!$          call outfld('MET_TAUX',srfflx_state2d(c)%wsx , pcols   ,c   )
!!$          call outfld('MET_TAUY',srfflx_state2d(c)%wsy , pcols   ,c   )
!!$          call outfld('MET_SHFX',srfflx_state2d(c)%shf , pcols   ,c   )
!!$          call outfld('MET_QFLX',srfflx_state2d(c)%cflx(:,1) , pcols   ,c   )
!!$       enddo
!!$    endif

    call t_stopf('MET__GET_SRF_FLDS')

  end subroutine get_srf_flds

!-------------------------------------------------------------------
! allows access to physics state fields 
!   q  : water vapor
!   ps : surface pressure
!   t  : temperature
!-------------------------------------------------------------------
  subroutine get_dyn_flds( state, tend, dt )

    use shr_kind_mod,   only: r8 => shr_kind_r8
    use physics_types,  only: physics_state, physics_tend, physics_dme_adjust
    use ppgrid,         only: pcols, pver, begchunk, endchunk
    use pmgrid,         only: plon, beglat, endlat, beglev, endlev
    use phys_grid,      only: get_ncols_p
    use dynamics_vars,  only: ng_d, ng_s
!    use history,        only: outfld 

    implicit none

    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: tend
    real(r8),            intent(in   ) :: dt                  ! model physics timestep

    integer :: lats(pcols)           ! array of latitude indices
    integer :: lons(pcols)           ! array of longitude indices
    integer :: c, ncol, i,j,k
    real(r8):: qini(pcols,pver)   ! initial specific humidity

    real(r8) :: tmp(pcols,pver)

    call t_startf('MET__GET_DYN2')

    
    do c = begchunk, endchunk
       ncol = get_ncols_p(c)
       do i=1,ncol

          do k=1,pver

             state(c)%t(i,k)   = met_t(i,c,k)

             qini(i,k) = state(c)%q(i,k,1)

             ! at this point tracer mixing ratios have already been
             ! converted from dry to moist
!!$             if (  moist_q_mmr .and. (.not. online_test)) then
                state(c)%q(i,k,1) = alpha*state(c)%q(i,k,1) + &
                     (1-alpha)*met_q(i,c,k)
!!$             else 
!!$                ! dry-to-moist conversion
!!$                state(c)%q(i,k,1) = alpha*state(c)%q(i,k,1) + &
!!$                     (1.-alpha)*met_q(i,c,k) &
!!$                     * state(c)%pdeldry(i,k)/state(c)%pdel(i,k)
!!$             endif

             if ((state(c)%q(i,k,1) < 0.0).and. (alpha .ne. 1. )) state(c)%q(i,k,1) = 0.0

          end do

       end do

       ! now adjust mass of each layer now that water vapor has changed
       if (( .not. online_test ) .and. (alpha .ne. 1. )) then
          call physics_dme_adjust(state(c), tend(c), qini, dt)
       endif

    end do

    if (debug) then
      print*,'METDATA maxval(met_t),minval(met_t): ', maxval(met_t),minval(met_t)
      print*,'METDATA maxval(met_ps_next),minval(met_ps_next): ',  maxval(met_ps_next),minval(met_ps_next)
    endif

!!$    if ( debug ) then
!!$       do c = begchunk, endchunk
!!$          call outfld('MET_T  ',state(c)%t , pcols   ,c   )
!!$       enddo
!!$       do j = beglat, endlat
!!$          call outfld('MET_PS ',met_ps_next(1,j), plon   ,j   )
!!$       enddo
!!$    endif
    call t_stopf('MET__GET_DYN2')

  end subroutine get_dyn_flds

!------------------------------------------------------------------------
! get the meteorological winds on the grid cell centers (A-grid)
!------------------------------------------------------------------------
  subroutine get_uv_centered( u, v )

    use shr_kind_mod,   only: r8 => shr_kind_r8
    use pmgrid,         only: plon, beglat, endlat, beglev, endlev
    use dynamics_vars,  only: ng_d, ng_s
!    use history,        only: outfld

    implicit none

    real(r8), intent(out) :: u(plon, beglat-ng_d:endlat+ng_d, beglev:endlev)  ! u-wind on A-grid
    real(r8), intent(out) :: v(plon, beglat-ng_s:endlat+ng_d, beglev:endlev)  ! v-wind on A-grid

    integer :: i,j,k

    real(r8) u3s_tmp(plon,beglev:endlev), v3s_tmp(plon,beglev:endlev)

    u(:,:,:) = 0.
    v(:,:,:) = 0.

    u(         :,  max(1,beglat-ng_d):min(plat,endlat+ng_d), beglev:endlev ) = &
         met_u(:,  max(1,beglat-ng_d):min(plat,endlat+ng_d), beglev:endlev )

    v(         :,  max(1,beglat-ng_s):min(plat,endlat+ng_d), beglev:endlev ) = &
         met_v(:,  max(1,beglat-ng_s):min(plat,endlat+ng_d), beglev:endlev )

    if (debug) print*,'METDATA maxval(u),minval(u),maxval(v),minval(v) : ',&
         maxval(u(:,  max(1,beglat-ng_d):min(plat,endlat+ng_d), beglev:endlev )),&
         minval(u(:,  max(1,beglat-ng_d):min(plat,endlat+ng_d), beglev:endlev )),&
         maxval(v(:,  max(1,beglat-ng_s):min(plat,endlat+ng_d), beglev:endlev )),&
         minval(v(:,  max(1,beglat-ng_s):min(plat,endlat+ng_d), beglev:endlev ))

!!$    if ( debug ) then
!!$       do j = beglat, endlat
!!$          do k = beglev, endlev
!!$             do i = 1, plon
!!$                u3s_tmp(i,k) = u(i,j,k)
!!$                v3s_tmp(i,k) = v(i,j,k)
!!$             enddo
!!$          enddo
!!$          call outfld ('MET_U ', u3s_tmp, plon, j )
!!$          call outfld ('MET_V ', v3s_tmp, plon, j )
!!$       enddo
!!$    endif

  end subroutine get_uv_centered

!------------------------------------------------------------------------
! get the meteorological winds on the grid cell walls (vorticity winds)
!   us : staggered zonal wind
!   vs : staggered meridional wind
!------------------------------------------------------------------------
  subroutine get_us_vs( us, vs )

    use shr_kind_mod,   only: r8 => shr_kind_r8
    use pmgrid,         only: plon, beglat, endlat, beglev, endlev
    use dynamics_vars,  only: ng_d, ng_s
!    use history,        only: outfld

    implicit none

    real(r8), intent(inout) :: us(plon, beglat-ng_d:endlat+ng_s, beglev:endlev)    ! u-wind on d-grid
    real(r8), intent(inout) :: vs(plon, beglat-ng_s:endlat+ng_d, beglev:endlev)       ! v-wind on d-grid

    real(r8) u3s_tmp(plon,beglev:endlev), v3s_tmp(plon,beglev:endlev)

    integer :: i,j,k

    call t_startf('MET__get_us_vs')

    us(          :, max(2,beglat-ng_d):  min(plat,endlat+ng_s), beglev:endlev) =      &
         met_us( :, max(2,beglat-ng_d):  min(plat,endlat+ng_s), beglev:endlev)

    vs(          :, max(1,beglat-ng_s):  min(plat,endlat+ng_d),    beglev:endlev) =      &
         met_vs( :, max(1,beglat-ng_s):  min(plat,endlat+ng_d),    beglev:endlev)

    if (debug) print*,iam,': METDATA maxval(us),minval(us),maxval(vs),minval(vs) : ',&
         maxval(us(          :, max(2,beglat-ng_d):  min(plat,endlat+ng_s),    beglev:endlev)),&
         minval(us(          :, max(2,beglat-ng_d):  min(plat,endlat+ng_s),    beglev:endlev)),&
         maxval(vs(          :, max(1,beglat-ng_s):  min(plat,endlat+ng_d),    beglev:endlev)),&
         minval(vs(          :, max(1,beglat-ng_s):  min(plat,endlat+ng_d),    beglev:endlev))

!!$    if (debug) then
!!$        u3s_tmp = 1.e36
!!$       do j = beglat, endlat
!!$          do k = beglev, endlev
!!$             do i = 1, plon
!!$                if (j >= 2) u3s_tmp(i,k) = us(i,j,k)
!!$                v3s_tmp(i,k) = vs(i,j,k)
!!$             enddo
!!$          enddo
!!$          call outfld ('MET_US ', u3s_tmp, plon, j )
!!$          call outfld ('MET_VS ', v3s_tmp, plon, j )
!!$       enddo
!!$    endif
!!$
    call t_stopf('MET__get_us_vs')

  end subroutine get_us_vs

!-------------------------------------------------------------------------
! writes file names to restart file
!-------------------------------------------------------------------------
  subroutine write_met_restart( nrg )
    implicit none
    integer,intent(in) :: nrg     ! Unit number
    integer :: ioerr   ! error status

    if (masterproc) then
       write(nrg, iostat=ioerr) curr_filename
       if (ioerr /= 0 ) then
          write (6,*) 'WRITE ioerror ',ioerr,' on i/o unit = ',nrg
          call endrun ('WRITE_RESTART_DYNAMICS')
       end if
       write(nrg, iostat=ioerr) next_filename
       if (ioerr /= 0 ) then
          write (6,*) 'WRITE ioerror ',ioerr,' on i/o unit = ',nrg
          call endrun ('WRITE_RESTART_DYNAMICS')
       end if
    end if
  end subroutine write_met_restart

!-------------------------------------------------------------------------
! reads file names from restart file
!-------------------------------------------------------------------------
  subroutine read_met_restart( nrg )
    implicit none
    integer,intent(in) :: nrg     ! Unit number
    integer :: ioerr   ! error status

    if (masterproc) then
       read(nrg, iostat=ioerr) curr_filename
       if (ioerr /= 0 ) then
          write (6,*) 'READ ioerror ',ioerr,' on i/o unit = ',nrg
          call endrun ('READ_RESTART_DYNAMICS')
       end if
       read(nrg, iostat=ioerr) next_filename
       if (ioerr /= 0 ) then
          write (6,*) 'READ ioerror ',ioerr,' on i/o unit = ',nrg
          call endrun ('READ_RESTART_DYNAMICS')
       end if
    end if

#if ( defined SPMD )
    call mpibcast ( curr_filename ,len(curr_filename) ,mpichar,0,mpicom)    
    call mpibcast ( next_filename ,len(next_filename) ,mpichar,0,mpicom)    
#endif
  end subroutine read_met_restart

!-------------------------------------------------------------------------
! returns true if the met winds are defined on cell walls
!-------------------------------------------------------------------------
  function met_winds_on_walls()
    logical :: met_winds_on_walls

    met_winds_on_walls = cell_wall_winds
  end function met_winds_on_walls

! internal methods :

!-------------------------------------------------------------------------
! transfers cell-centered winds to cell walls
!-------------------------------------------------------------------------
  subroutine transfer_windsToWalls()

    implicit none

    integer :: i,j,k

    call t_startf('MET__transfer_windsToWalls')

    !$omp parallel do private (i, j, k)
    do k = beglev, endlev

       do j = beglat+1,endlat
          do i = 1,plon
             met_us(i,j,k) = ( met_u(i,j,k) + met_u(i,j-1,k) )*.5
          end do
       end do

#if defined( SPMD )
       if ( beglat .gt. 1 ) then
          do i = 1, plon
             ! met_u is alread ghosted at this point
             met_us(i,beglat,k) = ( met_u(i,beglat,k) + met_u(i,beglat-1,k) )*.5
          enddo
       endif
#endif

       do j = beglat,endlat
          met_vs(1,j,k) = ( met_v(1,j,k) + met_v(plon,j,k) )*.5
          do i = 2,plon
             met_vs(i,j,k) = ( met_v(i,j,k) + met_v(i-1,j,k) )*.5
          end do
       end do
    end do

    call t_stopf('MET__transfer_windsToWalls')

  end subroutine transfer_windsToWalls

  subroutine get_model_time()
    implicit none
    integer yr, mon, day, ncsec  ! components of a date

    call t_startf('MET__get_model_time')

    call get_curr_date(yr, mon, day, ncsec)

    curr_mod_time = get_time_float( yr, mon, day, ncsec )
    next_mod_time = curr_mod_time + get_step_size()/86400._r8

    call t_stopf('MET__get_model_time')

  end subroutine get_model_time

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine check_files( )

    use shr_sys_mod, only: shr_sys_system
    implicit none

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
    character(len=128) :: ctmp
    integer ::  istat


    if (next_mod_time > curr_data_times(size(curr_data_times))) then
       if ( .not. associated(next_data_times) ) then
          ! open next file...
          next_filename = incr_filename( curr_filename )
          call open_met_datafile( next_filename, next_fileid, next_data_times )
       endif
    endif

    if ( associated(next_data_times) ) then
       if (curr_mod_time >= next_data_times(1)) then

          ! close current file ...
          if (masterproc) then
             call wrap_close( curr_fileid )
             ! remove if requested
             if( remove_met_file ) then
                write(*,*) 'check_files: removing file = ',trim(curr_filename) 
                ctmp = 'rm -f ' // trim(curr_filename) 
                write(*,*) 'check_files: fsystem issuing command - '
                write(*,*) trim(ctmp)
                !!call shell_cmd( ctmp, istat )
                call shr_sys_system( ctmp, istat )
             end if
          endif

          curr_filename = next_filename
          curr_fileid = next_fileid

          deallocate( curr_data_times )
          allocate( curr_data_times( size( next_data_times ) ) )
          curr_data_times(:) = next_data_times(:)

          next_fileid = -1
          next_filename = ''

          deallocate( next_data_times )
          nullify(  next_data_times )

       endif
    endif

  end subroutine check_files

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  function incr_filename( filename )

    !-----------------------------------------------------------------------
    ! 	... Increment or decrement a date string withing a filename
    !           the filename date section is assumed to be of the form
    !           yyyy-dd-mm
    !-----------------------------------------------------------------------

    use string_utils,  only : incstr

    implicit none


    character(len=*), intent(in) :: filename ! present dynamical dataset filename
    character(len=256) :: incr_filename      ! next filename in the sequence

    ! set new next_filename ...

    !-----------------------------------------------------------------------
    !	... local variables
    !-----------------------------------------------------------------------
    integer :: pos, pos1, istat
    character(len=256) :: fn_new
    character(len=6)   :: seconds
    character(len=5)   :: num

    !-----------------------------------------------------------------------
    !	... ccm type filename
    !-----------------------------------------------------------------------
    pos = len_trim( filename )
    fn_new = filename(:pos)
    write(*,*) 'incr_flnm: old filename = ',trim(fn_new)
    if( fn_new(pos-2:) == '.nc' ) then
       pos = pos - 3
    end if
    istat = incstr( fn_new(:pos), 1 )
    if( istat /= 0 ) then
       write(*,*) 'incr_flnm: incstr returned ', istat
       write(*,*) '           while trying to decrement ',trim( fn_new )
       call endrun
    end if

    incr_filename = trim(fn_new)
    write(*,*) 'incr_flnm: new filename = ',incr_filename

  end function incr_filename

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine find_times( itms, fids, datatm, datatp, time )

    implicit none

    integer, intent(out) :: itms(2) ! record numbers that bracket time
    integer, intent(out) :: fids(2) ! ids of files that contains these recs
    real(r8), intent(in) :: time    ! time of interest
    real(r8), intent(out):: datatm, datatp

    integer np1        ! current forward time index of dataset
    integer n,i      ! 
    integer :: curr_tsize, next_tsize, all_tsize

    real(r8), allocatable, dimension(:):: all_data_times

    curr_tsize = size(curr_data_times)
    next_tsize = 0
    if ( associated(next_data_times)) next_tsize = size(next_data_times)

    all_tsize = curr_tsize + next_tsize

    allocate( all_data_times( all_tsize ) )

    all_data_times(:curr_tsize) = curr_data_times(:)
    if (next_tsize > 0) all_data_times(curr_tsize+1:all_tsize) = next_data_times(:)

    ! find bracketing times 
    do n=1, all_tsize-1
       np1 = n + 1
       datatm = all_data_times(n)
       datatp = all_data_times(np1)
       if ( (time .ge. datatm) .and. (time .le. datatp) ) then
          goto 20
       endif
    enddo

    write(6,*)'FIND_TIMES: Failed to find dates bracketing desired time =', time
    write(6,*)' datatm = ',datatm
    write(6,*)' datatp = ',datatp
    write(6,*)' all_data_times = ',all_data_times

    call endrun

20  continue

    deallocate( all_data_times )
  
    itms(1) = n
    itms(2) = np1
    fids(:) = curr_fileid
  
    do i=1,2
       if ( itms(i) > curr_tsize ) then 
          itms(i) = itms(i) - curr_tsize 
          fids(i) = next_fileid
       endif
    enddo

  end subroutine find_times

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine read_next_ps()
    implicit none

    integer :: recnos(2),fids(2)       
    integer :: cnt3(3)            ! array of counts for each dimension
    integer :: strt3(3)           ! array of starting indices

    !
    ! Set up hyperslab counters
    !
    strt3(1:2) = 1

    cnt3(1)  = plon
    cnt3(2)  = plat
    cnt3(3)  = 1

    call find_times( recnos, fids, datatimemn, datatimepn, next_mod_time )

    strt3(3) = recnos(1)

    call read_2d_met( fids(1), ps_id, met_psi_next(nm)%data, strt3, cnt3, 'dyn')

    strt3(3) = recnos(2)

    call read_2d_met( fids(2), ps_id, met_psi_next(np)%data, strt3, cnt3, 'dyn')

    write(6,*)'READ_NEXT_PS: Read meteorological data  '

  end subroutine read_next_ps

!------------------------------------------------------------------------
!------------------------------------------------------------------------
  subroutine read_next_metdata()
    implicit none

    integer recnos(2),fids(2),i      ! 
    integer cnt4(4)            ! array of counts for each dimension
    integer strt4(4)           ! array of starting indices
    integer cnt4s(4)            ! array of counts for each dimension
    integer cnt3(3)            ! array of counts for each dimension
    integer strt3(3)           ! array of starting indices

    call t_startf('MET__read_next_metdata')

    call find_times( recnos, fids, datatimem, datatimep, curr_mod_time )
    !
    ! Set up hyperslab corners
    !
    strt3(1:2) = 1
    strt4(1:3) = 1

    cnt3(1)  = plon
    cnt3(2)  = plat
    cnt3(3)  = 1

    cnt4(1)  = plon
    cnt4(2)  = plat
    cnt4(3)  = plev
    cnt4(4)  = 1

    cnt4s(1)  = plon
    cnt4s(2)  = plat-1
    cnt4s(3)  = plev
    cnt4s(4)  = 1

    do i=1,2

       strt4(4) = recnos(i)
       strt3(3) = recnos(i)

       call read_3d_met( fids(i), t_id, met_ti(i)%data, strt4, cnt4, 'phys')

       if (cell_wall_winds) then 
          call read_3d_met( fids(i), us_id, met_usi(i)%data, strt4, cnt4s, 'dyn'  , staggered=.true.)
          call read_3d_met( fids(i), vs_id, met_vsi(i)%data, strt4, cnt4, 'dyn' )
       else
          call read_3d_met( fids(i), u_id, met_ui(i)%data, strt4, cnt4, 'dyn' )
          call read_3d_met( fids(i), v_id, met_vi(i)%data, strt4, cnt4, 'dyn' )
       endif

       call read_3d_met( fids(i), q_id, met_qi(i)%data, strt4, cnt4, 'phys')

       call read_2d_met( fids(i), shflx_id,met_shflxi(i)%data,strt3, cnt3, 'phys')
       call read_2d_met( fids(i), qflx_id, met_qflxi(i)%data, strt3, cnt3, 'phys')
       call read_2d_met( fids(i), taux_id, met_tauxi(i)%data, strt3, cnt3, 'phys')
       call read_2d_met( fids(i), tauy_id, met_tauyi(i)%data, strt3, cnt3, 'phys')

    enddo

    write(6,*)'READ_NEXT_METDATA: Read meteorological data '

    call t_stopf('MET__read_next_metdata')

  end subroutine read_next_metdata


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine interpolate_metdata()

#if defined( SPMD )
    use mod_comm, only : mp_send4d_ns, mp_recv4d_ns
#endif

    implicit none

    real(r4) fact1, fact2
    real(r4) nfact1, nfact2
    real(r8) deltat,deltatn
    integer i,j,k

    call t_startf('MET__interpolate_metdata')

    deltat = datatimep - datatimem
    deltatn = datatimepn - datatimemn

    fact1 = (datatimep - curr_mod_time)/deltat
!    fact2 = (curr_mod_time - datatimem)/deltat
    fact2 = 1.-fact1

    nfact1 = (datatimepn - next_mod_time)/deltatn
!    nfact2 = (next_mod_time - datatimemn)/deltatn
    nfact2 = 1.-nfact1

    met_t(:,:,:) = fact1*met_ti(nm)%data(:,:,:) + fact2*met_ti(np)%data(:,:,:)
    met_q(:,:,:) = fact1*met_qi(nm)%data(:,:,:) + fact2*met_qi(np)%data(:,:,:)

    if (.not. online_test) where (met_q .lt. 0.0) met_q = 0.0


    met_ps_next(:,:) =  nfact1*met_psi_next(nm)%data(:,:) + nfact2*met_psi_next(np)%data(:,:)

    met_shflx(:,:) = fact1*met_shflxi(nm)%data(:,:) + fact2*met_shflxi(np)%data(:,:)
    met_qflx(:,:) = fact1*met_qflxi(nm)%data(:,:) + fact2*met_qflxi(np)%data(:,:)
    met_taux(:,:) = fact1*met_tauxi(nm)%data(:,:) + fact2*met_tauxi(np)%data(:,:)
    met_tauy(:,:) = fact1*met_tauyi(nm)%data(:,:) + fact2*met_tauyi(np)%data(:,:)

    if ( .not. cell_wall_winds ) then

       met_u(1:plon,beglat:endlat,beglev:endlev) = fact1*met_ui(nm)%data(1:plon,beglat:endlat,beglev:endlev) &
                                                 + fact2*met_ui(np)%data(1:plon,beglat:endlat,beglev:endlev)
       met_v(1:plon,beglat:endlat,beglev:endlev) = fact1*met_vi(nm)%data(1:plon,beglat:endlat,beglev:endlev) &
                                                 + fact2*met_vi(np)%data(1:plon,beglat:endlat,beglev:endlev)

       ! ghost u,v

#if defined( SPMD )
       call mp_send4d_ns( plon, plat, plev, 1, beglat, endlat,       &
            beglev, endlev, ng_d, ng_d, met_u )
       call mp_send4d_ns( plon, plat, plev, 1, beglat, endlat,       &
            beglev, endlev, ng_d, ng_s, met_v )
       call mp_recv4d_ns( plon, plat, plev, 1, beglat, endlat,       &
            beglev, endlev, ng_d, ng_d, met_u )
       call mp_recv4d_ns( plon, plat, plev, 1, beglat, endlat,       &
            beglev, endlev, ng_d, ng_s, met_v )
#endif

       ! average to cell walls (vorticity winds)
       call transfer_windsToWalls()
    else
       met_us(:,beglat:endlat,beglev:endlev) = fact1*met_usi(nm)%data(:,beglat:endlat,beglev:endlev) + &
            fact2*met_usi(np)%data(:,beglat:endlat,beglev:endlev)
       met_vs(:,beglat:endlat,beglev:endlev) = fact1*met_vsi(nm)%data(:,beglat:endlat,beglev:endlev) + &
            fact2*met_vsi(np)%data(:,beglat:endlat,beglev:endlev)

    endif

    ! ghost staggered u,v
#if defined( SPMD )
    call mp_send4d_ns( plon, plat, plev, 1, beglat, endlat,       &
         beglev, endlev, ng_s, ng_d, met_us )
    call mp_send4d_ns( plon, plat, plev, 1, beglat, endlat,       &
         beglev, endlev, ng_d, ng_s, met_vs )
    call mp_recv4d_ns( plon, plat, plev, 1, beglat, endlat,       &
         beglev, endlev, ng_s, ng_d, met_us )
    call mp_recv4d_ns( plon, plat, plev, 1, beglat, endlat,       &
         beglev, endlev, ng_d, ng_s, met_vs )
#endif

!    write(6,*)'INTERPOLATE_METDATA: complete.'

    call t_stopf('MET__interpolate_metdata')

  end subroutine interpolate_metdata

!------------------------------------------------------------------------
!------------------------------------------------------------------------
  subroutine read_3d_met(  fid, vid, loc_arr, strt, cnt, grid_map , staggered )

    implicit none

    integer, intent(in) :: fid, vid, strt(:), cnt(:)
    character(len=*), intent(in) :: grid_map
    real(r8),intent(out)  :: loc_arr(:,:,:)
    logical, optional, intent(in) :: staggered

    real(r8), allocatable :: global(:,:,:)
    integer :: i,j,k

    integer :: jstart,jend

    jstart = 1
    jend = cnt(2)

    if ( present( staggered ) ) then
       if ( staggered ) then
          jend=jend+1
          jstart=2
       endif
    endif

    allocate( global(cnt(1),1:jend,cnt(3)) )

    if (masterproc) then
       call wrap_get_vara_realx( fid, vid, strt, cnt, global(:,jstart:jend,:) )
    endif

    if (grid_map .eq. 'phys') then
       call scatter_3d_phys(global,loc_arr)
    else
       call scatter_3d_dyn(global,loc_arr)
    endif

    deallocate( global )

  end subroutine read_3d_met

!------------------------------------------------------------------------
!------------------------------------------------------------------------
  subroutine read_2d_met( fid, vid, loc_arr, strt, cnt, grid_map )

    implicit none

    integer, intent(in) :: fid, vid, strt(:), cnt(:)
    character(len=*), intent(in) :: grid_map

    real(r8) ,intent(out) :: loc_arr(:,:)

    real(r8), allocatable :: global(:,:)

    allocate( global(cnt(1),cnt(2)) )

    if (masterproc) then
       call wrap_get_vara_realx( fid, vid, strt, cnt, global )
    endif

    if (grid_map .eq. 'phys') then
       call scatter_2d_phys(global,loc_arr)
    else
       call scatter_2d_dyn(global,loc_arr)
    endif

    deallocate( global )

  end subroutine read_2d_met

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine scatter_2d_phys( glob_arr, local_arr )

    implicit none

    real(r8), intent(in) :: glob_arr(:,:)
    real(r8), intent(out) :: local_arr(pcols,begchunk:endchunk)

    call scatter_field_to_chunk(1,1,1,plon,glob_arr,local_arr)

  endsubroutine scatter_2d_phys

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine scatter_3d_phys( glob_arr, local_arr )

    implicit none

    real(r8), intent(in) :: glob_arr(:,:,:)
    real(r8), intent(out) :: local_arr(pcols,begchunk:endchunk,plev)

    local_arr(:,:,:) = -9999.

    call scatter_field_to_chunk(1,1,plev,plon,glob_arr,local_arr)

  end subroutine scatter_3d_phys

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine scatter_3d_dyn( glob_arr, local_arr )

    implicit none

    real(r8), intent(in) :: glob_arr(:,:,:)
    real(r8), intent(out) :: local_arr(plon, beglat:endlat, beglev:endlev)

#if ( defined SPMD )
    call scatter( glob_arr, strip3dxyz, local_arr, mpicom )
!    call parscatter( mpicom, 0, glob_arr, strip3dxyz, local_arr )
#else      
    local_arr(:,:,:) = glob_arr(:,:,:)
#endif

  endsubroutine scatter_3d_dyn

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine scatter_2d_dyn( glob_arr, local_arr )

    implicit none

    real(r8), intent(in) :: glob_arr(:,:)
    real(r8), intent(out) :: local_arr(plon, beglat:endlat)

#if ( defined SPMD )
    if (myid_z .eq. 0) then
       call scatter( glob_arr, strip2d, local_arr, comm_y )
!       call parscatter(  comm_y, 0, glob_arr, strip2d, local_arr )
    endif
    if (twod_decomp .eq. 1) then
       call parcollective2d( comm_z, BCSTOP, plon, endlat-beglat+1, local_arr ) 
    endif
#else      
    local_arr(:,:) = glob_arr(:,:)
#endif

  end subroutine scatter_2d_dyn

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine get_dimension( fid, dname, dsize )
    implicit none
    integer, intent(in) :: fid
    character(*), intent(in) :: dname
    integer, intent(out) :: dsize

    integer :: dimid

    if (masterproc) then 
       call wrap_inq_dimid( fid, dname, dimid )
       call wrap_inq_dimlen( fid, dimid, dsize )
    endif ! masterproc

#if (defined SPMD )
    call mpibcast( dsize, 1, mpiint, 0, mpicom )
#endif

  end subroutine get_dimension

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine open_met_datafile( fname, fileid, times, check_dims )

    use ioFileMod, only: getfil

    implicit none

    character(*), intent(in) :: fname
    integer, intent(inout) :: fileid
    real(r8), pointer :: times(:)
    logical, optional, intent(in) :: check_dims

    character(len=256) :: filen   
    integer :: year, month, day, dsize, i, timesize
    integer :: dateid,secid
    integer, allocatable , dimension(:) :: dates, datesecs

    !
    ! open file and get fileid
    !
    if (masterproc) then 
       call getfil( fname, filen, 0 )
       call wrap_open( filen, 0, fileid )
       write(6,*)'open_met_datafile: ',trim(filen)
    endif

    call get_dimension( fileid, 'time', timesize )

    if ( associated(times) ) deallocate(times)
    allocate( times(timesize) )

    if (masterproc) then

       allocate( dates(timesize) )
       allocate( datesecs(timesize) )

       call wrap_inq_varid( fileid, 'date',    dateid  )
       call wrap_inq_varid( fileid, 'datesec', secid  )

       call wrap_get_var_int( fileid, dateid, dates )
       call wrap_get_var_int( fileid, secid,  datesecs  )

       do i=1,timesize
          year = dates(i) / 10000
          month = mod(dates(i),10000)/100
          day = mod(dates(i),100)
          times(i) = get_time_float( year, month, day, datesecs(i) )
       enddo

       deallocate( dates )
       deallocate( datesecs )       

    endif ! masterproc

#if (defined SPMD )
    call mpibcast( times, timesize, mpir8, 0, mpicom )
#endif

!
! check that the data dim sizes match models dimensions
!
    if (present(check_dims)) then
       if (check_dims) then

          call get_dimension( fileid, 'lon', dsize )
          if (dsize /= plon) then
             write(6,*)'open_met_datafile: lonsiz=',dsize,' must = ',plon
             call endrun
          endif
          call get_dimension( fileid, 'lat', dsize )
          if (dsize /= plat) then
             write(6,*)'open_met_datafile: latsiz=',dsize,' must = ',plat
             call endrun
          endif
          call get_dimension( fileid, 'lev', dsize )
          if (dsize /= plev) then
             write(6,*)'open_met_datafile: levsiz=',dsize,' must = ',plev
             call endrun
          endif

       endif
    endif

  end subroutine open_met_datafile

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  function get_time_float( year, month, day, sec )

! returns float representation of time -- number of days
! since 1 jan 0001 00:00:00.000

    implicit none

    integer, intent(in) :: year, month, day
    integer, intent(in) :: sec
    real(r8) :: get_time_float

! ref date is 1 jan 0001

    integer :: refyr, refmn, refdy
    real(r8) :: refsc, fltdy
    integer :: doy(12)

!              jan feb mar apr may jun jul aug sep oct nov dec
!              31  28  31  30  31  30  31  31  31  31  30  31
    data doy /  1, 32, 60, 91,121,152,182,213,244,274,305,335 /

    refyr = 1
    refmn = 1
    refdy = 1
    refsc = 0._r8

    if ( calendar == 'GREGORIAN' ) then
       fltdy = greg2jday(year, month, day) - greg2jday(refyr,refmn,refdy)
    else ! assume no_leap (all years are 365 days)
       fltdy = (year - refyr)*365. + &
               (doy(month)-doy(refmn)) + &
               (day-refdy) 
    endif

    get_time_float = fltdy + ((sec-refsc)/86400._r8)

  endfunction get_time_float

!-----------------------------------------------------------------------
! 	... Return Julian day number given Gregorian date.
!
! Algorithm from Hatcher,D.A., Simple Formulae for Julian Day Numbers
! and Calendar Dates, Q.Jl.R.astr.Soc. (1984) v25, pp 53-55.
!-----------------------------------------------------------------------
  function greg2jday( year, month, day )

    implicit none

    integer, intent(in) :: year, month, day
    integer :: greg2jday

    !-----------------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------------
    integer :: ap, mp
    integer :: y, d, n, g

    !-----------------------------------------------------------------------
    !     	... Modify year and month numbers
    !-----------------------------------------------------------------------
    ap = year - (12 - month)/10
    mp = MOD( month-3,12 )
    if( mp < 0 ) then
       mp = mp + 12
    end if

    !-----------------------------------------------------------------------
    !     	... Julian day
    !-----------------------------------------------------------------------
    y = INT( 365.25*( ap + 4712 ) )
    d = INT( 30.6*mp + .5 )
    n = y + d + day  + 59
    g = INT( .75*INT( ap/100 + 49 ) ) - 38
    greg2jday = n - g

  end function greg2jday

end module metdata
