#include <misc.h>
#include <params.h>

subroutine physpkg(phys_state, gw, ztodt, phys_tend, pbuf)


!----------------------------------------------------------------------- 
! 
! Purpose: 
! Loop over time, calling driving routines for physics
! 
! Method: 
! COUP_CSM and must be checked in order to invoke the proper calling
! sequence for running the CSM model
! 
! Author: 
! Original version:  CCM3
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plon, plat, masterproc
   use ppgrid, only: pcols
   use buffer, only: pblht, tpert, qpert
   use check_energy, only: check_energy_gmean
   use comsrf
   use comsrfdiag
#ifdef COUP_CSM
   use ccsm_msg, only: ccsmave, dorecv, dosend, ccsmsnd, ccsmrcv
#else
   use atm_lndMod, only: atmlnd_drv
#endif
#ifdef SPMD
   use mpishorthand
#endif
   use phys_buffer,    only: pbuf_size_max, pbuf_fld, pbuf_allocate, pbuf_deallocate, &
                             pbuf_update_tim_idx
   use phys_grid,      only: get_ncols_p, get_lat_all_p, get_lon_all_p
   use physics_types,  only: physics_state, physics_tend
   use diagnostics,    only: diag_allocate, diag_deallocate, diag_surf, &
                             diag_physvar_ic
   use time_manager,   only: get_nstep, is_first_step, is_first_restart_step, &
                             is_end_curr_month, get_curr_date, is_end_curr_day
   use physconst,      only: stebol, latvap

#if ( defined WACCM_MOZART )
   use mo_lightning,   only: lightning_no_prod
#endif

#if (!defined COUP_CSM)
   use ice_constants, only: TfrezK
#endif
#if (defined BFB_CAM_SCAM_IOP )
   use history, only: outfld
#endif
#if (!defined COUP_CSM)
#if (!defined COUP_SOM)
   use sst_data, only: sstint
   use ice_data, only: iceint
#endif
#endif
#if ( defined SCAM )
#include <max.h>
   use scamMod, only: switch,use_srfprop,lhflxobs,shflxobs,tground,have_lhflx,have_shflx,have_tg
#endif
#if ( defined OFFLINE_DYN )
   use metdata, only: get_met_fields
#endif
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comsol.h>
!-----------------------------------------------------------------------
!
! Input arguments
!
   real(r8), intent(in) :: gw(plat)                    ! Gaussian weights
   real(r8), intent(in) :: ztodt                       ! physics time step unless nstep=0
!
! Input/Output arguments
!
   type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
   type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
   type(pbuf_fld),      intent(inout), dimension(pbuf_size_max)     :: pbuf
!
!---------------------------Local workspace-----------------------------
!
   integer :: i,c,k                             ! indices
   integer :: ncol                              ! number of columns
   integer :: nstep                             ! current timestep number
   integer :: yr, mon, day                      ! year, month, and day components of a date

   real(r8) fsds(pcols,begchunk:endchunk)        ! Surface solar down flux
   real(r8) :: tmp(pcols,begchunk:endchunk)
!-----------------------------------------------------------------------

   call t_startf ('physpkg_st')
   nstep = get_nstep()

   call pbuf_allocate('physpkg')
   call diag_allocate()

! Compute total energy of input state and previous output state
   call t_startf ('chk_en_gmean')
   call check_energy_gmean(phys_state, pbuf, ztodt, nstep)
   call t_stopf ('chk_en_gmean')

!-----------------------------------------------------------------------
! Advance time information
!-----------------------------------------------------------------------

   call advnce( phys_state )
   call t_stopf ('physpkg_st')

#ifdef TRACER_CHECK
   call gavglook ('before tphysbc DRY', phys_state, gw)
#endif


!-----------------------------------------------------------------------
! Tendency physics before flux coupler invocation
!-----------------------------------------------------------------------
!

#if (defined BFB_CAM_SCAM_IOP )
   do c=begchunk, endchunk
      call outfld('Tg',srfflx_state2d(c)%ts,pcols   ,c     )
   end do
#endif
   call t_startf ('bc_physics')

!$OMP PARALLEL DO PRIVATE (C)

   do c=begchunk, endchunk
!
! Output physics terms to IC file
!
      call diag_physvar_ic ( c, pbuf )

      call t_startf ('tphysbc')
      call tphysbc (ztodt, pblht(1,c), tpert(1,c), qpert(1,1,c), snowhland(1,c),      &
                    fsns(1,c), fsnt(1,c), flns(1,c), flnt(1,c), phys_state(c),        &
                    phys_tend(c), pbuf,  fsds(1,c), landm(1,c),           &
                    landfrac(1,c), ocnfrac(1,c), icefrac(1,c), surface_state2d(c), srfflx_state2d(c))

      call t_stopf ('tphysbc')

      if (dosw .or. dolw) then
	call output_flns_fsns_fluxes(surface_state2d(c),c)
      end if	

#if ( ! defined COUP_CSM )
!
! zero surface fluxes at beginning of each time step.  Land Ocean and Ice
! processes will will write into process specific flux variables
! at the end of the time step these separate fluxes will be combined over the
! entire grid
!
      call srfflx_state_reset (srfflx_state2d(c))
#endif

   end do
   call t_stopf ('bc_physics')

#if ( defined SCAM )
   ! Don't call the rest in CRM mode
   if(switch(CRM_SW+1)) return
#endif

#ifdef TRACER_CHECK
   call gavglook ('between DRY', phys_state, gw)
#endif

#if ( ! defined COUP_CSM )
!
!-----------------------------------------------------------------------
! Determine surface quantities - no flux coupler
!-----------------------------------------------------------------------
!
   if (.not. aqua_planet) then
!
! Call land model driving routine
!
#ifdef TIMING_BARRIERS
      call t_startf ('sync_tphysbc_lnd')
      call mpibarrier (mpicom)
      call t_stopf ('sync_tphysbc_lnd')
#endif
      call t_startf ('atmlnd_drv')

#if ( defined SCAM )
       if (landfrac(1,begchunk).gt.0) &
#endif

      call atmlnd_drv(nstep, iradsw, eccen, obliqr, lambm0,&
                      mvelpp,surface_state2d,srfflx_parm2d)

      call t_stopf ('atmlnd_drv')
#ifdef TIMING_BARRIERS
      call t_startf ('sync_after_lnd')
      call mpibarrier (mpicom)
      call t_stopf ('sync_after_lnd')
#endif
!
! save off albedos and longwave for som offline vars
!
!$OMP PARALLEL DO PRIVATE (C,NCOL,I)
      do c=begchunk,endchunk
         ncol = get_ncols_p(c)
         do i=1,ncol
            if (landfrac(i,c) > 0.) then
               asdirlnd(i,c) = srfflx_parm2d(c)%asdir(i)
               asdiflnd(i,c) = srfflx_parm2d(c)%asdif(i)
               aldirlnd(i,c) = srfflx_parm2d(c)%aldir(i)
               aldiflnd(i,c) = srfflx_parm2d(c)%aldif(i)
               lwuplnd(i,c)  = srfflx_parm2d(c)%lwup(i)
            else
               asdirlnd(i,c) = 0. 
               asdiflnd(i,c) = 0. 
               aldirlnd(i,c) = 0. 
               aldiflnd(i,c) = 0. 
               lwuplnd(i,c)  = 0. 
            end if
         end do
!
!output shf/lhf fluxes for land model
!
         call output_shf_lhf_fluxes(srfflx_parm2d(c), c, ncol, landfrac(1,c), 'LND')
         call update_srf_fluxes (srfflx_state2d(c), srfflx_parm2d(c), landfrac(1,c), ncol)
      end do
   end if                    ! end of not aqua_planet if block

#if (defined COUP_SOM)
!
! Set ocean surface quantities - ocn model internal to atm
!
   if (is_end_curr_day ()) then
      call print_coverage ('icefrac', ' million km^2', icefrac, 1.d-12)
      do c=begchunk,endchunk
         ncol = get_ncols_p(c)
         do i=1,ncol
            tmp(i,c) = icefrac(i,c)*sicthk(i,c)
         end do
      end do
      call print_coverage ('icevol ', ' 10^13m^3', tmp, 1.d-13)

      do c=begchunk,endchunk
         ncol = get_ncols_p(c)
         do i=1,ncol
            tmp(i,c) = icefrac(i,c)*snowhice(i,c)
         end do
      end do
      call print_coverage ('snowvol', ' 10^13m^3', tmp, 1.d-13)
   end if

   call t_startf ('somint')
   call somint ()
   call t_stopf ('somint')

   call t_startf ('somoce')
   call somoce (surface_state2d, srfflx_parm2d_ocn)
   call t_stopf ('somoce')

#else

   call t_startf ('sstint')
   call sstint ()
   call t_stopf ('sstint')
!
! iceint may change ocean fraction, so call it before camoce
!
   call t_startf ('iceint')
   call iceint ()
   call t_stopf ('iceint')

   call t_startf ('camoce')
   call camoce (surface_state2d, srfflx_parm2d_ocn)
   call t_stopf ('camoce')
#endif
!
! Set ice surface quantities - icn model internal to atm
!
   call t_startf('camice')
   call camice (surface_state2d, srfflx_parm2d)
   call t_stopf('camice')
!
! output shf/lhf fluxes for ice/ocn/som_offline 
!
!$OMP PARALLEL DO PRIVATE (C, NCOL, I)
   do c=begchunk,endchunk
      ncol = get_ncols_p(c)
      do i=1,ncol
         if(icefrac(i,c) > 0.) then
            tsice_rad(i,c) = sqrt(sqrt(srfflx_parm2d(c)%lwup(i)/stebol))
         else
            tsice_rad(i,c) = TfrezK
         endif
      end do
      call output_shf_lhf_fluxes (srfflx_parm2d(c), c, ncol, icefrac(1,c), 'ICE')
      call output_shf_lhf_fluxes (srfflx_parm2d_ocn(c), c, ncol, ocnfrac(1,c), 'OCN')
      call output_shfoi_lhfoi_fluxes (srfflx_parm2d_ocn(c), srfflx_parm2d(c), c)

!JR SOM case: Have to wait to call update routine till after both ocean and ice have
!JR operated, since the fractions can change internal to the parameterization
      do i = 1, ncol
         srfflx_state2d(c)%sst(i) = srfflx_parm2d_ocn(c)%ts(i)
      enddo
      call update_srf_fluxes (srfflx_state2d(c), srfflx_parm2d_ocn(c), ocnfrac(1,c), ncol)
      call update_srf_fluxes (srfflx_state2d(c), srfflx_parm2d(c), icefrac(1,c), ncol)
   end do
#endif

#if ( defined COUP_CSM )
!
!-----------------------------------------------------------------------
! Determine surface quantities using csm flux coupler
!-----------------------------------------------------------------------
!
! If send data to flux coupler only on radiation time steps:
!
   if (flxave) then
!
! Average the precipitation input to lsm between radiation calls.
!
      call ccsmave(iradsw, nstep, dosw)
!
! Use solar radiation flag to determine data exchange steps 
! with flux coupler. This processes are not independent since 
! instantaneous radiative fluxes are passed, valid over the 
! interval to the next radiation calculation. The same 
! considerations apply to the long and shortwave fluxes, so 
! the intervals must be the same. Data is received from the 
! coupler one step after it is sent.
!
      if (nstep == 0) then
         dorecv = .true.
         dosend = .true.
      else if (nstep == 1) then
         dorecv = .false.
         dosend = .false.
      else if ( (nstep == 2) .and. (iradsw == 1) ) then
         dorecv = .true.
         dosend = dosw
      else
         dorecv = dosend
         dosend = dosw
      end if
   endif
!
! If send data to flux coupler on every time step
!
   if (.not. flxave) then
      if (nstep /= 1) then
         dorecv = .true.
         dosend = .true.
      else 
         dorecv = .false.
         dosend = .false.
      endif
   endif
!
! Send/recv data to/from the csm flux coupler.
!
   if (dosend) call ccsmsnd ( )
   if (dorecv) call ccsmrcv ( )
#endif
#if ( defined OFFLINE_DYN )
!
! if offline mode set SHFLX QFLX TAUX TAUY for vert diffusion
!
   call get_met_fields( srfflx_state2d )
#endif
!
!-----------------------------------------------------------------------
! Tendency physics after coupler 
! Not necessary at terminal timestep.
!-----------------------------------------------------------------------
!
#if ( defined WACCM_MOZART )
! Set lightning no production
   call lightning_no_prod( phys_state, pbuf )
#endif
   call t_startf ('ac_physics')

!$OMP PARALLEL DO PRIVATE (C, NCOL)

   do c=begchunk,endchunk
      ncol = get_ncols_p(c)
!
! replace surface fluxes with observed values for IOP forcing if
! requested by switch settings in the GUI
!
#if ( defined SCAM )
      if (use_srfprop) then
         if(have_lhflx) then
            srfflx_state2d(c)%lhf(1) = lhflxobs(1)
            srfflx_state2d(c)%cflx(1,1) = lhflxobs(1)/latvap
         endif
         if(have_shflx) srfflx_state2d(c)%shf(1) = shflxobs(1)
         if(have_tg) then
            srfflx_state2d(c)%ts(1) = tground(1)
            srfflx_state2d(c)%lwup(1) = stebol * tground(1)**4
            do k=1,plevmx
               surface_state2d(c)%tssub(1,k) = tground(1)
            end do
         endif
      endif
#endif
!
! surface diagnostics for history files
!
      call diag_surf (srfflx_state2d(c), surface_state2d(c), icefrac(1,c), ocnfrac(1,c), landfrac(1,c), &
                      sicthk(1,c), snowhland(1,c), snowhice(1,c), tsice(1,c), trefmxav(1,c), &
                      trefmnav(1,c) )

      call t_startf ('tphysac')

      call tphysac (ztodt, pblht(1,c), qpert(1,1,c), tpert(1,c), srfflx_state2d(c)%shf,        &
                    srfflx_state2d(c)%wsx,srfflx_state2d(c)%wsy, srfflx_state2d(c)%cflx, &
                    sgh(1,c), sgh30(1,c), srfflx_state2d(c)%lhf, landfrac(1,c), snowhland(1,c), &
                    srfflx_state2d(c)%tref, surface_state2d(c)%precc, surface_state2d(c)%precl,    &
                    surface_state2d(c)%precsc, surface_state2d(c)%precsl, phys_state(c), phys_tend(c), pbuf, &
                    ocnfrac(1,c), fsds(1,c), icefrac(1,c), fv(1,c), ram1(1,c))
      call t_stopf ('tphysac')
   end do                    ! Chunk loop

   call t_stopf('ac_physics')

#ifdef TRACER_CHECK
   call gavglook ('after tphysac FV:WET)', phys_state, gw )
#endif

   call pbuf_deallocate('physpkg')
   call pbuf_update_tim_idx()
   call diag_deallocate()

end subroutine physpkg



subroutine gavglook (title, state, gw)

  !
  ! process info from state vectors. Only useful when data in all chunks are in sync
  ! e.g. before and after tphysac and tphysbc
  !

  use physics_types,  only: physics_state, physics_tend
!  use comsrf,         only: srfflx_state
  use phys_grid,      only: get_ncols_p, get_lat_all_p, get_lon_all_p
  use ppgrid,         only: begchunk, endchunk
  use ppgrid,         only: pcols, pver
  use constituents,   only: ppcnst, cnst_name, cnst_need_pdeldry
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use pmgrid,         only: plon, plat, masterproc
  use physconst,      only: gravit
  use time_manager,   only: dtime
#ifdef SPMD
  use mpishorthand
#endif

  implicit none

  ! arguments
  character(len=*), intent(in) :: title
  type(physics_state), intent(in), dimension(begchunk:endchunk) :: state
  real(r8), intent(in) :: gw(plat)                    ! Gaussian weights
  
  ! local
  integer i, lat, c, lon, k
  integer :: lats(pcols)                       ! array of latitude indices
  integer :: lons(pcols)                       ! array of longitude indices
  integer m
  integer :: ncol                              ! number of columns
  real(r8) twodfld(plon,plat,ppcnst)          ! summed at each grid point
  real(r8) twodfle(plon,plat,ppcnst)          ! summed at each grid point
  real(r8) twodflx(plon,plat,ppcnst)          ! summed at each grid point
  real(r8) twodfly(plon,plat,ppcnst)          ! summed at each grid point
#ifdef SPMD                                     
  real(r8) :: twodfld_glob(plon,plat,ppcnst)   ! global summed at each grid point
  real(r8) :: twodfle_glob(plon,plat,ppcnst)   ! global summed at each grid point
  real(r8) :: twodflx_glob(plon,plat,ppcnst)   ! global summed at each grid point
  real(r8) :: twodfly_glob(plon,plat,ppcnst)   ! global summed at each grid point
#endif                                          
  real(r8) :: zonal(plat), zonalw(plat)              !  summed along each latitude
  real(r8) gavg, gavgw
  real(r8) col, wmin, wmax, colw

  !!--------------------------------------------------------


  ! operations on each processor
  twodfld(:,:,:) = 0.
  twodfle(:,:,:) = 0.
  twodflx(:,:,:) = 0.
  twodfly(:,:,:) = 0.
  do c=begchunk, endchunk
     ncol = get_ncols_p(c)
     call get_lat_all_p(c, ncol, lats)
     call get_lon_all_p(c, ncol, lons)
     do m = 1,ppcnst
        do i=1,ncol
           lat = lats(i)
           lon = lons(i)
           col = 0.
           colw = 0.
!           fluxcol = 0.
           wmax = -1.e36
           wmin = 1.e36
           do k = 1,pver
              if ( cnst_need_pdeldry)  col  = col + state(c)%pdeldry(i,k)*state(c)%q(i,k,m)*gw(lats(i))
              colw = colw + state(c)%pdel(i,k)  *state(c)%q(i,k,m)*gw(lats(i))
              wmax = max(wmax,state(c)%q(i,k,m))
              wmin = min(wmin,state(c)%q(i,k,m))
           end do ! k
           if ( cnst_need_pdeldry)  col = col/gravit
           colw = colw/gravit
           if ( cnst_need_pdeldry) twodfld(lons(i),lats(i),m) = twodfld(lons(i),lats(i),m) + col
           twodfle(lons(i),lats(i),m) = twodfle(lons(i),lats(i),m) + colw
           twodflx(lons(i),lats(i),m) = twodflx(lons(i),lats(i),m) + wmin
           twodfly(lons(i),lats(i),m) = twodfly(lons(i),lats(i),m) + wmax
        enddo ! i
     enddo ! m
  end do ! c

  ! move data to masterproc
#ifdef SPMD
#ifdef TIMING_BARRIERS
  call mpibarrier (mpicom)
#endif
  if ( cnst_need_pdeldry ) call mpisum(twodfld, twodfld_glob, plon*plat*ppcnst, mpir8, 0, mpicom)
  call mpisum(twodfle, twodfle_glob, plon*plat*ppcnst, mpir8, 0, mpicom)
  call mpisum(twodflx, twodflx_glob, plon*plat*ppcnst, mpir8, 0, mpicom)
  call mpisum(twodfly, twodfly_glob, plon*plat*ppcnst, mpir8, 0, mpicom)
  if (masterproc) then
     if ( cnst_need_pdeldry ) twodfld(:,:,:) = twodfld_glob(:,:,:) 
     twodfle(:,:,:) = twodfle_glob(:,:,:) 
     twodflx(:,:,:) = twodflx_glob(:,:,:) 
     twodfly(:,:,:) = twodfly_glob(:,:,:) 
  endif
#endif

  ! process the data
  if (masterproc) then
     do m = 1,ppcnst
        wmax = -1.e36
        wmin = 1.e36
        do lat=1,plat
           if ( cnst_need_pdeldry ) zonal(lat) = 0.
           zonalw(lat) = 0.
           do i=1,plon
              if ( cnst_need_pdeldry ) zonal(lat) = zonal(lat) + twodfld(i,lat,m)
              zonalw(lat) = zonalw(lat) + twodfle(i,lat,m)
              wmax = max(wmax,twodfly(i,lat,m))
              wmin = min(wmin,twodflx(i,lat,m))
           end do
        end do
        if ( cnst_need_pdeldry ) gavg = 0.
        gavgw = 0.
        do lat=1,plat
           if ( cnst_need_pdeldry ) gavg = gavg + zonal(lat)
           gavgw = gavgw + zonalw(lat)
        end do
        if ( cnst_need_pdeldry )  gavg = gavg/(2.*plon)
        gavgw = gavgw/(2.*plon)

        if ( .not. cnst_need_pdeldry ) then
             write (6,67) trim(title)//' m=',m,'name='//trim(cnst_name(m))//' gavg wet, min, max ' &
                  , gavgw,wmin,wmax
67           format (a24,i2,a36,1p,4g25.14)
          else
             write (6,66) trim(title)//' m=',m,'name='//trim(cnst_name(m))//' gavg dry, wet, min, max ' &
                  , gavg, gavgw,wmin,wmax
!66           format (a24,i2,a36,1p,4g25.14)
66           format (a24,i2,a36,1p,4e25.13)
     endif

     end do
  endif

end subroutine gavglook
