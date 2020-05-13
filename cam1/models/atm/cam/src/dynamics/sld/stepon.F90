#include <misc.h>
#include <params.h>

subroutine stepon
!-----------------------------------------------------------------------
!
! Purpose:
! Loop over time, calling driving routines for physics, dynamics, transport
!
! Original version:  CCM1
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use history, only: wshist, wrapup
   use pmgrid,  only: plev, plevp, beglat, endlat, masterproc, plat, plon
   use rgrid,   only: nlon
   use prognostics, only: u3, v3, t3, q3, n3, div, dpsl, dpsm, ps, omga, &
                          phis
   use commap,  only: clat
   use restart, only: write_restart
#if (defined COUP_CSM)
   use ccsm_msg, only: csmstop, ccsmfin
#endif
   use ppgrid   ,      only: begchunk, endchunk
   use physics_types,  only: physics_state, physics_tend
   use phys_buffer,    only: pbuf
   use dp_coupling,    only: d_p_coupling, p_d_coupling
   use time_manager, only: advance_timestep, get_step_size, get_nstep, &
                           is_first_step, is_first_restart_step, &
                           is_last_step, is_end_curr_day
   use scanslt,        only: slt_initial, slt_final, advection_state

   implicit none

#include <comctl.h>
#include <comhyb.h>
#include <comlun.h>

   type(physics_state), allocatable :: phys_state(:)
   type(physics_tend ), allocatable :: phys_tend(:)

   real(r8) gw     (plat)             ! Gaussian weights
   real(r8) detam  (plev)             ! intervals between vert full levs.
   real(r8) cwava  (plat)             ! weight applied to global integrals
   real(r8) etamid (plev)             ! vertical coords at midpoints 
   real(r8), allocatable :: t2(:,:,:) ! temp tendency
   real(r8), allocatable :: fu(:,:,:) ! u wind tendency
   real(r8), allocatable :: fv(:,:,:) ! v wind tendency
   real(r8) flx_net(plon,beglat:endlat)       ! net flux from physics
   real(r8) coslat(plon)              ! cosine of latitude
   real(r8) rcoslat(plon)             ! 1 over cosine of latitude
   real(r8) rpmid(plon,plev)          ! inverse of pressure at mid points
   real(r8) pdel(plon,plev)           ! Difference in pressure at level
   real(r8) pint(plon,plevp)          ! Interface pressure
   real(r8) pmid(plon,plev)           ! Mid-point pressure
   real(r8) ztodt                     ! time step
   real(r8) :: wcstart, wcend   ! wallclock timestamp at start, end of timestep
   real(r8) :: usrstart, usrend ! user timestamp at start, end of timestep
   real(r8) :: sysstart, sysend ! sys timestamp at start, end of timestep

   integer i,k,lat                    ! longitude,level,latitude indices
   
   type(advection_state) :: adv_state ! advection type
!
! Externals
!
   logical, external :: rstwr  ! whether or not to write restart files
!
!-----------------------------------------------------------------------
!
! Define eta coordinates on midpoints
!
   call t_startf ('stepon_startup')

   do k=1,plev
      etamid(k) = hyam(k) + hybm(k)
   end do

   if (is_first_step()) then
      do lat=beglat,endlat
         do i=1,nlon(lat)
            coslat(i) = cos(clat(lat))
            rcoslat(i) = 1./coslat(i)
         end do
!     
! Set current time pressure arrays for model levels etc.
!
         call plevs0(nlon(lat), plon, plev, ps(1,lat,n3), pint, pmid, pdel)
!
         do k=1,plev
            do i=1,nlon(lat)
               rpmid(i,k) = 1./pmid(i,k)
            end do
         end do
!
! Calculate vertical motion field
!
         call omcalc (rcoslat, div(1,1,lat,n3), u3(1,1,lat,n3), v3(1,1,lat,n3), dpsl(1,lat), &
                      dpsm(1,lat), pmid, pdel, rpmid, pint(1,plevp), &
                      omga(1,1,lat), nlon(lat))
      end do
   end if

   call slt_initial( detam, gw, cwava, etamid, adv_state )

   allocate(phys_state(begchunk:endchunk))
   allocate(phys_tend(begchunk:endchunk))
   allocate(t2(plon,plev,beglat:endlat))
   allocate(fu(plon,plev,beglat:endlat))
   allocate(fv(plon,plev,beglat:endlat))
!
! Beginning of basic time step loop
!
   call t_stopf ('stepon_startup')


   do

      if (masterproc .and. print_step_cost) then
         call t_stampf (wcstart, usrstart, sysstart)
      end if

      ztodt = get_step_size()

!
! Dump state variables to IC file
!
      call diag_dynvar_ic (phis, ps(1,beglat,n3), t3(1,1,beglat,n3), u3(1,1,beglat,n3), &
                           v3(1,1,beglat,n3), q3(1,1,1,beglat,n3) )
!
!----------------------------------------------------------
! PHYSPKG  Call the Physics package
!----------------------------------------------------------
!
      call t_startf('d_p_coupling')
      call d_p_coupling (ps(1,beglat,n3), t3(1,1,beglat,n3), u3(1,1,beglat,n3), &
                         v3(1,1,beglat,n3), q3(1,1,1,beglat,n3), &
                         omga, phis, phys_state, phys_tend, pbuf)
      call t_stopf('d_p_coupling')

      call t_startf('phys_driver')
      if (ideal_phys) then
         call phys_idealized(phys_state, phys_tend, ztodt, etamid)
      else if (adiabatic) then
         call phys_adiabatic(phys_state, phys_tend)
      else
         call physpkg (phys_state, gw, ztodt, phys_tend, pbuf)
      end if
      call t_stopf('phys_driver')

      call t_startf('p_d_coupling')
      call p_d_coupling (phys_state, phys_tend, t2, fu, fv, flx_net,q3(1,1,1,beglat,n3))
      call t_stopf('p_d_coupling')

      if (is_first_step() .or. is_first_restart_step()) then
         call print_memusage ('stepon after p_d_coupling')
      end if
!----------------------------------------------------------
! DYNPKG Call the Dynamics Package
!----------------------------------------------------------

      call t_startf ('dynpkg')
      call dynpkg(t2      ,fu      ,fv      ,etamid  , &
                  cwava   ,detam   ,flx_net , ztodt  ,adv_state)
      call t_stopf ('dynpkg')

      if (is_first_restart_step()) then
         call print_memusage ('stepon after dynpkg')
      end if

      call t_startf('stepon_single')

! Set end of run flag.

#if ( ! defined COUP_CSM )
      if (is_last_step()) nlend = .true.
#else
      if (csmstop) then
         if ( masterproc ) write(6,*)'atm: Stopping at the end of this day'
         if (is_end_curr_day()) nlend = .true.
      end if
#endif
!
!----------------------------------------------------------
! History and restart logic
!----------------------------------------------------------
!
! Write and/or dispose history tapes if required
!
      call t_startf ('wshist')
      call wshist ()
      call t_stopf ('wshist')
!
! Write restart file
!
      if (rstwr() .and. nrefrq /= 0) then
         call t_startf ('write_restart')
         call write_restart
         call t_stopf ('write_restart')
      end if
!
! Dispose necessary files
!
      call t_startf ('wrapup')
      call wrapup
      call t_stopf ('wrapup')
   
      if (masterproc .and. print_step_cost) then
         call t_stampf (wcend, usrend, sysend)
         write(6,'(a,3f8.3,a)')'Prv timestep wallclock, usr, sys=', &
                               wcend-wcstart, usrend-usrstart, sysend-sysstart, ' seconds'
      end if
!
! Increment nstep before returning to top of loop
!
      call advance_timestep()

      call t_stopf('stepon_single')
!
! Check for end of run
!
      if (nlend) then
         call print_memusage ('End stepon')
         deallocate(phys_state)
         deallocate(phys_tend)
         deallocate(t2)
         deallocate(fu)
         deallocate(fv)
#ifdef COUP_CSM
         call ccsmfin
#endif
         return
      end if

   end do  ! End of timestep loop
   call slt_final( adv_state )

end subroutine stepon
