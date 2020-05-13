#include <misc.h>
#include <params.h>

subroutine advnce( phys_state )
!-----------------------------------------------------------------------
!
! Purpose: 
! Advance time information
!
! Method: 
!
! Author: CCM1, CMS Contact: J. Truesdale
!
!-----------------------------------------------------------------------
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use chemistry,           only: chem_timestep_init
  use chem_surfvals,       only: chem_surfvals_set
  use ppgrid,              only: begchunk, endchunk
  use physics_types,       only: physics_state
  use aerosol_intr,        only: aerosol_time_interp
  use ozone_data,          only: ozone_data_timestep_init
  use tracers,             only: tracers_timestep_init
  use time_manager,        only: get_nstep
  use volcanicmass,        only: read_volcanic_mass
  use prescribed_aerosols, only: aerint, strat_volcanic
  use ramp_scon,           only: ramp_sconst
  use vertical_diffusion,  only: vertical_diffusion_ts_init
  use radheat,             only: radheat_timestep_init

  implicit none

#include <comctl.h>
#include <comlun.h>
!-----------------------------------------------------------------------
! dummy arguments
!-----------------------------------------------------------------------
   type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
!-----------------------------------------------------------------------
!
! Local workspace
!
  integer :: nstep             ! current timestep number
!-----------------------------------------------------------------------

  nstep = get_nstep()
!
! Determine whether it is time for a shortwave or longwave radiation 
! calculation
!
  dosw = nstep.eq.0 .or. iradsw.eq.1 .or. (mod(nstep-1,iradsw).eq.0 .and. nstep.ne.1) .or. nstep <= irad_always
  dolw = nstep.eq.0 .or. iradlw.eq.1 .or. (mod(nstep-1,iradlw).eq.0 .and. nstep.ne.1) .or. nstep <= irad_always
!
! Determine whether it is time for an absorptivity/emissivity calculation
!
  doabsems = nstep.eq.0 .or. iradae.eq.1 .or. (mod(nstep-1,iradae).eq.0 .and. nstep.ne.1)
  aeres = (mod(nstep,iradae).ne.0)
!
! Update aerosol data on shortwave or longwave time step. 
!
  if (dosw .or. dolw) then
     call aerint ()
     if(strat_volcanic) then
       call read_volcanic_mass
     end if
  end if
!
! Ramping ghg if appropriate
!
  call chem_surfvals_set( phys_state )
!
! Ramp solar constant if appropraite
!
  if (doRamp_scon) call ramp_sconst
!
! Time interpolate for chemistry.
!
  call chem_timestep_init

! Time interpolate ozone data
  call ozone_data_timestep_init()

  ! Upper atmosphere radiative processes
  call radheat_timestep_init
 
! Time interpolate for vertical diffusion upper boundary condition
  call vertical_diffusion_ts_init

! get the aerosol surface fluxes for this time step
  call aerosol_time_interp
!
! Time interpolate for tracers, if appropriate
!
  call tracers_timestep_init(phys_state)


  return
end subroutine advnce
