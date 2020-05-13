#include <misc.h>
#include <params.h>

subroutine inti ()
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Set constants and call initialization procedures for time independent
! physics routines
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Rosinski
! 
!-----------------------------------------------------------------------
   use shr_kind_mod,       only: r8 => shr_kind_r8
   use pmgrid,             only: plev, plevp
   use chemistry,          only: chem_init
   use ghg_defaults,       only: ghg_defaults_init
   use ozone_data,         only: ozone_data_init
   use tracers,            only: tracers_init
   use aerosol_intr,       only: prognostic_aerosol_initialize
   use gw_drag,            only: gw_inti
   use vertical_diffusion, only: vertical_diffusion_init
   use convect_shallow,    only: convect_shallow_init
   use cloud_fraction,     only: cldfrc_init
   use stratiform,         only: stratiform_init
   use param_cldoptics,    only: param_cldoptics_init
   use shr_const_mod,      only: shr_const_zvir, shr_const_cpwv, shr_const_rwv
   use physconst,          only: rair, cpair, cpwv, gravit, stebol, epsilo, tmelt, &
                                 latvap, latice, rh2o, zvir, cpvir, rhoh2o, pstd,  &
                                 karman, rhodair, mwo3, mwdry, avogad, boltz
   use prescribed_aerosols,only: aerosol_initialize
   use aer_optics,         only: aer_optics_initialize
   use check_energy,       only: check_energy_init
   use convect_deep,       only: convect_deep_init
   use cloudsimulator,     only: doisccp, cloudsimulator_init
   use radiation,          only: radiation_init
   use radheat,            only: radheat_init
   use diagnostics,        only: diag_init

#if ( defined WACCM_GHG || defined WACCM_MOZART )
   use ctem,                    only: ctem_inti
   use ion_drag_calc,           only: iondragi
#endif
#if ( defined OFFLINE_DYN )
   use metdata, only:  init_met
#endif
 
   implicit none

#include <comctl.h>
#include <comhyb.h>

!
! Initialize radiation data for aerosol forcing calculation
! Initialize aerosol fields from files
!
   call aer_optics_initialize()

   call aerosol_initialize()
   
   call prognostic_aerosol_initialize()
!
!-----------------------------------------------------------------------
!
! Initialize physconst variables
! In adiabatic case, set zvir and cpvir explicitly to zero instead of 
! computing as (rh2o/rair - 1.) and (cpwv/cpair - 1.) respectively, in order 
! to guarantee an identical zero.
!
   if (adiabatic .or. ideal_phys) then
      rh2o  = rair
      zvir  = 0.
      cpwv  = cpair
      cpvir = 0.
   else
      rh2o  = shr_const_rwv
      zvir  = shr_const_zvir
      cpwv  = shr_const_cpwv
      cpvir = cpwv/cpair - 1.
   end if
!
! Call time independent initialization routines for parameterizations.
!
   ! Prognostic chemistry.
   call chem_init()

   ! Default distributions for CH4, N2O, CFC11 and CFC12.
   call ghg_defaults_init()

   ! Initialize reading of ozone dataset.
   call ozone_data_init(mwo3)

   call tracers_init

   call gw_inti (cpair   ,cpwv    ,gravit  ,rair    ,hypi    )

   call vertical_diffusion_init

   call tsinti  (tmelt   ,latvap  ,rair    ,stebol  ,latice  )

   call radiation_init

   call radheat_init(hypm)

   call esinti  (epsilo  ,latvap  ,latice  ,rh2o    ,cpair  , &
                 tmelt   )
   call convect_shallow_init(hypi)

   call cldfrc_init

   call convect_deep_init(hypi)

   call cldinti ()

   call stratiform_init

   call param_cldoptics_init

   call check_energy_init

   if (doisccp) call cloudsimulator_init

   call diag_init()

#if ( defined WACCM_GHG || defined WACCM_MOZART )
   call iondragi (cpair,hypm)

   call ctem_inti
#endif

#if ( defined OFFLINE_DYN )
!
! Initialize the offline meteorology data
!
   call init_met()
#endif
 
   return
end subroutine inti
