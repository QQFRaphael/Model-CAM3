#include <misc.h>
subroutine initindx
!----------------------------------------------------------------------- 
! 
! Purpose: Register constituents and physics buffer fields.
! 
! Author:    CSM Contact: M. Vertenstein, Aug. 1997
!            B.A. Boville, Oct 2001
! 
!-----------------------------------------------------------------------
  use shr_kind_mod,       only: r8 => shr_kind_r8
  use constituents,       only: pcnst, ppcnst, cnst_add, advected, cnst_chk_dim, cnst_name, &
                                dcconnam, sflxnam, tendnam, tottnam, fixcnam, hadvnam, vadvnam
  use phys_buffer,        only: pbuf_init
  use chemistry,          only: chem_register
  use stratiform,         only: stratiform_register
  use physconst,          only: mwdry, cpair, mwh2o, cph2o
  use tracers,            only: tracers_register
  use check_energy,       only: check_energy_register
  use aerosol_intr,       only: aerosol_register_cnst
  use vertical_diffusion, only: vd_register
  use convect_deep,       only: convect_deep_register
  use convect_shallow,    only: convect_shallow_register
  use radiation,          only: radiation_register

#if ( defined BFB_CAM_SCAM_IOP )
  use iop
  use string_utils, only: to_lower
#endif

  implicit none
!-----------------------------------------------------------------------
#include <comctl.h>
!---------------------------Local variables-----------------------------
!
  integer m            ! loop index
  integer mm           ! constituent index 
!-----------------------------------------------------------------------

  ! Initialize physics buffer
  call pbuf_init()

  ! Register water vapor.
  ! ***** N.B. ***** This must be the first call to cnst_add so that
  !                  water vapor is constituent 1.
  call cnst_add('Q', advected, mwh2o, cph2o, 1.E-12_r8, mm, &
                longname='Specific humidity', readiv=.true.)

  ! cloud water
  call stratiform_register()

  ! chemical constituents
  call chem_register()

  ! aerosols
  call aerosol_register_cnst()

  ! test tracers
  call tracers_register()

  ! All tracers registered, check that the dimensions are correct
  call cnst_chk_dim()

  ! ***NOTE*** No registering constituents below this point.  Only registration
  !            of physics buffer fields.

  ! Set default names for non-water advected and non-advected tracers
  ! Set names of advected and non-advected tracer diagnostics
  do m=1,ppcnst
     sflxnam(m)  = 'SF'//cnst_name(m)
  end do
  do m=1,pcnst
     hadvnam(m)  = 'HA'//cnst_name(m)
     vadvnam(m)  = 'VA'//cnst_name(m)
     fixcnam(m)  = 'DF'//cnst_name(m)
     tendnam(m)  = 'TE'//cnst_name(m)
     tottnam(m)  = 'TA'//cnst_name(m)
  end do


  call convect_deep_register

  call convect_shallow_register

  call check_energy_register

  call radiation_register

  ! vertical diffusion
  call vd_register()

#if ( defined BFB_CAM_SCAM_IOP )
  do m=1,pcnst
     alphanam(m) = 'AFIX'//cnst_name(m)
     alphanam(m)=to_lower(alphanam(m))
     dqfxnam(m) = 'DQFX'//cnst_name(m)
     dqfxnam(m) = to_lower(dqfxnam(m))
  end do
#endif

end subroutine initindx
