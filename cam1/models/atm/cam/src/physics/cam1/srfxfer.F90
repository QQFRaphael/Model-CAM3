#include <misc.h>
#include <params.h>
subroutine srfxfer(state,surface_state2d,prec_zmc,snow_zmc, &
	   prec_cmf,snow_cmf,prec_sed,snow_sed, &
	   prec_pcw,snow_pcw,fsns)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Transfer atmospheric fields into common block /comsrf/
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: L. Bath  CMS Contact: M. Vertenstein
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use physics_types,   only: physics_state
   use ppgrid, only : pver, pcols
   use history, only: outfld
   use comsrf, only: psm1, srfrpdel, prcsnw
   use comsrf, only: surface_state   ! user defined type. The field
                                     ! surface_state2d is now in the argument list

#if ( defined COUP_CSM )
   use ccsm_msg, only: rho, netsw
#endif
   implicit none

#include <comtsc.h>
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   type(physics_state), intent(in) :: state
   type (surface_state), intent(inout) :: surface_state2d
! convective precipitation variables
   real(r8), intent(in) :: prec_zmc(pcols)                ! total precipitation   from ZM convection
   real(r8), intent(in) :: snow_zmc(pcols)                ! snow from ZM   convection
   real(r8), intent(in) :: prec_cmf(pcols)                ! total precipitation   from Hack convection
   real(r8), intent(in) :: snow_cmf(pcols)                ! snow from   Hack   convection
   real(r8), intent(in) :: prec_sed(pcols)                ! total precipitation   from ZM convection
   real(r8), intent(in) :: snow_sed(pcols)                ! snow from ZM   convection
   real(r8), intent(in) :: prec_pcw(pcols)                ! total precipitation   from Hack convection
   real(r8), intent(in) :: snow_pcw(pcols)                ! snow from Hack   convection
   real(r8), intent(in) :: fsns(pcols)                    ! Surface solar absorbed flux

!
!---------------------------Local variables-----------------------------
!
   integer i                 ! Longitude index
   integer :: lchnk                 ! Chunk index
   integer :: ncol
   real(r8):: prect(pcols)                ! total (conv+large scale) precip rate



!
!-----------------------------------------------------------------------

   lchnk = state%lchnk
   ncol = state%ncol

!
! Stuff global fluxes and state variables into common
!
   do i=1,ncol
      surface_state2d%tbot(i) = state%t(i,pver)
      surface_state2d%thbot(i)= state%t(i,pver) * state%exner(i,pver)
      surface_state2d%zbot(i) = state%zm(i,pver)
      surface_state2d%ubot(i) = state%u(i,pver)
      surface_state2d%vbot(i) = state%v(i,pver)
      surface_state2d%qbot(i) = state%q(i,pver,1) 
      surface_state2d%pbot(i) = state%pmid(i,pver)
      psm1(i,lchnk) = state%ps(i)
      srfrpdel(i,lchnk) = state%rpdel(i,pver)
#if ( defined COUP_CSM )
      rho(i,lchnk)   = state%pmid(i,pver)/(rair*surface_state2d%tbot(i))
      netsw(i,lchnk) = surface_state2d%srfrad(i) - &
	surface_state2d%flwds(i)
#endif
   end do

! Compute net surface radiative flux for use by surface temperature code.
! Note that units have already been converted to mks in RADCTL.  Since
! fsns and flwds are in the buffer, array values will be carried across
! timesteps when the radiation code is not invoked.
!
   surface_state2d%srfrad(:ncol) = fsns(:ncol) + surface_state2d%flwds(:ncol)
   call outfld('SRFRAD  ',surface_state2d%srfrad,pcols,lchnk)

!
! Compile and write out precipation and snow rates from shallow convection, 
! deep convection and stratiform processes.
!

! Compute total convective and stratiform precipitation and snow rates
   do i=1,ncol
      surface_state2d%precc (i) = prec_zmc(i) + prec_cmf(i)
      surface_state2d%precl (i) = prec_sed(i) + prec_pcw(i)
      surface_state2d%precsc(i) = snow_zmc(i) + snow_cmf(i)
      surface_state2d%precsl(i) = snow_sed(i) + snow_pcw(i)
! jrm These checks should not be necessary if they exist in the parameterizations
      if(surface_state2d%precc(i).lt.0.) surface_state2d%precc(i)=0.
      if(surface_state2d%precl(i).lt.0.) surface_state2d%precl(i)=0.
      if(surface_state2d%precsc(i).lt.0.) surface_state2d%precsc(i)=0.
      if(surface_state2d%precsl(i).lt.0.) surface_state2d%precsl(i)=0.
      if(surface_state2d%precsc(i).gt.surface_state2d%precc(i)) surface_state2d%precsc(i)=surface_state2d%precc(i)
      if(surface_state2d%precsl(i).gt.surface_state2d%precl(i)) surface_state2d%precsl(i)=surface_state2d%precl(i)
! end jrm
   end do
   prcsnw(:ncol,lchnk) = surface_state2d%precsc(:ncol) + surface_state2d%precsl(:ncol)   ! total snowfall rate: needed by slab ocean model

   call outfld('PRECZ   ',                prec_zmc,pcols   ,lchnk       )
   call outfld('PRECL   ',surface_state2d%precl   ,pcols   ,lchnk       )
   call outfld('PRECC   ',surface_state2d%precc   ,pcols   ,lchnk       )
   call outfld('PRECSL  ',surface_state2d%precsl  ,pcols   ,lchnk       )
   call outfld('PRECSC  ',surface_state2d%precsc  ,pcols   ,lchnk       )
   
   prect(:ncol) = surface_state2d%precc(:ncol) + surface_state2d%precl(:ncol)
   call outfld('PRECT   ',prect   ,pcols   ,lchnk       )
   call outfld('PRECTMX ',prect   ,pcols   ,lchnk       )

#if ( defined COUP_CSM )
   call outfld('PRECLav ',surface_state2d%precl   ,pcols   ,lchnk   )
   call outfld('PRECCav ',surface_state2d%precc   ,pcols   ,lchnk   )
#endif

#if ( defined BFB_CAM_SCAM_IOP )
   call outfld('Prec   ',prect   ,pcols   ,lchnk       )
#endif

   return
end subroutine srfxfer







