#include <misc.h>
#include <params.h>

module restart_dynamics

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plat, beglat, plnlv
   use prognostics
   use constituents, only: pnats, pcnst, ppcnst
   use ppgrid, only: pcols, pver
   use binary_io
   use scanslt,      only: scanslt_alloc, lammp, phimp, sigmp, qfcst
   use massfix,        only: alpha, hw1, hw2, hw3
#if ( defined BFB_CAM_SCAM_IOP )
   use iop,          only: alphasav, dqfx3savm1, divq3dsav, divt3dsav, t3sav, u3sav, &
                           v3sav, t2sav, q3sav, pssav, tssav, fixmassav, betasav,    &
                           init_iop_fields
#endif
   use abortutils, only: endrun
   use constituents, only: cnst_need_pdeldry

   implicit none

CONTAINS

   subroutine write_restart_dynamics (nrg)

#include <comqfl.h>

!
! Input arguments
!
      integer, intent(in) :: nrg     ! Unit number
!
! Local workspace
!
      integer :: ioerr   ! error status
!
      call wrtout_r8 (nrg,vort(1,1,beglat,n3m1), plnlv)
      call wrtout_r8 (nrg,vort(1,1,beglat,n3m2), plnlv)

      call wrtout_r8 (nrg,div(1,1,beglat,n3m1) , plnlv)
      call wrtout_r8 (nrg,div(1,1,beglat,n3m2) , plnlv)

      call wrtout_r8 (nrg,dpsl  ,plon )
      call wrtout_r8 (nrg,dpsm  ,plon )
      call wrtout_r8 (nrg,dps   ,plon )
      call wrtout_r8 (nrg,phis  ,plon )
      call wrtout_r8 (nrg,omga  ,plnlv)
!
! Write fields u3,v3,t3,q3,ps at time indices n3 and n3m1
!
      call wrtout_r8 (nrg,u3(1,1,beglat,n3m1)  ,plnlv)
      call wrtout_r8 (nrg,v3(1,1,beglat,n3m1)  ,plnlv)
      call wrtout_r8 (nrg,t3(1,1,beglat,n3m1)  ,plnlv)
      call wrtout_r8 (nrg,ps(1,beglat,n3m1)  ,plon)

      call wrtout_r8 (nrg,u3(1,1,beglat,n3m2)  ,plnlv)
      call wrtout_r8 (nrg,v3(1,1,beglat,n3m2)  ,plnlv)
      call wrtout_r8 (nrg,t3(1,1,beglat,n3m2)  ,plnlv)
      call wrtout_r8 (nrg,ps(1,beglat,n3m2)  ,plon)
      
      call wrtout_r8 (nrg,q3(1,1,1,beglat,n3m1),plnlv*(pcnst+pnats))
      call wrtout_r8 (nrg,q3(1,1,1,beglat,n3m2),plnlv*(pcnst+pnats))
      if (cnst_need_pdeldry) then
         call wrtout_r8 (nrg,pdeld(1,1,beglat,n3)  ,plnlv)
         call wrtout_r8 (nrg,pdeld(1,1,beglat,n3m1)  ,plnlv)
         call wrtout_r8 (nrg,pdeld(1,1,beglat,n3m2)  ,plnlv)
      endif!
! Write slt arrays (trajectory mid-point coordinates and 
! slt forcast of moisture and constituents
!
      call wrtout_r8 (nrg,lammp,plnlv)
      call wrtout_r8 (nrg,phimp,plnlv)
      call wrtout_r8 (nrg,sigmp,plnlv)
      call wrtout_r8 (nrg,qfcst,plnlv*pcnst)
!
! Write global integrals
!
      if (masterproc) then
         write(nrg, iostat=ioerr) tmass0, fixmas, hw1,    hw2,  &
                                  hw3, alpha
         if (ioerr /= 0 ) then
            write (6,*) 'WRITE ioerror ',ioerr,' on i/o unit = ',nrg
            call endrun ('WRITE_RESTART_DYNAMICS')
         end if
      end if

#if ( defined BFB_CAM_SCAM_IOP )
!
! Write scam values
!
     call wrtout_r8 (nrg,alphasav(1,beglat),pcnst)
     call wrtout_r8 (nrg,dqfx3savm1(1,1,1,beglat),plnlv*pcnst)       
     call wrtout_r8 (nrg,divq3dsav(1,1,1,beglat),plnlv*ppcnst)
     call wrtout_r8 (nrg,divt3dsav(1,1,beglat),plnlv)       
     call wrtout_r8 (nrg,t3sav(1,1,beglat),plnlv)       
     call wrtout_r8 (nrg,u3sav(1,1,beglat),plnlv)
     call wrtout_r8 (nrg,v3sav(1,1,beglat),plnlv)
     call wrtout_r8 (nrg,t2sav(1,1,beglat),plnlv)
     call wrtout_r8 (nrg,q3sav(1,1,1,beglat),plnlv*ppcnst)
     call wrtout_r8 (nrg,pssav(1,beglat),plon)
     call wrtout_r8 (nrg,tssav(1,beglat),plon)
     call wrtout_r8 (nrg,fixmassav(beglat),1)
     call wrtout_r8 (nrg,betasav(beglat),1)
#endif
      return
   end subroutine write_restart_dynamics

!#######################################################################

   subroutine read_restart_dynamics (nrg)

#if ( defined SPMD )
      use mpishorthand
#endif

#include <comqfl.h>
!
! Input arguments
!
      integer, intent(in) :: nrg     ! Unit number
!
! Local workspace
!
      integer :: ioerr   ! error status
!
      call initialize_prognostics
      call readin_r8 (nrg,vort(1,1,beglat,n3m1), plnlv)
      call readin_r8 (nrg,vort(1,1,beglat,n3m2), plnlv)

      call readin_r8 (nrg,div(1,1,beglat,n3m1) , plnlv)
      call readin_r8 (nrg,div(1,1,beglat,n3m2) , plnlv)

      call readin_r8 (nrg,dpsl  ,plon )
      call readin_r8 (nrg,dpsm  ,plon )
      call readin_r8 (nrg,dps   ,plon )
      call readin_r8 (nrg,phis  ,plon )
      call readin_r8 (nrg,omga  ,plnlv)
!
! Write fields u3,v3,t3,q3,ps at time indices n3 and n3m1
!
      call readin_r8 (nrg,u3(1,1,beglat,n3m1)  ,plnlv)
      call readin_r8 (nrg,v3(1,1,beglat,n3m1)  ,plnlv)
      call readin_r8 (nrg,t3(1,1,beglat,n3m1)  ,plnlv)
      call readin_r8 (nrg,ps(1,beglat,n3m1)  ,plon)

      call readin_r8 (nrg,u3(1,1,beglat,n3m2)  ,plnlv)
      call readin_r8 (nrg,v3(1,1,beglat,n3m2)  ,plnlv)
      call readin_r8 (nrg,t3(1,1,beglat,n3m2)  ,plnlv)
      call readin_r8 (nrg,ps(1,beglat,n3m2)  ,plon)
      
      call readin_r8 (nrg,q3(1,1,1,beglat,n3m1),plnlv*(pcnst+pnats))
      call readin_r8 (nrg,q3(1,1,1,beglat,n3m2),plnlv*(pcnst+pnats))
      if (cnst_need_pdeldry) then
         call readin_r8 (nrg,pdeld(1,1,beglat,n3)  ,plnlv)
         call readin_r8 (nrg,pdeld(1,1,beglat,n3m1)  ,plnlv)
         call readin_r8 (nrg,pdeld(1,1,beglat,n3m2)  ,plnlv)
      endif!
! Write slt arrays (trajectory mid-point coordinates and 
! slt forcast of moisture and constituents
!
      call scanslt_alloc()
      call readin_r8 (nrg,lammp,plnlv)
      call readin_r8 (nrg,phimp,plnlv)
      call readin_r8 (nrg,sigmp,plnlv)
      call readin_r8 (nrg,qfcst,plnlv*pcnst)
!
! Read global integrals
!
      if (masterproc) then
         read (nrg, iostat=ioerr) tmass0, fixmas, hw1,    hw2,  &
                                  hw3, alpha
         if (ioerr /= 0 ) then
            write (6,*) 'WRITE ioerror ',ioerr,' on i/o unit = ',nrg
            call endrun ('READ_RESTART_DYNAMICS')
         end if
      end if

#if ( defined SPMD )
   call mpibcast (tmass0,1         ,mpir8  ,0,mpicom)      
   call mpibcast (fixmas,1         ,mpir8  ,0,mpicom)
   call mpibcast (hw1   ,pcnst     ,mpir8  ,0,mpicom)
   call mpibcast (hw2   ,pcnst     ,mpir8  ,0,mpicom)
   call mpibcast (hw3   ,pcnst     ,mpir8  ,0,mpicom)   
   call mpibcast (alpha ,pcnst     ,mpir8  ,0,mpicom)
#endif
#if ( defined BFB_CAM_SCAM_IOP )
!
! Read scam values
!
     call init_iop_fields(ps, t3, u3, v3, q3, nocopy=.true. )

     call readin_r8 (nrg,alphasav(1,beglat),pcnst)
     call readin_r8 (nrg,dqfx3savm1(1,1,1,beglat),plnlv*pcnst)       
     call readin_r8 (nrg,divq3dsav(1,1,1,beglat),plnlv*ppcnst)
     call readin_r8 (nrg,divt3dsav(1,1,beglat),plnlv)       
     call readin_r8 (nrg,t3sav(1,1,beglat),plnlv)       
     call readin_r8 (nrg,u3sav(1,1,beglat),plnlv)
     call readin_r8 (nrg,v3sav(1,1,beglat),plnlv)
     call readin_r8 (nrg,t2sav(1,1,beglat),plnlv)
     call readin_r8 (nrg,q3sav(1,1,1,beglat),plnlv*ppcnst)
     call readin_r8 (nrg,pssav(1,beglat),plon)
     call readin_r8 (nrg,tssav(1,beglat),plon)
     call readin_r8 (nrg,fixmassav(beglat),1)
     call readin_r8 (nrg,betasav(beglat),1)
#endif
      return

   end subroutine read_restart_dynamics

end module restart_dynamics
