#include <misc.h>
#include <params.h>

module restart_dynamics

   use shr_kind_mod, only: r8 => shr_kind_r8
   use constituents, only: pcnst, pnats
   use prognostics, only: u3, v3, q3, qm1, div, dps, dpsl, dpsm, phis, phisl, &
                          phism, ps, urhs, vrhs, trhs, prhs, ql, qm, &
                          tl, tm, omga, n3, n3m1, parrsld, t3, etadot, &
                          initialize_prognostics
   use pmgrid, only: plon, plat, plnlv, plevp, masterproc, plev, beglat, endlat
   use ppgrid, only: pcols, pver
   use scanslt,only: slt_alloc, qfcst, hw1, hw2, hw3, alpha
   use binary_io, only: wrtout_r8, readin_r8
#ifdef SPMD
   use pspect, only: psp
#endif
   use comspe, only: lnpstar
   use abortutils, only: endrun

   implicit none

   private

!
! Public interfaces
!
   public write_restart_dynamics   ! Write out restart information
   public read_restart_dynamics    ! Read in restart information

CONTAINS

   subroutine write_restart_dynamics (nrg)

#include <comqfl.h>

!
! Input arguments
!
      integer :: nrg     ! Unit number
!
! Local workspace
!
      integer :: ioerr   ! error status
      integer :: num     ! number of values
!
      call wrtout_r8 (nrg,div(1,1,beglat,n3  ) , plnlv)
      call wrtout_r8 (nrg,div(1,1,beglat,n3m1) , plnlv)

      call wrtout_r8 (nrg,urhs,plnlv)
      call wrtout_r8 (nrg,vrhs,plnlv)
      call wrtout_r8 (nrg,trhs,plnlv)
      call wrtout_r8 (nrg,prhs,plnlv)

      call wrtout_r8 (nrg,dpsl  ,plon )
      call wrtout_r8 (nrg,dpsm  ,plon )
      call wrtout_r8 (nrg,ql    ,plnlv)
      call wrtout_r8 (nrg,qm    ,plnlv)
      call wrtout_r8 (nrg,dps   ,plon )
      call wrtout_r8 (nrg,phis  ,plon )
      call wrtout_r8 (nrg,phisl ,plon )
      call wrtout_r8 (nrg,phism ,plon )
      call wrtout_r8 (nrg,tl    ,plnlv)
      call wrtout_r8 (nrg,tm    ,plnlv)
      call wrtout_r8 (nrg,omga  ,plnlv)
!
! Write fields u3,v3,t3,q3,qm1,ps at time indices n3 and n3m1
!
      call wrtout_r8 (nrg,u3(1,1,beglat,n3  )  ,plnlv)
      call wrtout_r8 (nrg,v3(1,1,beglat,n3  )  ,plnlv)
      call wrtout_r8 (nrg,t3(1,1,beglat,n3  )  ,plnlv)
      call wrtout_r8 (nrg,ps(1,beglat,n3  )  ,plon)

      call wrtout_r8 (nrg,u3(1,1,beglat,n3m1)  ,plnlv)
      call wrtout_r8 (nrg,v3(1,1,beglat,n3m1)  ,plnlv)
      call wrtout_r8 (nrg,t3(1,1,beglat,n3m1)  ,plnlv)
      call wrtout_r8 (nrg,ps(1,beglat,n3m1)  ,plon)
      
      call wrtout_r8 (nrg,q3 (1,1,1,beglat,n3  ),plnlv*(pcnst+pnats))
      call wrtout_r8 (nrg,q3 (1,1,1,beglat,n3m1),plnlv*(pcnst+pnats))
      call wrtout_r8 (nrg,qm1(1,1,1,beglat     ),plnlv*(pcnst+pnats))

      num = plnlv
      call wrtout_r8 (nrg,parrsld(1,1,beglat),num)

      num = plon*plevp
      call wrtout_r8 (nrg,etadot(1,1,beglat),num)

      num = plnlv*pcnst
      call wrtout_r8  (nrg,qfcst(1,1,1,beglat),num)
!
! Write global integrals
!
      if (masterproc) then
         write(nrg, iostat=ioerr) tmass0, hw1, hw2, hw3, alpha, lnpstar
         if (ioerr /= 0 ) then
            write (6,*) 'WRITE ioerror ',ioerr,' on i/o unit = ',nrg
            call endrun ('WRITE_RESTART_DYNAMICS')
         end if
      end if

      return
   end subroutine write_restart_dynamics

   subroutine read_restart_dynamics (nrg)

#if ( defined SPMD )
      use mpishorthand
#endif

#include <comqfl.h>
!
! Input arguments
!
      integer :: nrg     ! Unit number
!
! Local workspace
!
      integer :: ioerr   ! error status
      integer :: num     ! number of values
!
      call initialize_prognostics
      call readin_r8 (nrg,div(1,1,beglat,n3  ) , plnlv)
      call readin_r8 (nrg,div(1,1,beglat,n3m1) , plnlv)

      call readin_r8 (nrg,urhs,plnlv)
      call readin_r8 (nrg,vrhs,plnlv)
      call readin_r8 (nrg,trhs,plnlv)
      call readin_r8 (nrg,prhs,plnlv)

      call readin_r8 (nrg,dpsl  ,plon )
      call readin_r8 (nrg,dpsm  ,plon )
      call readin_r8 (nrg,ql    ,plnlv)
      call readin_r8 (nrg,qm    ,plnlv)
      call readin_r8 (nrg,dps   ,plon )
      call readin_r8 (nrg,phis  ,plon )
      call readin_r8 (nrg,phisl ,plon )
      call readin_r8 (nrg,phism ,plon )
      call readin_r8 (nrg,tl    ,plnlv)
      call readin_r8 (nrg,tm    ,plnlv)
      call readin_r8 (nrg,omga  ,plnlv)
!
! Read fields u3,v3,t3,q3,qm1,ps at time indices n3 and n3m1
!
      call readin_r8 (nrg,u3(1,1,beglat,n3  )  ,plnlv)
      call readin_r8 (nrg,v3(1,1,beglat,n3  )  ,plnlv)
      call readin_r8 (nrg,t3(1,1,beglat,n3  )  ,plnlv)
      call readin_r8 (nrg,ps(1,beglat,n3  )  ,plon)

      call readin_r8 (nrg,u3(1,1,beglat,n3m1)  ,plnlv)
      call readin_r8 (nrg,v3(1,1,beglat,n3m1)  ,plnlv)
      call readin_r8 (nrg,t3(1,1,beglat,n3m1)  ,plnlv)
      call readin_r8 (nrg,ps(1,beglat,n3m1)  ,plon)
      
      call readin_r8 (nrg,q3 (1,1,1,beglat,n3  ),plnlv*(pcnst+pnats))
      call readin_r8 (nrg,q3 (1,1,1,beglat,n3m1),plnlv*(pcnst+pnats))
      call readin_r8 (nrg,qm1(1,1,1,beglat     ),plnlv*(pcnst+pnats))

      num = plnlv
      call readin_r8 (nrg,parrsld(1,1,beglat),num)

      num = plon*plevp
      call readin_r8 (nrg,etadot(1,1,beglat),num)

!
! Read slt forcast of moisture and constituents
!
      call slt_alloc()
      num = plnlv*pcnst
      call readin_r8 (nrg,qfcst(1,1,1,beglat),num)
!
! Read global integrals
!
      if (masterproc) then
         read(nrg, iostat=ioerr) tmass0, hw1, hw2, hw3, alpha, lnpstar
         if (ioerr /= 0 ) then
            write (6,*) 'READ ioerror ',ioerr,' on i/o unit = ',nrg
            call endrun ('READ_RESTART_DYNAMICS')
         end if
      end if

#if ( defined SPMD )
      call mpibcast (tmass0,1         ,mpir8  ,0,mpicom)      
      call mpibcast (hw1   ,pcnst     ,mpir8  ,0,mpicom)
      call mpibcast (hw2   ,pcnst     ,mpir8  ,0,mpicom)
      call mpibcast (hw3   ,pcnst     ,mpir8  ,0,mpicom)   
      call mpibcast (alpha ,pcnst     ,mpir8  ,0,mpicom)
      call mpibcast (lnpstar, psp     ,mpir8  ,0,mpicom)
#endif

      return
   end subroutine read_restart_dynamics

end module restart_dynamics
