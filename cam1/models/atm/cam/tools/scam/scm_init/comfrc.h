!   
!    File:   /u1/hoswell/home/sccm/sccm/com/comfrc.com  (376 bytes) 
!    Date:   Tue Apr 18 15:51:55 1995
!    Author: Mike Hoswell <hoswell@francke.cgd.ucar.edu>
!   
#ifndef _COMFRC_H
#define _COMFRC_H
!
! Forcing variables
!
      common/comfrc/ wfld    ,wfldh   ,divq    ,divt    ,divu    ,divv
!
      real(r8) wfld(plev)         ! Vertical motion (slt)
      real(r8) wfldh(plevp)         ! Vertical motion (slt)
      real(r8) divq(plev,pcnst)   ! Divergence of moisture
      real(r8) divt(plev)         ! Divergence of temperature
      real(r8) divu(plev)         ! Horiz Divergence of E/W
      real(r8) divv(plev)          ! Horiz Divergence of N/S
!
#endif

