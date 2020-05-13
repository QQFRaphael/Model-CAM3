#include <misc.h>
#include <params.h>
#include <max.h>

module scamMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: scamMod
! 
! !DESCRIPTION: 
! scam specific routines and data
!
! !USES:
!
  use shr_kind_mod, only: r8 => shr_kind_r8
  use constituents, only: pcnst
  use pmgrid,       only: plon, plev
!
! !PUBLIC TYPES:
  implicit none

   public   ! By default all data is public to this module

  real(r8) pressure_levels(PLEV)
  real(r8) columnLat   ! The closest lat. in the input dataset
  real(r8) columnLon   ! The closest lon. in the input dataset
  integer  initLatIdx  ! The relative position of the data in the dataset.
  integer  initLonIdx  ! The relative position of the data in the dataset.

  logical use_iop
  logical use_analysis
  logical use_saveinit
  logical use_pert_init            ! perturb initial values
  logical use_pert_frc            ! perturb forcing 
  logical fix_div3dfrc            ! fix divt3d and divq3d for cam bfb tests
  logical use_diurnal_avg
  logical use_userdata
  logical isrestart
  logical switch(NUM_SWITCHES)
  logical l_uvphys ! If true, update u/v after TPHYS
  logical l_uvadvect! If true, T, U & V will be passed to SLT
  logical l_conv    ! use flux divergence terms for T and q?     
  logical l_divtr   ! use flux divergence terms for constituents?
  logical l_diag    ! do we want available diagnostics?

  integer initTimeIdx
  integer seedval

  character*(MAX_PATH_LEN) modelfile
  character*(MAX_PATH_LEN) analysisfile
  character*(MAX_PATH_LEN) sicfile
  character*(MAX_PATH_LEN) userfile
  character*(MAX_PATH_LEN) sstfile
  character*(MAX_PATH_LEN) lsmpftfile
  character*(MAX_PATH_LEN) pressfile
  character*(MAX_PATH_LEN) ozonefile
  character*(MAX_PATH_LEN) iopfile
  character*(MAX_PATH_LEN) absemsfile
  character*(MAX_PATH_LEN) aermassfile
  character*(MAX_PATH_LEN) aeropticsfile
  character*(MAX_PATH_LEN) timeinvfile
  character*(MAX_PATH_LEN) lsmsurffile
  character*(MAX_PATH_LEN) lsminifile

  real(r8) fixmascam(plon)
  real(r8) betacam(plon)
  real(r8) alphacam(pcnst)
  real(r8) dqfxcam(plon,plev,pcnst)

  real(r8)     divq3d(plev,pcnst)  ! 3D q advection
  real(r8)     divt3d(plev)        ! 3D T advection
  real(r8)     vertdivq(plev,pcnst)! vertical q advection
  real(r8)     vertdivt(plev)      ! vertical T advection
  real(r8)     ptend               ! surface pressure tendency
  real(r8)     qdiff(plev)         ! model minus observed humidity
  real(r8)     qobs(plev)          ! actual W.V. Mixing ratio
  real(r8)     precobs(1)             ! observed precipitation 
  real(r8)     lhflxobs(1)            ! observed surface latent heat flux 
  real(r8)     shflxobs(1)            ! observed surface sensible heat flux
  real(r8)     q1obs(plev)         ! observed apparent heat source
  real(r8)     q2obs(plev)         ! observed apparent heat sink
  real(r8)     tdiff(plev)         ! model minus observed temp 
  real(r8)     tground(1)             ! ground temperature
  real(r8)     tobs(plev)          ! actual temperature
  real(r8)     tsair(1)               ! air temperature at the surface
  real(r8)     udiff(plev)         ! model minus observed uwind
  real(r8)     uobs(plev)          ! actual u wind
  real(r8)     vdiff(plev)         ! model minus observed vwind
  real(r8)     vobs(plev)          ! actual v wind
  real(r8)     cldobs(plev)        ! observed cld
  real(r8)     clwpobs(plev)        ! observed clwp
  real(r8)     aldirobs(1)        ! observed aldir
  real(r8)     aldifobs(1)        ! observed aldif
  real(r8)     asdirobs(1)        ! observed asdir
  real(r8)     asdifobs(1)        ! observed asdif

  integer iopTimeIdx           ! index into iop dataset

  logical*4 doiopupdate   ! do we need to read next iop timepoint
  logical*4 have_divq     ! dataset contains divq 
  logical*4 have_divt     ! dataset contains divt
  logical*4 have_divq3d   ! dataset contains divq3d 
  logical*4 have_vertdivt ! dataset contains vertdivt
  logical*4 have_vertdivq ! dataset contains vertdivq 
  logical*4 have_divt3d   ! dataset contains divt3d
  logical*4 have_divu     ! dataset contains divu
  logical*4 have_divv     ! dataset contains divv 
  logical*4 have_omega    ! dataset contains omega
  logical*4 have_phis     ! dataset contains phis
  logical*4 have_ptend    ! dataset contains ptend
  logical*4 have_ps       ! dataset contains ps
  logical*4 have_q        ! dataset contains q
  logical*4 have_q1       ! dataset contains Q1
  logical*4 have_q2       ! dataset contains Q2
  logical*4 have_prec     ! dataset contains prec 
  logical*4 have_lhflx    ! dataset contains lhflx 
  logical*4 have_shflx    ! dataset contains shflx
  logical*4 have_t        ! dataset contains t
  logical*4 have_tg       ! dataset contains tg
  logical*4 have_tsair    ! dataset contains tsair
  logical*4 have_u        ! dataset contains u 
  logical*4 have_v        ! dataset contains v 
  logical*4 have_cld      ! dataset contains cld
  logical*4 have_clwp     ! dataset contains clwp
  logical*4 have_aldir     ! dataset contains aldir
  logical*4 have_aldif     ! dataset contains aldif
  logical*4 have_asdir     ! dataset contains asdir
  logical*4 have_asdif     ! dataset contains asdif
  logical*4 use_srfprop   ! use the specified surface properties
  logical*4 use_relax     ! use relaxation
  logical*4 use_ccmiop    ! use ccm generated forcing 
  logical*4 use_3dfrc     ! use 3d forcing

end module scamMod
