#include <misc.h>
#include <params.h>

module constituents

!----------------------------------------------------------------------------------------------
! 
! Purpose: Contains data and functions for manipulating advected and non-advected constituents.
!
! Revision history:
!             B.A. Boville    Original version
! June 2003   P. Rasch        Add wet/dry m.r. specifier
! 2004-08-28  B. Eaton        Add query function to allow turning off the default CAM output of
!                             constituents so that chemistry module can make the outfld calls.
!                             Allow cnst_get_ind to return without aborting when constituent not
!                             found.
!----------------------------------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst,    only: r_universal
  use pmgrid,       only: masterproc
  use abortutils,   only: endrun

  implicit none
  private
  save
!
! Public interfaces
!
  public cnst_add          ! add a constituent to the list of advected (or nonadvected) constituents
  public cnst_get_ind      ! get the index of a constituent
  public cnst_get_type_byind ! get the type of a constituent
  public cnst_get_type_byname ! get the type of a constituent
  public cnst_read_iv         ! query whether constituent initial values are read from initial file
  public cnst_chk_dim         ! check that number of constituents added equals dimensions (pcnst, pnats)
  public cnst_cam_outfld      ! Returns true if default CAM output was specified in the cnst_add calls.

! Public data

  integer, parameter, public :: pcnst  = PCNST      ! number of advected constituents (including water vapor)
  integer, parameter, public :: pnats  = PNATS      ! number of non-advected constituents
  integer, parameter, public :: ppcnst = pcnst+pnats! total number of constituents

  integer, parameter, public :: advected = 0        ! type value for constituents which are advected
  integer, parameter, public :: nonadvec = 1        ! type value for constituents which are not advected

  character(len=8),public :: cnst_name(ppcnst)      ! constituent names (including any non-advected)
  character(len=128),public :: cnst_longname(ppcnst)! long name of constituents
  logical, public ::  cnst_need_pdeldry = .false.   ! true if any dry constituents

! Namelist variables
  logical, public :: readtrace = .true.             ! true => obtain initial tracer data from IC file

!
! Constants for each tracer
  real(r8), public :: cnst_cp  (ppcnst)             ! specific heat at constant pressure (J/kg/K)
  real(r8), public :: cnst_cv  (ppcnst)             ! specific heat at constant volume (J/kg/K)
  real(r8), public :: cnst_mw  (ppcnst)             ! molecular weight (kg/kmole)
  character*3, public :: cnst_type(ppcnst)          ! wet or dry mixing ratio
  real(r8), public :: cnst_rgas(ppcnst)             ! gas constant ()
  real(r8), public :: qmin     (ppcnst)             ! minimum permitted constituent concentration (kg/kg)
  real(r8), public :: qmincg   (ppcnst)             ! for backward compatibility only
  logical,  public :: cnst_fixed_ubc(ppcnst) = .false.  ! upper bndy condition = fixed ?

!++bee - temporary...
! Lists of tracer names and diagnostics
   character(len=8), public :: hadvnam(pcnst)        ! names of horizontal advection tendencies
   character(len=8), public :: vadvnam(pcnst)        ! names of vertical advection tendencies
   character(len=8), public :: dcconnam(ppcnst)       ! names of convection tendencies
   character(len=8), public :: fixcnam(pcnst)         ! names of species slt fixer tendencies
   character(len=8), public :: tendnam(pcnst)         ! names of total tendencies of species
   character(len=8), public :: sflxnam(ppcnst)        ! names of surface fluxes of species
   character(len=8), public :: tottnam(pcnst)         ! names for horz + vert + fixer tendencies

!--bee

! Private data

  integer :: padv = 0                      ! index pointer to last advected tracer
  integer :: pnad = pcnst                  ! index pointer to last non-advected tracer
  logical :: read_init_vals(ppcnst)        ! true => read initial values from initial file
  logical :: cam_outfld_(ppcnst)           ! true  => default CAM output of constituents in kg/kg
                                           ! false => chemistry is responsible for making outfld
                                           !          calls for constituents

!==============================================================================================
CONTAINS
!==============================================================================================

  subroutine cnst_add (name, type, mwc, cpc, qminc, &
                       ind, longname, readiv, mixtype, cam_outfld, fixed_ubc)
!----------------------------------------------------------------------- 
! 
! Purpose: Register a constituent for treament by physical parameterizations and 
!          transport (if type=advected)
!
!---------------------------------------------------------------------------------
!
    character(len=*), intent(in) :: &
       name      ! constituent name used as variable name in history file output (8 char max)
    integer, intent(in)    :: type   ! flag indicating advected or nonadvected
    real(r8),intent(in)    :: mwc    ! constituent molecular weight (kg/kmol)
    real(r8),intent(in)    :: cpc    ! constituent specific heat at constant pressure (J/kg/K)
    real(r8),intent(in)    :: qminc  ! minimum value of mass mixing ratio (kg/kg)
                                     ! normally 0., except water 1.E-12, for radiation.
    integer, intent(out)   :: ind    ! global constituent index (in q array)

    character(len=*), intent(in), optional :: &
       longname    ! value for long_name attribute in netcdf output (128 char max, defaults to name)
    logical,          intent(in), optional :: &
       readiv      ! true => read initial values from initial file (default: true)
    character(len=*), intent(in), optional :: &
       mixtype     ! mixing ratio type (dry, wet)
    logical,          intent(in), optional :: &
       cam_outfld  ! true => default CAM output of constituent in kg/kg
    logical,          intent(in), optional :: &
       fixed_ubc ! true => const has a fixed upper bndy condition

!-----------------------------------------------------------------------

! set tracer index and check validity, advected tracer
    if (type == advected) then
       padv = padv+1
       ind  = padv
       if (padv > pcnst) then
          write(6,*) 'CNST_ADD: advected tracer index greater than pcnst = ', pcnst
          call endrun
       end if

! set tracer index and check validity, non-advected tracer
    else if (type == nonadvec) then
       pnad = pnad+1
       ind  = pnad
       if (pnad > ppcnst) then
          write(6,*) 'CNST_ADD: non-advected tracer index greater than pcnst+pnats = ', ppcnst
          call endrun
       end if

! unrecognized type value
    else
       write(6,*) 'CNST_ADD, input type flag not valid, type=', type
       call endrun
    end if

! set tracer name and constants
    cnst_name(ind) = name
    if ( present(longname) )then
       cnst_longname(ind) = longname
    else
       cnst_longname(ind) = name
    end if

! set whether to read initial values from initial file
    if ( present(readiv) ) then
       read_init_vals(ind) = readiv
    else
       read_init_vals(ind) = readtrace
    end if

! set constituent mixing ratio type
    if ( present(mixtype) )then
       cnst_type(ind) = mixtype
    else
       cnst_type(ind) = 'wet'
    end if

! set outfld type 
! (false: the module declaring the constituent is responsible for outfld calls)
    if ( present(cam_outfld) ) then
       cam_outfld_(ind) = cam_outfld
    else
       cam_outfld_(ind) = .true.
    end if

! set upper boundary condition type
    if ( present(fixed_ubc) ) then
       cnst_fixed_ubc(ind) = fixed_ubc
    else
       cnst_fixed_ubc(ind) = .false.
    end if

    cnst_cp  (ind) = cpc
    cnst_mw  (ind) = mwc
    qmin     (ind) = qminc
    qmincg   (ind) = qminc
    if (ind == 1) qmincg = 0.  ! This crap is replicate what was there before ****

    cnst_rgas(ind) = r_universal * mwc
    cnst_cv  (ind) = cpc - cnst_rgas(ind)

    return
  end subroutine cnst_add

!==============================================================================

  subroutine cnst_get_ind (name, ind, abort)
!----------------------------------------------------------------------- 
! 
! Purpose: Get the index of a constituent 
! 
! Author:  B.A. Boville
! 
!-----------------------------Arguments---------------------------------
!
    character(len=*),  intent(in)  :: name  ! constituent name
    integer,           intent(out) :: ind   ! global constituent index (in q array)
    logical, optional, intent(in)  :: abort ! optional flag controlling abort

!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index
    logical :: abort_on_error
!-----------------------------------------------------------------------

! Find tracer name in list
    do m = 1, ppcnst
       if (name == cnst_name(m)) then
          ind  = m
          return
       end if
    end do

! Unrecognized name
    abort_on_error = .true.
    if ( present(abort) ) abort_on_error = abort

    if ( abort_on_error ) then
       write(6,*) 'CNST_GET_IND, name:', name,  ' not found in list:', cnst_name(:)
       call endrun('CNST_GET_IND: name not found')
    end if

! error return
    ind = -1

  end subroutine cnst_get_ind

!==============================================================================================

  character*3 function cnst_get_type_byind (ind)
!----------------------------------------------------------------------- 
! 
! Purpose: Get the type of a constituent 
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author:  P. J. Rasch
! 
!-----------------------------Arguments---------------------------------
!
    integer, intent(in)   :: ind    ! global constituent index (in q array)

!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index

!-----------------------------------------------------------------------

    if (ind.le.ppcnst) then
       cnst_get_type_byind = cnst_type(ind)
    else
       ! Unrecognized name
       write(6,*) 'CNST_GET_TYPE_BYIND, ind:', ind
       call endrun
    endif


  end function cnst_get_type_byind

!==============================================================================================

  character*3 function cnst_get_type_byname (name)
!----------------------------------------------------------------------- 
! 
! Purpose: Get the type of a constituent 
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author:  P. J. Rasch
! 
!-----------------------------Arguments---------------------------------
!
    character(len=*), intent(in) :: name ! constituent name

!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index

!-----------------------------------------------------------------------

    do m = 1, ppcnst
       if (name == cnst_name(m)) then
          cnst_get_type_byname = cnst_type(m)
          return
       end if
    end do

! Unrecognized name
    write(6,*) 'CNST_GET_TYPE_BYNAME, name:', name,  ' not found in list:', cnst_name(:)
    call endrun

  end function cnst_get_type_byname

!==============================================================================
  function cnst_read_iv(m)
!----------------------------------------------------------------------- 
! 
! Purpose: Query whether constituent initial values are read from initial file.
! 
! Author:  B. Eaton
! 
!-----------------------------Arguments---------------------------------
!
    integer, intent(in) :: m    ! constituent index

    logical :: cnst_read_iv     ! true => read initial values from inital file
!-----------------------------------------------------------------------

    cnst_read_iv = read_init_vals(m)
 end function cnst_read_iv

!==============================================================================
  subroutine cnst_chk_dim
!----------------------------------------------------------------------- 
! 
! Purpose: Check that the number of registered constituents of each type is the
!          same as the dimension
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author:  B.A. Boville
! 
    integer i,m
!-----------------------------------------------------------------------
!
    if (padv /= pcnst) then
       write(6,*)'CNST_CHK_DIM: number of advected tracer ',padv, ' not equal to pcnst = ',pcnst
       call endrun ()
    endif
    if (pnad /= ppcnst) then
       write(6,*)'CNST_CHK_DIM: number of non-advected tracers ',pnad, ' not equal to pcnst+pnats = ', &
            ppcnst
       call endrun ()
    endif

    if (masterproc) then
       write (6,*) 'Advected constituent list:'
       do i = 1, pcnst
          write(6,'(i4,2x,a8,2x,a128,2x,a3)') i, cnst_name(i), cnst_longname(i), cnst_type(i)
       end do
       if (ppcnst > pcnst) then
          write (6,*) 'Nonadvected constituent list:'
          do i = pcnst+1, ppcnst
             write(6,'(i4,2x,a8,2x,a128)') i, cnst_name(i), cnst_longname(i)
          end do
       else
          write (6,*) 'No nonadvected constituents'
       end if
    end if

    !set cnst_need_pdeldry
    do m = 1,pcnst
       if (cnst_type(m).eq.'dry') then
          cnst_need_pdeldry = .true.
       endif
    end do

  end subroutine cnst_chk_dim

!==============================================================================

function cnst_cam_outfld(m)
!----------------------------------------------------------------------- 
! 
! Purpose:
! Query whether default CAM outfld calls should be made.
! 
!----------------------------------------------------------------------- 
   integer, intent(in) :: m                ! constituent index
   logical             :: cnst_cam_outfld  ! true => use default CAM outfld calls
!-----------------------------------------------------------------------

   cnst_cam_outfld = cam_outfld_(m)

end function cnst_cam_outfld

!==============================================================================

end module constituents
