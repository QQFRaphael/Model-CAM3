#include <misc.h>
#include <preproc.h>

module accumulMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: accumulMod
!
! !DESCRIPTION:
! This module contains generic subroutines that can be used to
! define, accumulate and extract  user-specified fields over
! user-defined intervals. Each interval  and accumulation type is
! unique to each field processed.
! Subroutine [init_accumulator] defines the values of the accumulated
! field data structure. Subroutine [update_accum_field] does
! the actual accumulation for a given field.
! Four types of accumulations are possible:
! - Average over time interval. Time average fields are only
!   valid at the end of the averaging interval.
! - Running mean over time interval. Running means are valid once the
!   length of the simulation exceeds the
! - Running accumulation over time interval. Accumulated fields are
!   continuously accumulated. The trigger value "-99999." resets
!   the accumulation to zero.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use abortutils,   only: endrun
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: init_accum_field     ! Initialize an accumulator field
  public :: print_accum_fields   ! Print info about accumulator fields
  public :: extract_accum_field  ! Extracts the current value of an accumulator field
  interface extract_accum_field
     module procedure extract_accum_field_sl ! Extract current val of single-level accumulator field
     module procedure extract_accum_field_ml ! Extract current val of multi-level accumulator field
  end interface
  public :: update_accum_field
  interface update_accum_field               ! Updates the current value of an accumulator field
     module procedure update_accum_field_sl  ! Update single-level accumulator field
     module procedure update_accum_field_ml  ! Update multi-level accumulator field
  end interface
  public :: restart_accum                    ! Write/read restart of accumulation fields
!
! !REVISION HISTORY:
! Created by Sam Levis
! Updated to clm2.1 data structures by Mariana Vertenstein
! Updated to include all subgrid type and multilevel fields, M. Vertenstein 03/2003
!
!EOP
!
  private
!
! PRIVATE TYPES:
!
  type accum_field
     character(len=  8) :: name     !field name
     character(len=128) :: desc     !field description
     character(len=  8) :: units    !field units
     character(len=  8) :: acctype  !accumulation type: ["timeavg","runmean","runaccum"]
     integer :: period              !field accumulation period (in model time steps)
     character(len= 8) :: type1d    !subgrid type: ["gridcell","landunit","column" or "pft"]
     integer :: beg1d               !subgrid type beginning index
     integer :: end1d               !subgrid type ending index
     integer :: num1d               !total subgrid points
     integer :: numlev              !number of vertical levels in field
     real(r8), pointer :: val(:,:)  !accumulated field
     real(r8) :: initval            !initial value of accumulated field
  end type accum_field

  integer, parameter :: max_accum = 100    !maximum number of accumulated fields
  type (accum_field) :: accum(max_accum)   !array accumulated fields
  integer :: naccflds = 0                  !accumulator field counter
!------------------------------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_accum_field
!
! !INTERFACE:
  subroutine init_accum_field (name, units, desc, &
       accum_type, accum_period, numlev, subgrid_type, init_value)
!
! !DESCRIPTION:
! Initialize accumulation fields. This subroutine sets:
! o name  of accumulated field
! o units of accumulated field
! o accumulation type of accumulated field
! o description of accumulated fields: accdes
! o accumulation period for accumulated field (in iterations)
! o initial value of accumulated field
!
! !USES:
    use shr_const_mod, only: SHR_CONST_CDAY
    use time_manager, only : get_step_size
    use decompMod, only : get_proc_bounds, get_proc_global
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: name           !field name
    character(len=*), intent(in) :: units          !field units
    character(len=*), intent(in) :: desc           !field description
    character(len=*), intent(in) :: accum_type     !field type: tavg, runm, runa, ins
    integer , intent(in) :: accum_period           !field accumulation period
    character(len=*), intent(in)   :: subgrid_type !["gridcell","landunit","column" or "pft"]
    integer , intent(in) :: numlev                 !number of vertical levels
    real(r8), intent(in) :: init_value             !field initial or reset value
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 03/2003
!
!EOP
!
! LOCAL VARIABLES:
    integer :: nf           ! field index
    integer :: beg1d,end1d  ! beggining and end subgrid indices
    integer :: num1d        ! total number subgrid indices
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    integer :: numg         ! total number of gridcells across all processors
    integer :: numl         ! total number of landunits across all processors
    integer :: numc         ! total number of columns across all processors
    integer :: nump         ! total number of pfts across all processors
!------------------------------------------------------------------------

    ! Determine necessary indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)

    ! update field index
    ! Consistency check that number of accumulated does not exceed maximum.

    naccflds = naccflds + 1
    if (naccflds > max_accum) then
       write (6,*) 'INIT_ACCUM_FIELD error: user-defined accumulation fields ', &
            'equal to ',naccflds,' exceeds max_accum'
       call endrun
    end if
    nf = naccflds

    ! Note accumulation period must be converted from days
    ! to number of iterations

    accum(nf)%name = trim(name)
    accum(nf)%units = trim(units)
    accum(nf)%desc = trim(desc)
    accum(nf)%acctype = trim(accum_type)
    accum(nf)%initval = init_value
    accum(nf)%period = accum_period
    if (accum(nf)%period < 0) then
       accum(nf)%period = -accum(nf)%period * nint(SHR_CONST_CDAY) / get_step_size()
    end if

    select case (trim(subgrid_type))
    case ('gridcell')
       beg1d = begg
       end1d = endg
       num1d = numg
    case ('landunit')
       beg1d = begl
       end1d = endl
       num1d = numl
    case ('column')
       beg1d = begc
       end1d = endc
       num1d = numc
    case ('pft')
       beg1d = begp
       end1d = endp
       num1d = nump
    case default
       write(6,*)'INIT_ACCUM_FIELD: unknown subgrid type ',subgrid_type
       call endrun ()
    end select

    accum(nf)%type1d = trim(subgrid_type)
    accum(nf)%beg1d = beg1d
    accum(nf)%end1d = end1d
    accum(nf)%num1d = num1d
    accum(nf)%numlev = numlev

    ! Allocate and initialize accumulation field

    allocate(accum(nf)%val(beg1d:end1d,numlev))
    accum(nf)%val(beg1d:end1d,numlev) = init_value

  end subroutine init_accum_field

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: print_accum_fields
!
! !INTERFACE:
  subroutine print_accum_fields()
!
! !DESCRIPTION:
! Diagnostic printout of accumulated fields
!
! !USES:
    use spmdMod, only : masterproc
    use nanMod, only : bigint
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 03/2003
!
!EOP
!
    integer :: i,nf   !indices
!------------------------------------------------------------------------

    if (masterproc) then
       write(6,*)
       write(6,*) 'Initializing variables for time accumulation .....'
       write(6,'(72a1)') ("-",i=1,60)
       write(6,*) 'Accumulated fields'
       write(6,1002)
       write(6,'(72a1)') ("_",i=1,71)
       do nf = 1, naccflds
          if (accum(nf)%period /= bigint) then
             write (6,1003) nf,accum(nf)%name,accum(nf)%units,&
                  accum(nf)%acctype, accum(nf)%period, accum(nf)%initval, &
                  accum(nf)%desc
          else
             write (6,1004) nf,accum(nf)%name,accum(nf)%units,&
                  accum(nf)%acctype, accum(nf)%initval, accum(nf)%desc
          endif
       end do
       write(6,'(72a1)') ("_",i=1,71)
       write(6,*)
       write(6,'(72a1)') ("-",i=1,60)
       write(6,*) 'Successfully initialized variables for accumulation'
       write(6,*)
    endif

1002 format(' No',' Name    ',' Units   ',' Type    ','Period',' Inival',' Description')
1003 format((1x,i2),(1x,a8),(1x,a8),(1x,a8), (1x,i5),(1x,f4.0),(1x,a40))
1004 format((1x,i2),(1x,a8),(1x,a8),(1x,a8),'  N.A.',(1x,f4.0),(1x,a40))

  end subroutine print_accum_fields

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: extract_accum_field_sl
!
! !INTERFACE:
  subroutine extract_accum_field_sl (name, field, nstep)
!
! !DESCRIPTION:
! Extract single-level accumulated field.
! This routine extracts the field values from the multi-level
! accumulation field. It extracts the current value except if
! the field type is a time average. In this case, an absurd value
! is assigned to  indicate the time average is not yet valid.
!
! !USES:
    use clm_varcon, only : spval
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: name     !field name
    real(r8), pointer, dimension(:) :: field !field values for current time step
    integer , intent(in) :: nstep            !timestep index
!
! !REVISION HISTORY:
! Created by Sam Levis
! Updated to clm2.1 data structures by Mariana Vertenstein
! Updated to include all subgrid type and multilevel fields, Mariana Vertenstein 03-2003
!
!EOP
!
! LOCAL VARIABLES:
    integer :: i,k,nf        !indices
    integer :: beg,end         !subgrid beginning,ending indices
!------------------------------------------------------------------------

    ! find field index. return if "name" is not on list

    nf = 0
!dir$ concurrent
!cdir nodep
    do i = 1, naccflds
       if (name == accum(i)%name) nf = i
    end do
    if (nf == 0) then
       write(6,*) 'EXTRACT_ACCUM_FIELD_SL error: field name ',name,' not found'
       call endrun
    endif

    ! error check

    beg = accum(nf)%beg1d
    end = accum(nf)%end1d
    if (size(field,dim=1) /= end-beg+1) then
       write(6,*)'ERROR in extract_accum_field for field ',accum(nf)%name
       write(6,*)'size of first dimension of field is ',&
            size(field,dim=1),' and should be ',end-beg+1
       call endrun
    endif

    ! extract field

    if (accum(nf)%acctype == 'timeavg' .and. &
         mod(nstep,accum(nf)%period) /= 0) then
!dir$ concurrent
!cdir nodep
       do k = beg,end
          field(k) = spval  !assign absurd value when avg not ready
       end do
    else
!dir$ concurrent
!cdir nodep
       do k = beg,end
          field(k) = accum(nf)%val(k,1)
       end do
    end if

  end subroutine extract_accum_field_sl

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: extract_accum_field_ml
!
! !INTERFACE:
  subroutine extract_accum_field_ml (name, field, nstep)
!
! !DESCRIPTION:
! Extract mutli-level accumulated field.
! This routine extracts the field values from the multi-level
! accumulation field. It extracts the current value except if
! the field type is a time average. In this case, an absurd value
! is assigned to  indicate the time average is not yet valid.
!
! !USES:
    use clm_varcon, only : spval
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: name       !field name
    real(r8), pointer, dimension(:,:) :: field !field values for current time step
    integer, intent(in) :: nstep               !timestep index
!
! !REVISION HISTORY:
! Created by Sam Levis
! Updated to clm2.1 data structures by Mariana Vertenstein
! Updated to include all subgrid type and multilevel fields, M. Vertenstein 03/2003
!
!EOP
!
! LOCAL VARIABLES:
    integer :: i,j,k,nf        !indices
    integer :: beg,end         !subgrid beginning,ending indices
    integer :: numlev          !number of vertical levels
!------------------------------------------------------------------------

    ! find field index. return if "name" is not on list

    nf = 0
    do i = 1, naccflds
       if (name == accum(i)%name) nf = i
    end do
    if (nf == 0) then
       write(6,*) 'EXTRACT_ACCUM_FIELD_ML error: field name ',name,' not found'
       call endrun
    endif

    ! error check

    numlev = accum(nf)%numlev
    beg = accum(nf)%beg1d
    end = accum(nf)%end1d
    if (size(field,dim=1) /= end-beg+1) then
       write(6,*)'ERROR in extract_accum_field for field ',accum(nf)%name
       write(6,*)'size of first dimension of field is ',&
            size(field,dim=1),' and should be ',end-beg+1
       call endrun
    else if (size(field,dim=2) /= numlev) then
       write(6,*)'ERROR in extract_accum_field for field ',accum(nf)%name
       write(6,*)'size of second dimension of field iis ',&
            size(field,dim=2),' and should be ',numlev
       call endrun
    endif

    !extract field

    if (accum(nf)%acctype == 'timeavg' .and. &
         mod(nstep,accum(nf)%period) /= 0) then
       do j = 1,numlev
!dir$ concurrent
!cdir nodep
          do k = beg,end
             field(k,j) = spval  !assign absurd value when avg not ready
          end do
       end do
    else
       do j = 1,numlev
!dir$ concurrent
!cdir nodep
          do k = beg,end
             field(k,j) = accum(nf)%val(k,j)
          end do
       end do
    end if

  end subroutine extract_accum_field_ml

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: update_accum_field_sl
!
! !INTERFACE:
  subroutine update_accum_field_sl (name, field, nstep)
!
! !DESCRIPTION:
! Accumulate single level field over specified time interval.
! The appropriate field is accumulated in the array [accval].
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: name     !field name
    real(r8), pointer, dimension(:) :: field !field values for current time step
    integer , intent(in) :: nstep            !time step index
!
! !REVISION HISTORY:
! Created by Sam Levis
! Updated to clm2.1 data structures by Mariana Vertenstein
! Updated to include all subgrid type and multilevel fields by M. Vertenstein 03/2003
!
!EOP
!
! LOCAL VARIABLES:
    integer :: i,k,nf              !indices
    integer :: accper              !temporary accumulation period
    integer :: beg,end             !subgrid beginning,ending indices
!------------------------------------------------------------------------

    ! find field index. return if "name" is not on list

    nf = 0
    do i = 1, naccflds
       if (name == accum(i)%name) nf = i
    end do
    if (nf == 0) then
       write(6,*) 'UPDATE_ACCUM_FIELD_SL error: field name ',name,' not found'
       call endrun
    endif

    ! error check

    beg = accum(nf)%beg1d
    end = accum(nf)%end1d
    if (size(field,dim=1) /= end-beg+1) then
       write(6,*)'ERROR in UPDATE_ACCUM_FIELD_SL for field ',accum(nf)%name
       write(6,*)'size of first dimension of field is ',size(field,dim=1),&
            ' and should be ',end-beg+1
       call endrun
    endif

    ! accumulate field

    if (accum(nf)%acctype /= 'timeavg' .AND. &
        accum(nf)%acctype /= 'runmean' .AND. &
        accum(nf)%acctype /= 'runaccum') then
       write(6,*) 'UPDATE_ACCUM_FIELD_SL error: incorrect accumulation type'
       write(6,*) ' was specified for field ',name
       write(6,*)' accumulation type specified is ',accum(nf)%acctype
       write(6,*)' only [timeavg, runmean, runaccum] are currently acceptable'
       call endrun()
    end if


    ! reset accumulated field value if necessary and  update
    ! accumulation field
    ! running mean never reset

    if (accum(nf)%acctype == 'timeavg') then

       !time average field reset every accumulation period
       !normalize at end of accumulation period

       if ((mod(nstep,accum(nf)%period) == 1) .and. (nstep /= 0)) then
          accum(nf)%val(beg:end,1) = 0._r8
       end if
       accum(nf)%val(beg:end,1) =  accum(nf)%val(beg:end,1) + field(beg:end)
       if (mod(nstep,accum(nf)%period) == 0) then
          accum(nf)%val(beg:end,1) = accum(nf)%val(beg:end,1) / accum(nf)%period
       endif

    else if (accum(nf)%acctype == 'runmean') then

       !running mean - reset accumulation period until greater than nstep

       accper = min (nstep,accum(nf)%period)
       accum(nf)%val(beg:end,1) = ((accper-1)*accum(nf)%val(beg:end,1) + field(beg:end)) / accper

    else if (accum(nf)%acctype == 'runaccum') then

       !running accumulation field reset at trigger -99999

!dir$ concurrent
!cdir nodep
       do k = beg,end
          if (nint(field(k)) == -99999) then
             accum(nf)%val(k,1) = 0._r8
          end if
       end do
       accum(nf)%val(beg:end,1) = min(max(accum(nf)%val(beg:end,1) + field(beg:end), 0.), 99999.)

    end if

  end subroutine update_accum_field_sl

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: update_accum_field_ml
!
! !INTERFACE:
  subroutine update_accum_field_ml (name, field, nstep)
!
! !DESCRIPTION:
! Accumulate multi level field over specified time interval.
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: name       !field name
    real(r8), pointer, dimension(:,:) :: field !field values for current time step
    integer , intent(in) :: nstep              !time step index
!
! !REVISION HISTORY:
! Created by Sam Levis
! Updated to clm2.1 data structures by Mariana Vertenstein
! Updated to include all subgrid type and multilevel fields by M. Vertenstein 03/2003
!
!EOP
!
! LOCAL VARIABLES:
    integer :: i,j,k,nf            !indices
    integer :: accper              !temporary accumulation period
    integer :: beg,end             !subgrid beginning,ending indices
    integer :: numlev              !number of vertical levels
!------------------------------------------------------------------------

    ! find field index. return if "name" is not on list

    nf = 0
    do i = 1, naccflds
       if (name == accum(i)%name) nf = i
    end do
    if (nf == 0) then
       write(6,*) 'UPDATE_ACCUM_FIELD_ML error: field name ',name,' not found'
       call endrun
    endif

    ! error check

    numlev = accum(nf)%numlev
    beg = accum(nf)%beg1d
    end = accum(nf)%end1d
    if (size(field,dim=1) /= end-beg+1) then
       write(6,*)'ERROR in UPDATE_ACCUM_FIELD_ML for field ',accum(nf)%name
       write(6,*)'size of first dimension of field is ',size(field,dim=1),&
            ' and should be ',end-beg+1
       call endrun
    else if (size(field,dim=2) /= numlev) then
       write(6,*)'ERROR in UPDATE_ACCUM_FIELD_ML for field ',accum(nf)%name
       write(6,*)'size of second dimension of field is ',size(field,dim=2),&
            ' and should be ',numlev
       call endrun
    endif

    ! accumulate field

    if (accum(nf)%acctype /= 'timeavg' .AND. &
        accum(nf)%acctype /= 'runmean' .AND. &
        accum(nf)%acctype /= 'runaccum') then
       write(6,*) 'UPDATE_ACCUM_FIELD_ML error: incorrect accumulation type'
       write(6,*) ' was specified for field ',name
       write(6,*)' accumulation type specified is ',accum(nf)%acctype
       write(6,*)' only [timeavg, runmean, runaccum] are currently acceptable'
       call endrun()
    end if

    ! accumulate field

    ! reset accumulated field value if necessary and  update
    ! accumulation field
    ! running mean never reset

    if (accum(nf)%acctype == 'timeavg') then

       !time average field reset every accumulation period
       !normalize at end of accumulation period

       if ((mod(nstep,accum(nf)%period) == 1) .and. (nstep /= 0)) then
          accum(nf)%val(beg:end,1:numlev) = 0._r8
       endif
       accum(nf)%val(beg:end,1:numlev) =  accum(nf)%val(beg:end,1:numlev) + field(beg:end,1:numlev)
       if (mod(nstep,accum(nf)%period) == 0) then
          accum(nf)%val(beg:end,1:numlev) = accum(nf)%val(beg:end,1:numlev) / accum(nf)%period
       endif

    else if (accum(nf)%acctype == 'runmean') then

       !running mean - reset accumulation period until greater than nstep

       accper = min (nstep,accum(nf)%period)
       accum(nf)%val(beg:end,1:numlev) = &
            ((accper-1)*accum(nf)%val(beg:end,1:numlev) + field(beg:end,1:numlev)) / accper

    else if (accum(nf)%acctype == 'runaccum') then

       !running accumulation field reset at trigger -99999

       do j = 1,numlev
!dir$ concurrent
!cdir nodep
          do k = beg,end
             if (nint(field(k,j)) == -99999) then
                accum(nf)%val(k,j) = 0._r8
             end if
          end do
       end do
       accum(nf)%val(beg:end,1:numlev) = &
            min(max(accum(nf)%val(beg:end,1:numlev) + field(beg:end,1:numlev), 0.), 99999.)

    end if

  end subroutine update_accum_field_ml

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: restart_accum
!
! !INTERFACE:
  subroutine restart_accum (nio, flag)
!
! !DESCRIPTION:
! Read/write accumulation restart data
!
! !USES:
    use decompMod, only : get_proc_bounds, get_proc_global
    use iobinary
#if (defined SPMD)
    use spmdMod, only : mpicom, MPI_CHARACTER, MPI_INTEGER
#endif
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: nio             !restart unit
    character(len=*), intent(in) :: flag   !'read' or 'write'
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer :: nf,k,j            ! indices
    integer :: ier               ! error status
    integer :: beg1d,end1d,num1d ! starting,ending 1d indices
    integer :: numlev            ! number of vertical leveles
    integer :: begp, endp        ! per-proc beginning and ending pft indices
    integer :: begc, endc        ! per-proc beginning and ending column indices
    integer :: begl, endl        ! per-proc beginning and ending landunit indices
    integer :: begg, endg        ! per-proc gridcell ending gridcell indices
    integer :: numg              ! total number of gridcells across all processors
    integer :: numl              ! total number of landunits across all processors
    integer :: numc              ! total number of columns across all processors
    integer :: nump              ! total number of pfts across all processors
    character(len=8) :: type1d   ! subgrid type ("gridcell","landunit","column","pft")
    real(r8), pointer :: rbuf2d(:,:) ! temporary 2d buf
!------------------------------------------------------------------------

    ! Determine necessary indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)

    ! data for accumulation variables

    if (masterproc) then
       if (flag=='read' ) then
          read  (nio, iostat=ier) naccflds
       else if (flag == 'write') then
          write (nio, iostat=ier) naccflds
       endif
       if (ier /= 0) then
          write(6,*)'RESTART_ACCUM read/write error 1: ',ier
          call endrun()
       end if
       do nf = 1,naccflds
          if (flag == 'write') then
             write (nio,iostat=ier) &
                  accum(nf)%name, &
                  accum(nf)%desc,   &
                  accum(nf)%units,  &
                  accum(nf)%acctype,&
                  accum(nf)%period, &
                  accum(nf)%type1d, &
                  accum(nf)%numlev
          else if (flag == 'read') then
             read  (nio,iostat=ier) &
                  accum(nf)%name, &
                  accum(nf)%desc,   &
                  accum(nf)%units,  &
                  accum(nf)%acctype,&
                  accum(nf)%period, &
                  accum(nf)%type1d, &
                  accum(nf)%numlev
          endif
          if (ier /= 0) then
             write(6,*)'RESTART_ACCUM read/write error 2:',ier
             call endrun()
          end if
       end do
    endif  ! end of if-masterproc block

#if ( defined SPMD )
    if (flag == 'read') then
       call mpi_bcast (naccflds, 1, MPI_INTEGER, 0, mpicom, ier)
       do nf = 1,naccflds
          call mpi_bcast (accum(nf)%name   , len(accum(nf)%name)   , MPI_CHARACTER, 0, mpicom, ier)
          call mpi_bcast (accum(nf)%desc   , len(accum(nf)%desc)   , MPI_CHARACTER, 0, mpicom, ier)
          call mpi_bcast (accum(nf)%units  , len(accum(nf)%units)  , MPI_CHARACTER, 0, mpicom, ier)
          call mpi_bcast (accum(nf)%acctype, len(accum(nf)%acctype), MPI_CHARACTER, 0, mpicom, ier)
          call mpi_bcast (accum(nf)%period , 1                     , MPI_INTEGER  , 0, mpicom, ier)
          call mpi_bcast (accum(nf)%type1d , len(accum(nf)%type1d) , MPI_CHARACTER, 0, mpicom, ier)
          call mpi_bcast (accum(nf)%numlev , 1                     , MPI_INTEGER  , 0, mpicom, ier)
       end do
     endif
#endif

     ! Read/write accumulated fields

     do nf = 1,naccflds
        if (flag == 'read') then
           numlev = accum(nf)%numlev
           type1d = accum(nf)%type1d
           select case (type1d)
           case ('gridcell')
              num1d = numg
              beg1d = begg
              end1d = endg
           case ('landunit')
              num1d = numl
              beg1d = begl
              end1d = endl
           case ('column')
              num1d = numc
              beg1d = begc
              end1d = endc
           case ('pft')
              num1d = nump
              beg1d = begp
              end1d = endp
           case default
              write(6,*)'RESTART_HISTORY: read unknown 1d vertical type=',type1d
              call endrun ()
           end select
           accum(nf)%num1d = num1d
           accum(nf)%beg1d = beg1d
           accum(nf)%end1d = end1d

           allocate (accum(nf)%val(beg1d:end1d,numlev))
           if (ier /= 0) then
              write (6,*) 'restart_accum(): allocation error for field nf ',nf
              call endrun
           end if
           accum(nf)%val(beg1d:end1d,numlev) = 0.
        endif

        ! Read/write accumulated value

        type1d = accum(nf)%type1d
        numlev = accum(nf)%numlev
        num1d = accum(nf)%num1d
        beg1d = accum(nf)%beg1d
        end1d = accum(nf)%end1d
        allocate(rbuf2d(numlev,num1d))
        if (flag == 'read') then
           call readin (nio, rbuf2d, clmlevel=type1d)
           do j = 1,numlev
!dir$ concurrent
!cdir nodep
              do k = beg1d,end1d
                 accum(nf)%val(k,j) = rbuf2d(j,k)
              end do
           end do
        else if (flag == 'write') then
           do j = 1,numlev
!dir$ concurrent
!cdir nodep
              do k = beg1d,end1d
                 rbuf2d(j,k) = accum(nf)%val(k,j)
              end do
           end do
           call wrtout (nio, rbuf2d, clmlevel=type1d)
        endif
        deallocate(rbuf2d)
     end do

  end subroutine restart_accum

end module accumulMod
