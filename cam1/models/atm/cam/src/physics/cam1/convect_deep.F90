
#include <misc.h>
#include <params.h>

module convect_deep
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to the Zhang-McFarlane deep convection scheme
!
! Author: D.B. Coleman
!
!---------------------------------------------------------------------------------
   use shr_kind_mod, only: r8=>shr_kind_r8
   use physconst,    only: cpair                              
   use ppgrid,       only: pver, pcols, pverp, begchunk, endchunk
   use zm_conv,         only: zm_conv_evap, zm_convr, convtran
   use history,         only: outfld, addfld, add_default, phys_decomp

   implicit none

   save
   private                         ! Make default type private to the module

! Public methods

   public ::&
      convect_deep_register,           &! register fields in physics buffer
      convect_deep_init,               &! initialize donner_deep module
      convect_deep_tend,               &! return tendencies
      convect_deep_tend_2               ! return tendencies


!
! Private module data
!

   real(r8), allocatable, dimension(:,:,:) :: mu  !(pcols,pver,begchunk:endchunk)
   real(r8), allocatable, dimension(:,:,:) :: eu  !(pcols,pver,begchunk:endchunk)
   real(r8), allocatable, dimension(:,:,:) :: du  !(pcols,pver,begchunk:endchunk)
   real(r8), allocatable, dimension(:,:,:) :: md  !(pcols,pver,begchunk:endchunk)
   real(r8), allocatable, dimension(:,:,:) :: ed  !(pcols,pver,begchunk:endchunk)
   real(r8), allocatable, dimension(:,:,:) :: dp  !(pcols,pver,begchunk:endchunk) 
	! wg layer thickness in mbs (between upper/lower interface).
   real(r8), allocatable, dimension(:,:)   :: dsubcld  !(pcols,begchunk:endchunk)
	! wg layer thickness in mbs between lcl and maxi.

   integer, allocatable, dimension(:,:) :: jt   !(pcols,begchunk:endchunk)
        ! wg top  level index of deep cumulus convection.
   integer, allocatable, dimension(:,:) :: maxg !(pcols,begchunk:endchunk)
        ! wg gathered values of maxi.
   integer, allocatable, dimension(:,:) :: ideep !(pcols,begchunk:endchunk)               
	! w holds position of gathered points vs longitude index

   integer, allocatable, dimension(:) :: lengath !(begchunk:endchunk)


! integer ::& ! indices for fields in the physics buffer
!      omint_idx,    &

!=========================================================================================
  contains

!=========================================================================================
subroutine convect_deep_register

!----------------------------------------
! Purpose: register fields with the physics buffer
!----------------------------------------

  use phys_buffer, only: pbuf_times, pbuf_add

  implicit none

  integer idx

   
  call pbuf_add('ICWMRDP' , 'physpkg', 1,pver,      1, idx)
  call pbuf_add('RPRDDP' , 'physpkg', 1,pver,      1, idx)


end subroutine convect_deep_register

!=========================================================================================

subroutine convect_deep_init(hypi)

!----------------------------------------
! Purpose:  declare output fields, initialize variables needed by convection
!----------------------------------------

  use history,         only: outfld, addfld, add_default, phys_decomp
  use phys_buffer, only: pbuf_times, pbuf_add
  use ppgrid,          only: pcols, pver
  use zm_conv,        only: zm_convi
  use pmgrid,       only: masterproc, plev,plevp
  use error_messages, only: alloc_err	


  implicit none

  real(r8),intent(in) :: hypi(plevp)        ! reference pressures at interfaces
      
  integer  limcnv       ! top interface level limit for convection
  integer k, istat


!
! Allocate space for arrays private to this module
!
     allocate( mu(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'convect_deep_tend', 'mu', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( eu(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'convect_deep_tend', 'eu', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( du(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'convect_deep_tend', 'du', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( md(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'convect_deep_tend', 'md', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( ed(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'convect_deep_tend', 'ed', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( dp(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'convect_deep_tend', 'dp', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( dsubcld(pcols,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'convect_deep_tend', 'dsubcld', &
                      pcols*((endchunk-begchunk)+1) )
     allocate( jt(pcols,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'convect_deep_tend', 'jt', &
                      pcols*((endchunk-begchunk)+1) )
     allocate( maxg(pcols,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'convect_deep_tend', 'maxg', &
                      pcols*((endchunk-begchunk)+1) )
     allocate( ideep(pcols,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'convect_deep_tend', 'ideep', &
                      pcols*((endchunk-begchunk)+1) )
     allocate( lengath(begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'convect_deep_tend', 'lengath', &
                      ((endchunk-begchunk)+1) )


! 
! Register fields with the output buffer
!


    call addfld ('ZMDT    ','K/s     ',pver, 'A','T tendency - Zhang-McFarlane moist convection',phys_decomp)
    call addfld ('ZMDQ    ','kg/kg/s ',pver, 'A','Q tendency - Zhang-McFarlane moist convection',phys_decomp)
    call addfld ('ZMDICE ','kg/kg/s ',pver, 'A','Cloud ice tendency - Zhang-McFarlane convection',phys_decomp)
    call addfld ('ZMDLIQ ','kg/kg/s ',pver, 'A','Cloud liq tendency - Zhang-McFarlane convection',phys_decomp)
    call addfld ('EVAPTZM ','K/s     ',pver, 'A','T tendency - Evaporation from Zhang-McFarlane moist convection',phys_decomp)
    call addfld ('EVAPQZM ','kg/kg/s ',pver, 'A','Q tendency - Evaporation from Zhang-McFarlane moist convection',phys_decomp)
    
    call addfld ('ZMFLXPRC','kg/m2/s ',pverp, 'A','Flux of precipitation from ZM convection'       ,phys_decomp)
    call addfld ('ZMFLXSNW','kg/m2/s ',pverp, 'A','Flux of snow from ZM convection'                ,phys_decomp)
    call addfld ('ZMNTPRPD','kg/kg/s ',pver , 'A','Net precipitation production from ZM convection',phys_decomp)
    call addfld ('ZMNTSNPD','kg/kg/s ',pver , 'A','Net snow production from ZM convection'         ,phys_decomp)
    call addfld ('ZMEIHEAT','W/kg'    ,pver , 'A','Heating by ice and evaporation in ZM convection',phys_decomp)
    
    call addfld ('HKFLXPRC','kg/m2/s ',pverp, 'A','Flux of precipitation from HK convection'       ,phys_decomp)
    call addfld ('HKFLXSNW','kg/m2/s ',pverp, 'A','Flux of snow from HK convection'                ,phys_decomp)
    call addfld ('HKNTPRPD','kg/kg/s ',pver , 'A','Net precipitation production from HK convection',phys_decomp)
    call addfld ('HKNTSNPD','kg/kg/s ',pver , 'A','Net snow production from HK convection'         ,phys_decomp)
    call addfld ('HKEIHEAT','W/kg'    ,pver , 'A','Heating by ice and evaporation in HK convection',phys_decomp)
    
    
!
! Limit deep convection to regions below 40 mb
! Note this calculation is repeated in the shallow convection interface
!
    limcnv = 0   ! null value to check against below
    if (hypi(1) >= 4.e3) then
       limcnv = 1
    else
       do k=1,plev
          if (hypi(k) < 4.e3 .and. hypi(k+1) >= 4.e3) then
             limcnv = k
             exit
          end if
       end do
       if ( limcnv == 0 ) limcnv = plevp
    end if
    
    if (masterproc) then
       write(6,*)'CONVECT_DEEP_INIT: Deep convection will be capped at intfc ',limcnv, &
            ' which is ',hypi(limcnv),' pascals'
    end if
        
    call zm_convi(limcnv)


end subroutine convect_deep_init
!=========================================================================================
!subroutine convect_deep_tend(state, ptend, tdt, pbuf)

subroutine convect_deep_tend(prec    , &
     pblh    ,mcon    ,cme     ,          &
     tpert   ,dlf     ,pflx    ,zdu      , &
     rliq    , &
     ztodt   ,snow    ,&
     state   ,ptend_all   ,pbuf  )
  

   use history,       only: outfld
   use physics_types, only: physics_state, physics_ptend, physics_tend
   use physics_types, only: physics_ptend_init,  physics_tend_init,physics_update
   use physics_types, only: physics_state_copy
   use physics_types, only: physics_ptend_sum

   use phys_grid,     only: get_lat_p, get_lon_p
   use time_manager,  only: get_nstep, is_first_step
   use phys_buffer,   only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
   use constituents, only: pcnst, pnats,ppcnst, cnst_get_ind
   use check_energy,    only: check_energy_chng
   use physconst,       only: gravit

! Arguments
   type(physics_state), intent(in ) :: state          ! Physics state variables
   type(physics_ptend), intent(out) :: ptend_all          ! indivdual parameterization tendencies
   type(pbuf_fld), intent(inout), dimension(pbuf_size_max) :: pbuf  ! physics buffer

   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)
   real(r8), intent(in) :: pblh(pcols)
   real(r8), intent(in) :: tpert(pcols)

   real(r8), intent(out) :: mcon(pcols,pverp)
   real(r8), intent(out) :: dlf(pcols,pver)    ! scattrd version of the detraining cld h2o tend
   real(r8), intent(out) :: pflx(pcols,pverp)  ! scattered precip flux at each level
   real(r8), intent(out) :: cme(pcols,pver)
   real(r8), intent(out) :: zdu(pcols,pver)

   real(r8), intent(out) :: prec(pcols)   ! total precipitation
   real(r8), intent(out) :: snow(pcols)   ! snow from ZM convection 
   real(r8), intent(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals


! Local variables


   integer :: i,k,m
   integer :: ilon                      ! global longitude index of a column
   integer :: ilat                      ! global latitude index of a column
   integer :: nstep
   integer :: ixcldice, ixcldliq              ! constituent indices for cloud liquid and ice water.
   integer :: lchnk                   ! chunk identifier
   integer :: ncol                    ! number of atmospheric columns
   integer itim, ifld  ! for physics buffer fields 

   real(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables
   real(r8) ntprprd(pcols,pver)    ! evap outfld: net precip production in layer
   real(r8) ntsnprd(pcols,pver)    ! evap outfld: net snow production in layer
   real(r8) flxprec(pcols,pverp)   ! evap outfld: Convective-scale flux of precip at interfaces (kg/m2/s)
   real(r8) flxsnow(pcols,pverp)   ! evap outfld: Convective-scale flux of snow   at interfaces (kg/m2/s)

! physics types
    type(physics_state) :: state1        ! locally modify for evaporation to use, not returned
    type(physics_tend ) :: tend          ! Physics tendencies (empty, needed for physics_update call)
    type(physics_ptend)  :: ptend_loc   ! package tendencies

! physics buffer fields 

   real(r8), pointer, dimension(:,:) :: cld
   real(r8), pointer, dimension(:,:) :: ql           ! wg grid slice of cloud liquid water.
   real(r8), pointer, dimension(:,:) :: rprd         ! rain production rate
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble
 
   real(r8), pointer, dimension(:) :: jctop
   real(r8), pointer, dimension(:) :: jcbot

!----------------------------------------------------------------------

   call t_startf ('convect_deep')

!
! initialize 
!
   ftem = 0.   

  call physics_state_copy(state,state1)   ! copy state to local state1.
  call physics_ptend_init(ptend_loc)  ! initialize local ptend type
  call physics_ptend_init(ptend_all)  ! initialize output ptend type
  call physics_tend_init(tend)        ! tend type here is a null place holder


   lchnk = state%lchnk
   ncol = state%ncol

!
! Associate pointers with physics buffer fields
!
   itim = pbuf_old_tim_idx()
   ifld = pbuf_get_fld_idx('CLD')
   cld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

   ifld = pbuf_get_fld_idx('ICWMRDP')
   ql => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
   ifld = pbuf_get_fld_idx('RPRDDP')
   rprd => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
   ifld = pbuf_get_fld_idx('FRACIS')
   fracis  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1:ppcnst)

   ifld = pbuf_get_fld_idx('CLDTOP')
   jctop => pbuf(ifld)%fld_ptr(1,1:pcols,1,lchnk,1)
   ifld = pbuf_get_fld_idx('CLDBOT')
   jcbot => pbuf(ifld)%fld_ptr(1,1:pcols,1,lchnk,1)
!
! Begin with Zhang-McFarlane (1996) convection parameterization
!
   call t_startf ('zm_convr')

   call zm_convr(   lchnk   ,ncol    , &
                    state%t       ,state%q     ,prec    ,jctop   ,jcbot   , &
                    pblh    ,state%zm      ,state%phis    ,state%zi      ,ptend_loc%q(:,:,1)    , &
                    ptend_loc%s    ,state%pmid     ,state%pint    ,state%pdel     , &
                    .5*ztodt    ,mcon    ,cme     ,          &
                    tpert   ,dlf     ,pflx    ,zdu     ,rprd    , &
		    mu(:,:,lchnk),md(:,:,lchnk),du(:,:,lchnk),eu(:,:,lchnk),ed(:,:,lchnk)      , &
                    dp(:,:,lchnk) ,dsubcld(:,lchnk) ,jt(:,lchnk),maxg(:,lchnk),ideep(:,lchnk)   , &
                    lengath(lchnk) ,ql      ,rliq    )


!
! Convert mass flux from reported mb/s to kg/m^2/s
!
   mcon(:ncol,:pver) = mcon(:ncol,:pver) * 100./gravit

   ptend_loc%name  = 'zm_convr'
   ptend_loc%ls    = .TRUE.
   ptend_loc%lq(1) = .TRUE.

   ftem(:ncol,:pver) = ptend_loc%s(:ncol,:pver)/cpair
   call outfld('ZMDT    ',ftem           ,pcols   ,lchnk   )
   call outfld('ZMDQ    ',ptend_loc%q(1,1,1) ,pcols   ,lchnk   )
   call t_stopf ('zm_convr')

  ! add tendency from this process to tendencies from other processes
  call physics_ptend_sum(ptend_loc,ptend_all, state)

  ! update physics state type state1 with ptend_loc 
  call physics_update(state1, tend, ptend_loc, ztodt)

  ! initialize ptend for next process
  call physics_ptend_init(ptend_loc)

   call t_startf ('zm_conv_evap')
!
! Determine the phase of the precipitation produced and add latent heat of fusion
! Evaporate some of the precip directly into the environment (Sundqvist)
! Allow this to use the updated state1 and the fresh ptend_loc type
! heating and specific humidity tendencies produced
!
    ptend_loc%name  = 'zm_conv_evap'
    ptend_loc%ls    = .TRUE.
    ptend_loc%lq(1) = .TRUE.


    call zm_conv_evap(state1%ncol,state1%lchnk, &
         state1%t,state1%pmid,state1%pdel,state1%q(:pcols,:pver,1), &
         ptend_loc%s, ptend_loc%q(:pcols,:pver,1), &
         rprd, cld, ztodt, &
         prec, snow, ntprprd, ntsnprd , flxprec, flxsnow)
!
! Write out variables from zm_conv_evap
!
   ftem(:ncol,:pver) = ptend_loc%s(:ncol,:pver)/cpair
   call outfld('EVAPTZM ',ftem           ,pcols   ,lchnk   )
   call outfld('EVAPQZM ',ptend_loc%q(1,1,1) ,pcols   ,lchnk   )
   call outfld('ZMFLXPRC', flxprec, pcols, lchnk)
   call outfld('ZMFLXSNW', flxsnow, pcols, lchnk)
   call outfld('ZMNTPRPD', ntprprd, pcols, lchnk)
   call outfld('ZMNTSNPD', ntsnprd, pcols, lchnk)
   call outfld('ZMEIHEAT', ptend_loc%s, pcols, lchnk)
   call t_stopf ('zm_conv_evap')

  ! add tendency from this process to tend from other processes here
  call physics_ptend_sum(ptend_loc,ptend_all, state)

  ! update physics state type state1 with ptend_loc 
  call physics_update(state1, tend, ptend_loc, ztodt)

  ! initialize ptend for next process
  call physics_ptend_init(ptend_loc)

   call t_startf ('convtran1')
!
! Transport cloud water and ice only
!


   nstep = get_nstep()

   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   ptend_loc%name = 'convtran1'
   ptend_loc%lq(ixcldice) = .true.
   ptend_loc%lq(ixcldliq) = .true.

   call convtran (lchnk,                                        &
                  ptend_loc%lq(1),state1%q(1,1,1), ppcnst,  mu(1,1,lchnk), md(1,1,lchnk),   &
                  du(1,1,lchnk), eu(1,1,lchnk), ed(1,1,lchnk), dp(1,1,lchnk), dsubcld(1,lchnk),  &
                  jt(1,lchnk),maxg(1,lchnk), ideep(1,lchnk), 1, lengath(lchnk),  &
                  nstep,   fracis,  ptend_loc%q(1,1,1) )  
   call t_stopf ('convtran1')

   call outfld('ZMDICE ',ptend_loc%q(1,1,ixcldice) ,pcols   ,lchnk   )
   call outfld('ZMDLIQ ',ptend_loc%q(1,1,ixcldliq) ,pcols   ,lchnk   )

  ! add tendency from this process to tend from other processes here
  call physics_ptend_sum(ptend_loc,ptend_all, state)

  ! ptend_all will be applied to original state on return to tphysbc
  ! This name triggers a special case in physics_types.F90:physics_update()
  ptend_all%name = 'convect_deep'


  call t_stopf('convect_deep')


end subroutine convect_deep_tend
!=========================================================================================


subroutine convect_deep_tend_2( state,  ptend,  ztodt, pbuf  )

   use physics_types, only: physics_state, physics_ptend, physics_ptend_init
   use time_manager,  only: get_nstep
   use phys_buffer,   only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
   use constituents, only: pcnst, pnats,ppcnst, cnst_get_ind,  cnst_need_pdeldry
   use error_messages, only: alloc_err	
 
! Arguments
   type(physics_state), intent(in ) :: state          ! Physics state variables
   type(physics_ptend), intent(out) :: ptend          ! indivdual parameterization tendencies
   type(pbuf_fld), intent(inout), dimension(pbuf_size_max) :: pbuf  ! physics buffer
   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)

! Local variables
   integer :: i, lchnk, istat
   integer :: nstep
   integer :: ixcldice, ixcldliq              ! constituent indices for cloud liquid and ice water.
   real(r8), dimension(pcols,pver) :: dpdry

! physics buffer fields 
   integer itim, ifld
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble

!
! Initialize
!
  call physics_ptend_init(ptend)

!
! Associate pointers with physics buffer fields
!
   ifld = pbuf_get_fld_idx('FRACIS')
   fracis  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,state%lchnk,1:ppcnst)

!
! Transport all constituents except cloud water and ice
!

  lchnk = state%lchnk

   nstep = get_nstep()

!
!     Convective transport of all trace species except cloud liquid 
!     and cloud ice done here because we need to do the scavenging first
!     to determine the interstitial fraction.
!
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)


   ptend%name  = 'convtran2'
   ptend%lq(:) = .true.
   ptend%lq(ixcldice) = .false.
   ptend%lq(ixcldliq) = .false.

!
! initialize dpdry for call to convtran 
! it is only used for dry mixing ratio tracers
!

   dpdry = 0.
   if (cnst_need_pdeldry ) then
      do i = 1,lengath(lchnk)
         dpdry(i,:) = state%pdeldry(ideep(i,lchnk),:)/100.
      end do
   endif


   call convtran (lchnk,                                        &
                  ptend%lq(1),state%q(1,1,1), ppcnst,  mu(1,1,lchnk), md(1,1,lchnk),   &
                  du(1,1,lchnk), eu(1,1,lchnk), ed(1,1,lchnk), dp(1,1,lchnk), dsubcld(1,lchnk),  &
                  jt(1,lchnk),maxg(1,lchnk),ideep(1,lchnk), 1, lengath(lchnk),  &
                  nstep,   fracis,  ptend%q(1,1,1), dpdry(1,1)  )


end subroutine convect_deep_tend_2

!=========================================================================================



end module
