
#include <misc.h>
#include <params.h>

module convect_shallow
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to the Hack shallow convection scheme
!
! Author: D.B. Coleman
!
!---------------------------------------------------------------------------------
   use shr_kind_mod, only: r8=>shr_kind_r8
   use physconst,    only: cpair                              
   use ppgrid,       only: pver, pcols, pverp
   use zm_conv,      only: zm_conv_evap
   use history,       only: outfld, addfld, add_default, phys_decomp

   implicit none
   private
   save

! Public methods

   public ::&
        convect_shallow_register,           &! register fields in physics buffer
        convect_shallow_init,               &! initialize donner_shallow module
        convect_shallow_tend                 ! return tendencies

!=========================================================================================
contains
!=========================================================================================

subroutine convect_shallow_register

!--------------------------------------------------
! Purpose:  register fields with the physics buffer
!--------------------------------------------------

  use phys_buffer, only: pbuf_add

  integer idx

  call pbuf_add('ICWMRSH', 'physpkg', 1, pver, 1, idx)
  call pbuf_add('RPRDSH',  'physpkg', 1, pver, 1, idx)

  ! CMFDQR is the total rain production from the deep and shallow schemes
  call pbuf_add('RPRDTOT', 'physpkg', 1, pver, 1, idx)
  ! CLDTOP and CLDBOT are combined values from deep and shallow convection.
  call pbuf_add('CLDTOP',  'physpkg', 1,    1, 1, idx)
  call pbuf_add('CLDBOT',  'physpkg', 1,    1, 1, idx)
  
end subroutine convect_shallow_register

!=========================================================================================

subroutine convect_shallow_init(hypi)

!----------------------------------------
! Purpose:  declare output fields, initialize variables needed by convection
!----------------------------------------

  use history,         only: addfld, add_default, phys_decomp
  use ppgrid,          only: pcols, pver
  use hk_conv,        only: mfinti
  use physconst,      only: rair    ,gravit  ,latvap  ,rhoh2o
  use pmgrid,       only: masterproc, plev,plevp

  implicit none

  real(r8),intent(in) :: hypi(plevp)        ! reference pressures at interfaces
    

  integer limcnv          ! top interface level limit for convection
  integer k


    call addfld ('CMFDT   ','K/s     ',pver, 'A','T tendency - Hack convection',phys_decomp)
    call addfld ('CMFDQ   ','kg/kg/s ',pver, 'A','Q tendency - Hack convection',phys_decomp)
    call addfld ('CMFDICE ','kg/kg/s ',pver, 'A','Cloud ice tendency - Hack convection',phys_decomp)
    call addfld ('CMFDLIQ ','kg/kg/s ',pver, 'A','Cloud liq tendency - Hack convection',phys_decomp)
    call addfld ('CMFDQR  ','kg/kg/s ',pver, 'A','Q tendency - shallow convection rainout',phys_decomp)
    call addfld ('EVAPTCM ','K/s     ',pver, 'A','T tendency - Evaporation from Hack convection',phys_decomp)
    call addfld ('EVAPQCM ','kg/kg/s ',pver, 'A','Q tendency - Evaporation from Hack convection',phys_decomp)
    call addfld ('QC      ','kg/kg/s ',pver, 'A','Q tendency - shallow convection LW export',phys_decomp)
    call addfld ('PRECSH  ','m/s     ',1,    'A','Shallow Convection precipitation rate',phys_decomp)
    call addfld ('CMFMC   ','kg/m2/s ',pverp,'A','Moist convection mass flux',phys_decomp)
    call addfld ('CMFSL   ','W/m2    ',pverp,'A','Moist convection liquid water static energy flux',phys_decomp)
    call addfld ('CMFLQ   ','W/m2    ',pverp,'A','Moist convection total water flux',phys_decomp)

    call add_default ('CMFDT   ', 1, ' ')
    call add_default ('CMFDQ   ', 1, ' ')
    call add_default ('CMFDQR  ', 1, ' ')
    call add_default ('QC      ', 1, ' ')
    call add_default ('PRECSH  ', 1, ' ')
    call add_default ('CMFMC   ', 1, ' ')
    

!
! Limit shallow convection to regions below 40 mb
! Note this calculation is repeated in the deep convection interface
!
   if (hypi(1) >= 4.e3) then
      limcnv = 1
   else
      do k=1,plev
         if (hypi(k) < 4.e3 .and. hypi(k+1) >= 4.e3) then
            limcnv = k
            goto 10
         end if
      end do
      limcnv = plevp
   end if

10 continue

   if (masterproc) then
      write(6,*)'MFINTI: Convection will be capped at intfc ',limcnv, &
                ' which is ',hypi(limcnv),' pascals'
   end if



    call mfinti(rair    ,cpair   ,gravit  ,latvap  ,rhoh2o, limcnv) ! get args from inti.F90
    

end subroutine convect_shallow_init
!=========================================================================================

subroutine convect_shallow_tend( ztodt     , &
     qpert   ,&
     pblht   ,&
     cmfmc   ,cmfmc2  ,precc   , &
     qc      ,rliq    ,rliq2   , & 
     snow    ,state, ptend_all    , pbuf    )

   use history,       only: outfld
   use physics_types, only: physics_state, physics_ptend, physics_tend
   use physics_types, only: physics_ptend_init,  physics_tend_init,physics_update
   use physics_types, only: physics_state_copy
   use physics_types, only: physics_ptend_sum

   use phys_buffer,   only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
   use constituents, only: pcnst, pnats, cnst_get_ind
   use hk_conv,       only: cmfmca
   use time_manager,  only: get_nstep
    

! Arguments
   type(physics_state), intent(in ) :: state          ! Physics state variables
   type(physics_ptend), intent(out) :: ptend_all         ! indivdual parameterization tendencies
   type(pbuf_fld), intent(inout), dimension(pbuf_size_max) :: pbuf  ! physics buffer
   real(r8), intent(in) :: ztodt               ! 2 delta-t (seconds)

!
! Input arguments
!
   real(r8), intent(in) :: pblht(pcols)        ! PBL height (provided by PBL routine)
   real(r8), intent(inout) :: qpert(pcols,pcnst+pnats)  ! PBL perturbation specific humidity

!
! fields which combine deep convection with shallow
!
   real(r8), intent(inout) :: cmfmc(pcols,pverp)  ! moist deep + shallow convection cloud mass flux
   real(r8), intent(inout) :: qc(pcols,pver)      ! dq/dt due to export of cloud water
   real(r8), intent(inout) :: rliq(pcols) 

!
! Output arguments
!

   real(r8), intent(out) :: cmfmc2(pcols,pverp)  ! moist shallow convection cloud mass flux
   real(r8) cmfsl(pcols,pver )  ! convective lw static energy flux
   real(r8) :: cmflq(pcols,pver )  ! convective total water flux
   real(r8), intent(out) :: precc(pcols)        ! convective precipitation rate
   real(r8), intent(out)  :: rliq2(pcols) 
   real(r8), intent(out) :: snow(pcols)  ! snow from this convection 


! Local variables

   integer :: i,k,m
   integer :: n,x
   integer :: ilon                      ! global longitude index of a column
   integer :: ilat                      ! global latitude index of a column
   integer :: lchnk                ! chunk identifier
   integer :: ncol                 ! number of atmospheric columns
   integer :: nstep                ! current time step index
   integer :: ixcldice, ixcldliq   ! constituent indices for cloud liquid and ice water.
   real(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables
   real(r8) :: qc2(pcols,pver)      ! dq/dt due to export of cloud water
   real(r8) :: cnt2(pcols)          ! top level of convective activity
   real(r8) :: cnb2(pcols)          ! bottom level of convective activity
   real(r8) tpert(pcols)        ! PBL perturbation theta
   real(r8) ntprprd(pcols,pver)    ! evap outfld: net precip production in layer
   real(r8) ntsnprd(pcols,pver)    ! evap outfld: net snow production in layer
   real(r8) flxprec(pcols,pverp)   ! evap outfld: Convective-scale flux of precip at interfaces (kg/m2/s)
   real(r8) flxsnow(pcols,pverp)   ! evap outfld: Convective-scale flux of snow   at interfaces (kg/m2/s)

! physics types
   type(physics_state) :: state1        ! locally modify for evaporation to use, not returned
   type(physics_ptend) :: ptend_loc     ! local tend from processes, added up to return as ptend_all
   type(physics_tend ) :: tend          ! Physics tendencies (empty, needed for physics_update call)

! physics buffer fields 
   integer itim, ifld
   real(r8), pointer, dimension(:,:) :: cld
   real(r8), pointer, dimension(:,:) :: icwmr    ! in cloud water mixing ratio
   real(r8), pointer, dimension(:,:) :: rprddp   ! dq/dt due to deep convective rainout
   real(r8), pointer, dimension(:,:) :: rprdsh   ! dq/dt due to deep and shallow convective rainout
   real(r8), pointer, dimension(:)   :: cnt
   real(r8), pointer, dimension(:)   :: cnb

!----------------------------------------------------------------------

   call t_startf('convect_shallow')


   nstep = get_nstep()

!
! initialize 
!
   lchnk = state%lchnk
   ncol = state%ncol
  

  call physics_state_copy(state,state1)   ! copy state to local state1.
  call physics_ptend_init(ptend_loc)  ! initialize local ptend type
  call physics_ptend_init(ptend_all)  ! initialize output ptend type
  call physics_tend_init(tend)        ! tend type here is a null place holder

!
! Associate pointers with physics buffer fields
!

   itim = pbuf_old_tim_idx()
   ifld = pbuf_get_fld_idx('CLD')
   cld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

   ifld = pbuf_get_fld_idx('ICWMRSH')
   icwmr => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)

   ifld = pbuf_get_fld_idx('RPRDDP')
   rprddp => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
   ifld = pbuf_get_fld_idx('RPRDSH')
   rprdsh => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)

   ifld = pbuf_get_fld_idx('CLDTOP')
   cnt => pbuf(ifld)%fld_ptr(1,1:pcols,1,lchnk,1)
   ifld = pbuf_get_fld_idx('CLDBOT')
   cnb => pbuf(ifld)%fld_ptr(1,1:pcols,1,lchnk,1)
 
! Initialize

   tpert(:ncol  ) =0.
   qpert(:ncol,2:pcnst+pnats) = 0.0

!
!  Call Hack convection 
!
   call cmfmca (lchnk,   ncol, &
                nstep,   ztodt,   state%pmid,  state%pdel,   &
                state%rpdel,   state%zm,      tpert,  qpert,  state%phis,     &
                pblht,   state%t,   state%q,   ptend_loc%s,   ptend_loc%q,      &
                cmfmc2,  rprdsh, cmfsl,  cmflq,  precc,   &
                qc2,     cnt2,    cnb2,    icwmr   , rliq2, & 
                state%pmiddry, state%pdeldry, state%rpdeldry)
   
   ptend_loc%name  = 'cmfmca'
   ptend_loc%ls    = .TRUE.
   ptend_loc%lq(:) = .TRUE.

   !
   ! Merge shallow/mid-level output with prior results from deep
   !

   !combine cmfmc2 (from shallow)  with cmfmc (from deep)
   cmfmc (:ncol,:pver) = cmfmc (:ncol,:pver) + cmfmc2 (:ncol,:pver)

   ! cnt2, cnb2 are from shallow, cnt & cnb are from deep
   do i=1,ncol
      if (cnt2(i) < cnt(i)) cnt(i) = cnt2(i)
      if (cnb2(i) > cnb(i)) cnb(i) = cnb2(i)
   end do
 
   ! This quantity was previously known as CMFDQR.  Now CMFDQR is the shallow rain production only.
   ifld = pbuf_get_fld_idx( 'RPRDTOT' )
   pbuf(ifld)%fld_ptr(1,1:ncol,1:pver,lchnk,1) = rprdsh(:ncol,:pver) + rprddp(:ncol,:pver)
 
   ! Add shallow cloud water detrainment to cloud water detrained from deep convection
   qc(:ncol,:pver) = qc(:ncol,:pver) + qc2(:ncol,:pver)
   rliq(:ncol) = rliq(:ncol) + rliq2(:ncol)    

   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   ftem(:ncol,:pver) = ptend_loc%s(:ncol,:pver)/cpair
   call outfld('CMFDT   ',ftem          ,pcols   ,lchnk   )
   call outfld('CMFDQ   ',ptend_loc%q(1,1,1),pcols   ,lchnk   )
   call outfld('CMFDICE ',ptend_loc%q(1,1,ixcldice),pcols   ,lchnk   )
   call outfld('CMFDLIQ ',ptend_loc%q(1,1,ixcldliq),pcols   ,lchnk   )
   call outfld('CMFMC' , cmfmc  , pcols, lchnk)
!  output new partition of cloud condensate variables, as well as precipitation 
   call outfld('QC      ',qc2    ,pcols, lchnk)
   call outfld('CMFDQR', rprdsh,  pcols, lchnk)
   call outfld('CMFSL' , cmfsl  , pcols, lchnk)
   call outfld('CMFLQ' , cmflq  , pcols, lchnk)
   call outfld('DQP'   , qc2    , pcols, lchnk)


 ! add tendency from this process to tend from other processes here
  call physics_ptend_sum(ptend_loc,ptend_all, state)

  ! update physics state type state1 with ptend_loc 
  call physics_update(state1, tend, ptend_loc, ztodt)

  ! initialize ptend for next process
  call physics_ptend_init(ptend_loc)

!
! Determine the phase of the precipitation produced and add latent heat of fusion
! Evaporate some of the precip directly into the environment (Sundqvist)
! Allow this to use the updated state1 and a fresh ptend_loc type
! Heating and specific humidity tendencies produced
!
    ptend_loc%name  = 'zm_conv_evap'
    ptend_loc%ls    = .TRUE.
    ptend_loc%lq(1) = .TRUE.
   call zm_conv_evap(state1%ncol,state1%lchnk, &
        state1%t,state1%pmid,state1%pdel,state1%q(:pcols,:pver,1), &
        ptend_loc%s, ptend_loc%q(:pcols,:pver,1), &
        rprdsh, cld, ztodt, &
        precc, snow, ntprprd, ntsnprd , flxprec, flxsnow )


! record history variables from zm_conv_evap
   ftem(:ncol,:pver) = ptend_loc%s(:ncol,:pver)/cpair
   call outfld('EVAPTCM ',ftem           ,pcols   ,lchnk   )
   call outfld('EVAPQCM ',ptend_loc%q(1,1,1),pcols   ,lchnk   )
   call outfld('PRECSH  ',precc      ,pcols   ,lchnk       )
   call outfld('HKFLXPRC', flxprec, pcols, lchnk)
   call outfld('HKFLXSNW', flxsnow, pcols, lchnk)
   call outfld('HKNTPRPD', ntprprd, pcols, lchnk)
   call outfld('HKNTSNPD', ntsnprd, pcols, lchnk)
   call outfld('HKEIHEAT', ptend_loc%s, pcols, lchnk)

      
 ! add tendency from this process to tend from other processes here
  call physics_ptend_sum(ptend_loc,ptend_all, state)

 ! update name of parameterization tendencies to send to tphysbc
  ptend_all%name = 'convect_shallow'

  call t_stopf('convect_shallow')


end subroutine convect_shallow_tend
!=========================================================================================


end module
