#include <params.h>
#include <max.h>
!------------------------------------------------------------------------
! File: readiopdata.F 
! Author: John Truesdale (jet@ucar.edu) 
! $Id: readiopdata.F90,v 1.1.6.3 2004/09/23 17:46:59 erik Exp $
!
!------------------------------------------------------------------------
subroutine readiopdata( error_code )


!-----------------------------------------------------------------------
!     
!     Open and read netCDF file containing initial IOP  conditions
!     
!---------------------------Code history--------------------------------
!     
!     Written by J.  Truesdale    August, 1996, revised January, 1998
!     
!-----------------------------------------------------------------------
	use pmgrid
	use prognostics
	use buffer
	use comsrf
	use phys_grid, only: clat_p
	use commap
        use time_manager, only : get_nstep
        use constituents, only : readtrace,cnst_get_ind
        use string_utils, only: to_lower
        use scamMod
        use iop
        use getnetcdfdata
!-----------------------------------------------------------------------
   implicit none
#if ( defined RS6000 )
   implicit automatic ( a-z )
#endif
#include <runtype.h>
!------------------------------Commons----------------------------------
#include <comhyb.h>
!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comfrc.h>
!-----------------------------------------------------------------------
#include <netcdf.inc>

!------------------------------Inputs-----------------------------------

   integer error_code        ! returns netcdf errors

!------------------------------Locals-----------------------------------
!     
   integer NCID, STATUS
   integer time_dimID, lev_dimID,  lev_varID
   integer tsec_varID, bdate_varID
   integer i,j
   integer nlev
   integer total_levs

   integer bdate, tsec( MAX_TIME_DIM ), ntime
   integer k, m
   integer icldliq,icldice

   logical have_srf              ! value at surface is available
   logical use_nf_real           ! nctype for 4byte real
   logical fill_ends             ! 
   logical have_dcldliq,have_dcldice
   real(r8) dplevs( MAX_DATASET_LEVS +1 )
   real(r8) dummy
   real(r8) lat,xlat
   real(r8) srf(1)                  ! value at surface
   real(r8) ptendarr(1)                  ! value at surface
   real(r8) pmid(plev)  ! pressure at model levels (time n)
   real(r8) pint(plevp) ! pressure at model interfaces (n  )
   real(r8) pdel(plev)  ! pdel(k)   = pint  (k+1)-pint  (k)
   real(r8) weight
   real(r8) tmpdata(1)
   real(r8) coldata(plev)

   character(len=8) lowername

#if ( defined sun )
!     Trap ieee exceptions on SUN for debugging purposes
   external myhandler
   integer iexcept, ieee_handler, myhandler
   iexcept = ieee_handler( 'set', 'common', myhandler )
   if ( iexcept .ne. 0 ) write( 6,* )'ieee trapping not supported here'
#endif

#if USE_4BYTE_REAL
   use_nf_real = .true.
#else
   use_nf_real = .false.
#endif

   fill_ends= .false.
   if ( use_userdata ) then 
      error_code = USER
   else
      error_code = IOP
   endif

!     
!     Open IOP dataset
!     
   STATUS = NF_OPEN( iopfile, NF_NOWRITE, NCID )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readiopdata.F:Cant open iop dataset: ' ,iopfile
      return
   end if

!
!     if the dataset is a CCM generated dataset set use_ccmiop to true
!       CCM IOP datasets have a global attribute called CCM_GENERATED_IOP      
!
   if ( nf_inq_attid( NCID, NF_GLOBAL, 'CAM_GENERATED_FORCING', i ).EQ. NF_NOERR ) then
      use_ccmiop = .true.
   else
      use_ccmiop = .false.
   endif

!=====================================================================
!     
!     Read time variables


   STATUS = NF_INQ_DIMID( NCID, 'time', time_dimID )
   if ( STATUS .NE. NF_NOERR )  then
     STATUS = NF_INQ_DIMID( NCID, 'tsec', time_dimID )
     if ( STATUS .NE. NF_NOERR )  then
       write( 6,* )'ERROR - readiopdata.F:Could not find variable dim ID for time'
       STATUS = NF_CLOSE ( NCID )
       return
     end if
   end if

   STATUS = NF_INQ_DIMLEN( NCID, time_dimID, ntime )
   if ( STATUS .NE. NF_NOERR )  then
     write( 6,* )'ERROR - readiopdata.F:Could not find variable dim len for time'
     STATUS = NF_CLOSE ( NCID )
     return
   end if


   STATUS = NF_INQ_VARID( NCID, 'tsec', tsec_varID )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readiopdata.F:Could not find variable ID for tsec'
      STATUS = NF_CLOSE ( NCID )
      return
   end if


   STATUS = NF_GET_VAR_INT( NCID, tsec_varID, tsec )
   if ( STATUS .NE. NF_NOERR ) then
     write( 6,* )'ERROR - readiopdata.F:Could not read variable tsec'
     STATUS = NF_CLOSE ( NCID )
     return
   end if

   STATUS = NF_GET_VAR_INT( NCID, tsec_varID, tsec )
   if ( STATUS .NE. NF_NOERR ) then
     write( 6,* )'ERROR - readiopdata.F:Could not read variable tsec'
     STATUS = NF_CLOSE ( NCID )
     return
   end if

   STATUS = NF_INQ_VARID( NCID, 'nbdate', bdate_varID )
   if ( STATUS .NE. NF_NOERR ) then
      STATUS = NF_INQ_VARID( NCID, 'bdate', bdate_varID )
      if ( STATUS .NE. NF_NOERR ) then
         write( 6,* )'ERROR - readiopdata.F:Could not find variable ID for bdate'
         STATUS = NF_CLOSE ( NCID )
         return
      end if
   end if

   STATUS = NF_GET_VAR_INT( NCID, bdate_varID, bdate )
   if ( STATUS .NE. NF_NOERR )then
     write( 6,* )'ERROR - readiopdata.F:Could not find variable bdate'
     STATUS = NF_CLOSE ( NCID )
     return
   end if

!     
!======================================================
!     read level data
!     
   STATUS = NF_INQ_DIMID( NCID, 'lev', lev_dimID )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readiopdata.F:Could not find variable dim ID  for lev'
      STATUS = NF_CLOSE ( NCID )
      return
   end if

   STATUS = NF_INQ_DIMLEN( NCID, lev_dimID, nlev )

   STATUS = NF_INQ_VARID( NCID, 'lev', lev_varID )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readiopdata.F:Could not find variable ID for lev'
      STATUS = NF_CLOSE ( NCID )
      return
   end if

   if (use_nf_real) then
      STATUS = NF_GET_VAR_REAL( NCID, lev_varID, dplevs )
   else
      STATUS = NF_GET_VAR_DOUBLE( NCID, lev_varID, dplevs )
   endif

   if (STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - readiopdata.F:Could not find variable lev'
      STATUS = NF_CLOSE( NCID )
      return
   endif
!
!CAM generated forcing already has pressure on millibars
!
   if (.not. use_ccmiop) then
!
!     convert pressure to millibars ( lev is expressed in pascals in iop datasets )
!
      do i=1,nlev
         dplevs( i ) = dplevs( i )/100.
      end do
   endif

   call getncdata ( NCID, 1, 1, iopTimeIdx, 'Ps', tmpdata, STATUS )
   ps(1,1,n3) = tmpdata(1)
   if ( STATUS .NE. NF_NOERR ) then
      have_ps = .false.
      if (get_nstep() .eq. 0 ) write(6,*)'Could not find variable Ps'
      if ( .not. use_userdata ) then
         STATUS = NF_CLOSE( NCID )
         return
      else
         if ( get_nstep() .eq. 0 ) write(*,*) 'Using value from Analysis Dataset'
      endif
   else
      have_ps = .true.
   endif


!$$$JP  for reproducing CCM output don't do interpolation.
!$$$JP  the most expedient way of doing this is to set      
!$$$JP  the dataset pressure levels to the current
!$$$JP  scam model levels
	
   if ( use_ccmiop ) then
      do i = 1, plev
         dplevs( i ) = 1000.0 * hyam( i ) + ps(1,1,n3) * hybm( i ) / 100.0
      end do
   endif

!     add the surface pressure to the pressure level data, so that
!     surface boundary condition will be set properly,
!     making sure that it is the highest pressure in the array.
!
   total_levs = nlev+1
   dplevs(nlev+1) = ps(1,1,n3)/100 ! ps is expressed in pascals
   do i= nlev, 1, -1
      if ( dplevs(i) .GT. ps(1,1,n3)/100.0) then
         total_levs = i
         dplevs(i) = ps(1,1,n3)/100
      end if
   end do
   if (.not. use_ccmiop ) then
      nlev = total_levs
   endif
   if ( nlev .eq. 1 ) then
      write(*,*) 'Error - Readiopdata.F: Ps too low!'
      return
   endif


!     ===================================================================
!     
!     Get variable IDs/Data
!     
!     
!     Extract the data for each of the fields.
!     
   if ( get_nstep() .eq. 0  ) then
      call getncdata ( NCID, 1, 1, ioptimeidx,'phis', phis, STATUS)
      if ( STATUS .NE. NF_NOERR ) then
         have_phis = .false.
         if ( get_nstep() .eq. 0 ) write(6,*)'Could not find variable phis'
         if ( .not. use_userdata ) then
            STATUS = NF_CLOSE( NCID )
            return
         else
            if ( get_nstep() .eq. 0 ) write(*,*) 'Using value from Analysis Dataset'
         endif
      else
         have_phis = .true.
      endif
   endif

!=====================================================================
!     

   call getncdata( NCID, 1, 1, ioptimeidx,'Tsair', tsair, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_tsair = .false.
   else
      have_tsair = .true.
   endif
!
!      read in Tobs  For cam generated iop readin small t to avoid confusion
!      with capital T defined in cam
!
!  We want to overlay the iop data on top of the analysis data already read
!  in.  Since IOP data can be on only a subset of model levels we need to
!  fill in missing levels with reasonable data.  T3 already contains data
!  read in from read_inidat.
!  read tobs t3(:,:,:,1) initialized on start from read_inidat
!  otherwise t3(:,:,:,n3) is the time level we want when this
!  routine is called from within the timestep

   if ( get_nstep() .eq. 0 ) then
     tobs(:)=t3(1,:,1,1)
   else	
     tobs(:)=t3(1,:,1,n3)
   endif


   if ( use_ccmiop ) then
     call getinterpncdata( NCID, 1, 1, ioptimeidx,'t', have_tsair, &
          tsair(1), fill_ends, &
          dplevs, nlev, tobs, STATUS )
   else
     call getinterpncdata( NCID, 1, 1, ioptimeidx,'T', have_tsair, &
          tsair(1), fill_ends, &
          dplevs, nlev, tobs, STATUS )
   endif
   if ( STATUS .NE. NF_NOERR ) then
      have_t = .false.
      if ( get_nstep() .eq. 0 ) write(6,*)'Could not find variable T'
      if ( .not. use_userdata ) then
         STATUS = NF_CLOSE( NCID )
         return
      else
         if ( get_nstep() .eq. 0  ) write(*,*) 'Using value from Analysis Dataset'
      endif
!     
!     set T3 to Tobs on first time step
!     
   else
      have_t = .true.
      if ( get_nstep() .eq. 0 .and. .not.isrestart ) then
         do i=1, plev
            t3(1,i,1,n3) = tobs(i)
         end do
      endif
   endif

   call getncdata (NCID, 1, 1, ioptimeidx,'Tg', tground, STATUS)
   if (STATUS .NE. NF_NOERR) then
      if ( get_nstep() .eq. 0  ) then
         write(6,*)'Could not find variable Tg'
      endif
      if ( have_tsair ) then
         if ( get_nstep() .eq. 0  ) then
            write(6,*) 'Using Tsair'
         endif
         tground = tsair     ! use surface value from T field
      else
         if ( get_nstep() .eq. 0 ) then
            write(6,*) 'Using T at lowest level'
         endif
         tground = t3(1,nlev,1,n3)
      endif
   else
      have_tg = .true.
   endif

   if ( get_nstep() .eq. 0  ) then
      srfflx_state2d(begchunk)%ts(1) = tground(1)
      surface_state2d(begchunk)%tssub(1,1) = tground(1)
      surface_state2d(begchunk)%tssub(1,2) = tground(1)
      surface_state2d(begchunk)%tssub(1,3) = tground(1)
      surface_state2d(begchunk)%tssub(1,4) = tground(1)
   endif

   call getncdata( NCID, 1, 1, ioptimeidx, 'qsrf', srf(1), STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_srf = .false.
   else
      have_srf = .true.
   endif
!
!     read Qobs q3(:,:,:,1) initialized on start from read_inidat
!     otherwise q3(:,:,:,n3) is the time level we want when this
!     routine is called from within the timestep
!
   if ( get_nstep() .le. 1  ) then
     qobs(:)= q3(1,:,1,1,1)
   else
     qobs(:)= q3(1,:,1,1,n3)
   endif
   call getinterpncdata( NCID, 1, 1, ioptimeidx,  'q', have_srf, &
      srf(1), fill_ends, &
      dplevs, nlev, qobs, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_q = .false.
      if ( get_nstep() .eq. 0 )  write(6,*)'Could not find variable q'
      if ( .not. use_userdata ) then
         STATUS = NF_CLOSE( NCID )
         return
      else
         if ( get_nstep() .eq. 0 )  write(*,*) 'Using values from Analysis Dataset'
      endif
   else
      have_q = .true.
!
!     set Q3 to Qobs at first time step
!         
      if ( get_nstep() .eq. 0  .and. .not.isrestart) then
         do i=1, PLEV
            q3(1,i,1,1,n3) = qobs(i)
         end do
      endif
   endif

   call getinterpncdata( NCID, 1, 1, ioptimeidx,  'cld', .false., &
      dummy, fill_ends, dplevs, nlev, cldobs, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      !      write(6,*)'Could not read cld from iop data - leaving as is.'
      have_cld = .false.
   else
      have_cld = .true.
   endif
   
   call getinterpncdata( NCID, 1, 1, ioptimeidx,  'clwp', .false., &
      dummy, fill_ends, dplevs, nlev, clwpobs, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      !      write(6,*)'Could not read clwp from iop data - leaving as is.'
      have_clwp = .false.
   else
      have_clwp = .true.
   endif

!
!	read divq (horizontal advection)
!      
   call getncdata( NCID, 1, 1, ioptimeidx, 'divqsrf', srf(1), STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_srf = .false.
   else
      have_srf = .true.
   endif

   call getinterpncdata( NCID, 1, 1, ioptimeidx, &
      'divq', have_srf, srf(1), fill_ends, &
      dplevs, nlev, divq, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_divq = .false.
      if ( get_nstep() .eq. 0 ) then
         write(6,*)'Could not find variable divq'
         if ( use_userdata ) then
            write(*,*) 'Using value from Analysis Dataset'
         endif
      endif
   else
      have_divq = .true.
   endif

!
!     read vertdivq if available
!
   call getncdata( NCID, 1, 1, ioptimeidx, 'vertdivqsrf', &
      srf, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_srf = .false.
   else
      have_srf = .true.
   endif

   call getinterpncdata( NCID, 1, 1, ioptimeidx, 'vertdivq', &
      have_srf, srf(1), fill_ends, &
      dplevs, nlev, vertdivq, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_vertdivq = .false.
   else
      have_vertdivq = .true.
   endif
!
!     read divq3d (combined vertical/horizontal advection)
!     if available

   call getncdata( NCID, 1, 1, ioptimeidx, 'divq3dsrf', &
      srf(1), STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_srf = .false.
   else
      have_srf = .true.
   endif

   call getinterpncdata( NCID, 1, 1, ioptimeidx, 'divq3d', &
      have_srf, srf(1), fill_ends, &
      dplevs, nlev, divq3d, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_divq3d = .false.
      divq3d(1:,1)=0.
   else
      have_divq3d = .true.
   endif


   call cnst_get_ind('CLDLIQ', icldliq)
   have_srf = .false.
   call getinterpncdata( NCID, 1, 1, ioptimeidx, 'dcldliq', &
      have_srf, srf(1), fill_ends, &
      dplevs, nlev, divq3d(1,icldliq), STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_dcldliq = .false.
      divq3d(1:,icldliq)=0.
   else
      have_dcldliq = .true.
   endif

   call cnst_get_ind('CLDICE', icldice)
   have_srf = .false.
   call getinterpncdata( NCID, 1, 1, ioptimeidx, 'dcldice', &
      have_srf, srf(1), fill_ends, &
      dplevs, nlev, divq3d(1,icldice), STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_dcldice = .false.
      divq3d(1:,icldice)=0.
   else
      have_dcldice = .true.
   endif

!      if ( get_nstep() .eq. 0  .and. .not.isrestart) then
      if ( get_nstep() .eq. 0) then
         
         do i=1, PLEV
            q3(1,i,icldice,1,n3) = divq3d(i,icldice)
            q3(1,i,icldliq,1,n3) = divq3d(i,icldliq)
         end do
      endif
!jt    write(6,*)'q3 after filled with divq3d=',q3
!
!	read divu (optional field)
!      
   call getncdata( NCID, 1, 1, ioptimeidx, 'divusrf', &
      srf(1), STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_srf = .false.
   else
      have_srf = .true.
   endif

   call getinterpncdata( NCID, 1, 1, ioptimeidx, 'divu', &
      have_srf, srf(1), fill_ends, &
      dplevs, nlev, divu, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_divu = .false.
   else
      have_divu = .true.
   endif
!
!	read divv (optional field)
!      
   call getncdata( NCID, 1, 1, ioptimeidx, 'divvsrf', &
      srf, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_srf = .false.
   else
      have_srf = .true.
   endif

   call getinterpncdata( NCID, 1, 1, ioptimeidx, 'divv', &
      have_srf, srf(1), fill_ends, &
      dplevs, nlev, divv, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_divv = .false.
   else
      have_divv = .true.
   endif
!     
!     read divt
!      
   call getncdata( NCID, 1, 1, ioptimeidx, 'divTsrf', &
      srf, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_srf = .false.
   else
      have_srf = .true.
   endif

   call getinterpncdata( NCID, 1, 1, ioptimeidx, &
      'divT', have_srf, srf(1), fill_ends, &
      dplevs, nlev, divt, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_divt = .false.
      if ( get_nstep() .eq. 0 ) then
         write(6,*)'Could not find variable divT'
         if ( use_userdata ) write(*,*) 'Using value from Analysis Dataset'
      endif
   else
      have_divt = .true.
   endif

!
!     read vertdivt if available
!
   call getncdata( NCID, 1, 1, ioptimeidx, 'vertdivTsrf', &
      srf, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_srf = .false.
   else
      have_srf = .true.
   endif

   call getinterpncdata( NCID, 1, 1, ioptimeidx, 'vertdivT', &
      have_srf, srf(1), fill_ends, &
      dplevs, nlev, vertdivt, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_vertdivt = .false.
   else
      have_vertdivt = .true.
   endif
!
!	read divt3d (combined vertical/horizontal advection)
!      (optional field)




   call getncdata( NCID, 1, 1, ioptimeidx, 'divT3dsrf',  &
      srf, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_srf = .false.
   else
      have_srf = .true.
   endif

   call getinterpncdata( NCID, 1, 1, ioptimeidx, 'divT3d', &
      have_srf, srf(1), fill_ends, &
      dplevs, nlev, divt3d, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_divt3d = .false.
   else
      have_divt3d = .true.
   endif


!     zero out omega field before interpolation if at first time step
   if ( get_nstep() .eq. 0 )then
      do i=1, plev
         wfld( i ) = 0
      end do
   endif

   call getncdata( NCID, 1, 1, ioptimeidx, 'Ptend', &
      ptendarr, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_ptend = .false.
      if ( get_nstep() .eq. 0 ) &
         write(6,*)'Could not find variable Ptend. Setting to zero'
      ptend = 0.0
   else
      have_ptend = .true.
      ptend= ptendarr(1)
   endif

   call getinterpncdata( NCID, 1, 1, ioptimeidx, &
      'omega', .true., ptend, fill_ends, &
      dplevs, nlev, wfld, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_omega = .false.
      if ( get_nstep() .eq. 0 )  write(6,*)'Could not find variable omega'
      if ( .not. use_userdata ) then
         STATUS = NF_CLOSE( NCID )
         return
      else
         if ( get_nstep() .eq. 0 ) write(*,*) 'Using value from Analysis Dataset'
      endif
   else
      have_omega = .true.
   endif
   call plevs0(1    ,plon   ,plev    ,ps(1,1,n3)   ,pint,pmid ,pdel)
!
! Build interface vector for the specified omega profile
! (weighted average in pressure of specified level values)
!
   wfldh(1) = 0.0

   do k=2,plev
      weight = (pint(k) - pmid(k-1))/(pmid(k) - pmid(k-1))
      wfldh(k) = (1.0 - weight)*wfld(k-1) + weight*wfld(k)
   end do

   wfldh(plevp) = 0.0


   call getncdata( NCID, 1, 1, ioptimeidx, 'usrf',  srf, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_srf = .false.
   else
      have_srf = .true.
   endif

   call getinterpncdata( NCID, 1, 1, ioptimeidx, &
      'u', have_srf, srf(1), fill_ends, &
      dplevs, nlev, uobs, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_u = .false.
   else
      have_u = .true.
      if ( get_nstep() .eq. 0 ) then
         do i=1, PLEV
            u3(1,i,1,n3) = uobs(i)  !     set u to uobs at first time step
         end do
      endif
   endif

   call getncdata( NCID, 1, 1, ioptimeidx, 'vsrf', srf, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_srf = .false.
   else
      have_srf = .true.
   endif

   call getinterpncdata( NCID, 1, 1, ioptimeidx, &
      'v', have_srf, srf(1), fill_ends, &
      dplevs, nlev, vobs, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_v = .false.
   else
      have_v = .true.
      if ( get_nstep() .eq. 0 ) then
         do i=1, PLEV
            v3(1,i,1,n3) = vobs(i)  !     set v to vobs at first time step
         end do
      endif
   endif


   call getncdata( NCID, 1, 1, ioptimeidx, 'Prec', precobs, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_prec = .false.
   else
      have_prec = .true.
   endif

   call getinterpncdata( NCID, 1, 1, ioptimeidx, 'Q1', &
      .false., dummy, fill_ends, & ! datasets don't contain Q1 at surface
      dplevs, nlev, q1obs, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_q1 = .false.
   else
      have_q1 = .true.
   endif

   call getinterpncdata( NCID, 1, 1, ioptimeidx, 'Q2', &
      .false., dummy, fill_ends, & ! datasets don't contain Q2 at surface
      dplevs, nlev, q1obs, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_q2 = .false.
   else
      have_q2 = .true.
   endif

   call getncdata( NCID, 1, 1, ioptimeidx, 'lhflx', lhflxobs, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_lhflx = .false.
   else
      have_lhflx = .true.
   endif

   call getncdata( NCID, 1, 1, ioptimeidx, 'shflx', shflxobs, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      have_shflx = .false.
   else
      have_shflx = .true.
   endif
!
!     old names - backwards compatibility
!
   call getncdata( NCID, 1, 1, ioptimeidx,  'sh', shflxobs, STATUS )
   if ( STATUS .EQ. NF_NOERR ) then
      have_shflx = .true.
   endif


   call getncdata( NCID, 1, 1, ioptimeidx,  'lh', lhflxobs, STATUS )
   if ( STATUS .EQ. NF_NOERR ) then
      have_lhflx = .true.
   endif

!
!     fill in 3d forcing variables if we have both horizontal
!     and vertical components, but not the 3d
!
   if ( .not. have_divq3d .and. have_divq .and. have_vertdivq ) then
      do k=1,plev
         do m=1,pcnst
            divq3d(k,m) = divq(k,m) + vertdivq(k,m)
         enddo
      enddo
      have_divq3d = .true.
   endif

   if ( .not. have_divt3d .and. have_divt .and. have_vertdivt ) then
      print *, 'Don''t have divt3d - using divt and vertdivt'
      do k=1,plev
         divt3d(k) = divt(k) + vertdivt(k)
      enddo
      have_divt3d = .true.
   endif
!
!     make sure that use_3dfrc flag is set to true if we only have
!     3d forcing available
!
   if ( .not. have_divt .or. .not. have_divq ) then
      use_3dfrc = .true.
   endif

   call getncdata( NCID, 1, 1, ioptimeidx, 'CLAT', tmpdata, STATUS )
   if ( STATUS .EQ. NF_NOERR ) then
      clat(1)=tmpdata(1)
      clat_p(1)=tmpdata(1)
      latdeg(1) = clat(1)*45._r8/atan(1._r8)
   endif

   call getncdata( NCID, 1, 1, ioptimeidx,  'beta', tmpdata, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      betacam = 0.
   else
      betacam = tmpdata(1)
   endif
   call getncdata( NCID, 1, 1, ioptimeidx,  'fixmas', tmpdata, STATUS )
   if ( STATUS .NE. NF_NOERR ) then
      fixmascam=1.0
   else
      fixmascam=tmpdata(1)
   endif
   do m=1,pcnst
      lowername=to_lower(alphanam(m))
      call getncdata( NCID,1,1,ioptimeidx,lowername,tmpdata,STATUS )
      alphacam(m)=tmpdata(1)
      if ( STATUS .NE. NF_NOERR ) then
         alphacam(:)=0.
      endif
      lowername=to_lower(dqfxnam(m))
      call getncdata( NCID,1,1,ioptimeidx,lowername,coldata,STATUS )
      if ( STATUS .NE. NF_NOERR ) then
         dqfxcam=0.
      else
         dqfxcam(1,:,m)=coldata(:)
      endif
   enddo

   STATUS = NF_CLOSE( NCID )
   error_code = 0

   return
end subroutine readiopdata

