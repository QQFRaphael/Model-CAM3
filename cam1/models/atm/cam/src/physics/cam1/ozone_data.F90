#include <misc.h>
#include <params.h>

module ozone_data
!----------------------------------------------------------------------- 
! Purpose:
!
! Interpolation of ozone datasets.
! 
! Revision history:
! 2004-07-31  B. Eaton       Assemble module from comozp.F90, oznini.F90, oznint.F90, radozn.F90
! 2004-08-19  B. Eaton       Modify ozone_data_vert_interp to return mass mixing ratio.
! 2004-08-30  B. Eaton       Add ozone_data_get_cnst method.
!-----------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use pmgrid,         only: plon, plat, masterproc
use ppgrid,         only: begchunk, endchunk, pcols, pver
use phys_grid,      only: scatter_field_to_chunk, get_ncols_p, get_lat_all_p, get_lon_all_p
use physics_types,  only: physics_state
use rgrid,          only: nlon
use commap,         only: londeg, latdeg
use time_manager,   only: get_curr_date, get_perp_date, get_curr_calday, &
                          is_perpetual
use abortutils,     only: endrun
use error_messages, only: handle_err
#if ( defined SPMD )
use mpishorthand
#endif

implicit none
private
save

! Public methods
public ::&
   ozone_data_defaultopts,   &! set default values of namelist variables
   ozone_data_setopts,       &! get namelist input
   ozone_data_init,          &! open dataset and spatially interpolate data bounding initial time
   ozone_data_timestep_init, &! interpolate to current time
   ozone_data_get_cnst,      &! return pointer to o3 interpolated in time and vertically
   ozone_data_final           ! deallocate module data and close dataset

! Public data


! Private data

logical            :: ozncyc = .true.  ! .true. => assume annual cycle ozone data
character(len=256) :: bndtvo = ' '     ! full pathname for time-variant ozone dataset

real(r8) :: cdayozm  ! dataset calendar day previous month
real(r8) :: cdayozp  ! dataset calendar day next month

integer :: nm        ! Array indices for previous month ozone data
integer :: np        ! Array indices for next month ozone data
integer :: ncid_oz   ! netcdf ID for ozone dataset 
integer :: oznid     ! netcdf id for ozone variable
integer :: lonsiz    ! size of longitude dimension on ozone dataset
integer :: levsiz    ! size of level dimension on ozone dataset
integer :: latsiz    ! size of latitude dimension on ozone dataset
integer :: timesiz   ! size of time dimension on ozone dataset
integer :: np1       ! current forward time index of ozone dataset

type ozmixm_pters
   real(r8), dimension(:,:,:), pointer :: val  ! (pcols,levsiz,begchunk:endchunk)
end type ozmixm_pters
type (ozmixm_pters) :: ozmixm(2)          ! mixing ratios for lower and upper bounds
real(r8), allocatable :: ozmix(:,:,:)     ! mixing ratio
                                          ! (pcols,levsiz,begchunk:endchunk)
real(r8), allocatable :: pin(:)           ! ozone pressure level (levsiz)
real(r8), allocatable :: ozlon(:)         ! Longitudes of bdy dataset (lonsiz)
real(r8), allocatable :: ozlat(:)         ! Latitudes of bdy dataset (latsiz)

integer, allocatable :: date_oz(:)        ! Date on ozone dataset (YYYYMMDD)
                                          ! (timesiz)
integer, allocatable :: sec_oz(:)         ! seconds of date (0-86399)
                                          ! (timesiz)
! storage for ozone data interpolated in time and vertically
real(r8), allocatable, target, dimension(:,:,:) :: o3_data   ! mass mixing ratio

! *** N.B. this hardwired mw of dry air needs to be changed to the share value
real(r8), parameter :: mwdry = 28.9644    ! Effective molecular weight of dry air (g/mol)
real(r8)            :: mwo3               ! Molecular weight of ozone (g/mol)

!===============================================================================
contains
!===============================================================================

subroutine ozone_data_defaultopts(ozncyc_out, bndtvo_out)
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
!-----------------------------------------------------------------------

   logical,          intent(out), optional :: ozncyc_out ! annual cycle ozone dataset 
   character(len=*), intent(out), optional :: bndtvo_out ! full pathname for ozone dataset

!-----------------------------------------------------------------------
   if ( present(ozncyc_out) ) then
      ozncyc_out = ozncyc
   endif
   if ( present(bndtvo_out) ) then
      bndtvo_out = bndtvo
   endif
end subroutine ozone_data_defaultopts

!===============================================================================

subroutine ozone_data_setopts(ozncyc_in, bndtvo_in)
!----------------------------------------------------------------------- 
! Purpose: Set runtime options
!-----------------------------------------------------------------------

   logical,          intent(in), optional :: ozncyc_in ! annual cycle ozone dataset 
   character(len=*), intent(in), optional :: bndtvo_in ! full pathname for ozone dataset

!-----------------------------------------------------------------------
   if ( present(ozncyc_in) ) then
      ozncyc = ozncyc_in
   end if
   if ( present(bndtvo_in) ) then
      bndtvo = bndtvo_in
   endif
end subroutine ozone_data_setopts

!===============================================================================

subroutine ozone_data_init(mwo3x)
!----------------------------------------------------------------------- 
! 
! Purpose: Do initial read of time-variant ozone boundary dataset, containing
!          ozone mixing ratios as a function of latitude and pressure.  Read two
!          consecutive months between which the current date lies.  Routine
!          RADOZ2 then evaluates the two path length integrals (with and without
!          pressure weighting) from zero to the interfaces between the input
!          levels.  It also stores the contribution to the integral from each
!          layer.
! 
! Method: Call appropriate netcdf wrapper routines and interpolate to model grid
! 
! Author: CCM Core Group
! Modified: P. Worley, August 2003, for chunking and performance optimization
! 
!-----------------------------------------------------------------------
   use iofilemod,    only: getfil
   use history,      only: addfld, phys_decomp

   real(r8), intent(in) :: mwo3x           ! Molecular weight of ozone
!
! Local workspace
!
   integer dateid                          ! netcdf id for date variable
   integer secid                           ! netcdf id for seconds variable
   integer londimid                        ! netcdf id for longitude dimension
   integer latdimid                        ! netcdf id for latitude dimension
   integer levdimid                        ! netcdf id for level dimension
   integer lonid                           ! netcdf id for longitude variable
   integer latid                           ! netcdf id for latitude variable
   integer levid                           ! netcdf id for level variable
   integer timeid                          ! netcdf id for time variable
   integer dimids(7)                       ! variable shape
   integer cnt4(4)                         ! array of counts for each dimension
   integer strt4(4)                        ! array of starting indices
   integer i, k, lat, n                    ! longitude, level, latitude, time indices
   integer  :: yr, mon, day                ! components of a date
   integer  :: ncdate                      ! current date in integer format [yyyymmdd]
   integer  :: ncsec                       ! current time of day [seconds]
   integer  :: istat
   real(r8) :: calday                      ! current calendar day
   real(r8) caldayloc                      ! calendar day (includes yr if no cycling)
   real(r8), allocatable :: ozmix2D(:,:)   ! temporary ozmix arrays
   real(r8), allocatable :: ozmix3D(:,:,:)
   real(r8), allocatable :: oznbdym(:,:,:) ! ozone data previous time sample
   real(r8), allocatable :: oznbdyp(:,:,:) ! ozone data next time sample
   character(len=256) :: locfn    ! netcdf local filename to open 
!-----------------------------------------------------------------------

   ! Initialize module data
   mwo3 = mwo3x

   call addfld ('O3VMR', 'm3/m3', pver, 'A', 'Ozone volume mixing ratio', phys_decomp, sampling_seq='rad_lwsw')

   nm = 1
   np = 2
!
! SPMD: Master does all the work.  Sends needed info to slaves
!
   if (masterproc) then

! need declaration for locfn

      call getfil (bndtvo, locfn, iflag=0)
      call wrap_open (locfn, 0, ncid_oz)
      write(6,*)'ozone_data_init: successfully opened ',trim(locfn)
!
! Use year information only if not cycling ozone dataset
!
      calday = get_curr_calday()
      if ( is_perpetual() ) then
         call get_perp_date(yr, mon, day, ncsec)
      else
         call get_curr_date(yr, mon, day, ncsec)
      end if
      ncdate = yr*10000 + mon*100 + day
      if (ozncyc) then
         caldayloc = calday
      else
         caldayloc = calday + yr*365.
      end if
!
! Get and check dimension info
!
      CALL WRAP_INQ_DIMID( ncid_oz, 'lon', londimid   )
      CALL WRAP_INQ_DIMID( ncid_oz, 'lev', levdimid   )
      CALL WRAP_INQ_DIMID( ncid_oz, 'time', timeid  )
      CALL WRAP_INQ_DIMID( ncid_oz, 'lat', latdimid   )

      CALL WRAP_INQ_DIMLEN( ncid_oz, londimid, lonsiz   )
      CALL WRAP_INQ_DIMLEN( ncid_oz, levdimid, levsiz   )
      CALL WRAP_INQ_DIMLEN( ncid_oz, latdimid, latsiz   )
      CALL WRAP_INQ_DIMLEN( ncid_oz, timeid, timesiz   )

      CALL WRAP_INQ_VARID( ncid_oz, 'date', dateid   )
      CALL WRAP_INQ_VARID( ncid_oz, 'datesec', secid   )
      CALL WRAP_INQ_VARID( ncid_oz, 'OZONE', oznid   )
      CALL WRAP_INQ_VARID( ncid_oz, 'lon', lonid   )
      CALL WRAP_INQ_VARID( ncid_oz, 'lat', latid   )
      CALL WRAP_INQ_VARID( ncid_oz, 'lev', levid   )

      CALL WRAP_INQ_VARDIMID (ncid_oz, oznid, dimids)
      if (dimids(1) /= londimid .or. dimids(2) /= levdimid .or. dimids(3) /= latdimid) then
         call endrun ('OZONE_DATA_INIT: Data must be ordered lon, lev, lat, time')
      end if
   end if

#if (defined SPMD )
   call mpibcast( lonsiz, 1, mpiint, 0, mpicom )
   call mpibcast( latsiz, 1, mpiint, 0, mpicom )
   call mpibcast( levsiz, 1, mpiint, 0, mpicom )
   call mpibcast( timesiz, 1, mpiint, 0, mpicom )
#endif

!
! Dynamically allocated memory for module data
!
   allocate( date_oz(timesiz), &
             sec_oz(timesiz),  &
             pin(levsiz),      &
             ozmixm(nm)%val(pcols,levsiz,begchunk:endchunk), &
             ozmixm(np)%val(pcols,levsiz,begchunk:endchunk), &
             ozmix(pcols,levsiz,begchunk:endchunk),          &
             o3_data(pcols,pver,begchunk:endchunk), stat=istat )
   call handle_err(istat, 'ozone_data_init: ERROR allocating module data')
!
! Locally dynamic that will be deallocated before "return"
!
   allocate (ozmix3D(plon,levsiz,plat))

   if (masterproc) then
!
! More dynamically allocated memory for module comozp 
! (for masterproc only)
!
      allocate (ozlon(lonsiz))
      allocate (ozlat(latsiz))
!
! More locally dynamic that will be deallocated before "return"
!
      allocate (oznbdym(lonsiz,levsiz,latsiz))
      allocate (oznbdyp(lonsiz,levsiz,latsiz))
      allocate (ozmix2D(levsiz,plat))
!
! Retrieve longitude, latitude and level arrays for interpolation.
!
      CALL WRAP_GET_VAR_REALX (NCID_OZ, lonid,ozlon)
      CALL WRAP_GET_VAR_REALX (NCID_OZ, latid,ozlat)
      CALL WRAP_GET_VAR_REALX (NCID_OZ, levid,pin)
!
! Convert from millibars to pascals
!
      do k=1,levsiz
         pin(k) = pin(k)*100.
      end do

!
! Retrieve entire date and sec variables.
!
      CALL WRAP_GET_VAR_INT (ncid_oz,dateid,date_oz)
      CALL WRAP_GET_VAR_INT (ncid_oz,secid,sec_oz)
      if (ozncyc) then
         if (timesiz < 12) then 
            write(6,*)'OZONE_DATA_INIT: When cycling ozone, dataset must have 12 consecutive ', &
                      'months of data starting with Jan'
            write(6,*)'Current dataset has only ',timesiz,' months'
            call endrun
         end if
         do n = 1,12
            if (mod(date_oz(n),10000)/100 /= n) then
               write(6,*)'OZONE_DATA_INIT: When cycling ozone, dataset must have 12 consecutive ', &
                         'months of data starting with Jan'
               write(6,*)'Month ',n,' of dataset says date=',date_oz(n)
               call endrun
            end if
         end do
      end if

      strt4(1) = 1
      strt4(2) = 1
      strt4(3) = 1
      cnt4(1)  = lonsiz
      cnt4(2)  = levsiz
      cnt4(3)  = latsiz
      cnt4(4)  = 1
!
! Special code for interpolation between December and January
!
      if (ozncyc) then
         n = 12
         np1 = 1
         call bnddyi(date_oz(n  ), sec_oz(n  ), cdayozm)
         call bnddyi(date_oz(np1), sec_oz(np1), cdayozp)
         if (caldayloc <= cdayozp .or. caldayloc > cdayozm) then
            strt4(4) = n
            call wrap_get_vara_realx (ncid_oz,oznid,strt4,cnt4,oznbdym)

            strt4(4) = np1
            call wrap_get_vara_realx (ncid_oz,oznid,strt4,cnt4,oznbdyp)
            goto 10
         end if
      end if
!
! Normal interpolation between consecutive time slices.
!
      do n=1,timesiz-1
         np1 = n + 1
         call bnddyi(date_oz(n  ), sec_oz(n  ), cdayozm)
         call bnddyi(date_oz(np1), sec_oz(np1), cdayozp)
         if (.not.ozncyc) then
            yr = date_oz(n)/10000
            cdayozm = cdayozm + yr*365.
            yr = date_oz(np1)/10000
            cdayozp = cdayozp + yr*365.
         end if
         if (caldayloc > cdayozm .and. caldayloc <= cdayozp) then
            strt4(4) = n
            call wrap_get_vara_realx (ncid_oz,oznid,strt4,cnt4,oznbdym)

            strt4(4) = np1
            call wrap_get_vara_realx (ncid_oz,oznid,strt4,cnt4,oznbdyp)
            goto 10
         end if
      end do
      write(6,*)'OZONE_DATA_INIT: Failed to find dates bracketing ncdate, ncsec=', ncdate, ncsec
      call endrun
10    continue
      write(6,*)'OZONE_DATA_INIT: Read ozone data for dates ',date_oz(n), &
                sec_oz(n),' and ',date_oz(np1),sec_oz(np1)
!
! Spatial interpolation.  If ozone dataset is 2-d (i.e. lonsiz = 1) and 
! thus only latitude interpolation is necessary, expand to 3-d after 
! interpolation.
!
      if (lonsiz == 1) then
         call lininterp (oznbdym ,ozlat   ,levsiz  ,latsiz  ,ozmix2D, &
                         latdeg  ,plat    )
         do lat=1,plat
            do k=1,levsiz
               do i=1,nlon(lat)
                  ozmix3D(i,k,lat) = ozmix2D(k,lat)
               end do
            end do
         end do
 
         call scatter_field_to_chunk(1,levsiz,1,plon,ozmix3D,ozmixm(nm)%val)

         call lininterp (oznbdyp ,ozlat   ,levsiz  ,latsiz  ,ozmix2D, &
                         latdeg  ,plat)
         do lat=1,plat
            do k=1,levsiz
               do i=1,nlon(lat)
                  ozmix3D(i,k,lat) = ozmix2D(k,lat)
               end do
            end do
         end do

        call scatter_field_to_chunk(1,levsiz,1,plon,ozmix3D,ozmixm(np)%val)

      else

         call bilin (oznbdym, ozlon, ozlat, lonsiz, lonsiz, &
                     levsiz, levsiz, latsiz, ozmix3D, londeg, &
                     latdeg, plon, nlon, levsiz, plat)

         call scatter_field_to_chunk(1,levsiz,1,plon,ozmix3D,ozmixm(nm)%val)

         call bilin (oznbdyp, ozlon, ozlat, lonsiz, lonsiz, &
                     levsiz, levsiz, latsiz, ozmix3D, londeg, &
                     latdeg, plon, nlon, levsiz, plat)

        call scatter_field_to_chunk(1,levsiz,1,plon,ozmix3D,ozmixm(np)%val)
      end if
!
! Deallocate dynamic memory for local workspace.  NOT for pointers in common.
!
      deallocate (oznbdym)
      deallocate (oznbdyp)
      deallocate (ozmix2D)
      deallocate (ozmix3D)

#if (defined SPMD )
   else

      call scatter_field_to_chunk(1,levsiz,1,plon,ozmix3D,ozmixm(nm)%val)
      call scatter_field_to_chunk(1,levsiz,1,plon,ozmix3D,ozmixm(np)%val)
!
! Deallocate dynamic memory for local workspace.  NOT for pointers in common.
!
      deallocate (ozmix3D)

#endif
   end if

#if (defined SPMD )
   call mpibcast (np1, 1, mpiint, 0, mpicom)
   call mpibcast (pin, levsiz, mpir8, 0, mpicom)
   call mpibcast (date_oz, timesiz, mpiint, 0, mpicom )
   call mpibcast (sec_oz, timesiz, mpiint, 0, mpicom )
   call mpibcast (cdayozm, 1, mpir8, 0, mpicom)
   call mpibcast (cdayozp, 1, mpir8, 0, mpicom)
#endif

end subroutine ozone_data_init

!===============================================================================

subroutine ozone_data_timestep_init
!----------------------------------------------------------------------- 
! 
! Purpose: Interpolate ozone mixing ratios to current time, reading in new monthly
!          data if necessary, and spatially interpolating it.
! 
! Method: Find next month of ozone data to interpolate.  Linearly interpolate 
!         vertically and horizontally
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
   use timeinterp, only: getfactors
!
! Local workspace
!
   integer cnt4(4)                ! array of counts for each dimension
   integer strt4(4)               ! array of starting indices
   integer i, k, lat              ! longitude/column, level, latitude indices
   integer c, ncol                ! chunk index, number of columns in chunk
   integer ntmp                   ! temporary
   integer :: yr, mon, day        ! components of a date
   integer :: ncdate              ! current date in integer format [yyyymmdd]
   integer :: ncsec               ! current time of day [seconds]

   real(r8) fact1, fact2          ! time interpolation factors
   real(r8) :: calday             ! current calendar day
   real(r8) caldayloc             ! calendar day (includes yr if no cycling)
   real(r8) deltat                ! time (days) between interpolating ozone data

   real(r8), allocatable :: oznbdyp(:,:,:)  ! ozone data from next time sample
   real(r8), allocatable :: ozmix2D(:,:)    ! temporary ozmix arrays
   real(r8), allocatable :: ozmix3D(:,:,:)
!
! Use year information only if a multiyear dataset
!
   calday = get_curr_calday()
   if ( is_perpetual() ) then
      call get_perp_date(yr, mon, day, ncsec)
   else
      call get_curr_date(yr, mon, day, ncsec)
   end if
   ncdate = yr*10000 + mon*100 + day
   if (ozncyc) then
      caldayloc = calday
   else
      caldayloc = calday + yr*365.
   end if
!
! If model time is past current forward ozone timeslice, then
! masterproc reads in the next timeslice for time interpolation.  Messy logic is 
! for ozncyc = .true. interpolation between December and January (np1 == 1).  Note 
! that np1 is never 1 when ozncyc is .false.
!
   if (caldayloc > cdayozp .and. .not. (np1 == 1 .and. caldayloc > cdayozm)) then
!
      if (ozncyc) then
         np1 = mod(np1,12) + 1
      else
         np1 = np1 + 1
      end if
      if (np1 > timesiz) then
         call endrun ('OZONE_DATA_TIMESTEP_INIT: Attempt to read past end of O3 dataset')
      end if
      cdayozm = cdayozp
      call bnddyi(date_oz(np1), sec_oz(np1), cdayozp)
      if (.not.ozncyc) then
         yr = date_oz(np1)/10000
         cdayozp = cdayozp + yr*365.
      end if

      if (.not. (np1 == 1 .or. caldayloc <= cdayozp)) then
         if (masterproc) then
            write(6,*)'OZONE_DATA_TIMESTEP_INIT: Input ozone for date',date_oz(np1),' sec ',sec_oz(np1), &
                      'does not exceed model date',ncdate,' sec ',ncsec,' Stopping.'
         endif
         call endrun ()
      end if

      ntmp = nm
      nm = np
      np = ntmp
!
! Allocate memory for dynamic local workspace
!
      allocate (ozmix3D(plon,levsiz,plat))

      if (masterproc) then
         strt4(1) = 1
         strt4(2) = 1
         strt4(3) = 1
         strt4(4) = np1
         cnt4(1)  = lonsiz
         cnt4(2)  = levsiz
         cnt4(3)  = latsiz
         cnt4(4)  = 1
!
! Allocate memory for more dynamic local workspace
!
         allocate (oznbdyp(lonsiz,levsiz,latsiz))
         allocate (ozmix2D(levsiz,plat))

         call wrap_get_vara_realx (ncid_oz,oznid,strt4,cnt4,oznbdyp)
         write(6,*)'OZONE_DATA_TIMESTEP_INIT: Read ozone for date (yyyymmdd) ', date_oz(np1),' sec ',sec_oz(np1)
!
! Spatial interpolation.  If ozone dataset is only 2-d (i.e. lonsiz = 1) and 
! thus only latitude interpolation is necessary, expand to 3-d after 
! interpolation.
!
         if (lonsiz == 1) then
            call lininterp (oznbdyp, ozlat, levsiz, latsiz, ozmix2D, &
                            latdeg, plat)
            do lat=1,plat
               do k=1,levsiz
                  do i=1,nlon(lat)
                     ozmix3D(i,k,lat) = ozmix2D(k,lat)
                  end do
               end do
            end do
         else
            call bilin (oznbdyp ,ozlon   ,ozlat   ,lonsiz  ,lonsiz  , &
                        levsiz  ,levsiz  ,latsiz  ,ozmix3D , londeg  , &
                        latdeg  ,plon   ,nlon    ,levsiz  ,plat    )
         end if
      end if

      call scatter_field_to_chunk(1,levsiz,1,plon,ozmix3D,ozmixm(np)%val)
!
! Deallocate dynamic memory for local workspace.
!
      deallocate (ozmix3D)
      if (masterproc) then
         deallocate (oznbdyp)
         deallocate (ozmix2D)
      end if
   end if
!
! Determine time interpolation factors.
!
   call getfactors (ozncyc, np1, cdayozm, cdayozp, caldayloc, &
                    fact1, fact2, 'OZONE_DATA_TIMESTEP_INIT:')
!
! Time interpolation.
!
   do c=begchunk,endchunk
      ncol = get_ncols_p(c)
      do k=1,levsiz
         do i=1,ncol
            ozmix(i,k,c) = ozmixm(nm)%val(i,k,c)*fact1 + ozmixm(np)%val(i,k,c)*fact2
         end do
      end do
   end do

end subroutine ozone_data_timestep_init

!================================================================================================

subroutine ozone_data_get_cnst(state, q)
!-------------------------------------------------------------------------------
! 
! Purpose: 
! Return pointer to constituent concentrations.
!
!-------------------------------------------------------------------------------
   type(physics_state),               intent(in) :: state
   real(r8), pointer, dimension(:,:)             :: q     ! constituent mass mixing ratio

   ! local variables
   integer :: lchnk            ! chunk identifier
!-------------------------------------------------------------------------------

   lchnk = state%lchnk
   call ozone_data_vert_interp(lchnk, state%ncol, state%pmid, o3_data(:,:,lchnk))
   q => o3_data(:,:,lchnk)

end subroutine ozone_data_get_cnst

!================================================================================================

subroutine ozone_data_vert_interp(lchnk, ncol, pmid, o3)
!----------------------------------------------------------------------- 
! 
! Purpose: Interpolate ozone from current time-interpolated values to model levels
! 
! Method: Use pressure values to determine interpolation levels
! 
! Author: Bruce Briegleb
! 
!--------------------------------------------------------------------------
   use history, only: outfld

! Arguments
!
   integer, intent(in) :: lchnk               ! chunk identifier
   integer, intent(in) :: ncol                ! number of atmospheric columns

   real(r8), intent(in) :: pmid(pcols,pver)   ! level pressures (mks)

   real(r8), intent(out) :: o3(pcols,pver)    ! ozone mass mixing ratio
!
! local storage
!
   integer i                   ! longitude index
   integer k, kk, kkstart      ! level indices
   integer kupper(pcols)       ! Level indices for interpolation
   integer kount               ! Counter
   integer lats(pcols)         ! latitude indices
   integer lons(pcols)         ! latitude indices

   real(r8) dpu                ! upper level pressure difference
   real(r8) dpl                ! lower level pressure difference
   real(r8) mwr                ! molecular weight ratio
!--------------------------------------------------------------------------
!
! Initialize latitude indices
!
   call get_lat_all_p(lchnk, ncol, lats)
   call get_lon_all_p(lchnk, ncol, lons)
!
! Initialize index array
!
   do i=1,ncol
      kupper(i) = 1
   end do

   do k=1,pver
!
! Top level we need to start looking is the top level for the previous k
! for all longitude points
!
      kkstart = levsiz
      do i=1,ncol
         kkstart = min0(kkstart,kupper(i))
      end do
      kount = 0
!
! Store level indices for interpolation
!
      do kk=kkstart,levsiz-1
         do i=1,ncol
            if (pin(kk).lt.pmid(i,k) .and. pmid(i,k).le.pin(kk+1)) then
               kupper(i) = kk
               kount = kount + 1
            end if
         end do
!
! If all indices for this level have been found, do the interpolation and
! go to the next level
!
         if (kount.eq.ncol) then
            do i=1,ncol
               dpu = pmid(i,k) - pin(kupper(i))
               dpl = pin(kupper(i)+1) - pmid(i,k)
               o3(i,k) = (ozmix(i,kupper(i)  ,lchnk)*dpl + &
                             ozmix(i,kupper(i)+1,lchnk)*dpu)/(dpl + dpu)
            end do
            goto 35
         end if
      end do
!
! If we've fallen through the kk=1,levsiz-1 loop, we cannot interpolate and
! must extrapolate from the bottom or top ozone data level for at least some
! of the longitude points.
!
      do i=1,ncol
         if (pmid(i,k) .lt. pin(1)) then
            o3(i,k) = ozmix(i,1,lchnk)*pmid(i,k)/pin(1)
         else if (pmid(i,k) .gt. pin(levsiz)) then
            o3(i,k) = ozmix(i,levsiz,lchnk)
         else
            dpu = pmid(i,k) - pin(kupper(i))
            dpl = pin(kupper(i)+1) - pmid(i,k)
            o3(i,k) = (ozmix(i,kupper(i)  ,lchnk)*dpl + &
                          ozmix(i,kupper(i)+1,lchnk)*dpu)/(dpl + dpu)
         end if
      end do

      if (kount.gt.ncol) then
         call endrun ('ozone_data_vert_interp: Bad ozone data: non-monotonicity suspected')
      end if
35    continue
   end do

   call outfld('O3VMR', o3, pcols, lchnk)

   ! convert from the dataset values of vmr to mmr
   mwr = mwo3/mwdry
   do k=1,pver
      do i=1,ncol
         o3(i,k) = mwr*o3(i,k)
      end do
   end do

end subroutine ozone_data_vert_interp

!===============================================================================

subroutine ozone_data_final()

! Deallocate dynamically allocated memory for module comozp and close dataset.

   deallocate (date_oz)
   deallocate (sec_oz)
   deallocate (pin)
   deallocate (ozmixm(1)%val)
   deallocate (ozmixm(2)%val)
   deallocate (ozmix)
   deallocate (o3_data)
   if (masterproc) then
      deallocate (ozlon)
      deallocate (ozlat)
      call wrap_close (ncid_oz)
   end if

end subroutine ozone_data_final

!===============================================================================

end module ozone_data
