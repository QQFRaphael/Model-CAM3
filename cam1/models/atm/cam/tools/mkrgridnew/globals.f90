module globals
!
! data module containing global data for convenience of minimizing size
! of argument lists.
!
   use prec, only: r8
!
! def_fillvalue is the value to use for _FillValue if it does not already
! exist for the input variable (it might already exist for e.g. ISCCP fields).  
! This value will be placed in output fields outside the bounds of the reduced 
! grid (since the full rectangular grid is used to store the data).
!
   real(r8), parameter :: def_fillvalue  = 1.d36
!
! nlon and wnummax are namelist variables.  Thus they cannot be dynamically
! allocated.  So have to live with a "max" length
!
   integer, parameter  :: maxnlat = 10000  ! max number of latitudes
   integer :: nlon(maxnlat) = -1           ! number of longitudes
   integer :: wnummax(maxnlat) = -1        ! max fourier wave number to keep
!
! netcdf file and dimension ids: init to bad value
!
   integer :: ncidi = -1       ! input file
   integer :: ncido = -1       ! output file
   integer :: londimid = -1    ! longitude
   integer :: latdimid = -1    ! latitude
   integer :: levdimid = -1    ! vertical level
   integer :: ilevdimid = -1   ! vertical interfaces
   integer :: unlimdimid = -1  ! unlimited (time) dimension
end module globals
