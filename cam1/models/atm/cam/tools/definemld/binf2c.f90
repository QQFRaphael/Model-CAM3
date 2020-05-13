subroutine binf2c (plonin,  nlatin,  lonin,   latin, mldin, fillvaluein, &
                   plonout, nlatout, nlon, rlon,  latout, mldout, fillvalueout, &
                   verbose)
!-----------------------------------------------------------------------
!
! Purpose: bin from a finer grid to a coarser one, taking account of missing or
! filled data points.  Configured to work correctly on a reduced grid.  Algorithm: 
! For each point on the fine grid with valid data, 1st find the nearest latitude 
! on the coarse mesh.  Then for that latitude find the nearest coarse longitude and 
! put the fine grid point into that coarse "bin".
!
!-----------------------------------------------------------------------
   use precision

   implicit none
!
! Arguments
!
   integer, intent(in) :: plonin                    ! longitude dimension of input
   integer, intent(in) :: nlatin                    ! latitude dimension of input
   real(r8), intent(in) :: lonin(plonin)            ! input longitudes
   real(r8), intent(in) :: latin(nlatin)            ! input latitudes
   real(r8), intent(in) :: mldin(plonin,nlatin)     ! input mixed layer depths
   real(r8), intent(in) :: fillvaluein              ! land flag (input grid)

   integer, intent(in) :: plonout                   ! longitude dimension of output
   integer, intent(in) :: nlatout                   ! latitude dimension of output
   integer, intent(in) :: nlon(nlatout)             ! lons (deg.) at each lat (maybe on reduced grid)
   real(r8), intent(in) :: rlon(plonout,nlatout)    ! longitudes (deg.) at each latitude (maybe on reduced grid)
   real(r8), intent(in) :: latout(nlatout)          ! output latitudes
   real(r8), intent(in) :: fillvalueout             ! land flag (output grid)
   real(r8), intent(out) :: mldout(plonout,nlatout) ! output mixed layer depths

   logical, intent(in) :: verbose                   ! added printout
!
! Local workspace
!  
   integer :: i,j
   integer :: ii,jj
   integer :: iiarr(1), jjarr(1)
   integer :: num
   integer :: bincount(plonout,nlatout)

   real(r8) :: deltay(nlatout)
   real(r8) :: deltax(plonout)

   if (nlatin < nlatout) then
      write(6,*)'Warning: input grid coarser than output grid'
   end if

   bincount(:,:) = 0
   mldout(:,:) = 0.
!
! Find closest output grid lon index (ii) and lat index (jj) to input grid lon (i) and lat (j)
! Then sum mld into appropriate bin
!
   do j=1,nlatin
      deltay(:) = abs (latin(j) - latout(:))
      jjarr = minloc (deltay(:))                  ! lat index of closest output grid point
      jj = jjarr(1)
      do i=1,plonin
         if (mldin(i,j) /= fillvaluein) then
            num = nlon(jj)
            deltax(:num) = abs (lonin(i) - rlon(:num,jj))
            iiarr = minloc (deltax(:num))         ! lon index of closest output grid point
            ii = iiarr(1)
            mldout(ii,jj) = mldout(ii,jj) + mldin(i,j)
            bincount(ii,jj) = bincount(ii,jj) + 1
         end if
      end do
   end do
!
! Normalize by bin count and put fillvalue where count is zero
!
   do jj=1,nlatout
      do ii=1,nlon(jj)
         if (bincount(ii,jj) > 0) then
            mldout(ii,jj) = mldout(ii,jj)/bincount(ii,jj)
         else
            mldout(ii,jj) = fillvalueout
         end if
      end do
   end do

   if (verbose) then
      write(6,*)'bincount:'
      do jj=1,nlatout
         write(6,'(1000i1)') (bincount(ii,jj),ii=1,nlon(jj))
      end do
   end if

   return
end subroutine binf2c
