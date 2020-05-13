subroutine neighborfill (plonout, nlatout, nlon, landfrac, mldout, &
                         fillpts, fillvalueout)
!----------------------------------------------------------------------------------
!
! Purpose: Perform nearest neighbor fill for those points which need to have a value
! defined (i.e. LANDFRAC < 1), but for which previous binning procedure did not produce an
! entry.  Coded as an iterative procedure which only looks north, east, south, and west
! a distance of one grid point.  Iteration is complete when when all points have data
! where LANDFRAC says non-land exists.
!
!----------------------------------------------------------------------------------
   use precision

   implicit none

   include 'netcdf.inc'
!
! Arguments
!
   integer, intent(in) :: plonout                     ! longitude dimension of output
   integer, intent(in) :: nlatout                     ! latitude dimension of output
   integer, intent(in) :: nlon(nlatout)               ! lons (deg.) at each lat (maybe on reduced grid)
   real(r8), intent(in) :: fillvalueout               ! land flag
   real(r8), intent(in) :: landfrac(plonout,nlatout)  ! land fraction
   real(r8), intent(out) :: fillpts(plonout,nlatout)  ! flag indicates how filling was done
   real(r8), intent(inout) :: mldout(plonout,nlatout) ! output mixed layer depths
!
! Local workspace
!  
   integer :: i,j        ! lon, lat indices
   integer :: im, ip     ! previous, next lon indices
   integer :: iter       ! iteration counter
   integer :: nfill      ! number of points which got filled each iteration
   integer :: npts       ! number of non-land points surrounding a land point

   real(r8) :: avg       ! average of nearest neighbors
!
! Fill with nearest neighbor until all points dictated by LANDFRAC field are accounted for
!
   fillpts(:,:) = 0
   iter = 0
10 continue
   iter = iter + 1
   nfill = 0
   do j=1,nlatout
      do i=1,nlon(j)
         if (mldout(i,j) == fillvalueout .and. nint (landfrac(i,j)) < 1.) then
            ip = mod (i,nlon(j)) + 1
            im = i - 1
            if (im == 0) im = nlon(j)

            avg = 0.
            npts = 0

            if (mldout(im,j) /= fillvalueout) then
               avg = avg + mldout(im,j)
               npts = npts + 1
            end if

            if (mldout(ip,j) /= fillvalueout) then
               avg = avg + mldout(ip,j)
               npts = npts + 1
            end if

            if (j > 1) then
               if (mldout(i,j-1) /= fillvalueout) then
                  avg = avg + mldout(i,j-1)
                  npts = npts + 1
               end if
            end if

            if (j < nlatout) then
               if (mldout(i,j+1) /= fillvalueout) then
                  avg = avg + mldout(i,j+1)
                  npts = npts + 1
               end if
            end if

            if (npts > 0) then
               mldout(i,j) = avg/npts
               nfill = nfill + 1
            end if

         end if
      end do
   end do

   write(6,*)'NEIGHBORFILL: iter ', iter,' filled ', nfill, ' points'
   if (nfill > 0) goto 10

   return
end subroutine neighborfill
