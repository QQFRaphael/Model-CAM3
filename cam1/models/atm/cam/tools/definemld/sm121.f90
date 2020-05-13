subroutine sm121 (plonout, nlatout, nlon, mldout, fillvalueout)
   use precision

   implicit none
!
! Arguments
!
   integer, intent(in) :: plonout                   ! longitude dimension of output
   integer, intent(in) :: nlatout                   ! latitude dimension of output
   integer, intent(in) :: nlon(nlatout)             ! lons (deg.) at each lat (maybe on reduced grid)
   real(r8), intent(in) :: fillvalueout             ! land flag
   real(r8), intent(inout) :: mldout(plonout,nlatout) ! output mixed layer depths
!
! Local workspace
!
   real(r8) :: tmp(plonout,nlatout)
   real(r8) :: sum1
   real(r8) :: sum2
   integer :: i, j
   integer :: im, ip

   do j=1,nlatout
      if (nlon(j) /= plonout) then
         write(6,*) 'sm121 not yet ready for reduced grid'
         stop 999
      end if

      do i=1,plonout
         if (mldout(i,j) == fillvalueout) then
            tmp(i,j) = fillvalueout
         else
            sum1 = 4.*mldout(i,j)
            sum2 = 4.
            ip = mod (i,plonout) + 1
            im = i - 1
            if (im == 0) im = plonout
            if (mldout(im,j) /= fillvalueout) then
               sum1 = sum1 + mldout(im,j)
               sum2 = sum2 + 1.
            end if
            if (mldout(ip,j) /= fillvalueout) then
               sum1 = sum1 + mldout(ip,j)
               sum2 = sum2 + 1.
            end if
            if (j < nlatout) then 
               if (mldout(i,j+1) /= fillvalueout) then
                  sum1 = sum1 + mldout(i,j+1)
                  sum2 = sum2 + 1.
               end if
            end if
            if (j > 1) then
               if (mldout(i,j-1) /= fillvalueout) then
                  sum1 = sum1 + mldout(i,j-1)
                  sum2 = sum2 + 1.
               end if
            end if
            tmp(i,j) = sum1/sum2
         end if
      end do
   end do

   mldout(:,:) = tmp(:,:)

   return
end subroutine sm121
