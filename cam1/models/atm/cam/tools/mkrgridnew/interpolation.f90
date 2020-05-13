module interpolation
!
! Routines to interpolate between full and reduced grids
!
   use prec,    only: r8
   use globals, only: wnummax, nlon

PRIVATE

   public :: lininterp, cubinterp, fourier

CONTAINS

   subroutine lininterp (plon, nlat, numlev, reverse, arrin, &
                         arrot, fillvaluein, fillvalueot)
!
! Purpose: linear interpolator (in one dimension only)
!
! Method:  Assume latitudes match, as do 1st longitude points.  Assume equal
!          spacing in longitude.  Use "nlon" values to determine reduced
!          grid, then interpolate in longitude.
!
      implicit none
!
! Arguments
!
      integer, intent(in) :: plon     ! longitude dimension
      integer, intent(in) :: nlat     ! latitude dimension
      integer, intent(in) :: numlev   ! size of 3rd dimension (normally vertical levels)

      logical, intent(in) :: reverse  ! true => convert reduced->full

      real(r8), intent(in) :: arrin(0:plon-1,nlat,numlev)  ! input array
      real(r8), intent(out) :: arrot(0:plon-1,nlat,numlev) ! output array
!
! If fillvaluein is non-zero, this means it was found as an attribute on the
! input field.  If non-zero, check for it when interpolating.
!
      real(r8), intent(in) :: fillvaluein   ! _FillValue to be checked for on input
                                            !   Ignore if zero.
      real(r8), intent(in) :: fillvalueot   ! _FillValue to use on output
!
! Local workspace
!
      integer :: i, ii, j, k  ! spatial indices
      integer :: iiright      ! input right-hand interp index
      integer :: nlonin       ! input number of longitudes
      integer :: nlonot       ! output number of longitudes

      real(r8) :: dxin        ! delta-x input
      real(r8) :: dxot        ! delta-x output
      real(r8) :: xlocot      ! x-location on output grid
      real(r8) :: left        ! left-hand input interp point
      real(r8) :: right       ! right-hand input interp point
      real(r8) :: lfact       ! left-hand interp factor
      real(r8) :: rfact       ! right-hand interp factor
!
! Do the interpolation on a latitude by latitude basis
!
      do j=1,nlat
         if (reverse) then
            nlonin = nlon(j)
            nlonot = plon
         else
            nlonin = plon
            nlonot = nlon(j)
         end if

         if (nlonin > plon .or. nlonot > plon) then
            write(6,*)'lininterp: number of longitudes cannot exceed dimension size'
            stop 999
         end if
         
         if (nlonin < 2 .or. nlonot < 1) then
            write(6,*)'lininterp: must have at least 2 points for input, 1 for output'
            stop 999
         end if

         if (fillvaluein /= 0. .and. fillvaluein /= fillvalueot) then
            write(6,*)'lininterp: fillvaluein must = fillvalueot if non-zero'
            stop 999
         end if
!
! 1st longitude point: *assume* overlapping
!
         arrot(0,j,:) = arrin(0,j,:)
         dxin = 1./nlonin
         dxot = 1./nlonot
!
! Loop over longitude of output grid, finding interpolation points
!
         do i=1,nlonot-1
            xlocot = i*dxot
            right = 0.
            do ii=1,nlonin
               left = right
               right = ii*dxin
               if (right >= xlocot) then
                  lfact = (right-xlocot)/(right-left)
                  rfact = (xlocot-left)/(right-left)
!
! Wrap rhs index of input grid if interpolating reduced -> full, and we are 
! at the rhs of the input grid
!
                  iiright = ii
                  if (iiright == nlonin) iiright = 0
!
! If either interp point is fillvaluein, put fillvalueot in the output field.
! This can happen in ISCCP fields, for example.
!
                  if (fillvaluein /= 0.) then
                     do k=1,numlev
                        if (arrin(ii-1,j,k) == fillvaluein .or. &
                            arrin(iiright,j,k) == fillvaluein) then
                           arrot(i,j,k) = fillvalueot
                        else
                           arrot(i,j,k) = arrin(ii-1,j,k)*lfact + arrin(iiright,j,k)*rfact
                        end if
                     end do
                  else
                     do k=1,numlev
                        arrot(i,j,k) = arrin(ii-1,j,k)*lfact + arrin(iiright,j,k)*rfact
                     end do
                  end if
                  goto 10
               end if
            end do
            write(6,*)'Interpolation failed for i=',i
            stop 99
10          continue
         end do
!
! If interpolating to reduced grid, Fill remainder with fillvalue
!
         arrot(nlonot:,j,:) = fillvalueot
      end do

      return
   end subroutine lininterp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine cubinterp (plon, nlat, numlev, reverse, mono, &
                         arrin, arrot, fillvaluein, fillvalueot)
!
! Purpose: cubic interpolator (in one dimension only)
!
! Method:  Assume latitudes match, as do 1st longitude points.  Assume equal
!          spacing in longitude.  Use "nlon" values to determine reduced
!          grid, then interpolate in longitude.  The method follows that
!          defined in a paper by Williamson.
!
      implicit none
!
! Parameters
!
      real(r8), parameter :: epsmac = 1.e-13 ! epsilon guard
!
! Arguments
!
      integer, intent(in) :: plon     ! longitude dimension
      integer, intent(in) :: nlat     ! latitude dimension
      integer, intent(in) :: numlev   ! size of 3rd dimension

      logical, intent(in) :: reverse  ! true => convert reduced->full
      logical, intent(in) :: mono     ! whether or not do impose monotonicity

      real(r8), intent(in) :: arrin(0:plon-1,nlat,numlev)  ! input array
      real(r8), intent(out) :: arrot(0:plon-1,nlat,numlev) ! output array
!
! If fillvaluein is non-zero, this means it was found as an attribute on the
! input field.  If non-zero, check for it when interpolating.
!
      real(r8), intent(in) :: fillvaluein   ! _FillValue to be checked for on input
                                            !   Ignore if zero.
      real(r8), intent(in) :: fillvalueot   ! _FillValue to use on output
!
! Local workspace
!
      integer :: i, ii, j, k     ! spatial indices
      integer :: nlonin          ! input number of longitudes
      integer :: nlonot          ! output number of longitudes

      real(r8) :: dxin           ! delta-x on input grid
      real(r8) :: dxot           ! delta-x on output grid
      real(r8) :: xlocot         ! position on output grid (i * dxot)
      real(r8) :: left           ! left position on input grid (ii * dxin)
      real(r8) :: right          ! right position on input grid ((ii+1) * dxin)
      real(r8) :: beta           ! (Xi - Xl) / (Xr - Xl)
      real(r8) :: onemb          ! 1. - beta
      real(r8) :: onepb          ! 1. + beta
      real(r8) :: twomb          ! 2. - beta
      real(r8) :: delq(numlev)   ! delta-q (fm Williamson paper)
      real(r8) :: test(numlev)   ! compared with 0 and 3 (fm Williamson paper)
      real(r8) :: dl(numlev)     ! left derivative estimate (fm Williamson paper)
      real(r8) :: dr(numlev)     ! right derivative estimate (fm Williamson paper)
      real(r8) :: tmp(numlev)    ! temporary workspace

      logical :: found           ! whether fillvaluein found

      real(r8) :: arrin_wr(-1:plon+1,numlev) ! Input array with wrap points included
!
! Do the interpolation on a latitude by latitude basis
!
      do j=1,nlat
         if (reverse) then
            nlonin = nlon(j)
            nlonot = plon
         else
            nlonin = plon
            nlonot = nlon(j)
         end if

         if (nlonin > plon .or. nlonot > plon) then
            write(6,*)'cubinterp: number of longitudes cannot exceed dimension'
            stop 999
         end if

         if (nlonin < 2 .or. nlonot < 1) then
            write(6,*)'cubinterp: must have at least 2 points for input, 1 for output'
            stop 999
         end if

         if (fillvaluein /= 0. .and. fillvaluein /= fillvalueot) then
            write(6,*)'cubinterp: fillvaluein must = fillvalueot if non-zero'
            stop 999
         end if
!
! Build local copy of input array with wrap points for this latitude
!
         arrin_wr(-1,:)         = arrin(nlonin-1,j,:)
         arrin_wr(0:nlonin-1,:) = arrin(:,j,:)
         arrin_wr(nlonin,:)     = arrin(0,j,:)
         arrin_wr(nlonin+1,:)   = arrin(1,j,:)
!
! Put fillvalueot in input array positions that should NEVER be accessed.
! This should show up as unmistakable crap in the output field.
!
         if (nlonin < plon-1) then
            arrin_wr(nlonin+2:,:) = fillvalueot
         end if
!
! Check for fillvalue in input field (e.g. cloudsimulator fields).
! If found, put fillvalue into entire lon x lev section of output
! field and don't interpolate.
!
         found = .false.
         if (fillvaluein /= 0.) then
            do k=1,numlev
               do i=1,nlonin
                  if (arrin(i,j,k) == fillvaluein) then
                     found = .true.
                  end if
               end do
            end do
         end if

         if (found) then
            arrot(:,j,:) = fillvalueot
            cycle
         end if
!
! 1st longitude point: *assume* overlapping
!
         arrot(0,j,:) = arrin(0,j,:)
         dxin = 1.d0/nlonin
         dxot = 1.d0/nlonot
!
! Loop over longitude of output grid, finding interpolation points
!
         do i=1,nlonot-1
            xlocot = i*dxot
            right = 0.
            do ii=1,nlonin
               left = right
               right = ii*dxin
               if (right >= xlocot) then
                  beta = (xlocot-left)/(right-left)
                  onepb = 1. + beta
                  onemb = 1. - beta
                  twomb = 2. - beta

                  if (mono) then
                     dl(:) = -arrin_wr(ii-2,:)/(3.*dxin) - &
                              arrin_wr(ii-1,:)/(2.*dxin) + &
                              arrin_wr(ii  ,:)/dxin -      &
                              arrin_wr(ii+1,:)/(6.*dxin)

                     dr(:) =  arrin_wr(ii-2,:)/(6.*dxin) - &
                              arrin_wr(ii-1,:)/dxin +      &
                              arrin_wr(ii  ,:)/(2.*dxin) + &
                              arrin_wr(ii+1,:)/(3.*dxin)

                     delq(:) = (arrin_wr(ii,:) - arrin_wr(ii-1,:))/dxin
                     test(:) = abs(3.*delq(:))*(1. - epsmac)

                     where (dl(:)*delq(:) < 0.)
                        dl(:) = 0.
                     elsewhere
                        tmp(:) = abs(dl(:))
                        dl(:) = sign (min (tmp(:), test(:)), delq(:))
                     end where

                     where (dr(:)*delq(:) < 0.)
                        dr(:) = 0.
                     elsewhere
                        tmp(:) = abs(dr(:))
                        dr(:) = sign (min (tmp(:), test(:)), delq(:))
                     end where

                     arrot(i,j,:) = (3. - 2.*onemb)*onemb**2*arrin_wr(ii-1,:) + &
                          dxin*beta*onemb**2*dl(:) +                  &
                          (3. - 2.*beta)*beta**2*arrin_wr(ii,:) -     &
                          dxin*beta**2*onemb*dr(:)            
                  else

                     arrot(i,j,:) = -beta/6.*onemb*twomb*arrin_wr(ii-2,:) +     &
                          0.5*onepb*onemb*twomb*arrin_wr(ii-1,:) +    &
                          0.5*onepb*beta*twomb*arrin_wr(ii,:) -       &
                          onepb/6.*beta*onemb*arrin_wr(ii+1,:)
                  end if

                  goto 10

               end if
            end do
            write(6,*)'cubinterp: interpolation failed for i=',i
            stop 99
10          continue
         end do
!
! If interpolating to reduced grid, Fill remainder with fillvalue
!
         arrot(nlonot:,j,:) = fillvalueot
      end do

      return
   end subroutine cubinterp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine fourier (plon, nlat, numlev, reverse, arrin, &
                       arrot, fname, fillvaluein, fillvalueot)
!
! Purpose: fourier interpolator (in one dimension only)
!
! Method:  Assume latitudes match, as do 1st longitude points.  Assume equal
!          spacing in longitude.  Use "nlon" values to determine reduced
!          grid and interpolate in longitude.  Use "wnummax" to determine
!          max fourier wavenumber to keep.
!
      use fortfft

      implicit none
!
! Arguments
!
      integer, intent(in) :: plon     ! longitude dimension
      integer, intent(in) :: nlat     ! latitude dimension
      integer, intent(in) :: numlev   ! size of 3rd dimension
      logical, intent(in) :: reverse  ! true => convert reduced->full
      real(r8), intent(in) :: arrin(plon,nlat,numlev)  ! input array
      real(r8), intent(out) :: arrot(plon,nlat,numlev) ! output array
      character(len=*), intent(in) :: fname  ! field name (do something special for PS)
!
! If fillvaluein is non-zero, this means it was found as an attribute on the
! input field.  If non-zero, check for it when interpolating.
!
      real(r8), intent(in) :: fillvaluein   ! _FillValue to be checked for on input
                                            !   Ignore if zero.
      real(r8), intent(in) :: fillvalueot   ! _FillValue to use on output
!
! Local workspace
!
      integer :: i, j, k                ! spatial indices
      integer :: ifaxin(19)             ! needed by set99
      integer :: ifaxot(19)             ! needed by set99
      integer :: throwindx              ! cutoff fourier index 
      integer :: nlonin                 ! length of longitude index input grid
      integer :: nlonot                 ! length of longitude index output grid
!
! Need to use a temp array for transform due to length requirements of fft
!
      real(r8) :: x(plon+2,numlev)      ! workspace for transform
      real(r8) :: trigin(3*plon/2+1)    ! workspace for set99
      real(r8) :: trigot(3*plon/2+1)    ! workspace for set99
      real(r8) :: work((plon+1)*numlev) ! workspace for fft991

      logical :: found                  ! whether fillvaluein found
!
! Do the interpolation on a latitude by latitude basis
!
      do j=1,nlat
!
! Call fft setup routines for input and output grids
!
         if (reverse) then
            nlonin = nlon(j)
            nlonot = plon
         else
            nlonin = plon
            nlonot = nlon(j)
         end if

         call set99 (trigin, ifaxin, nlonin)
         call set99 (trigot, ifaxot, nlonot)

         if (nlonin > plon .or. nlonot > plon) then
            write(6,*)'fourier: number of longitudes cannot exceed dimension'
            stop 999
         end if

         if (nlonin < 2 .or. nlonot < 1) then
            write(6,*)'fourier: must have at least 2 points for input, 1 for output'
            stop 999
         end if

         if (fillvaluein /= 0. .and. fillvaluein /= fillvalueot) then
            write(6,*)'fourier: fillvaluein must = fillvalueot if non-zero'
            stop 999
         end if
!
! Check for fillvaluein in input field (e.g. cloudsimulator fields).
! If found, put fillvalueot into entire lon x lev section of output
! field and don't do transform.
!
         found = .false.
         if (fillvaluein /= 0.) then
            do k=1,numlev
               do i=1,nlonin
                  if (arrin(i,j,k) == fillvaluein) then
                     found = .true.
                  end if
               end do
            end do
         end if

         if (found) then
            arrot(:,j,:) = fillvalueot
            cycle
         end if
            
         if (fname == 'PS') then
            x(:nlonin,:) = log(arrin(:nlonin,j,:))
         else
            x(:nlonin,:) = arrin(:nlonin,j,:)
         end if
         x(nlonin+1:,:) = 0.

         call fft991 (x, work, trigin, ifaxin, 1, plon+2, nlonin, numlev, -1)
!
! Zero out waves which will not be included
! The 3 in throwindx counts 2 for the mean plus 1 for fortran indexing 
! starting at 1
!
         throwindx = 2*wnummax(j) + 3

         if (throwindx < nlonot) then
            x(throwindx:,:) = 0.
         end if

         call fft991 (x, work, trigot, ifaxot, 1, plon+2, nlonot, numlev, +1)
!
! If interpolating to reduced grid, Fill remainder with fillvalueot
!
         if (fname == 'PS') then
            arrot(:nlonot,j,:) = exp(x(:nlonot,:))
         else
            arrot(:nlonot,j,:) = x(:nlonot,:)
         end if
         arrot(nlonot+1:,j,:) = fillvalueot
      end do

      return
   end subroutine fourier
end module interpolation
