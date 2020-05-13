module commap

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plon, plat

   real(r8) w(plat)                  ! integration weights for physics grid
   real(r8) w_staggered(plat-1)      ! integration weights for the staggered wind arrays
!
   real(r8) clat(plat)             ! model latitudes (radians)
   real(r8) clat_staggered (plat-1)  ! model latitudes on staggered grid (radians)
   real(r8) clon(plon,plat)          ! model longitudes (radians)
   real(r8) latdeg(plat)           ! model latitudes (degrees)
   real(r8) londeg(plon,plat)        ! model longitudes (degrees)
end module commap
