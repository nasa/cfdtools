!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine earth_atmosphere (alt, rho, T, p, a)

!  For the given altitude (km), return the corresponding density, temperature,
!  pressure, and speed of sound, in SI units.  The US 1976 Standard Atmosphere
!  model is (evidently) employed, but no reference is provided in the original
!  source code from Grant Palmer (NASA ARC).  Seven atmosphere layers are
!  treated, with extrapolation allowed above the highest layer.
!
!  This is a minimally-modified form of the original main program, which read
!  a list of altitudes and printed altitude, density, temperature & pressure.
!
!  History:
!
!     ??/??/????  G.E.Palmer    Initial implementation as a main program.
!     01/27/2017  D.A.Saunders  Adaptation as a reusable subroutine, handling
!                               one altitude per call.
!
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real, intent (in)  :: alt       ! Geometric altitude, km
   real, intent (out) :: rho, &    ! Corresponding density, kg/m^3
                         T,   &    !   "   "   "   temperature, K
                         p,   &    !   "   "   "   pressure, Pa
                         a         !   "   "   "   speed of sound, m/s
!  Local constants:

   integer, parameter :: ntab = 8      ! # entries in the defining tables
   real, parameter :: gamma = 1.4      ! specific heat ratio
   real, parameter :: gmr = 34.163195  ! hydrostatic constant
   real, parameter :: rearth = 6369.0  ! radius of the Earth (km)
   real, parameter :: Rs = 287.058     ! specific gas constant, J/kg/K

!  Local variables:

   integer :: i, j, k

   real    :: delta         ! pressure/sea-level standard pressure
   real    :: deltah        ! height above base of this layer
   real    :: h             ! geopotential altitude (km)
   real    :: sigma         ! density/sea-level standard density
   real    :: tgrad, tbase  ! temperature gradient and base temp of this layer
   real    :: theta         ! temperature/sea-level standard temperature
   real    :: tlocal        ! local temperature

!  Data:

   real, dimension (ntab), parameter :: &
      htab = (/0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852/), &
      gtab = (/-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0/), &
      ttab = (/288.15, 216.65, 216.65, 228.65, 270.65,270.65,214.65,186.946/), &
      ptab = (/1.0, 2.233611E-1, 5.403295E-2, 8.5666784E-3, 1.0945601E-3, &
              6.6063531E-4, 3.9046834E-5, 3.68501E-6/)
!  Execution:

   h = alt*rearth/(alt+rearth)     ! convert geometric to geopotential altitude

   i = 1                           ! setting up for binary search
   j = ntab
   do  ! Until one interval
      k = (i+j)/2                  ! integer division
      if (h < htab(k)) then
         j = k
      else
         i = k
      end if
      if (j <= i+1) exit
   end do

   tgrad  = gtab(i)                ! i will be in [1, ntab-1]
   tbase  = ttab(i)
   deltah = h - htab(i)
   tlocal = tbase + tgrad*deltah
   theta  = tlocal/ttab(1)         ! temperature ratio

   if (tgrad == 0.0) then
       delta = ptab(i)*exp(-gmr*deltah/tbase)  ! pressure ratio
   else
       delta = ptab(i)*(tbase/tlocal)**(gmr/tgrad)
   end if

   sigma = delta/theta
   rho   = sigma*1.225             ! 1.225 kg/m^3 = sea level/15 deg C density
   T     = theta*288.15            ! 15 deg C = 288.15 K
   p     = delta*101325.           ! 101325 Pa = 1 atmosphere
   a     = sqrt (gamma*Rs*T)       ! speed of sound

   end subroutine earth_atmosphere
