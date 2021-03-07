!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program LoverD

!  Read (x,y,z) components of force (or force coefficients) and calculate lift/
!  drag forces or coefficients for a given angle of attack.  The body is assumed
!  to be aligned with the xyz axes, which are right-handed with Oy spanwise and
!  positive on the starboard side.  Deal with nonzero sideslip and with moments
!  when the need arises.  A few interactive prompts suffice.  SI units are used.
!
!  02/28/2019  D.A.Saunders  Initial implementation, for the Mars Sample Return
!                            Lander.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   integer, parameter :: lunkbd = 5, luncrt = 6
   real :: Aref, alpha, cosa, sina, CL, CD, Cx, Cz, Fx, Fz, LbyD, Q, rho, V

   write (luncrt, '(/, a)', advance='no') 'Fx & Fz or Cx & Cz: '
   read  (lunkbd, *) Fx, Fz
   write (luncrt, '(a)', advance='no') 'Angle of attack: '
   read  (lunkbd, *) alpha
   write (luncrt, '(a)', advance='no') 'Aref?  [0 => Cx & Cz input] '

   read  (lunkbd, *) Aref
   if (Aref == 0.) then  ! Assume they are force coefs.
      Cx = Fx
      Cz = Fz
   else                  ! Assume forces are input as from DPLR's PostFlow
      write (luncrt, '(a)', advance='no') 'Density & velocity: '
      read  (lunkbd, *) rho, V
      Q  = 0.5*rho*V**2
      Cx = Fx/(Q*Aref)
      Cz = Fz/(Q*Aref)
   end if

   cosa = cosd (alpha);  sina = sind (alpha)
   CL   = -sina*Cx + cosa*Cz
   CD   =  cosa*Cx + sina*Cz
   LbyD = CL/CD

   write (luncrt, '(3x, a, f10.6)') &
      'CL: ', CL, &
      'CD: ', CD, &
      'L/D:', LbyD

   end program LoverD
