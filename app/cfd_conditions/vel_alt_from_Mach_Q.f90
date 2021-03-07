!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine vel_alt_from_Mach_Q (M, Q, V, h, rho, T, p, a, ier)
!
!  Description:
!
!     For the given (Mach, dynamic pressure) pair, (M, Q), return the equivalent
!     (velocity, altitude) pair (V, h) along with corresponding density, temp-
!     erature, pressure, and speed of sound.  The atmosphere is assumed to be
!     air (U.S. 1976 Standard Atmosphere model).
!
!  Method (see program CFD_Conditions):
!
!     While M and Q can be derived from V and h explicitly, going the other way
!     involves solving a nonlinear equation in h as follows:
!
!            M = V/a                    ! a = speed of sound = sqrt (gamma Rs T)
!        =>  V = M a                    ! gamma = specific heat ratio, 1.4
!            Q = 0.5 rho V**2           ! Rs = specific gas constant
!              = 0.5 rho M**2 gamma Rs T
!        =>  rho T = 2Q/(M**2 gamma Rs)
!
!     From the atmosphere model, rho = rho(h) and T = T(h), so we need to find
!     the altitude h* at which rho(h*) T(h*) = rho T for the given M, Q.  Non-
!     derivative zero finder ZERORC is appropriate.  The solution is assumed to
!     be unique in a sensible altitude range.
!
!        Solve  rho(h) T(h) - rho T = 0  for h = h*
!        Then   V* = M a = M sqrt (gamma Rs T(h*))
!
!  History:
!     01/27/2017  D.A.Saunders  Initial implementation.
!     01/29/2017    "      "    Gamma was missing from the rho T expression.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real, intent (in)  :: &
      M, Q                  ! Input Mach number and dynamic pressure (Pa)

   real, intent (out) :: &
      V,                 &  ! Corresponding velocity, m/s
      h,                 &  !  "    "    "  altitude, km
      rho,               &  !  "    "    "  density, kg/m**3
      T,                 &  !  "    "    "  temperature, K
      p,                 &  !  "    "    "  pressure, Pa
      a                     !  "    "    "  speed of sound, m/s

   integer, intent (out) :: &
      ier                   ! ier /= 0 means the zero finder failed

!  Local constants:

   integer,   parameter :: maxfun = 30   ! Limit on the # function evals.
   integer,   parameter :: lunout = -6   ! Use lunout > 0 to print iterations
   real,      parameter :: gamma  = 1.4  ! Ratio of specific heats for air
   real,      parameter :: ha     = 2.0  ! Altitude range for the zero-finder
   real,      parameter :: hb     = 100. ! But 200. gives extrapolated T < 0.
   real,      parameter :: Rs     = 287.058  ! Specific gas constant, J/kg/K
   real,      parameter :: tol    = 0.0  ! Request full precision in the answer
   character (8), parameter :: subname = 'VhfromMQ'  ! If iterations are printed

!  Local variables:

   integer :: istat, lunerr, numfun
   real    :: f, hold(13), hopt, rholoc, rhoT, Tloc, ploc, aloc, Vloc

!  Execution:

   rhoT   = (Q + Q) / (M*M*gamma*Rs)  ! Target value of rho*temperature
   istat  = 2
   lunerr = abs (lunout) 
   numfun = maxfun

   do  ! Until convergence or an error

      call zerorc (ha, hb, hopt, f, tol, numfun, subname, &
                   lunout, hold, istat)

      if (istat < -1) then      ! Fatal error

         write (lunerr, '(2a, i3)') &
            subname, ' problem. Poor starting interval?  istat:', istat
         ier = istat
         exit

      else if (istat < 0) then  ! Iteration limit reached

         write (lunerr, '(2a)') &
            subname, ' hit iteration limit.  Proceeding.'
         ier = 0
         exit

      else if (istat == 0) then ! A zero has been found

         ier = 0
         exit

      else  ! istat > 0: Evaluate the function

         call earth_atmosphere (hopt, rholoc, Tloc, ploc, aloc)

         f = rholoc*Tloc - rhoT

      end if

   end do

   call earth_atmosphere (hopt, rholoc, Tloc, ploc, aloc)  ! Ensure best eval.

   V   = M*aloc
   h   = hopt
   rho = rholoc
   T   = Tloc
   p   = ploc
   a   = aloc

   end subroutine vel_alt_from_Mach_Q
