!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine catenary_parameters (semispan, deflection, a, ier)

!  Calculate the parameter a of the catenary curve  y  =  a cosh (x/a) - a
!  such that its maximum deflection below its height at x = semispan is the
!  given value.
!
!  If b = semispan and d = deflection, we have to solve the following equation:
!
!     0.5 a (exp (b/a) + exp (-b/a)) - a  =  d
!
!  Heuristically, a good starting guess for a is 0.5 b**2/d.
!
!  10/05/11  D. A. Saunders  Initial implementation, prompted by flexible
!                            thermal protection studies for Venus entry.
!  10/06/11    "       "     Avoiding two exponentials may be unwise: use COSH.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real,    intent (in)  :: semispan   ! Half width of catenary centered on Oy
   real,    intent (in)  :: deflection ! Catenary height y at x = +/-semispan
   real,    intent (out) :: a          ! Desired catenary parameter
   integer, intent (out) :: ier        ! ier /= 0 means the zero finder failed

!  Local constants:

   integer,   parameter :: maxfun = 30 ! Limit on the expected # function evals.
   integer,   parameter :: lunout = -6 ! Enter lunout < 0 to suppress iterations
   real,      parameter :: half = 0.5
   real,      parameter :: one5 = 1.5
   real,      parameter :: tol  = 0.0  ! Request full precision in the answer
   character, parameter :: subname*8 = 'CATENARY'  ! If iterations are printed

!  Local variables:

   integer :: istat, lunerr, numfun
   real    :: aleft, aright, f, hold(13)

!  Execution:

   istat  = 2                                ! Initialize the iteration
   aleft  = half * semispan**2 / deflection  ! Heuristic search interval
   aright = one5 * aleft
   aleft  = half * aleft
   lunerr = abs (lunout)
   numfun = maxfun

   do  ! Until convergence or an error

      call zerorc (aleft, aright, a, f, tol, numfun, subname, &
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

         f = a * cosh (semispan / a) - a - deflection

      end if

   end do

   end subroutine catenary_parameters
