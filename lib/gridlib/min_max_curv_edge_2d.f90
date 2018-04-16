!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine min_max_curv_edge_2d (n1, x1, y1, t1, i1, outward1,  &
                                    n2, x2, y2, t2, il, ir, outward2, topt, &
                                    n3, x3, y3, t3, ier)
!  Purpose:
!
!     Determine the curve orthogonal to two given discretized curves through a
!     specified point on the first curve (such as one end) and the presumably
!     unique point in a specified interval of the second curve that minimizes
!     the maximum curvature of the desired curve.  This third curve will be
!     defined by just 6 points (3 at each end as needed to force orthogonality)
!     but it will be returned in discrete form with the specified number of
!     points, n3, as used for the finite difference curvature calculations.
!
!     Such a functionality was prompted by program HB_GRID (hypersonic blunt
!     body grid generator) for its outflow boundaries between the shoulder and
!     the estimated bow shock curve.
!
!  Method:
!
!     Non-derivative minimizer FMINRC serves to solve what is cast as a 1-D
!     optimization problem: locate the arc length distance along curve 2 that
!     produces the minimum peak curvature on curve 3.
!
!     Determining which side of curves 1 and 2 to construct normals on is the
!     main awkwardness.  Current strategy:  define outward as meaning the
!     cross product t x n has positive z component (out of the page), where
!     t is a tangent vector in the direction of increasing arc length, and
!     n is the desired normal.  It can be shown that (nx ny) = (-ty tx).
!
!     Discretizing each 6-point curve via uniform chord length will not produce
!     uniform arc-length spacing, but that should not matter, since the intent
!     is for the points on the output third curve to be redistributed by the
!     application anyway.
!
!  History:
!
!     01/03/2011  D.A.Saunders  Initial implementation.
!     01/05/2011    "     "     Introduced outward arguments as described above.
!
!  Author:
!
!     David Saunders, ERC, Inc./NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)    :: n1        ! Number of points on input curve 1
   real,    intent (in)    :: x1(n1)    ! Corresponding abscissas
   real,    intent (in)    :: y1(n1)    ! Corresponding ordinates
   real,    intent (in)    :: t1(n1)    ! Corresponding chord lengths
   integer, intent (in)    :: i1        ! Index of fixed point on curve 1
   logical, intent (in)    :: outward1  ! T or F as described above, curve 1
   integer, intent (in)    :: n2        ! Number of points on input curve 2
   real,    intent (in)    :: x2(n2)    ! Corresponding abscissas
   real,    intent (in)    :: y2(n2)    ! Corresponding ordinates
   real,    intent (in)    :: t2(n2)    ! Corresponding chord lengths
   integer, intent (in)    :: il        ! Lower index of curve 2 search interval
   integer, intent (in)    :: ir        ! Upper index
   logical, intent (in)    :: outward2  ! T or F as described above, curve 2
   real,    intent (out)   :: topt      ! Optimized chord length along curve 2
   integer, intent (in)    :: n3        ! Number of points on output curve 3
   real,    intent (inout) :: x3(n3)    ! Corresponding abscissas
   real,    intent (inout) :: y3(n3)    ! Corresponding ordinates
   real,    intent (inout) :: t3(n3)    ! Corresponding chord lengths
   integer, intent (inout) :: ier       ! Input 0 suppresses iteration printing
                                        ! Output 0 => no error detected;
                                        !       -1 => 1-D minimizer fatal error
!  Constants:

   real,      parameter :: half   = 0.5
   real,      parameter :: zero   = 0.
   integer,   parameter :: nfmax  = 50   ! Limit on # function evaluations
   logical,   parameter :: false  = .false.
   logical,   parameter :: new    = .true.
   logical,   parameter :: true   = .true.
   character, parameter :: caller * 20 = 'min_max_curv_edge_2d'
   character, parameter :: method * 1  = 'B'  ! Loose "Bessel" spline fits

!  Variables:

   integer :: istat, lunerr, lunout, numfun
   real    :: dt, dx, dxdt, dy, dydt, peakcrv, tmag, tol, total, &
              t6(6), x6(6), y6(6), xm, xt, yt, ym
   real    :: curvature(n3), xp(n3), xpp(n3), yp(n3), ypp(n3)
!! logical :: first

!  Procedures:

   external :: chords2d ! Cumulative chord lengths for a 2-space curve
   external :: fminrc   ! Reverse-communication 1-D minimizer
   external :: lcsfit   ! Local cubic spline utility (including 1st derivative)

!  Execution:

!! first  = true
   lunout = -6
   if (ier /= 0) lunout = -lunout

!  A minimum can be found to within sqrt (machine epsilon), but avoid the sqrt:

   if (epsilon (tol) < 1.e-10) then
      tol = 1.e-8
   else
      tol = 1.e-4
   end if

!  Pick a reasonable arc length on curve 3 for the points forcing orthogonality:

   xm = half * (x2(il) + x2(ir))
   ym = half * (y2(il) + y2(ir))
   dt = sqrt ((x1(i1) - xm)**2 + (y1(i1) - ym)**2) / real (n3)
   dt = dt + dt

!  Find a gradient vector for curve 1 at point i1:

   call lcsfit (n1, t1, x1, new, method, 1, t1(i1), xt, dxdt)
   call lcsfit (n1, t1, y1, new, method, 1, t1(i1), yt, dydt)

   tmag = sqrt (dxdt**2 + dydt**2)  ! Magnitude of tangent vector
   dx   = -dydt * (dt/tmag)         ! Components of outward unit normal
   dy   =  dxdt * (dt/tmag)         ! are (-ty tx); scale them by dt

   if (.not. outward1) then
      dx = -dx
      dy = -dy
   end if

!  First 3 points of curve 3 orthogonal to curve 1 at i1:

   x6(1) = x1(i1);      y6(1) = y1(i1)
   x6(2) = x6(1) + dx;  y6(2) = y6(1) + dy
   x6(3) = x6(2) + dx;  y6(3) = y6(2) + dy

!! write (6, *) 'i1: ', i1
!! write (6, *) 'xt: ', xt
!! write (6, *) 'yt: ', yt
!! write (6, *) 'dxdt: ', dxdt
!! write (6, *) 'dydt: ', dydt
!! write (6, *) 'dt: ', dt
!! write (6, *) 'dx: ', dx
!! write (6, *) 'dy: ', dy
!! write (6, *) 'x6(1), y6(1): ', x6(1), y6(1)
!! write (6, *) 'x6(2), y6(2): ', x6(2), y6(2)
!! write (6, *) 'x6(3), y6(3): ', x6(3), y6(3)

!  Perform the optimization:

   lunerr = abs (lunout)
   numfun = nfmax      ! Limit; FMINRC typically takes 12-14 iterations
   istat  = 2          ! Initialize the minimization
   ier    = 0          ! Upon exit, normally

   do while (istat > 0)

      call fminrc (t2(il), t2(ir), topt, peakcrv, tol, numfun, caller, lunout, &
                   istat)

      if (istat < -1) then

         write (lunerr, '(/, 2a)') caller, ': FMINRC fatal error'
         ier = -1

      else if (istat < 0) then ! Iteration limit; may be usable

         write (lunerr, '(/, 2a)') caller, ': Iteration limit.'
         istat = 0

      else  ! istat >= 0       ! Evaluate the objective function at topt

         call max_curvature () ! Find the maximum curvature on curve 3 for topt

         if (istat == 0) exit  ! Success; the desired curve 3 is in x3, y3, t3
      end if

   end do

!! write (6, '(a, 1p, 3e15.7)') ' tl, t, tr: ', t2(il), topt, t2(ir)

99 return

   contains

!     Internal procedures for subroutine min_max_curv_edge_2d:

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine max_curvature ()

!     Complete the 6-point curve for this topt, discretize it more finely, and
!     return its peak curvature magnitude as the objective function.

!     Local variables:

      integer :: i
      real    :: du

!     Execution:

!     Find a gradient vector for curve 2 at point t = topt:

      call lcsfit (n2, t2, x2, new, method, 1, topt, xt, dxdt)
      call lcsfit (n2, t2, y2, new, method, 1, topt, yt, dydt)

      tmag = sqrt (dxdt**2 + dydt**2)  ! Magnitude of tangent vector
      dx   = -dydt * (dt/tmag)         ! Components of outward unit normal
      dy   =  dxdt * (dt/tmag)         ! are (-ty tx); scale them by dt

      if (.not. outward2) then
         dx = -dx
         dy = -dy
      end if

!     Last 3 points of the 6-point form of curve 3 normal to curve 2 at topt:

      x6(6) = xt;          y6(6) = yt
      x6(5) = x6(6) + dx;  y6(5) = y6(6) + dy
      x6(4) = x6(5) + dx;  y6(4) = y6(5) + dy

!!    if (first) then
!!       first = false
!!       write (6, *) 'x6(4), y6(4): ', x6(4), y6(4)
!!       write (6, *) 'x6(5), y6(5): ', x6(5), y6(5)
!!       write (6, *) 'x6(6), y6(6): ', x6(6), y6(6)
!!    end if

!     Discretize this 6-point form for higher resolution of its curvature:

      call chords2d (6, x6, y6, false, total, t6)

      du = total / real (n3 - 1);  t3(n3) = total

      do i = 1, n3 - 1
         t3(i) = du * real (i - 1)
      end do

      call lcsfit (6, t6, x6, new, method, n3, t3, x3, xp)  ! This xp and ...
      call fd12k  (n3, t3, x3, xp, xpp, curvature)          ! this curvature ...
      call lcsfit (6, t6, y6, new, method, n3, t3, y3, yp)  ! are not used;
      call fd12k  (n3, t3, y3, yp, ypp, curvature)          ! likewise for y
      call curv2d (n3, xp, xpp, yp, ypp, curvature)

      peakcrv = zero
      do i = 2, n3 - 1
         peakcrv = max (peakcrv, abs (curvature(i)))
      end do

      end subroutine max_curvature

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   end subroutine min_max_curv_edge_2d
