!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program redistribute_xy

!  Purpose:
!
!     Read an (x,y) line segment (2 columns, 2 or more points), assumed to be
!     geometric, and impose a specified number of points redistributed in terms
!     of arc length with specified first and last spacings. One-sided stretching
!     is also an option (spacing specified at one end or the other only).
!     Curvature-based redistribution is also an option now.
!
!  Motivation:
!
!     A capsule defined by a generatrix with corners that need to be captured
!     as well as possible needs those vertices well defined enough for splines
!     not to misbehave there.  Simply adding lots of points either side of a
!     vertex on its own may also misbehave because the curvature-based redis-
!     tribution in terms of arc length has more trouble converging as the number
!     of data points goes up.  This utility is intended to enable clustering of
!     the defining point towards a vertex with more moderate numbers of points.
!     A thruster nozzle profile prompted the curvature-based option.
!
!  Method:
!
!     Read the line segment and calculate its arc lengths.  Prompt for first
!     [and last?] desired arc lengths along with the output number of points.
!     A call to the expdis5 or vinokur utility gives the 1-sided or 2-sided
!     redistributed arc lengths from which the redistributed (x,y)s are
!     interpolated via local spline interpolation.  In the case of curvature-
!     based redistribution, curvdis2 performs that on [normalized] arc lengths
!     before doing the same sort of spline interpolation of x and y vs. s.
!
!  History:
!
!     02/03/2021  D.A.Saunders  Initial implementation to help CAPSULE_GRID
!                               (2-sided stretching only).
!     02/10/2021    "      "    Added the 1-sided stretching option that
!                               should have been there from the start.
!     03/03/2021    "      "    Added curvature-based redistribution option.
!     03/-5/2021    "      "    Arranged for suppressing or logging the
!                               (voluminous) diagnostic output from curvdis.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Constants:

   integer, parameter :: &
      lunin  = 1,        &  ! Input line segemt: 2+ (x,y)s in columns
      lunout = 2,        &  ! Output (x,y)s with desired number and stretching
      lunlog = 3,        &  ! Optional log file for curvdis output
      lunkbd = 5,        &  ! Keyboard inputs
      luncrt = 6            ! Screen prompts

   logical, parameter :: &
      false  = .false.,  &
      true   = .true.

!  Variables:

   integer :: &
      i, ier, ios, ismooth, itype, luninfo, nin, nout
   real :: &
      power, s1, s2, stotal
   real, allocatable, dimension (:) :: &
      sin, sout, unused, xin, xout, yin, yout
   logical :: &
      curvature_based, lowend, onesided
   character (1) :: &
      answer
   character (128) :: &
      filein, fileout

!  Execution:

   write (luncrt, '(/, a)', advance='no') 'Input line segment, two+ x/y pairs: '
   read  (lunkbd, '(a)') filein
   open  (lunin, file=filein, status='old', iostat=ios)
   if (ios /= 0) then
      write (luncrt, '(2a)') 'Unable to open ', trim (filein)
      go to 99
   end if
   nin = 0
   do  ! Until EOF
      read (lunin, *, iostat=ios)
      if (ios < 0) exit
      nin = nin + 1
   end do
   rewind (lunin)
   write (luncrt, '(a, i3)') 'Number of points found: ', nin
   allocate (xin(nin), yin(nin), sin(nin))
   do i = 1, nin
      read (lunin, *) xin(i), yin(i)  ! Allow for trailing comments
   end do
   close (lunin)

   write (luncrt, '(a)', advance='no') 'Desired number of points: '
   read  (lunkbd, *) nout
   allocate (xout(nout), yout(nout), sout(nout), unused(nout))
   write (luncrt, '(a)', advance='no') 'Output file name: '
   read  (lunkbd, '(a)') fileout
   open  (lunout, file=fileout, status='unknown')

   call chords2d (nin, xin, yin, .false., stotal, sin)  ! Input arc length(s)

   sout(1)    = sin(1)
   sout(nout) = stotal
   write (luncrt, '(a, es16.8)') 'Input total arc length:', stotal

   write (luncrt, '(a)', advance='no') &
      '1-sided or 2-sided stretching or curvature-based? [1|2|3]: '
   read  (lunkbd, *) itype
   onesided = itype == 1
   curvature_based = itype == 3

   if (curvature_based) then

      write (luncrt, '(a)', advance='no') &
         'Exponent in [0, 1] to try first; higher means more kappa effect: '
      read  (lunkbd, *) power
      write (luncrt, '(2a)', advance='no') &
         'Smoothing control: 0=none; 1=curvature; 2=redistributed arcs; ', &
         '3=both: '
      read  (lunkbd, *) ismooth
      write (luncrt, '(a)', advance='no') &
         'Save curvdis iteration history and diagnostic output? y|n: '
      read  (lunkbd, *) answer
      luninfo = -lunlog
      if (answer == 'y' .or. answer == 'Y') then
         luninfo = lunlog
         open (lunlog, file='redistribution.log', status='unknown')
      end if

      call curvdis2 (nin, xin, yin, nout, power, ismooth, luninfo, &
                     false, xout, yout, ier)

   else if (onesided) then

      write (luncrt, '(a)', advance='no') &
         'Stretch from first point or last point [f|l]: '
      read  (lunkbd, *) answer
      lowend = answer == 'f' .or. answer == 'F'
      write (luncrt, '(a)', advance='no') 'Initial arc length: '
      read  (lunkbd, *) s1
      call expdis5 (1, sin(1), stotal, s1, nout, sout, -luncrt)

      if (.not. lowend) then ! Reverse the clustering
         do i = 1, nout
            unused(i) = stotal - sout(i)
         end do
         sout(:) = unused(:)
      end if

   else

      write (luncrt, '(a)', advance='no') 'First and last arc lengths: '
      read  (lunkbd, *) s1, s2
      call vinokur (1, nout, s1, s2, sout, luncrt, ier)

   end if

   if (.not. curvature_based) then
      call lcsfit (nin, sin, xin, true, 'B', nout, sout, xout, unused) ! x vs. s
      call lcsfit (nin, sin, yin, true, 'B', nout, sout, yout, unused) ! y vs. s
   end if

   write (lunout, '(2es16.8)') (xout(i), yout(i), i = 1, nout)

99 continue

   end program redistribute_xy
