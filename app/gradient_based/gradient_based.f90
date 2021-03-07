!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program gradient_based
!
!     This is a driving program for subroutines GRADDIS2 and GRADDIS and their
!     3-space analogues, which are intended to redistribute grid points along
!     a 2- or 3-space line or curve based on the gradient magnitudes of the
!     function accompanying the (x,y[,z]) points.
!
!     It was prompted by a desire to thin the line-of-sight data used to perform
!     radiative heating calculations via NEQAIR, in a somewhat intelligent way.
!
!  Input data format:
!
!     x1  y1  [z1] f1
!     x2  y2  [z2] f2      Arc lengths along this curve will be redistributed
!     :   :   :    :
!     xn  yn  [zn] fn
!
!  Output:
!     X1  Y1  [z1] F1 
!     X2  Y2  [z2] F2      Here, m is probably less than n, but not necessarily
!     :   :   :    :
!     Xm  Ym  [zm] Fm
!
!     Chances are that multiple functions will need to be interpolated to the
!     new coordinates, but interpolating the driving function tests the new
!     variant of LCSFIT (LCSFIT2, which normalizes the data beforehand).
!
!  History:
!
!     11/21/2013  D.A.Saunders  Initial test program for GRADDIS[2] ...
!     11/22/2013    "      "    ... and for LCSFIT2 (interpln. w/ normalzn.).
!     12/03/2013    "      "    Added the 3-space option.
!
!  Author:  David Saunders, ERC Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Local constants:

   integer, parameter :: &
      lunin  = 1,        &     ! Input data
      lunout = 2,        &     ! Redistributed data
      lunkbd = 5,        &     ! Keyboard inputs
      luncrt = 6               ! Prompts and diagnostics
   logical, parameter :: &
      normalize = .false.      ! Don't normalize arc lengths here
   character (1), parameter :: &
      method = 'M'             ! Monotonic local cubic interpolations

!  Local variables:

   integer :: &
      i, ier, ios, ismooth, ndat, ndim, nnew
   real :: &
      power, total
   logical :: &
      twod
   character (64) :: &
      filename
   real, allocatable, dimension (:) ::&
      fdat, xdat, ydat, zdat, sdat, xnew, ynew, znew, fnew, snew

!  Procedures:

   external :: &
      graddis2,   &            ! Normalizes data; calls GRADDIS in a power loop
      graddis3d2, &            ! 3-space analogue of GRADDIS2; calls GRADDIS3D
      lcsfit2                  ! Normalizes data; calls LCSFIT

!  Execution:

   write (luncrt, '(a)', advance='no') '(x,y[,z],f) data file: '
   read  (lunkbd, *) filename
   open  (lunin, file=filename, status='old')

   i = 0
   do  ! Until EOF
      read (lunin, *, iostat=ios)
      if (ios < 0) exit
      i = i + 1
   end do
   ndat = i

!  Try to read 3-space data; interpret as 2-space if an error occurs:

   allocate (xdat(ndat), ydat(ndat), zdat(ndat), fdat(ndat), sdat(ndat))

   rewind (lunin)
   read   (lunin, *, iostat=ios) &
      (xdat(i), ydat(i), zdat(i), fdat(i), i = 1, ndat)

   twod = ios /= 0

   if (twod) then
      deallocate (zdat)
      rewind (lunin)
      read   (lunin, *, iostat=ios) (xdat(i), ydat(i), fdat(i), i = 1, ndat)
      if (ios /= 0) then
         write (luncrt, '(a)') 'Cannot read dataset as either 3 or 4 columns.'
         go to 99
      end if
   end if

   close (lunin)

   write (luncrt, '(a, i5, a)', advance='no') &
      '# data points found:', ndat, '; desired #: '
   read  (lunkbd, *) nnew

   write (luncrt, '(2a)', advance='no') 'Initial exponent to use in ', &
      'gradient-based shape fn. [0. <= e <= 2.; try 1.]: '
   read  (lunkbd, *) power

   write (luncrt, '(a)', advance='no') &
      'Smoothing? [0=none, 1=shape fn.; 2=revised arcs; 3=both; try 1]: '
   read  (lunkbd, *) ismooth

   allocate (xnew(nnew), ynew(nnew), fnew(nnew), snew(nnew))

   if (twod) then
      call chords2d (ndat, xdat, ydat, normalize, total, sdat)  ! Arcs, unnorm.
      call graddis2 (ndat, xdat, ydat, sdat, fdat, nnew, power, &
                     ismooth, luncrt, xnew, ynew, snew, ier)
   else
      allocate (znew(nnew))
      call chords3d (ndat, xdat, ydat, zdat, normalize, total, sdat)
      call graddis3d2 (ndat, xdat, ydat, zdat, sdat, fdat, nnew, power, &
                       ismooth, luncrt, xnew, ynew, znew, snew, ier)
   end if

!  Interpolate the function at the new arc lengths:

   call lcsfit2 (ndat, sdat, fdat, method, nnew, snew, fnew)

   open (lunout, file='gradient_based.dat', status='unknown')

   if (twod) then
      write (lunout, '(a)') '#  x  y  f  s'
      write (lunout, '(4es16.8)') &
         (xnew(i), ynew(i), fnew(i), snew(i), i = 1, nnew)
   else
      write (lunout, '(a)') '#  x  y  z  f  s'
      write (lunout, '(5es16.8)') &
         (xnew(i), ynew(i), znew(i), fnew(i), snew(i), i = 1, nnew)
      deallocate (zdat, znew)
   end if

   close (lunout)

   deallocate (xdat, ydat, fdat, xnew, ynew, snew, fnew)

99 continue

   end program gradient_based
