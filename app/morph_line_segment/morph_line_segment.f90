!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program morph_line_segment
!
!  Description:
!
!     For a given 3-space line segment, redistribute it to the desired number
!     of points by imposing the relative point distribution of a second 3-space
!     line segment via arc-length based local spline interpolation.
!     This capability facilitates deriving a new capsule generatrix with the
!     same point counts in the various segments as an existing generatrix, which
!     in turn can enable reusing computational flow field solutions as starting
!     guesses for the new geometry.
!
!     See also the earlier redistribute_xy for further line segment operations.
!
!  History:
!
!     10/03/2022  D.A.Saunders  Initial implementation, to simplify adapting the
!                               Mars SRL capsule DAC-2.1 as the similar DAC-2.
!
!  Aurhor:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Constants:

   integer, parameter :: &
      lundat1 = 1,       &  ! 1st line segment to be morphed
      lundat2 = 2,       &  ! 2nd line segment with desired relative arc lengths 
      lunout  = 3,       &  ! Output result
      lunkbd  = 5,       &  ! Keyboard inputs
      luncrt  = 6           ! Screen output

   logical, parameter :: &
      false  = .false.,  &
      true   = .true.

!  Variables:

   integer :: &
      i, ios, ndat1, ndat2, nin, nout

   real :: &
      a, b, ratio, stotal1, stotal2, unused

   real, allocatable, dimension (:) :: &
      sdat1, sdat2, smorph, stemp, xdat1, xdat2, xout, ydat1, ydat2, yout, &
      zdat1, zdat2, zout

   character (1) :: &
      answer
   character (128) :: &
      filein, fileout

!  Execution:

!  Read the line to morph:

   write (luncrt, '(/, a)', advance='no') &
      'Line segment to morph (two+ x/y/z triples): '
   read  (lunkbd, '(a)') filein
   open  (lundat1, file=filein, status='old', iostat=ios)
   if (ios /= 0) then
      write (luncrt, '(2a)') 'Unable to open ', trim (filein)
      go to 99
   end if
   nin = 0
   do  ! Until EOF
      read (lundat1, *, iostat=ios)
      if (ios < 0) exit
      nin = nin + 1
   end do
   rewind (lundat1)
   ndat1 = nin
   write (luncrt, '(a, i3)') 'Number of points found: ', ndat1
   allocate (xdat1(ndat1), ydat1(ndat1), zdat1(ndat1), sdat1(ndat1))

   do i = 1, ndat1
      read (lundat1, *) xdat1(i), ydat1(i), zdat1(i)  ! Allow trailing comments
   end do
   close (lundat1)

   write (luncrt, '(a)', advance='no') 'Desired number of points: '
   read  (lunkbd, *) nout

!  Read the line with the desired relative arc length distribution:

   write (luncrt, '(/, a)', advance='no') &
      'Line segment to morph (two+ x/y/z triples): '
   read  (lunkbd, '(a)') filein
   open  (lundat2, file=filein, status='old', iostat=ios)
   if (ios /= 0) then
      write (luncrt, '(2a)') 'Unable to open ', trim (filein)
      go to 99
   end if
   nin = 0
   do  ! Until EOF
      read (lundat2, *, iostat=ios)
      if (ios < 0) exit
      nin = nin + 1
   end do
   rewind (lundat2)
   ndat2 = nin
   write (luncrt, '(a, i3)') 'Number of points found: ', ndat2
   allocate (xdat2(ndat2), ydat2(ndat2), zdat2(ndat2), sdat2(ndat2))

   do i = 1, ndat2
      read (lundat2, *) xdat2(i), ydat2(i), zdat2(i)  ! Allow trailing comments
   end do
   close (lundat2)

   call chords3d (ndat1, xdat1, ydat1, zdat1, false, stotal1, sdat1)  ! 1 arcs
   call chords3d (ndat2, xdat2, ydat2, zdat2, false, stotal2, sdat2)  ! 2 arcs
!! write (51, '(i3, es16.8)') (i, sdat1(i), i = 1, ndat1)
!! write (52, '(i3, es16.8)') (i, sdat2(i), i = 1, ndat2)

   allocate (stemp(ndat2), smorph(nout))

!  Transform the line 2 arcs from [0, stotal2] to [0, stotal1]

   ratio = stotal1/stotal2
   stemp(:) = ratio * sdat2(:)
!! write (53, '(i3, es16.8)') (i, stemp(i), i = 1, ndat2)

   if (nout /= ndat2) then  ! Change the number of points:
      call changen1d (1, ndat2, stemp, 1, nout, smorph)
   else
      smorph(:) = stemp(:)
   end if

!! write (54, '(i3, es16.8)') (i, smorph(i), i = 1, nout)

!  Interpolate line 1 coordinates at the morphed arc lengths:

   allocate (xout(nout), yout(nout), zout(nout))

   call lcsfit (ndat1, sdat1, xdat1, true, 'B', nout, smorph, xout, unused)
   call lcsfit (ndat1, sdat1, ydat1, true, 'B', nout, smorph, yout, unused)
   call lcsfit (ndat1, sdat1, zdat1, true, 'B', nout, smorph, zout, unused)

!  Save the result:

   write (luncrt, '(a)', advance='no') 'Output file name: '
   read  (lunkbd, '(a)') fileout
   open  (lunout, file=fileout, status='unknown')

   write (lunout, '(3es16.8)') (xout(i), yout(i), zout(i), i = 1, nout)
   close (lunout)

99 continue

   end program morph_line_segment
