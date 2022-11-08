!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      program reshape
!
!  Purpose:
!
!        RESHAPE is a utility to transform (x,y) data in one or more of
!     a variety of ways such as scaling, rotating, reversing order, etc.
!     Since x and y may be in specified columns other than 1 and 2,
!     RESHAPE can also serve the purpose of extracting two columns of
!     data from many without other changes.  (But now, COLUMNEDIT can do
!     the same thing and more.)
!
!  Method:
!
!     > Deal with just one dataset at a time.  The formats chosen are
!       those of program SMOOTH:
!
!          [Title]             <Optional title, up to 80 characters>
!          [n]                 <Optional no. of pts - EOF also serves>
!                              <Blank lines are ignored>
!          x(1)  y(1)
!           :     :            <X and Y may be extracted from specified
!           :     :             columns if there are more than two>
!           :     :
!          ! x     y           <! suppresses points or adds comments>
!          x(n)  y(n)
!
!     > Transformations are done in-place (y and/or x).
!
!     > "Undo last" and "start over" operations are done with spare copies.
!
!  Procedures:
!
!     ALPHA      Distinguishes text and numeric data
!     CHANGEN2D  Change the number of points proportionally
!     GETLINE    Gets next significant line
!     OPENER     File opening utility
!     RDLIST     Gets an indefinite number of integers
!     RDXYZ      Gets 2 or 3 columns from many
!     READER     Prompting utility (entry pts. READI, READR, etc.)
!
!  History:
!
!     08/29/86   DAS   Initial implementation (in haste) -
!                      simple Y-translation or scaling.
!     09/03/86   DAS   Added rotation and X-translation/scaling.
!     10/21/86   DAS   Using EOF to determine number of points turns
!                      out to be inconvenient for other utilities.
!                      Expect N to be with the data.  (Changed later.)
!     08/19/88   DAS   Added option to reverse the order 1:N.
!     10/21/88   DAS   Bug fix for reverse-order option; added
!                      "start over", "switch X and Y", and
!                      "reflect" options.
!     10/13/92   DAS   Handled multi-column files with ! comments or
!                      blank lines ignored.  N and title are now optional.
!     01/03/98   DAS   Provided for rotating Y about (Zc=0, Yc).
!     05/20/99   DAS   Minor Fortran 90 changes.
!     08/30/10   DAS   The advance='no' prompting was misbehaving; added a
!                      "Done" item to the menu.
!     05/24/11   DAS   64-bit precision outputs now, not single precision.
!     07/10/21   DAS   In order to test the revised CHANGEN2D, install
!                      it as one more option here, as first done for testing
!                      CHANGEN via RESHAPE3D.  Add saving of before and
!                      after cell growth rates for this option.
!     12/08/21   DAS   Minimal Fortran 90 translation, prompted by a new option
!                      to tabulate angles between adjacent points (in turn,
!                      prompted by a capsule generatrix issue).
!
!  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
!           Later with ELORET, Inc. and AMA, Inc. at NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use trigd
      implicit none

!     Constants:

      integer, parameter :: &
         luncrt = 6,        &
         lunkbd = 5,        &
         lunin  = 7,        &
         lunout = 8,        &
         mxmenu = 14,       &
         mxpts  = 200000,   &
         halfm  = (mxmenu + 4) / 2 ! 4 here allows for -2, -1, 0, and mxmenu+1

!     Variables:

      integer :: &
         choice, columns(2), i, ier, ios, j, last, n, ncol, nnew
      real :: &
         angle, ci, dx, dy, growth, p, q, scale, shift, si, temp, total, xp, yp
      real, dimension (mxpts) :: &
         x, y, xlast, ylast, xorig, yorig
      real, allocatable, dimension (:) :: &
         arc, xnew, ynew
      logical :: &
         cr, eof, finis, yestitle
      character :: &
         dataset*128, menu(-2 : mxmenu + 1)*33, method*1, title*80

!     Procedures:

      logical :: &
         alpha
      external :: &
         alpha, changen2d, chords2d, getline, opener, rdlist, rdxyz, &
         readi, readr, reads

!     Storage:

      data menu / &
         '  -2: Start over',                  &
         '  -1: Undo last transformation',    &
         '   0: Review data',                 &
         '   1: Translate Y',                 &
         '   2: Scale Y',                     &
         '   3: Translate X',                 &
         '   4: Scale X',                     &
         '   5: Rotate (X,Y) about (p,q)',    &
         '   6: Reflect about the X-axis',    &
         '   7: Reflect about the Y-axis',    &
         '   8: Reverse the order (1:N)',     &
         '   9: Switch X and Y',              &
         '  10: Rotate Y about (Zc=0,Yc)',    &
         '  11: Scale X & Y the same way',    &
         '  12: Change N; same rel. sp.',     &
         '  13: Display adjacent pt. angles', &
         '  14: Done',                        &
         '                              '/    ! Last blank eases display
                                              ! of menu as two columns.
!     Execution:

      write (luncrt, 1000) &
         'RESHAPE applies simple transformations to an (X,Y) dataset.', &
         'Input and output files are in [program] SMOOTH format(s).', &
         ' '

!     !!!!!!!!!!!!!
!     Get the data:
!     !!!!!!!!!!!!!

      call opener (luncrt, 'Dataset name: ', lunkbd, dataset, lunin, 'old')

!     Look for the optional title:

  120 call getline (lunin, '!', title, last, ios)

      if (ios  /= 0) go to 900
      if (last == 0) go to 120

!     If we have strictly numeric data, assume there is no title:

      yestitle = alpha (title)
      if (.not. yestitle) then
         rewind lunin         ! RDXYZ should have an optional input buffer
      end if

!     Too awkward to avoid a prompt here (read/tokenize/backspace):

      columns(1) = 1
      columns(2) = 2
      ncol = 2
      call rdlist (luncrt, 'X and Y column numbers?  [1, 2]: ', &
                   lunkbd, ncol, columns)
      if (ncol < 0) go to 999           ! EOF = Quit

      write (luncrt, 1001)

!     Read the data proper:

      call rdxyz (2, lunin, luncrt, columns, mxpts, n, x, y, y, finis, ier)
      if (ier /= 0) go to 900  ! RDXYZ explains any errors
      if (n == 0)   go to 901  ! Doesn't make sense in this application

!     Keep some copies around:

      do i = 1, n
         xorig(i) = x(i)
         yorig(i) = y(i)
         xlast(i) = x(i)
         ylast(i) = y(i)
      end do


!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Loop over possibly several transformations per run:
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  200 continue

         write (luncrt, 1004) &                        ! Basically, I=1:MXMENU/2
            (menu(i), menu(i + halfm), i = -2, halfm - 3)

         if (mod (mxmenu, 2) == 0) write (luncrt, 1001)

  210    call readi (luncrt, 'Pick one. EOF (^Z or ^D) means no more. ', &
                     lunkbd, choice, cr, eof)
         if (eof) go to 800
         if (cr)  go to 210

         if (choice > 0 .or. choice == -2) then   ! Save current values:
            do i = 1, n
               xlast(i) = x(i)
               ylast(i) = y(i)
            end do
         end if

         if (choice == -2) then      ! "Start over" from scratch:
            do i = 1, n
               x(i) = xorig(i)
               y(i) = yorig(i)
            end do

         else if (choice == -1) then ! "Undo" previous operation:

            do i = 1, n
               x(i) = xlast(i)
               y(i) = ylast(i)
            end do

         else if (choice == 0) then  ! "Review": Display data pairs.

            write (luncrt, '(/, (i6, 2es16.7))') (i, x(i), y(i), i = 1, n)

         else if (choice == 1) then  ! "Translate Y":

            call readr (luncrt, '   Enter Y shift (+ or -): ', lunkbd, shift, &
                        cr, eof)
            if (cr .or. eof) go to 210

            do i = 1, n
               y(i) = y(i) + shift
            end do

         else if (choice == 2) then  ! "Scale Y":

            call readr (luncrt, '   Enter Y scale (+ or -): ', lunkbd, scale, &
                        cr, eof)
            if (cr .or. eof) go to 210

            do i = 1, n
               y(i) = y(i) * scale
            end do

         else if (choice == 3) then  ! "Translate X":

            call readr (luncrt, '   Enter X shift (+ or -): ', lunkbd, shift, &
                        cr, eof)
            if (cr .or. eof) go to 210

            do i = 1, n
               x(i) = x(i) + shift
            end do

         else if (choice == 4) then  ! "Scale X":

            call readr (luncrt, '   Enter X scale (+ or -): ', lunkbd, scale, &
                        cr, eof)
            if (cr .or. eof) go to 210

            do i = 1, n
               x(i) = x(i) * scale
            end do

         else if (choice == 5) then  ! "Rotate (x,y) about (p,q)":

            call readr (luncrt, &
                      '   Enter rotation in degrees (+ve is anticlockwise): ', &
                        lunkbd, angle, cr, eof)
            if (cr .or. eof) go to 210

  350       write (luncrt, 1001, advance='no') &
               '    Enter center of rotation (Xc,Yc) as two values: '
            read (lunkbd, *, err=350, end=210) p, q

            angle = angle * asin (1.e+0) / 90.e+0
            ci = cos (angle)
            si = sin (angle)

            do i = 1, n
               xp  = x(i) - p
               yp  = y(i) - q
               x(i) = xp * ci - yp * si + p
               y(i) = xp * si + yp * ci + q
            end do

         else if (choice == 6) then  ! "Reflect Y about the X-axis":

            do i = 1, n
               y(i) = -y(i)
            end do

         else if (choice == 7) then  ! "Reflect X about the Y-axis":

            do i = 1, n
               x(i) = -x(i)
            end do

         else if (choice == 8) then  ! "Reverse the order":

            do i = 1, (n + 1) / 2
               j = n + 1 - i
               temp = x(i)
               x(i) = x(j)
               x(j) = temp
               temp = y(i)
               y(i) = y(j)
               y(j) = temp
            end do

         else if (choice == 9) then  ! "Switch X and Y":

            do i = 1, n
               temp = x(i)
               x(i) = y(i)
               y(i) = temp
            end do

         else if (choice == 10) then  ! "Rotate Y about (Z=0,Yc)":

            call readr (luncrt, &
                      '   Enter rotation in degrees (+ve is anticlockwise): ', &
                        lunkbd, angle, cr, eof)
            if (cr .or. eof) go to 210

            call readr (luncrt, &
                        '   Enter Yc for center of rotation (Zc=0,Yc): ', &
                        lunkbd, q, cr, eof)
            if (cr .or. eof) go to 210

            angle = angle * asin (1.e+0) / 90.E+0
            ci = cos (angle)

            do i = 1, n
               y(i) = (y(i) - q) * ci + q
            end do

         else if (choice == 11) then  ! "Scale X & Y the same way":

            call readr (luncrt, '   Enter scale (+ or -): ', &
                        lunkbd, scale, cr, eof)
            if (cr .or. eof) go to 210

            do i = 1, n
               x(i) = x(i) * scale
               y(i) = y(i) * scale
            end do

         else if (choice == 12) then  ! "Change N; same relative spacing":

            write (luncrt, '(a, i7)') '   Current number of points: ', N
            call readi (luncrt,        '  Desired number of points: ', &
                        lunkbd, nnew, cr, eof)
            if (cr .or. eof) go to 210

            write (luncrt, '(a)') '   Interpolation methods:', &
                                  '   M (monotonic), B (loose), L (Linear)'
            call readc (luncrt, '  Interpolation choice: ', &
                        lunkbd, method, cr, eof)
            if (cr .or. eof) go to 210

!           Save the current growth rates:

            open (lunout, file='growth-rates-before.dat', status='unknown')

            allocate (arc(n))

            call chords2d (n, x, y, .false., total, arc)

            do i = 2, n-1
               growth = (arc(i+1) - arc(i)) / (arc(i) - arc(i-1))
               write (lunout, '(2es16.8)') arc(i), growth
            end do

            close (lunout)

            deallocate (arc)
            allocate (xnew(nnew), ynew(nnew))

            call changen2d (1, n, x, y, 1, nnew, xnew, ynew, method)

            n = nnew
            x(1:n) = xnew(:)
            y(1:n) = ynew(:)

            deallocate (xnew, ynew)

            open (lunout, file='growth-rates-after.dat', status='unknown')

            allocate (arc(n))

            call chords2d (n, x, y, .false., total, arc)

            do i = 2, n-1
               growth = (arc(i+1) - arc(i)) / (arc(i) - arc(i-1))
               write (lunout, '(2es16.8)') arc(i), growth
            end do
            close (lunout)
            deallocate (arc)

         else if (choice == 13) then  ! "Display angles between adjacent pts.":

            write (luncrt, '(/, a)') &
               '   i               x               y           angle'

            do i = 1, n-1
               dx = x(i+1) - x(i)
               dy = y(i+1) - y(i)
               angle = atan2d (dy, dx)
               write (luncrt, '(i4, 3es16.8)') i, x(i), y(i), angle
            end do

         else if (choice == 14) then  ! Done (works better than ^D)

            go to 800

         else

            go to 210

         end if

      go to 200

!     !!!!!!!!!!!!!
!     Save results:
!     !!!!!!!!!!!!!

  800 write (luncrt, 1001)

      call opener (luncrt, 'Output file name?  EOF = quit: ', &
         lunkbd, dataset, lunout, 'new')

      write (luncrt, 1001)

      if (yestitle) then
         call reads (luncrt, 'Output title line? <CR>=same: ', &
                     lunkbd, title, cr, eof)

         i = len_trim (title)
         write (lunout, 1001, err=902) title(1:i)
         write (lunout, 1003, err=902) n
      end if

      write (lunout, '(2es24.15)', err=902) (x(i), y(i), i = 1, n)

      go to 999


!     Error handling:

  900 write (luncrt, 1000) 'RESHAPE: Error reading the data.'
      go to 999
  901 write (luncrt, 1000) 'RESHAPE: No data points found.'
      go to 999
  902 write (luncrt, 1000) 'RESHAPE: Error saving results.'

  999 continue  ! Avoid system-dependencies with STOP


!     Formats:

 1000 format (/, (1x, a))
 1001 format (a)
 1002 format (a, i3)
 1003 format (i5)
 1004 format (/, (a, 10x, a))

      end program reshape
