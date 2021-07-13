!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      program reshape3d
!
!  PURPOSE:
!
!     Utility to transform XYZ data in one or more ways at a time.
!
!  METHOD:
!
!     > Deal with just one dataset at a time.  The simple format chosen is
!       compatible with other utilities ("SMOOTH2D" format).  Actually, if
!       the first line is found to be strictly numeric, the title and point
!       count are assumed to be omitted.
!
!          [TITLE
!          N]                  <No. of pts - trailing text here is ignored>
!          X (1)  Y (1)  Z (1)
!           :     :            <Read list-directed, one triple at a time>
!           :     :            <Additional columns will be ignored/lost>
!          X (N)  Y (N)  Z (N)
!
!     > Transformations are done in-place. 
!
!     > "Undo last" and "start over" operations are done with spare copies.
!
!  HISTORY:
!
!     08/29/86   David Saunders   Initial implementation of RESHAPE for
!                                 XY data. 
!
!     01/17/88   Michael Wong     RESHAPE3D developed from RESHAPE.
!                                 Added features are translation and
!                                 reflection about the XY, YZ, or ZX planes,
!                                 Z shifting and scaling, and YZ or ZX
!                                 data switching.  Deleted feature is
!                                 rotation about point (p,q).
!
!     12/24/97   DAS              Allowed for rotating Y and Z about (p,q),
!                                 as needed for pylon or fin sections.
!                                
!     11/25/98   DAS              Allowed for up to 10,000 points.
!
!     05/19/99   DAS              Minor Fortran 90 revisions.
!
!     02/18/00   DAS              Added rotations about axes parallel to
!                                 the Y and Z axes for completeness.
!
!     06/16/08   DAS              Added general rotation about an axis defined
!                                 by two points.  Free formatting now, and
!                                 dynamic allocation of work-space.  Guard
!                                 against displaying excessive amounts of data
!                                 for the "review" option.  Disallow starting
!                                 over and undoing if the number of points is
!                                 too large (as ~375000 seems to be).
!
!     07/24/08   DAS              Display the data range and center (prompted
!                                 by dealing with clouds of laser-scanned data
!                                 rather than 3-space curves).
!
!     08/13/08   DAS              Added splitting options.  Once a dataset is
!                                 split, only the upper portion of the repacked
!                                 coordinate arrays is operated on if further
!                                 changes are specified.
!
!     10/07/11   DAS              Added the option to apply a rigid transform-
!                                 ation defined by new curve end points (rigid
!                                 if they're the same distance apart as the
!                                 input end points, that is, otherwise some
!                                 other result is obtained that may nevertheless
!                                 be of interest).
!
!     08/21/12   DAS              Added NULINE3D option for comparison with
!                                 RIGID_TRANSFORM (which indeed gives useful
!                                 results even if the transformation isn't
!                                 really rigid because the size is changing).
!
!     04/06/18   DAS              Full precision output is long overdue.
!
!     07/10/21   DAS              In order to test the revised CHANGEN, install
!                                 it as one more option here.  Raise the "undo"
!                                 limit significantly.  Add saving of before and
!                                 after cell growth rates for this option.
!
!  AUTHORS (Original): David Saunders, Michael Wong,   Sterling Software/ARC, CA
!          (Later):    David Saunders, ELORET Corp/NASA Ames Research Center, CA
!                      Now with AMA, Inc. at NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Constants:

   integer, parameter :: &
      luncrt   = 6, lunkbd = 5, lunin = 7, lunout = 8, &
      mxmenu   = 24,     &
      ndisplay = 30,     &  ! Don't display huge numbers of data points
      nundo    = 999999, &  ! Allocating an extra 6 x N can be extravagant
      halfm    = (mxmenu + 4) / 2 ! 4 here allows for -2, -1, 0, and mxmenu + 1

   real, parameter :: &
      big  = 1.e+32,  &
      half = 0.5

!  Variables:

   integer :: &
      choice, i, ios, j, n, nnew, nread, nsplit
   real :: &
      angle, cutoff, growth, p, px, py, pz, q, qx, qy, qz, scale, shift, temp, &
      total, xp, yp, zp, xmax, xmin, xmid, ymax, ymin, ymid, zmax, zmin, zmid
   real, allocatable, dimension (:) :: &
      arc, x, y, z, xlast, ylast, zlast, xnew, ynew, znew, xorig, yorig, zorig
   logical :: &
      cr, datalost, eof, first, header, undo
   character :: &
      coord*1, dataset*48, menu (-2:mxmenu + 1)*35, method*1, test*2, title*80

!  Procedures:

   logical :: &
      alpha
   external &
      alpha  ! Distinguishes numeric text from alphanumeric

!  Data:

   data menu / &
      '  -2: Start over',                    &
      '  -1: Undo last transformation',      &
      '   0: Review data',                   &
      '   1: Translate X',                   &
      '   2: Scale X',                       &
      '   3: Translate Y',                   &
      '   4: Scale Y',                       &
      '   5: Translate Z',                   &
      '   6: Scale Z',                       &
      '   7: Reflect Z about the XY-plane',  &
      '   8: Reflect X about the YZ-plane',  &
      '   9: Reflect Y about the ZX-plane',  &
      '  10: Reverse the order (1:N)',       &
      '  11: Switch X & Y',                  &
      '  12: Switch Y & Z',                  &
      '  13: Switch Z & X',                  &
      '  14: Scale X & Y & Z the same way',  &
      '  15: Rotate (X,Y) about (Xc,Yc)',    &
      '  16: Rotate (Y,Z) about (Yc,Zc)',    &
      '  17: Rotate (Z,X) about (Zc,Xc)',    &
      '  18: Rotate (X,Y,Z) about line PQ',  &
      '  19: Display data range and center', &
      '  20: Split: X|Y|Z <|<=|==|>|>= cut', &
      '  21: RIGID_TRANSFORM (new end pts)', &
      '  22: NULINE3D morph  (new end pts)', &
      '  23: CHANGEN: new N; same rel. sp.', &
      '  24: Done',                          &
      '                                   '/ ! Last ' ' eases display
                                             ! of menu as two columns.

!  Execution:

   write (luncrt, 1001) ' ', &
      ' RESHAPE3D applies simple transformations to a 3D dataset.', &
      ' Input and output files are in Title/N/3-columns format.',   &
      ' The title and point count may be omitted, however.',        &
      ' '

!  :::::::::::::
!  Get the data:
!  :::::::::::::

   call opener (luncrt, 'Enter dataset name: ', &
      lunkbd, dataset, lunin, 'old')

   read (lunin, 1001) title
   header = alpha (title)  ! Strictly numeric first line means no header

   if (header) then
      read (lunin, *, err=900) nread
      n = nread
   else  ! Count the number of dataset lines
      rewind (lunin)
      n = -1;   ios =  0
      do while (ios == 0)
         n = n + 1
         read (lunin, *, iostat=ios)
      end do
      nread = n
      write  (luncrt, '(/, a, i8)') '# data points found:', n
      rewind (lunin)
   end if

   undo = n < nundo  ! n > 375000 seems to cause a lot of paging

   allocate (x (n), y (n), z (n))
   n = 0
   do i = 1, nread
      read (lunin, *, err=901, end=160) x (i), y (i), z (i)
      n = i
   end do

160 close (lunin)
   if (n /= nread) go to 902

   if (undo) then
      allocate (xorig (n), yorig (n), zorig (n), &
                xlast (n), ylast (n), zlast (n))
      do i = 1, n
         xorig (i) = x (i);  xlast = x (i)
         yorig (i) = y (i);  ylast = y (i)
         zorig (i) = z (i);  zlast = z (i)
      end do
   end if

   first    = .true.
   datalost = .false.

!  :::::::::::::::::::::::::::::::::::::::::::::::::::
!  Loop over possibly several transformations per run:
!  :::::::::::::::::::::::::::::::::::::::::::::::::::

200 continue

      if (first) then
         first = .false.
         write (luncrt, 1004) &                     ! Basically, i=1:mxmenu/2
            (menu (i), menu(i + halfm), i = -2, halfm - 3)
         if (mod (mxmenu, 2) == 0) write (luncrt, 1001)
      end if

210   call readi (luncrt, 'Pick one. EOF (^Z or ^D) means no more. ', &
                  lunkbd, choice, cr, eof)
      if (eof) go to 800
      if (cr) go to 210

      if (undo) then
         if (choice > 0 .or. choice == -2) then   ! Save current values:
            do i = 1, n
               xlast (i) = x (i)
               ylast (i) = y (i)
               zlast (i) = z (i)
            end do
         end if
      end if

      select case (choice)

         case (-2)  ! "Start over" from scratch:

            if (datalost) then
               write (luncrt, '(a)') &
                  ' Starting over disabled: some points were lost to splitting.'
            else
               if (undo) then
                  do i = 1, n
                     x (i) = xorig (i)
                     y (i) = yorig (i)
                     z (i) = zorig (i)
                  end do
               else
                  write (luncrt, '(a)') &
                     ' Starting over has been disabled: too many points.'
               end if
            end if

         case (-1)  ! "Undo" previous operation:

            if (datalost) then
               write (luncrt, '(a)') &
                  ' Undo disabled: some points were lost to splitting.'
            else
               if (undo) then
                  do i = 1, n
                     x (i) = xlast (i)
                     y (i) = ylast (i)
                     z (i) = zlast (i)
                  end do
               else
                  write (luncrt, '(a)') &
                     ' Undo has been disabled: too many points.'
               end if
            end if

         case (0)   ! "Review": Display the data.

            write (luncrt, '(/, (1x, i6, 3es24.15))') &
               (i, x(i), y(i), z(i), i = 1, min (n, ndisplay))

            if (n > ndisplay) then
               write (luncrt, '(a)') ' ::::::::::::::::::::::'
               write (luncrt, '(1x, i6, 3es24.15)') &
                  (i, x(i), y(i), z(i), i = max (ndisplay + 1, n - ndisplay), n)
            end if

            first = .true.  ! Allow the menu to be displayed again

         case (1)   ! "Translate X":

            call readr (luncrt, '   Enter X shift (+ or -): ', &
                        lunkbd, shift, cr, eof)
            if (cr .or. eof) go to 210

            do i = 1, n
               x (i) = x (i) + shift
            end do

         case (2)   ! "Scale X":

            call readr (luncrt, '   Enter X scale (+ or -): ', &
                        lunkbd, scale, cr, eof)
            if (cr .or. eof) go to 210

            do i = 1, n
               x (i) = x (i) * scale
            end do

         case (3)   ! "Translate Y":

            call readr (luncrt, '   Enter Y shift (+ or -): ', &
                        lunkbd, shift, cr, eof)
            if (cr .or. eof) go to 210

            do i = 1, n
               y (i) = y (i) + shift
            end do

         case (4)   ! "Scale Y":

            call readr (luncrt, '   Enter Y scale (+ or -): ', &
                        lunkbd, scale, cr, eof)
            if (cr .or. eof) go to 210

            do i = 1, n
               y (i) = y (i) * scale
            end do

         case (5)   ! "Translate Z"

            call readr (luncrt, '   Enter Z shift (+ or -): ', &
                        lunkbd, shift, cr, eof)
            if (cr .or. eof) go to 210

            do i = 1, n
               z (i) = z (i) + shift
            end do

         case (6)   ! "Scale Z"

            call readr (luncrt, '   Enter Z scale (+ or -): ', &
                        lunkbd, scale, cr, eof)
            if (cr .or. eof) GO TO 210

            do i = 1, n
               z (i) = z (i) * scale
            end do

         case (7)   ! "Reflect Z about the XY-plane":

            do i = 1, n
               z (i) = -z (i)
            end do

         case (8)   ! "Reflect X about the YZ-plane":

            do i = 1, n
               x (i) = -x (i)
            end do

         case (9)   ! "Reflect Y about the ZX-plane":

            do i = 1, n
               y (i) = -y (i)
            end do

         case (10)  ! "Reverse the order":

            do i = 1, (n + 1) / 2
               j = n + 1 - i
               temp  = x (i)
               x (i) = x (j)
               x (j) = temp
               temp  = y (i)
               y (i) = y (j)
               y (j) = temp
               temp  = z (i)
               z (i) = z (j)
               z (j) = temp
            end do

         case (11)  ! "Switch X and Y":

            do i = 1, n
               temp  = x (i)
               x (i) = y (i)
               y (i) = temp
            end do

         case (12)  ! "Switch Y and Z":

            do i = 1, n
               temp  = y (i)
               y (i) = z (i)
               z (i) = temp
            end do

         case (13)  ! "Switch Z and X":

            do i = 1, n
               temp  = x (i)
               x (i) = z (i)
               z (i) = temp
            end do
        
         case (14)  ! "Scale X & Y & Z"

            call readr (luncrt, '   Enter scale (+ or -): ', &
                        lunkbd, scale, cr, eof)
            if (cr .or. eof) go to 210

            do i = 1, n
               x (i) = x (i) * scale
               y (i) = y (i) * scale
               z (i) = z (i) * scale
            end do

         case (15)  ! "Rotate (X,Y)":

            call readr (luncrt, '   Degrees anticlockwise: ', &
                        lunkbd, angle, cr, eof)
            if (cr .or. eof) go to 210

            write (luncrt, 1001, advance='no') '  Center (Xc, Yc): '
            read  (lunkbd, *) p, q
            write (luncrt, 1001)

            call rotate2d (n, x, y, angle, p, q)

         case (16)  ! "Rotate (Y,Z)":

            call readr (luncrt, '   Degrees anticlockwise: ', &
                        lunkbd, angle, cr, eof)
            if (cr .or. eof) go to 210

            write (luncrt, 1001, advance='no') '  Center (Yc, Zc): '
            read  (lunkbd, *) p, q
            write (luncrt, 1001)

            call rotate2d (n, y, z, angle, p, q)

         case (17)  ! "Rotate (Z,X)":

            call readr (luncrt, '   Degrees anticlockwise: ', &
                        lunkbd, angle, cr, eof)
            if (cr .or. eof) go to 210

            write (luncrt, 1001, advance='no') '  Center (Zc, Xc): '
            read  (lunkbd, *) p, q
            write (luncrt, 1001)

            call rotate2d (n, z, x, angle, p, q)

         case (18)  ! "Rotate (X,Y,Z) about line joining P and Q":

            call readr (luncrt, '   Degrees (RH rule; thumb P -> Q): ', &
                        lunkbd, angle, cr, eof)
            if (cr .or. eof) go to 210

            write (luncrt, 1001, advance='no') '  (Px, Py, Pz): '
            read  (lunkbd, *) px, py, pz
            write (luncrt, 1001, advance='no') '  (Qx, Qy, Qz): '
            read  (lunkbd, *) qx, qy, qz
            write (luncrt, 1001)

            call rotate3d (n, x, y, z, angle, px, py, pz, qx, qy, qz)

         case (19)  ! "Display data range and center":

            xmax = -big;  xmin = big
            ymax = -big;  ymin = big
            zmax = -big;  zmin = big

            do i = 1, n
               xmax = max (x (i), xmax)
               xmin = min (x (i), xmin)
               ymax = max (y (i), ymax)
               ymin = min (y (i), ymin)
               zmax = max (z (i), zmax)
               zmin = min (z (i), zmin)
            end do

            xmid = (xmin + xmax) * half
            ymid = (ymin + ymax) * half
            zmid = (zmin + zmax) * half

            write (luncrt, '(/, 3(3x, a, 2es24.15, es26.15, /))') &
               'Xmin, Xmax, Xcenter:', xmin, xmax, xmid, &
               'Ymin, Ymax, Ycenter:', ymin, ymax, ymid, &
               'Zmin, Zmax, Zcenter:', zmin, zmax, zmid

         case (20)  ! "Split the dataset"

            call readc (luncrt, '   Split in which direction? [X|Y|Z]: ', &
                        lunkbd, coord, cr, eof)
            if (cr .or. eof) go to 210

            call readc (luncrt, '   Comparison test? [<|<=|==|>|>=]:   ', &
                        lunkbd, test, cr, eof)
            if (cr .or. eof) go to 210

            call readr (luncrt, '   Value to compare with:             ', &
                        lunkbd, cutoff, cr, eof)
            if (cr .or. eof) go to 210

            call splitxyz (coord, test, cutoff, n, x, y, z, nsplit)

            write (luncrt, '(a, i10)') &
               '    New # points:             ', nsplit

            if (.not. datalost) datalost = nsplit < n
            n = nsplit

         case (21, 22)  ! Transform the curve as defined by new end points

            allocate (xnew(n), ynew(n), znew(n))
            write (luncrt, 1001, advance='no') '  New x/y/z(1): '
            read  (lunkbd, *) xnew(1), ynew(1), znew(1)
            write (luncrt, 1001, advance='no') '  New x/y/z(n): '
            read  (lunkbd, *) xnew(n), ynew(n), znew(n)
            write (luncrt, 1001)

            if (choice == 21) then
               call rigid_transform (n, x, y, z, xnew, ynew, znew)
            else  ! choice = 22
               call nuline3d (1, n, x, y, z, xnew, ynew, znew)
            end if

            x(1:n) = xnew(1:n)
            y(1:n) = ynew(1:n)
            z(1:n) = znew(1:n)

            deallocate (xnew, ynew, znew)

         case (23)  ! Change the # points; preserve the relative spacing

            write (luncrt, '(a, i7)') '   Current number of points: ', n
            call readi (luncrt,        '  Desired number of points: ', &
                        lunkbd, nnew, cr, eof)
            if (cr .or. eof) go to 210

            write (luncrt, '(a)') '   Interpolation methods:', &
                                  '   M (monotonic), B (loose), L (linear);', &
                                  '   Uppercase => preserve relative spacing.',&
                                  '   Lowercase m|b|l => ~uniform spacing.'
            call reads (luncrt,    '  Interpolation choice: ', &
                        lunkbd, method, cr, eof)
            if (cr .or. eof) go to 210

!           Save the current growth rates:

            open (lunout, file='growth-rates-before.dat', status='unknown')

            allocate (arc(n))

            call chords3d (n, x, y, z, .false., total, arc)

            do i = 2, n-1
               growth = (arc(i+1) - arc(i)) / (arc(i) - arc(i-1))
               write (lunout, '(2es16.8)') arc(i), growth
            end do
            close (lunout)
            deallocate (arc)

            allocate (xnew(nnew), ynew(nnew), znew(nnew))

            call changen (1, n, x, y, z, 1, nnew, xnew, ynew, znew, method)

            deallocate (x, y, z);  allocate (x(nnew), y(nnew), z(nnew))

            n    = nnew
            x(:) = xnew(:)
            y(:) = ynew(:)
            z(:) = znew(:)

            deallocate (xnew, ynew, znew)

            open (lunout, file='growth-rates-after.dat', status='unknown')

            allocate (arc(n))

            call chords3d (n, x, y, z, .false., total, arc)

            do i = 2, n-1
               growth = (arc(i+1) - arc(i)) / (arc(i) - arc(i-1))
               write (lunout, '(2es16.8)') arc(i), growth
            end do
            close (lunout)
            deallocate (arc)

         case (24)  ! Done (may work better than ^D)

            go to 800

         case default

            go to 210

      end select

   go to 200  ! Next selection

!  :::::::::::::
!  Save results:
!  :::::::::::::

800  write (luncrt, 1001)
     call opener (luncrt, 'Output file name?  EOF = quit: ', &
        lunkbd, dataset, lunout, 'UNKNOWN')

     if (header) then
        call reads (luncrt, 'Output title line? <CR> = same: ', &
                    lunkbd, title, cr, eof)

        i = len_trim (title)
        write (lunout, 1001, err=903) title(1:i)
        write (lunout, 1003, err=903) n
     end if

     write (lunout, '(3es24.15)', err=903) (x(i), y(i), z(i), i = 1, n)

     deallocate (x, y, z)
     if (undo) deallocate (xorig, yorig, zorig, xlast, ylast, zlast)

     go to 999

!    Error handling:

900  write (luncrt, 1001) ' Error reading number of points.'
     go to 999
901  write (luncrt, 1001) ' Error reading (X,Y,Z) triple.'
     go to 999
902  write (luncrt, 1001) ' N on line 2 does not match # points found.'
     write (luncrt, '(a, 2i10)') ' N (file) and N (read): ', nread, n
     go to 999
903  write (luncrt, 1001) ' Error saving results.'
     go to 999

999  continue ! Avoid STOP machine dependencies

!    Formats:

1001 format (a)
1002 format (a, i3)
1003 format (i6)
1004 format (/, (a, 10x, a))

   end program reshape3d
