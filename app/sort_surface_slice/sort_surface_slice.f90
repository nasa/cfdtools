!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program sort_surface_slice
!
!  This utility reads a dataset resulting from Tecplot's option to slice a
!  surface dataset defined by a plane.  The dataset has been saved as a single
!  zone of nodes and 2-point line segment pointers.  Plotting the points as
!  symbols may not be appropriate: a continuous curve is commonly desired, and
!  sorting by x, y, or z may be all that is needed.  Plotting against arc length
!  along the curve may also be desirable, so the option to add cumulative arc
!  length as an additional variable is also provided.
!
!  Originally, the program just treated the common case where sorting the nodes
!  in x, y, or z suffices to produce a single meaningful continuous curve.
!  This version retains that option but can also allow for more than one
!  contiguous curve in the slice, as in the case of a cut through a fuselage
!  that also encounters the outboard portion of a swept wing.  Each such segment
!  appears as a distinct zone in the output file and may have arc lengths added.
!  This option makes use of the 2-point line segment pointers accompanying the
!  slice coordinates.
!
!  Further, a specialized option is provided for the common requirement of
!  plotting surface quantities versus run length from the stagnation point of
!  a hypersonic flow solution (capsule centerline, or possibly near the wing
!  leading edge of a wing section slice).  If a single contiguous curve is
!  found, and a flow variable name starts with p or P, and insertion of arc
!  lengths has been requested, then the peak pressure point along the curve is
!  taken to be the stag. pt., and arc lengths are adjusted to increase in both
!  directions away from it.  This form is appended as a second output zone,
!  preserving the initial ordered slice zone in case the stag. pt. option is
!  not intended, since no extra prompts have been added.
!
!  Sample input dataset (POINT order, single zone, one variable name per line):
!
!     TITLE     = "A"
!     VARIABLES = "x [m]"
!     "y [m]"
!     "z [m]"
!     "pw [Pa]"
!     "qw [W/m^2]"
!     ZONE T="Slc: Z=0"
!      STRANDID=0, SOLUTIONTIME=0
!      Nodes=118, Elements=109, ZONETYPE=FELineSeg
!      DATAPACKING=POINT
!      DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )
!      1.584258E-02 3.505522E-01 0.000000E+00 1.358903E+03 1.643338E+05
!      1.646541E-02 3.586402E-01 1.969463E-19 1.358848E+03 1.643071E+05
!       :            :            :            :            :
!      1.250000E+01 1.000000E+01 4.147692E-18 1.009208E+03 4.999492E+04
!      2 1
!      2 3
!      3 4
!      6 5
!      : :
!      117 118
!
!  Corresponding output (single segment case with arc lengths inserted):
!
!     TITLE     = "A"
!     VARIABLES = "x [m]"
!     "y [m]"
!     "z [m]"
!     "pw [Pa]"
!     "qw [W/m^2]"
!     "s [m]"
!     ZONE T="Slc: Z=0"
!      STRANDID=0, SOLUTIONTIME=0
!      I=118, J=1, K=1, ZONETYPE=Ordered
!      DATAPACKING=POINT
!      DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )
!      1.58425E-02 3.50552E-01 0.00000E+00 1.35890E+03 1.64333E+05 0.00000E+00
!      1.64654E-02 3.58640E-01 1.96946E-19 1.35884E+03 1.64307E+05 1.23456E-02
!       :           :           :           :           :           :
!
!     [ZONE T="Slc: Z=0, Stag. Pt. Run Length"
!      STRANDID=0, SOLUTIONTIME=0
!      I=118, J=1, K=1, ZONETYPE=Ordered
!      DATAPACKING=POINT
!      DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )
!      1.64654E-02 3.58640E-01 1.96946E-19 1.35884E+03 1.64307E+04 1.23456E+00
!      1.58425E-02 3.50552E-01 0.00000E+00 1.35890E+03 1.64333E+04 1.23400E+00
!       :           :           :           :           :           :
!      1.58425E-02 3.50552E-01 0.00000E+00 7.65432E+03 1.64333E+05 0.00000E+00
!       :           :           :           :           :           :
!      1.58425E-02 3.50552E-01 0.00000E+00 1.35890E+03 1.64333E+04 1.23400E+00
!      1.64654E-02 3.58640E-01 1.96946E-19 1.35884E+03 1.64307E+04 1.23456E+00]
!
!  I/O Strategy:
!
!     Keep the parsing that would elsewhere be de rigeur to a minimum.
!     This is a file written consistently by Tecplot, after all.
!     Scan the header lines as strings until a purely numeric string is found.
!     Rewind and reread the header lines as an array.  Find the line containing
!     NODES (after upcasing) and adjust as shown in the sample if simple sorting
!     is specified.  If arc length is to be appended as an extra column, insert
!     its name appropriately during writing of the output header lines.  Read
!     the indicated number of points, sort them as specified, and write them in
!     the desired order, possibly appending an arc length to each line.  If the
!     2-point line segments are to be pieced together as one or more contiguous
!     curves, a more awkward output strategy is required, but existing utility
!     CUTORDER does the hard part after extension to one or more functions.
!
!  History:
!
!     07/30/10  D. A. Saunders  Initial implementation of useful functionality
!                               prompted by John Theisinger's axisymmetric/3-D
!                               aero-heating correlation study.  Sorting should
!                               suffice for simple geometries, but the 2-point
!                               line segments should be pieced together properly
!                               for general cases where multiple curve segments
!                               should be written as separate zones.
!     08/13/10   "       "      Completed it with the option to process the
!                               2-point line segment information properly.
!     08/19/10   "       "      Data lines from BLAYER exceed 512 characters,
!                               so the apparent number of variables was wrong.
!     08/20/10   "       "      Plotting versus run length from stagnation point
!                               is desirable enough that it has been added as a
!                               somewhat specialized option for cases that find
!                               one contiguous curve and have a flow variable
!                               starting with p or P (since peak pressure is
!                               taken to be at the stag. pt.).
!     08/23/10   "       "      A slice that happens to be a closed curve and
!                               can start/end anywhere needs special treatment
!                               for the stag. pt. run-length option.  We assume
!                               that the stag. pt. is on y = 0, locate the
!                               aft-most point with y closest to 0, and shift
!                               the data circularly to make that point index 1,
!                               then proceed as for an open single curve.  Also,
!                               blunt bodies don't suit positive arc lengths on
!                               the lee side, so make those negative, and allow
!                               for the common case of an upside-down capsule
!                               by making the shorter of the two run lengths use
!                               negative arc lengths if aspect ratio < 1.5.
!     09/28/10   "       "      Slices from layouts with text boxes turn out to
!                               include the text boxes before the zone info., so
!                               inserting the arc length variable name had to be
!                               done more carefully.
!     10/06/10   "       "      Last_var_line wasn't being set if a pressure
!                               variable was not present.
!  Author:
!
!     David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Local constants:

   integer, parameter :: &
      lunin     = 1,     &  ! Input dataset
      lunout    = 2,     &  ! Output dataset
      lunkbd    = 5,     &  ! Keyboard inputs
      luncrt    = 6,     &  ! Screen outputs
      maxlen    = 999,   &  ! Limit on length of an input dataset line
      naddim    = 100       ! Max. # contiguous curves in the slice

   real, parameter :: &
      eps       = 1.e-12    ! Applied to data range -> point matching tolerance

   character, parameter :: &
      blank*1   = ' '

!  Local variables:

   integer :: &
      i, i1, i2, ier, index_pressure, ios, isortcol, l, last_var_line, line,   &
      nf, nheader, nodes, nodes_line, npts, nsegments, nvars, nzones,          &
      vars_line, zone_title_line

   integer :: &
      nad(naddim)                          ! nad(i) points to curve i of a slice
                                           ! within x/y/z/fslice(:)
   integer, allocatable :: &
      list(:), ndgcut(:,:)

   real :: &
      bbox(6), epsxyz, rnumber, total

   real, allocatable :: &
      arcs(:), columns(:,:), sortcol(:), x(:), y(:), z(:), f(:,:)

   logical :: &
      arc_lengths, run_lengths, sort_by_coordinate

   character :: &
      answer*1, buffer*(maxlen), filename*128, iformat*4

   character (len=maxlen), allocatable :: &
      header(:)

!  Procedures:

   logical  :: alpha
   external :: alpha         ! Identifies a non-numeric string
   external :: chords3d      ! Calculates cumulative arc lengths
   external :: decode_number ! Internal read from string containing integer/real
   external :: token_count   ! Counts the # items on a line
   external :: upcase        ! Uppercase utility
   external :: hsortri       ! Heap sort utility
   external :: cutorder_nf   ! Identifies contiguous curves among 2-pt. segments
   external :: circ_shift    ! Circular shift of a vector in-place

!  Execution:

   write (luncrt, '(/, a)', advance='no') 'Input Tecplot slice dataset: '
   read  (lunkbd, *) filename
   open  (lunin, file=filename, status='old')

   write (luncrt, '(a)', advance='no') 'Output file name: '
   read  (lunkbd, *) filename
   open  (lunout, file=filename, status='unknown')

   write (luncrt, '(a)', advance='no') &
      'Sort by x, y, or z? [y|n; no => allow for multiple curves]: '
   read  (lunkbd, *) answer

   sort_by_coordinate = answer == 'y' .or. answer == 'Y'

   if (sort_by_coordinate) then
      write (luncrt, '(a)', advance='no') 'Column number to sort on [1|2|3]: '
      read  (lunkbd, *) isortcol
   end if

   write (luncrt, '(a)', advance='no') &
      'Add arc lengths to output curve(s)? [y|n]: '
   read  (lunkbd, *) answer

   arc_lengths = answer == 'y' .or. answer == 'Y'

!  Determine the number of header lines for proper storage allocation.
!  Track the last variable name in case a text box is present ahead of the
!  first zone (because arc length may need to be inserted as an extra variable).

   nheader = 0;  vars_line = 0;  last_var_line = 0;  index_pressure = 0

   do  ! Until a purely numeric data line is found
      read (lunin, '(a)') buffer
      l = len_trim (buffer)
      if (.not. alpha (buffer(1:l))) exit

      nheader = nheader + 1

      if (vars_line == 0) then
         if (index (buffer(1:l), 'VARIABLES') > 0) vars_line = nheader
      else if (index_pressure == 0) then
         if (buffer(1:2) == '"p' .or. buffer(1:2) == '"P' .or. &  ! DPLR-ism
             buffer(1:4) == '"\"p' ) then
            index_pressure = nheader - vars_line - 2       ! Excludes x, y, z
         end if
      end if

      if (last_var_line == 0) then
         if (nheader > 2) then
            if (buffer(1:1) /= '"') last_var_line = nheader - 1
         end if
      end if
   end do

   run_lengths = arc_lengths .and. index_pressure > 0

!  Counting the tokens on the first numeric line gives the input # variables:

   call token_count (buffer(1:l), blank, nvars)

   nf = nvars - 3  ! >= 0

   write (luncrt, '(a, i7)') '# variables found:     ', nvars

!  Reread and store the header lines:

   rewind (lunin);  allocate (header(nheader))

   do line = 1, nheader
      read (lunin, '(a)') header(line)
      l = len_trim (header(line))
      if (index (header(line)(1:l), 'ZONE T') > 0) zone_title_line = line
   end do

   if (arc_lengths) header(nheader)(l:l+7) = 'DOUBLE )'

!  Transcribe the variable names, adding arc length if specified:

   write (lunout, '(a)') (trim (header(line)), line = 1, last_var_line)
   if (arc_lengths) write (lunout, '(a)') '"s, m"'
   if (zone_title_line > last_var_line + 1) then
      write (lunout, '(a)') (trim (header(line)), line = last_var_line + 1,    &
                             zone_title_line - 1)
   end if

!  Scan the remaining header lines for # points and # line segments, etc.:

   do line = zone_title_line + 1, nheader

      l = len_trim (header(line))
      call upcase  (header(line)(1:l))
      i1 =  index  (header(line)(1:l), 'NODES')

      if (i1 > 0) then

         nodes_line = line
         i1 = i1 + 6                              ! Start of node count
         i2 = index (header(line)(1:l), ',') - 1  ! End of node count

         call decode_number (header(line)(i1:i2), nodes, rnumber, ios) ! -> 2

         write (luncrt, '(a, i7)') '# data points found:   ', nodes

         if (sort_by_coordinate) then  ! Single output zone
            header(line)(1:3)       = ' I='
            header(line)(4:4+i2-i1) = header(line)(i1:i2)
            header(line)(5+i2-i1:)  = ', J=1, K=1, ZONETYPE=Ordered'
         else  ! We'll be reading the 2-point line segment pointers
            i1 = index (header(line)(1:l), 'ELEMENTS') + 9
            i  = i2 + 1              ! Kludge to allow indexing from 1
            header(line)(i:i) = ' '  ! by suppressing the first comma
            i2 = index (header(line)(1:l), ',') - 1
            call decode_number (header(line)(i1:i2), nsegments, rnumber, ios)
            write (luncrt, '(a, i7)') '# 2-pt. segments found:', nsegments
         end if

         exit

      end if

   end do

   allocate (columns(nvars,nodes))

   read (lunin, *, iostat=ios) columns

   if (ios /= 0) then
      write (luncrt, '(a, 2i7)') 'Trouble reading slice data. nvars, nodes: ', &
         nvars, nodes
      go to 99
   end if

   bbox(1) = minval (columns(1,:))
   bbox(2) = maxval (columns(1,:))
   bbox(3) = minval (columns(2,:))
   bbox(4) = maxval (columns(2,:))
   bbox(5) = minval (columns(3,:))
   bbox(6) = maxval (columns(3,:))

   epsxyz  = max (bbox(2) - bbox(1), &
                  bbox(4) - bbox(3), &
                  bbox(6) - bbox(5)) * eps

   if (arc_lengths .or. .not. sort_by_coordinate) then
      allocate (x(nodes), y(nodes), z(nodes))
      if (nf > 0) allocate (f(nf,nodes))
   end if

   if (arc_lengths) allocate (arcs(nodes))

   if (sort_by_coordinate) then

      close (lunin)  ! Don't need the pointers

      write (lunout, '(a)') &
         (trim (header(line)), line = zone_title_line, nheader)

      allocate (sortcol(nodes), list(nodes))

!     The heap sort utility reorders the column being sorted (inconvenient), so
!     make a copy of it for sorting:

      do i = 1, nodes
         list(i) = i
         sortcol(i) = columns(isortcol,i)
      end do

      call hsortri (sortcol, list, nodes)

      if (.not. arc_lengths) then

         do i = 1, nodes
            l = list(i)
            write (lunout, '(1p, 20e17.9)') columns(:,l)
         end do

      else

         do i = 1, nodes
            l = list(i)
            x(i) = columns(1,l)
            y(i) = columns(2,l)
            z(i) = columns(3,l)
            if (nf > 0) f(:,i) = columns(4:,l)
         end do

         call chords3d (nodes, x, y, z, .false., total, arcs)

         do i = 1, nodes
            l = list(i)
            write (lunout, '(1p, 20e17.9)') columns(:,l), arcs(i)
         end do

         npts = nodes  ! For stag. pt. run-length option

      end if

   else  ! Allow for more than one contiguous curve, one per output zone

      allocate (ndgcut(2,nsegments))

      read (lunin, *, iostat=ios) ndgcut

      if (ios /= 0) then
         write (luncrt, '(a)') 'Trouble reading the 2-point line pointers.'
         go to 99
      end if

      close (lunin)

!     Organize the 2-point line segments as contiguous curve(s):

      nad(1) = 0  ! # lines for this slice

      call cutorder_nf (nf, nodes, nsegments, nodes, naddim, nodes, nsegments, &
                        ndgcut, columns, epsxyz, x, y, z, f, nad, ier)

      if (ier /= 0) then
         write (luncrt, '(a, i3, a, i7)')   &
            'CUTORDER_NF error:', ier, '  # 2-pt. segments:', nsegments
         write (luncrt, '(a, i7)') '# sub-curves found so far:', nad(1)
      end if

!     Write each contiguous curve as a new zone:

      nzones = nad(1)

      nad(1) = 1  ! So the loop works for curve 1

      do line = 1, nzones
         i1   = nad(line)
         i2   = nad(line+1) - 1
         npts = i2 -i1 + 1

         write (lunout, '(a)') &
            (trim (header(l)), l = zone_title_line, nodes_line - 1)

         header(nodes_line)(1:3) = ' I='
         call ndigits (npts, l)
         iformat = '(i?)'
         write (iformat(3:3), '(i1)') l
         i = 4 + l - 1
         write (header(nodes_line)(4:i), iformat) npts
         header(nodes_line)(i+1:) = ', J=1, K=1, ZONETYPE=Ordered'
         write (lunout, '(a)') (trim (header(l)), l = nodes_line, nheader)

         if (.not. arc_lengths) then

            if (nf == 0) then
               write (lunout, '(1p, 3e17.9)') (x(i), y(i), z(i), i = i1, i2)
            else
               do i = i1, i2
                  write (lunout, '(1p, 20e17.9)') x(i), y(i), z(i), f(:,i)
               end do
            end if

         else

            call chords3d (npts, x(i1), y(i1), z(i1), .false., total, arcs(i1))

            if (nf == 0) then
               write (lunout, '(1p, 4e17.9)') &
                  (x(i), y(i), z(i), arcs(i), i = i1, i2)
            else
               do i = i1, i2
                  write (lunout, '(1p, 20e17.9)') &
                     x(i), y(i), z(i), f(:,i), arcs(i)
               end do
            end if

         end if

      end do  ! Next contiguous curve

      run_lengths = run_lengths .and. nzones == 1

   end if  ! Choice of sorting method

   if (run_lengths) call adjust_to_stag_pt ()  ! Local procedure below

99 continue

!  Local procedure for program sort_surface_slice:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine adjust_to_stag_pt ()

!     Adjust arc lengths to increase in both directions from the stag. pt. as
!     defined by the maximum pressure point.  Write a second zone accordingly.
!     Actually, this suffices only for a slender body (possibly a wing section).
!     Blunt forebodies and capsules with aft bodies that lead to closed curves
!     require more heuristics, all isolated here as part of a specialized option
!     whereas the rest of the program functionality is quite general-purpose.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local constants:

      real, parameter :: ar_limit = 1.5  ! Heuristic for bluntness

!     Local variables:

      integer :: i_peak, i_start, j, max_loc(1), n
      real    :: aspect_ratio, data_range(3), dsq, dsq_least, s_peak
      logical :: blunt, closed

      real, allocatable :: ftemp(:)

!     Execution:

!     For determining bluntness or high aspect ratio, we cannot avoid assuming
!     that X is streamwise.

      data_range(1) = bbox(2) - bbox(1)
      data_range(2) = bbox(4) - bbox(3)
      data_range(3) = bbox(6) - bbox(5)

      aspect_ratio  = data_range(1) / max (data_range(2), data_range(3))

      blunt  = aspect_ratio < ar_limit

!     If the Tecplot cut leads to a closed curve, index 1 could be anywhere.
!     We need to make the furthest aft point index 1, but there may be multiple
!     points with x = xmax, so pick the one nearest to (xmax, 0, 0):

      closed = (x(npts) - x(1))**2 + (y(npts) - y(1))**2 + &
               (z(npts) - z(1))**2 < epsxyz**2

      if (closed) then

         i_start = 1;  dsq_least = 1.e+10
         do i = 1, npts
            if (bbox(2) - x(i) < epsxyz) then
               dsq = y(i)**2 + z(i)**2
               if (dsq < dsq_least) then
                  dsq_least = dsq
                  i_start = i
               end if
            end if
         end do

         write (luncrt, '(a, i5, a, 1p, 2e15.6)') &
            'Adjusting closed-curve start index from 1 to', i_start, &
            '; x adjustment:', x(1), x(i_start)

         call circ_shift (npts, i_start, 1, x)
         call circ_shift (npts, i_start, 1, y)
         call circ_shift (npts, i_start, 1, z)

!        Now we have the original first and last point duplicated at indices
!        j and j + 1 where j = npts - i_start + 1.  Shift the data to remove
!        the duplicate point, then copy the new point 1 to the new end point:

         j = npts - i_start + 1

         x(j:npts-1) = x(j+1:npts);  x(npts) = x(1)
         y(j:npts-1) = y(j+1:npts);  y(npts) = y(1)
         z(j:npts-1) = z(j+1:npts);  z(npts) = z(1)

         call chords3d (npts, x, y, z, .false., total, arcs)  ! Recompute

         if (nf > 0) then
            allocate  (ftemp(npts))

            do n = 1, nf
               ftemp(:) = f(n,1:npts)
               call circ_shift (npts, i_start, 1, ftemp)
               f(n,1:npts)   = ftemp(:)
               f(n,j:npts-1) = f(n,j+1:npts);  f(n,npts) = f(n,1)
            end do

            deallocate (ftemp)
         end if

      end if

      max_loc = maxloc (f(index_pressure,1:npts))
      i_peak  = max_loc(1)
      s_peak  = arcs(i_peak)

!     Make arc length zero at the pressure peak:

      write (luncrt, '(a, 2i6, 2f16.5)') &
         'Pressure variable, and peak index, pressure & arc:', &
          index_pressure, i_peak, f(index_pressure,i_peak), s_peak

      arcs(1:i_peak)  = s_peak - arcs(1:i_peak)
      arcs(i_peak+1:npts) = arcs(i_peak+1:npts) - s_peak

      if (blunt) then  ! Make the wind-side (shorter) run lengths negative

         if (arcs(i_peak) < arcs(npts)) then
            arcs(1:i_peak)  = -arcs(1:i_peak)
         else
            arcs(i_peak+1:npts) = -arcs(i_peak+1:npts)
         end if

      end if

      l = len_trim (header(zone_title_line))
      header(zone_title_line)(l:) = ', Stag. Pt. Run Length"'

      write (lunout, '(a)') &
         (trim (header(line)), line = zone_title_line, nheader)

      do i = 1, npts
         write (lunout, '(1p, 20e17.9)') x(i), y(i), z(i), f(:,i), arcs(i)
      end do

      write (luncrt, '(a)') &
         'Use zone 2 for arc lengths starting from the stagnation point.'

      end subroutine adjust_to_stag_pt

   end program sort_surface_slice
