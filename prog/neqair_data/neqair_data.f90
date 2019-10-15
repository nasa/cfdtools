!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program neqair_data
!
!  Version Description:
!
!     This version of NEQAIR_DATA can write line-of-sight (LOS) data in both the
!  original rigid format read by NEQAIR versions up to and including v14 and the
!  simpler column-oriented format handled from NEQAIR v15 on.  Both data formats
!  involve a sample.LOS.dat file read here that indicates the species names.  If
!  a sample file is found with a NEQAIR version number on line 1, then the later
!  data format is indicated. The sample file is transcribed as the header of all
!  output LOS data files in either case.
!
!  Original Preamble:
!
!     Earlier utilities LINES_OF_SIGHT[_2D] and FLOW_INTERP[_2D] typically
!  produce multiple lines of sight and flow-field data interpolated to each
!  point of each line of sight, as a PLOT3D-type multiblock grid and function
!  file pair.  Note that the input coordinates should correspond to grid
!  vertices, not cell centers, because the lines of sight are discretized using
!  the relative distribution of the nearest radial line from the vertex-
!  centered grid.
!
!     Radiation program NEQAIR treats one line of sight at a time, and expects
!  the data in a rigid format containing arc lengths (measured in the direction
!  towards the body, and in cm), total number density, T, Tr, Tv, Te, and the
!  species number densities, at each grid point along the line.  The arc lengths
!  should correspond to the centers of cells at which the flow function data
!  are taken to apply (in NEQAIR), so they are adjusted accordingly here.  The
!  wall function values thus apply to a cell half the width of the grid cell
!  at the wall.
!
!     [Later:  This version has the option, prompted by meteor studies, to
!  output line-of-sight data in the direction away from the body.]
!
!     This program retrieves the PLOT3D-format files (e.g., los.g & los.f in
!  SI units) and writes a separate file in NEQAIR's LOS.dat format for each
!  line of sight it finds, into the working directory. These output file names
!  contain line numbers, such as LOS-1.dat.  A script can then run NEQAIR on
!  some range m:n of line numbers, by copying or linking to input files and
!  writing radiation results for each line one level down in directories
!  /LINE-m through /LINE-n (say).
!
!     The current version now offers options to redistribute the line-of-sight
!  data to a [presumably smaller] different number of points in order to speed
!  the expensive NEQAIR calculations.  This was initially done based on |dT/ds|,
!  where T is vibrational temperature and s is arc length along the line.  (Use
!  of the translational temperature, which is higher than Tv at the shock, tends
!  to produce spacing that is too coarse in the boundary layer.)  See below for
!  further redistribution options.
!
!  Method:
!
!     The grid and function file may be 2D or 3D (automatically determined).
!  The input grid coordinates are assumed to be in meters, emanating from the
!  body towards the shock.  Normal usage with the DPLR flow solver has block
!  dimensions (1,nj) or (1,1,nk), but these can safely be treated here as (nj,1)
!  or (nk,1,1) for more efficient Fortran indexing.  Points outside the shock
!  need to be suppressed to avoid NEQAIR trouble with low temperatures, and the
!  line data need to be output in the shock-to-body direction.
!
!               The following paragraph is NEQAIR V14-specific:
!
!     The function file should contain T, Tr, Tv, Te (K) and the species number
!  densities (m^-3).  Note that DPLR's postprocessor does NOT extract Tr unless
!  irot = 1.  BE SURE THAT THERE ARE 4 TEMPERATURES IN THE FUNCTION FILE.
!  Common usage with 2-temperature data (ivib = 1 in DPLR) means reading and
!  writing T, T, Tv, Tv for these 4 temperatures.
!
!               NEQAIR V15 differences:
!
!     NEQAIR V15[+] handles variables numbers of temperatures.  The number of
!  species names beginning with T determines the number of temperatures that
!  should be in the input function file in the indicated order.
!
!     One line of sight's worth of data (one input block) is read at a time.
!  The header of a file named sample.LOS.dat file in the working directory is
!  transcribed to each output file, down to the line above the first line of
!  the data proper for line 1.  This sample file should contain the species
!  names in the order matching the input function file.  This program is thus
!  independent of the number and order of the species, as long as the file
!  sample.LOS.dat lists them in the header in the correct order.
!
!  Minimum Contents of File sample.LOS.dat:
!
!     V14 or earlier:
!
!        >  A line of asterisks precedes the species section.
!        >  The species names may appear in any order.
!        >  Two lines of hyphens must appear.  The last line is all hyphens.
!        >  Replace the comments in column 1 below appropriately.
!
!*******************************************************************************
!        aaaaaaaa       aaaaaaaa       aaaaaaaa       aaaaaaaa    (2x,(7x,a8))
!        CO2            CO             N2             O2       :Species Symbols.
!        NO             C              N              O
!
!-------------------------------------------------------------------------------
!  no.   x,cm   total partcc       t        tr        tv        te  (i5,f8.3,
!iiii rrrrrrr rrrrrrrrrrrrrr rrrrrrrrr rrrrrrrrr rrrrrrrrr rrrrrrrrre15.6,4f10.1
!      rrrrrrrrrrrrrr rrrrrrrrrrrrrr rrrrrrrrrrrrrr rrrrrrrrrrrrr   (6x,4e15.6)
!      Include these 9 lines (from --- to --- lines) for first grid point only!!
!      End each grid point entry with a blank line.
!      End data file with a line of zero's as shown on the next line.
!   0     0.0            0.0       0.0       0.0       0.0       0.0
!-------------------------------------------------------------------------------
!
!     V15 sample.LOS.dat Format:
!
!        >  NEQAIR expects a version number on line 1.
!        >  Descriptive lines follow, with # comments in column 1.
!        >  The last line (no # comment) should contain species names preceded
!           by titles for data point number n, distance x (cm), and ntot (the
!           total number density, cm^-3, which is determined and written here
!           and is NOT expected in the input function file).
!
!  Point Redistribution Options:
!
!     If the option to reduce the number of points along each line of sight is
!  employed, this is done based on the magnitude of the distribution of the
!  gradient of Tv (3rd flow function) with respect to arc length.  All of the
!  function data profiles are interpolated nonlinearly to the redistributed
!  grid points here, as opposed to requiring the FLOW_INTERP[_2D] step to be
!  repeated before rerunning NEQAIR_DATA on redistributed line-of-sight data.
!  The typical function units are so large that these interpolations are
!  performed in normalized space.
!
!     It appears that 3/4ths of the input number of points per line may be the
!  best compromise, since halving the number tends to lose too much resolution
!  in the boundary layer of typical forebody lines.  For a similar reason, it
!  appears that the maximum effect of temperature gradient should be invoked
!  by using a power of 1.0 for the temperature-gradient-based shape function,
!  where lower power moderates the effect and power = 0. would give uniform
!  spacing at the other extreme (definitely NOT recommended).  Up to 2.0 is
!  permitted, but don't go that high without scrutinizing the consequences.
!
!     Extreme cases (high pressures/thin boundary layers) may not be able to
!  converge with high exponents anyway.  This can mean too much thinning in
!  the boundary layer, so a hybrid scheme is now an option: the shock region
!  (defined by the range from the outer boundary to the largest input spacing)
!  is redistributed by the gradient-base method to the SAME number of shock-
!  region points, and the rest of each profile employs a 2-sided Vinkour
!  distribution to obtain the specified total number of points, with the new
!  wall spacing determined by the ratio of the input to the output number of
!  points and a heuristic multiplier.  This preserves the character of the
!  boundary layer point distribution, where gradient-based everywhere may not.
!
!     Following the hybrid option, a curvature-based option has also been added.
!  This works with a normalized form of the Tv vs. LOS arc curve.  For the
!  extreme cases of high pressure, this method tends to improve the boundary
!  layer point distribution, but it does the same sort of thing in the shock
!  region and therefore tends to put too many points there.  Therefore, it
!  has now been changed to leave a heuristic fraction of the grid points
!  untouched in the shock region (0.2 currently), giving better resolution
!  of the boundary-layer-edge region, which appears to be where radiance
!  losses in NEQAIR occur compared with simple 2x thinning of the CFD points.
!  That thinning option is now (belatedly) a redistribution option as well.
!
!     Most of the above discussion applies to forebody lines of sight (only).
!  The hybrid and curvature-based methods both assume a shock is present, so do
!  NOT use either of them on lines in the wake where no shock is encountered.
!
!  History:
!
!     11/13/13  D.A.Saunders  Initial rewrite of a program by Aaron Brandis
!                             that expected LOS data in Tecplot format rather
!                             than PLOT3D multiblock format.
!     11/14/13    "     "     Low temperature data outside the shock cause
!                             trouble in NEQAIR.  Dinesh Prabhu recommends
!                             suppressing points with Tv < 500 K.
!     11/15/13    "     "     Guard against extremely low number densities,
!                             which can affect NEQAIR performance.
!     11/21/13    "     "     The arc length outputs now correspond to the
!                             centers of the cells at which the function data
!                             are taken to apply within NEQAIR.
!     11/23/13    "     "     An option has been provided to redistribute the
!                             line-of-sight coordinates based on |dTv/ds| and
!                             interpolate all the temperature and species
!                             profiles nonlinearly before writing the lines in
!                             NEQAIR format.
!     12/03/13    "     "     GRADDIS[3D]2 now return updated arc lengths
!                             because recovering them from the new x,y[,z]s
!                             is not identical to using the arcs that produced
!                             those new coordinates.
!     11/24/14    "     "     Extremely high pressure cases from asteroid entry
!                             studies prompted introduction of a hybrid scheme
!                             for redistributing the data points, as outlined
!                             above (appropriate for forebody data only).
!     11/25/14    "     "     Added a same-relative-spacing option that could
!                             be used on forebody lines of sight but is the
!                             likely choice for lines in the wake, where the
!                             hybrid scheme is definitely NOT recommended.
!     11/26/14    "     "     x/y/z were not being updated for the Vinokur
!                             part of the hybrid scheme (diagnostic only).
!     12/01/14    "     "     Provided a redistribution option that works with
!                             the curvature (in normalized space) of the
!                             (arc length, Tv) distribution.
!     12/02/14    "     "     The curvature-based method does what we want in
!                             the boundary layer, but puts too many points in
!                             the shock region.  Therefore, simply preserve
!                             the input data points in the shock region and
!                             confine the curvature-based method to the rest
!                             of the profile.  Lack of explicit blending of
!                             the two groups of points should not matter much.
!     12/03/14    "     "     Raised the heuristic shock-region fraction of the
!                             modified curvature-based method to 0.2 from 0.1.
!                             This puts more redistributed points in the
!                             boundary-layer-edge region, which is where the
!                             main radiance losses appear for NEQAIR compared
!                             with simple 2x thinning of the CFD grid.  Talking,
!                             of simple thinning, that is belatedly an option
!                             too now.
!     12/08/14    "     "     Averaging of sv(:) as sc(:) allowed an invalid
!                             index of 0 for a wake line of sight.
!     01/24/14    "     "     Added the option to output line-of-sight data in
!                             the order required for NEQAIR to integrate away
!                             from the body rather than towards it.  This is
!                             required for meteor studies, where the lines are
!                             probably parallel to the X axis, and a "light
!                             curve" is sought (requiring further processing
!                             of the NEQAIR results).
!     06/04/15    "     "     The redistribution prompt was missing the 'c'
!                             option; no real change.
!     03/04/19    "     "     Started option to output LOS data in the column-
!                             oriented format handled by NEQAIR v15.  The legacy
!                             format expected by V14 and its predecessors is
!                             retained as an option, and sample.LOS.dat is used
!                             to determine the kind of output produced.
!     03/22/19    "     "     Finished dual output format options and testing.
!     04/03/19    "     "     A two-temperature/8 species case encountered an
!                             ambiguity between variable names n and Na in the
!                             column-oriented data format of NEQAIR v15.  Also,
!                             the index for Tv is not necessarily 3 for this
!                             data format choice.
!
!  Author: David Saunders, ERC Inc./NASA Ames Research Center, Moffett Field/CA
!          Later with AMA, Inc. at NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure   ! Derived data type for one grid block
   use xyq_io_module          ! 2D I/O utilities
   use xyzq_io_module         ! 3D  "    "    "

   implicit none

!  Local constants:

   integer, parameter :: &
      lunf   = 1,        &    ! Input LOS function file matching input grid
      lung   = 2,        &    ! Input grid file (one LOS per block, 2D|3D)
      lunh   = 3,        &    ! File 'sample.LOS.dat' with header to transcribe
      lunout = 4,        &    ! Output LOSn.dat files
      lunkbd = 5,        &    ! Keyboard inputs
      luncrt = 6              ! Prompts and diagnostics

   real, parameter ::    &
      half     = 0.5,    &
      m_to_cm  = 100.,   &    ! Meters --> cms
      cm_to_cc = 1.e-6,  &    ! Per m^3 --> per cm^3
      Nmin     = 1000.,  &    ! Minimum number density allowed (any species)
      Tmin     = 500.         ! For suppressing points outside the shock via Tv

   logical, parameter ::   &
      false     = .false., &
      normalize = false,   &
      true      = .true.

   character, parameter :: &
      format*11 = 'unformatted'

!  Local variables:

   integer :: &
      i, i1, ilast, iline, ios, ir, is, is0, ismooth, ispecies1, iTv, length, &
      most_pts, nf, nlines, nnew, npts, nspecies, num_T

   real :: &
      arc, power, total, total_N

   logical :: &
      cr, curvature, eof, formatted_f, formatted_g, gradient, hybrid, &
      output_columns, output_n, output_ntot, redistribute, relative, thin2x, &
      towards_body, twod

   character (1) :: &
      answer
   character (32) :: &
      run_time_format
   character (128) :: &
      filename

   real, allocatable, dimension (:) :: &
      x, y, z, xnew, ynew, znew, sc, sv

   type (grid_type), pointer, dimension (:) :: &
      input_los

!  Execution:

   filename = 'sample.LOS.dat'
   open (lunh, file=filename, status='old', iostat=ios)
   if (ios /= 0) then
      write (luncrt, '(3a)') &
         'You need file ', trim (filename), ' with species names, etc., in it.'
      go to 99  ! Single stop preference
   end if

!  This sample file determines the output LOS format.  Check its first line:

   read (lunh, '(a)') filename  ! Convenient buffer
   rewind (lunh)
   output_columns = filename(1:4) /= '****' ! Easier than isolating a version #

!  Prompt for the input multiblock grid and function files of LOS data:

   ios = 1
   do while (ios /= 0)  ! Check existence and format
      filename = 'los.g'
      call reads (luncrt, &
                  'Input multiblock grid    [<CR> = los.g]: ',  &
                  lunkbd, filename, cr, eof)
      if (eof) go to 99  ! Abort, or whatever

      call determine_grid_form (filename, lung, formatted_g, ios)
   end do

   i1 = 1;  if (formatted_g) i1 = 3
   open (lung, file=filename, form=format(i1:11), status='OLD')

   ios = 1
   do while (ios /= 0)  ! Repeat for the function file
      filename = 'los.f'
      call reads (luncrt, &
                  'Associated function file [<CR> = los.f]: ',  &
                  lunkbd, filename, cr, eof)
      if (eof) go to 99

      call determine_grid_form (filename, lunf, formatted_f, ios)
   end do

   i1 = 1;  if (formatted_f) i1 = 3
   open (lunf, file=filename, form=format(i1:11), status='OLD')

!  Determine whether it's a 2D or 3D case by trying to read a 3D header:

   call xyz_header_io (1, lung, formatted_g, nlines, input_los, ios)
   twod = ios /= 0       ! Why isn't this reliable??

   if (.not. twod) then
      do i = 1, nlines
         if (input_los(i)%ni /= 1 .or. input_los(i)%nj /= 1 .or. &
             input_los(i)%nk <= 1) twod = true
      end do
   end if

!! write (*, '(a, /, (i4, 3i6))') ' LOS grid block dimensions:', &
!!    (i, input_los(i)%ni, input_los(i)%nj, input_los(i)%nk, i = 1, nlines)

   if (twod) then
      write (luncrt, '(a)') ' Flow solution appears to be 2D.'
      rewind (lung)
      call xy_header_io   (1, lung, formatted_g, nlines,     input_los, ios)
      call q_header_io_2d (1, lunf, formatted_f, nlines, nf, input_los, ios)
      npts = input_los(1)%nj
   else
      write (luncrt, '(a)') ' Flow solution appears to be 3D.'
      call q_header_io    (1, lunf, formatted_f, nlines, nf, input_los, ios)
      npts = input_los(1)%nk
   end if
   if (ios /= 0) go to 99

   write (luncrt, '(a, i7, i4)') ' # lines of sight and # functions found:', &
      nlines, nf
   write (luncrt, '(a, i5)') ' # points found on the first line of sight:', npts

   most_pts = 0  ! For picking # processors if NEQAIR is parallelized

   write (luncrt, '(a)') &
      ' Redistribute the line of sight points?', &
      '    n = no', &
      '    y = yes via dT/ds (gradient-based, normally satisfactory)', &
      '    h = hybrid dT/ds + Vinokur scheme (forebody only, extreme T, p)', &
      '    c = curvature-based along the (normalized) T profile', &
      '    r = same relative spacing'
   write (luncrt, '(a)', advance='no') &
      '    t = simple 2x thinning: [n|y|h|c|r|t]: '
   read  (lunkbd, *) answer

   redistribute = answer /= 'n'
   curvature    = answer == 'c'
   gradient     = answer == 'y'
   hybrid       = answer == 'h'
   relative     = answer == 'r'
   thin2x       = answer == 't'

   if (redistribute) then
      write (luncrt, '(2a)') ' Redistribution option chosen: ', answer
      if (.not. thin2x) then
         write (luncrt, '(a)', advance='no') ' Desired # pts. per output line: '
         read  (lunkbd, *) nnew
      end if
      if (.not. (relative .or. thin2x)) then
         write (luncrt, '(a)', advance='no') &
            ' Initial shape function exponent to use [0. -> 2.; try 1.]: '
         read  (lunkbd, *) power
         write (luncrt, '(a)', advance='no') &
            ' Smoothing control? [0|1|2|3; try 1 => shape fn. only]:     '
         read  (lunkbd, *) ismooth
         write (luncrt, '(a, f4.1)') ' Initial shape fn. exponent chosen:', &
            power
      end if
   end if

!  Option to have NEQAIR integrate away from the body instead of towards it:

   towards_body = true
   call ready (luncrt, &
              'Should NEQAIR integrate towards the body? [y|n; <CR> = yes]: ', &
               lunkbd, towards_body, cr, eof)
   write (luncrt, '(a, l2)') 'towards_body:', towards_body
   if (eof) go to 99

!  Process one line of sight at a time:

   do iline = 1, nlines

!     Construct the output LOS file name as LOS-n.dat:

      call numbered_name ('LOS-', iline, filename, length)

      length = length + 4
      filename(length-3:length) = '.dat'

      open (lunout, file=filename(1:length), status='unknown', iostat=ios)

      call transcribe_los_header ()  ! Transcribe sample LOS.dat header & rewind
      if (ios /= 0) go to 99         ! Species count mismatch?

!     Input lines of sight are dimensioned (1,nj) or (1,1,nk), but reversing
!     the order is preferable for the I/O of normal cases (towards the body):

      if (twod) then

         input_los(iline)%ni = input_los(iline)%nj
         input_los(iline)%nj = 1
         npts = input_los(iline)%ni
         input_los(iline)%mi = input_los(iline)%mj
         input_los(iline)%mj = 1

         call xy_allocate (input_los(iline), ios)
         call xy_block_io (1, lung, formatted_g, npts, &
                           input_los(iline)%x, input_los(iline)%y, ios)

         allocate (x(npts), y(npts), sv(npts), sc(npts))

!        Some redistributions require point indices increasing towards the body.
!        They may be reversed before output.

         x(:) = input_los(iline)%x(npts:1:-1,1,1)  ! k = 1 = body assumed
         y(:) = input_los(iline)%y(npts:1:-1,1,1)  ! for the input lines

         call chords2d (npts, x, y, normalize, total, sv) ! Vertex-centered arcs

         call q_allocate_2d (input_los(iline), nf, ios)
         call q_block_io_2d (1, lunf, formatted_f, nf, &
                             input_los(iline)%mi, 1,   &
                             input_los(iline)%q, ios)
      else  ! 3D

         input_los(iline)%ni = input_los(iline)%nk
         input_los(iline)%nk = 1
         npts = input_los(iline)%ni
         input_los(iline)%mi = input_los(iline)%mk
         input_los(iline)%mk = 1

         call xyz_allocate (input_los(iline), ios)
         call xyz_block_io (1, lung, formatted_g, npts, &
                            input_los(iline)%x, input_los(iline)%y, &
                            input_los(iline)%z, ios)

         allocate (x(npts), y(npts), z(npts), sv(npts), sc(npts))

!        Set up the lines with index increasing towards the body for some
!        redistribution options.  They may be reversed before output.

         x(:) = input_los(iline)%x(npts:1:-1,1,1)  ! k = 1 = body assumed
         y(:) = input_los(iline)%y(npts:1:-1,1,1)  ! for the input lines
         z(:) = input_los(iline)%z(npts:1:-1,1,1)

         call chords3d (npts, x, y, z, normalize, total, sv)  ! Vertex arcs

         call q_allocate (input_los(iline), nf, ios)
         call q_block_io (1, lunf, formatted_f, nf,  &
                          input_los(iline)%mi, 1, 1, &
                          input_los(iline)%q, ios)
      end if

!     Option to redistribute each line of sight to [presumably] fewer points:

      if (redistribute) then
         call redistribute_los ()  ! Update everything so that the rest
      end if                       ! works either way

      sv(:) = m_to_cm * sv(:)  ! NEQAIR expects arc lengths to be in cm.

!     Arrange to suppress points outside the shock.  Dinesh says use Tv, not T.

      do i = npts, 1, -1
         ilast = i
         if (input_los(iline)%q(iTv,i,1,1) >= Tmin) exit
      end do

      is0 = npts - ilast + 1
      most_pts = max (most_pts, ilast)

!     Make the arc lengths correspond to cell centers at which the function
!     data are taken to apply within NEQAIR:

      sc(is0) = sv(is0)
      do i = is0 + 1, npts
         sc(i) = (sv(i-1) + sv(i))*half
      end do
      sc(npts) = sv(npts)  ! Gives a half cell width for the wall function data

!     Output data in NEQAIR form, one grid point at a time:

      if (towards_body) then

!        The last arc length is taken by NEQAIR (14-?) to be total arc length,
!        so we subtract the shock arc length from all arc lengths upon output,
!        meaning they start at zero.

         if (output_columns) then  ! Optional outputs depend on header titles

            do i = is0, npts
               ir  = npts + 1 - i
               is  = i - is0 + 1
               input_los(iline)%q(ispecies1:,ir,1,1) = &
                  max (cm_to_cc*input_los(iline)%q(ispecies1:,ir,1,1), Nmin)
               total_N  =  sum (input_los(iline)%q(ispecies1:,ir,1,1))
               arc = sc(i) - sc(is0)
               if (output_n) then
                  if (output_ntot) then
                     write (lunout, run_time_format) &
                        is, arc, total_n, input_los(iline)%q(:,ir,1,1)
                  else
                     write (lunout, run_time_format) &
                        is, arc, input_los(iline)%q(:,ir,1,1)
                  end if
               else
                  if (output_ntot) then
                     write (lunout, run_time_format) &
                        arc, total_n, input_los(iline)%q(:,ir,1,1)
                  else
                     write (lunout, run_time_format) &
                        arc, input_los(iline)%q(:,ir,1,1)
                  end if
               end if
            end do

         else  ! Legacy line-oriented LOS format

            do i = is0, npts
               ir = npts + 1 - i
               is = i - is0 + 1
               input_los(iline)%q(ispecies1:,ir,1,1) = &
                  max (cm_to_cc*input_los(iline)%q(ispecies1:,ir,1,1), Nmin)
               total_N  =  sum (input_los(iline)%q(ispecies1:,ir,1,1))
               write (lunout, '(i5, 2es15.7, 4f10.1, /, (5x, 4es15.7))') &
                  is, sc(i) - sc(is0), total_N, input_los(iline)%q(:,ir,1,1)
               write (lunout, '(a)')
            end do
            write (lunout, '(a)') '    0  0.0  0.0  0.0  0.0  0.0  0.0'

         end if

      else  ! Meteor case?  We want NEQAIR to integrate away from the body

         total = sc(npts)
         do i = 1, npts
            sv(i) = total - sc(npts - i + 1)
         end do

         if (output_columns) then

            do i = 1, npts - is0 + 1
               input_los(iline)%q(ispecies1:,i,1,1) = &
                  max (cm_to_cc*input_los(iline)%q(ispecies1:,i,1,1), Nmin)
               total_N  =  sum (input_los(iline)%q(ispecies1:,i,1,1))
               if (output_n) then
                  if (output_ntot) then
                     write (lunout, run_time_format) &
                        i, sv(i), total_n, input_los(iline)%q(:,i,1,1)
                  else
                     write (lunout, run_time_format) &
                        i, sv(i), input_los(iline)%q(:,i,1,1)
                  end if
               else
                  if (output_ntot) then
                     write (lunout, run_time_format) &
                        sv(i), total_n, input_los(iline)%q(:,i,1,1)
                  else
                     write (lunout, run_time_format) &
                        sv(i), input_los(iline)%q(:,i,1,1)
                  end if
               end if
            end do

         else  ! Line-oriented output

            do i = 1, npts - is0 + 1
               input_los(iline)%q(ispecies1:,i,1,1) = &    ! 5 <-> 1st # density
                  max (cm_to_cc*input_los(iline)%q(ispecies1:,i,1,1), Nmin)
               total_N  =  sum (input_los(iline)%q(ispecies1:,i,1,1))
               write (lunout, '(i5, 2es15.7, 4f10.1, /, (5x, 4es15.7))') &
                  i, sv(i), total_N, input_los(iline)%q(:,i,1,1)
               write (lunout, '(a)')
            end do
            write (lunout, '(a)') '    0  0.0  0.0  0.0  0.0  0.0  0.0'

         end if

      end if

      close (lunout)

      deallocate (input_los(iline)%x, input_los(iline)%y, x, y, sv, sc)
      if (.not. twod) deallocate (input_los(iline)%z, z)
      deallocate (input_los(iline)%q)

   end do  ! Next line of sight

   close (lunh)
   close (lung)
   close (lunf)

   write (luncrt, '(a, i5)') ' Highest output LOS point count:', most_pts

99 continue  ! Avoid system-dependent STOP behavior

   contains  ! Internal procedures for program neqair_data

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine redistribute_los ()

!     Redistribute the current line of sight data by the indicated method.
!     Some methods are inappropriate for lines of sight in the wake that don't
!     encounter a shock.  See the program header for further discussion.
!
!     Reducing the number of data points is intended to speed the NEQAIR
!     calculations.
!
!     The LOS data are ordered outboard to inboard at this point.
!     Upon return, everything is made to look like the no-redistribute case
!     via deallocation and reallocation.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local constants:

      integer,       parameter :: lunplot = 8       ! For plottable output
      real,          parameter :: fudge   = 20.     ! x CFD wall ds x nold/nnew
      character (1), parameter :: method  = 'M'     ! Monotonic local cubics

!     Local variables:

      integer :: i, ier, n, n1, n2
      real    :: ds1, ds2, dsmax, Tshock
      real, allocatable, dimension (:)   :: f, xnew, ynew, znew, snew
      real, allocatable, dimension (:,:) :: fnew

!     Execution:

      if (thin2x) nnew = (npts + 1) / 2

      allocate (xnew(nnew), ynew(nnew), f(npts), snew(nnew))
      if (.not. twod) allocate (znew(nnew))

      f(:) = input_los(iline)%q(iTv,npts:1:-1,1,1)  ! Reversed to match x(:)

      if (thin2x) then  ! Simple 2x thinning

         i = 1
         do n = 1, npts, 2
            snew(i) = sv(n)
            xnew(i) =  x(n)
            ynew(i) =  y(n)
            if (.not. twod) znew(i) = z(n)
            i = i + 1
         end do

      else if (relative) then  ! Change the number of pts./same relative spacing

         call changen1d (1, npts, sv, 1, nnew, snew)
         call lcsfit2 (npts, sv, x, method, nnew, snew, xnew)
         call lcsfit2 (npts, sv, y, method, nnew, snew, ynew)
         if (.not. twod) then
            call lcsfit2 (npts, sv, z, method, nnew, snew, znew)
         end if

      else if (curvature) then

         if (twod) allocate (znew(1))  ! Can't pass an unallocated array

         call curvature_based (npts, sv, f, nnew, snew, xnew, ynew, znew)

         if (twod) deallocate (znew)

      else if (gradient) then

         if (twod) then
            call graddis2 (npts, x, y, sv, f, nnew, power, ismooth, luncrt, &
                           xnew, ynew, snew, ier)
         else
            call graddis3d2 (npts, x, y, z, sv, f, nnew, power,ismooth,luncrt, &
                             xnew, ynew, znew, snew, ier)
         end if

      else if (hybrid) then  ! Same # pts. around shock (dT/ds); Vinokur in b.l.
                             ! Inappropriate for wake lines of sight
         n1 = 2  ! Locate the largest (reversed) spacing after the shock
         dsmax  = sv(2) - sv(1)
         Tshock = f(1)*3.
         do n = 2, npts
            ds1 = sv(n) - sv(n-1)
            if (ds1 > dsmax .and. f(n) > Tshock) then
                dsmax = ds1
                n1 = n
            end if
         end do
         n1  = n1 - 1
         n2  = nnew - n1 + 1  ! Overlap the largest spacing

!        Redistribute the the first n1 LOS coordinates according to |dTv/ds|.
!        We don't change the number of points in the shock region up to the
!        largest spacing, because we don't want to give up resolution there.

         if (twod) then
            call graddis2 (n1, x, y, sv, f, n1, power, ismooth, luncrt, &
                           xnew, ynew, snew, ier)
         else  ! 3D
            call graddis3d2 (n1, x, y, z, sv, f, n1, power, ismooth, luncrt, &
                             xnew, ynew, znew, snew, ier)
         end if

!        Increase the new wall spacing considerably, because the temperatures
!        are much lower there.  We want to keep points where the temperatures
!        are high in the boundary layer.

         ds1 = snew(n1) - snew(n1-1)
         ds2 = fudge*(sv(npts) - sv(npts -1)) * nint (real (npts) / real (nnew))
         ds2 = min (ds1, ds2)  ! Unlikely adjustment
         snew(nnew) = sv(npts)

         call vinokur (n1, nnew, ds1, ds2, snew, luncrt, ier)  ! 2-sided stretch

         write (luncrt, '(a, 5i5, 2es16.8)') &
            ' Hybrid redistrib. details.  ier, npts, nnew, n1, n2, ds1, ds2:', &
            ier, npts, nnew, n1, n2, ds1, ds2

         call lcsfit2 (npts, sv, x, method, n2, snew(n1), xnew(n1))
         call lcsfit2 (npts, sv, y, method, n2, snew(n1), ynew(n1))
         if (.not. twod) then
            call lcsfit2 (npts, sv, z, method, n2, snew(n1), znew(n1))
         end if

      end if

!     For the first line, save plottable data to examine the behavior of the
!     redistribution method chosen:

      if (iline == 1) then
         open (lunplot, file='input_xsf.dat', status='unknown')
         if (twod) then
            write (lunplot, '(a)') '# Input x, y, s, f'
            write (lunplot, '(4es16.8)') &
               (x(n), y(n), sv(n), f(n), n = 1, npts)
         else
            write (lunplot, '(a)') '# Input x, y, z, s, f'
            write (lunplot, '(5es16.8)') &
               (x(n), y(n), z(n), sv(n), f(n), n = 1, npts)
         end if
         close (lunplot)
      end if

!     Interpolate all functions at the redistributed arc lengths:

      allocate (fnew(nnew,nf))

      do n = 1, nf
         f(:) = input_los(iline)%q(n,npts:1:-1,1,1)  ! Reverse to match x(:)

         call lcsfit2 (npts, sv, f, method, nnew, snew, fnew(:,n))
      end do

      if (iline == 1) then
         open  (lunplot, file='output_xsf.dat', status='unknown')
         if (twod) then
            write (lunplot, '(a)') '# Output x, y, s, f'
            write (lunplot, '(4es16.8)') &
               (xnew(n), ynew(n), snew(n), fnew(n,iTv), n = 1, nnew)
         else
            write (lunplot, '(a)') '# Output x, y, z, s, f'
            write (lunplot, '(5es16.8)') &
               (xnew(n), ynew(n), znew(n), snew(n), fnew(n,iTv), n = 1, nnew)
         end if
         close (lunplot)
      end if

      deallocate (f)

!     Make everything look like the no-distribute case upon return:

      deallocate (sv, sc, input_los(iline)%q)
      allocate (sv(nnew), sc(nnew), input_los(iline)%q(nf,nnew,1,1))

      sv(:) = snew(:)
      deallocate (snew)

      do n = 1, nf
         input_los(iline)%q(n,nnew:1:-1,1,1) = fnew(:,n)
      end do

      deallocate (xnew, ynew, fnew)
      if (.not. twod) deallocate (znew)

      npts = nnew

      end subroutine redistribute_los

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine curvature_based (npts, sv, f, nnew, snew, xnew, ynew, znew)

!     Redistribute LOS points via curvature with respect to a normalized form of
!     the Tv vs. LOS arc curve.  Since the plain curvature method puts too many
!     points in the shock region, that region is now left alone; the curvature
!     method should produce clustering near the boundary layer edge in place of
!     the unwanted clustering of the data points towards the wall.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)  :: npts                 ! Original LOS point count
      real,    intent (in)  :: sv(npts)             ! Original LOS arc lengths
      real,    intent (in)  :: f(npts)              ! Original LOS Tv data
      integer, intent (in)  :: nnew                 ! Desired point count
      real,    intent (out) :: snew(nnew)           ! Redistributed arc lengths
      real,    intent (out), dimension (nnew) :: &
                                 xnew, ynew, znew   ! For diagnostics only
!     Local constants:

      real, parameter :: one = 1.
      real, parameter :: shock_fraction = 0.2  ! Applied to npts; heuristic
                         ! because detecting the inboard edge of the shock is
                         ! not easy for extreme cases with no shock spike
!     Local variables:

      integer :: i, ic, icurv, id, ier, ioffset, ishock, left
      integer :: ncurv, nred, nshock
      real    :: r, rescale, scale(2), shift(2), Tshock, total
      real    :: arcs(npts), arcsnew(nnew), svnorm(npts), fnorm(npts)

!     Execution:

!     Estimate the shock location and avoid including the free-stream data
!     points.

      Tshock = f(1)*1.1
      do i = 1, npts
         if (f(i) > Tshock) then
            ishock = i - 1     ! First data point used by this method
            exit
         end if
      end do

      icurv  = ishock + nint (shock_fraction * real (npts))  ! Start of curva-
      nshock = icurv - ishock + 1                            ! ture-based zone
      snew(1:nshock) = sv(ishock:icurv)
      xnew(1:nshock) =  x(ishock:icurv)
      ynew(1:nshock) =  y(ishock:icurv)
      if (.not. twod) znew(1:nshock) = z(ishock:icurv)

      ioffset = icurv - 1          ! Used for interpolating s,x,y,z
      ncurv   = npts - ioffset     ! # data pts. to include in curvature calcs.
      nred    = nnew - nshock + 1  ! # redistributed points counting common pt.

      ier = -99  ! Kludge to normalize both "x" and "y" (i.e., sv and f)

      call curvdis2 (ncurv, sv(icurv), f(icurv), nred, power, ismooth, &
                     luncrt, true, snew(nshock), snew(nshock), ier)

!     These snew values are along the (sv,T) curve, not along the LOS arc.
!     Furthermore, they are badly scaled (having been denormalized).
!     We need to determine their equivalents within the sv(:) distribution.
!     Retrieve the normalized "arcs" that CURVDIS worked with in CURVDIS2:

      svnorm(1:ncurv) = sv(icurv:npts)
       fnorm(1:ncurv) =  f(icurv:npts)

      call getscale ('N', 2, ncurv, sv(icurv), f(icurv), f(icurv), &
                     scale, shift, ier)
      call usescale ('N', 2, ncurv, svnorm, fnorm, fnorm, scale, shift, ier)

      call chords2d (ncurv, svnorm, fnorm, false, total, arcs)

!     Retrieve the normalized form of the redistributed arcs:

      rescale = total / snew(nnew)
      arcsnew(1:nred) = snew(nshock:nnew) * rescale

!     Locate each new arc length within the original (normalized) arcs.
!     The new (x,y[,z])s can also be estimated (for diagnostic purposes only).

      do i = 1, nred
         call interval (ncurv, arcs, arcsnew(i), one, left)

         r  = (arcsnew(i) - arcs(left)) / (arcs(left+1) - arcs(left))
         ic = nshock - 1 + i
         id = ioffset + left
         snew(ic) = (one - r)*sv(id) + r*sv(id+1)
         xnew(ic) = (one - r)* x(id) + r* x(id+1)
         ynew(ic) = (one - r)* y(id) + r* y(id+1)
         if (.not. twod) znew(ic) = (one - r)*z(id) + r*z(id+1)
      end do

      end subroutine curvature_based

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine transcribe_los_header ()  ! From 'sample.LOS.dat'
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, parameter :: max_length = 512  ! Allow for many species
      integer, parameter :: num_titles = 10   ! Needs to include all possible Ts
      character (1), parameter :: blank = ' '

      integer :: i, l, ndashed, nheader, num_first
      character (3) :: seps
      character (4) :: first_titles(num_titles) ! Are n and/or ntot present?
      character (max_length) :: line

      if (output_columns) then  ! NEQAIR 15+

!        Transcribe down to and including the first non-comment (titles) line
         do
           read (lunh, '(a)') line;  l = len_trim (line)
           write (lunout, '(a)') line(1:l)
           if (line(1:1) /= '#') then
             if (l > 20) exit ! Handle a possible blank line preceding titles
           end if
         end do

!        For the first line of sight only, determine the data contents:

         if (iline == 1) then  ! Is n and/or ntot to be output? Also, count Ts.

           seps = ' , '; seps(3:3) = char (9)  ! Tab
           call token_count (line(1:l), seps, nheader)  ! Count all titles

           write (luncrt, '(a, i4)') '# variable names found: ', nheader

!          Tokenize leading titles to look for n, ntot and various temperatures:

           num_first = num_titles
           call tokens (line(1:l), num_first, first_titles)
           output_n = false;  output_ntot = false;  num_T = 0

           do i = 1, num_first  ! TOKENS converts to upper case
             l = len_trim (first_titles(i))
             if (l == 1 .and. first_titles(i)(1:1) == 'N')    output_n    = true
             if (l == 4 .and. first_titles(i)(1:4) == 'NTOT') output_ntot = true
             if (first_titles(i)(1:1) == 'T') num_T = num_T + 1
             if (l == 2 .and. first_titles(i)(1:2) == 'TV')   iTv = i
           end do

!          Handle a possible clash between 'n' (line number) and 'N' (nitrogen):

           if (output_n) then   ! This may have been because of N, not n
               output_n = false
               do i = 1, num_T + 1  ! n should precede a T if it's there
                 l = len_trim (first_titles(i))
                 if (l == 1 .and. first_titles(i)(1:1) == 'N') output_n = true
               end do
           end if

!          iTv needs to point into the function file, not the variable name list:

           iTv = iTv -1  ! For arc length
           if (output_n)    iTv = iTv - 1
           if (output_ntot) iTv = iTv - 1

!          Trap mismatched species counts:

           nspecies = nheader - num_T - 1  ! Don't count arc length
           if (output_n)    nspecies = nspecies - 1
           if (output_ntot) nspecies = nspecies - 1
           if (nspecies + num_T /= nf) then
             write (luncrt, '(a, 3i6)') &
               '*** Species mismatch? nf, # Ts, nspecies:', nf, num_T, nspecies
             l = len_trim (line)
             write (luncrt, '(a)') line(1:l)
             write (luncrt, '(a, l1)') 'output_n?    ', output_n, &
                                       'output_nTot? ', output_ntot
             ios = -1
           end if
           ispecies1 = num_T + 1  ! Index into input function data

           run_time_format = '(i4, 2es15.8, nes15.8, nnes15.8)'
!                             12345678901234567890123456789012
           if (.not. output_n)    run_time_format(2:4) = blank
           if (.not. output_ntot) run_time_format(6:6) = blank
           write (run_time_format(15:15), '(i1)') num_T
           write (run_time_format(24:25), '(i2)') nspecies
        end if

      else  ! Legacy data format, NEQAIR 14 and earlier

        iTv = 3   ! Because T, T, Tv, Tv are stipulated in the function file
        ios = 0;  ndashed = 0
        do while (ndashed < 2) ! Until the 2nd '-----' line has been transcribed
          read (lunh, '(a)') line
          write (lunout, '(a)') trim (line)
          if (line(1:1) == '-') ndashed = ndashed + 1
        end do
        if (ndashed /= 2) then
          write (luncrt, '(a)') '*** Missing asterisk lines in sample header?'
          ios = -1
        end if
        ispecies1 = 5  ! Index of first # density; 4 temperatures are assumed

      end if

      rewind (lunh)  ! For the next call

      end subroutine transcribe_los_header

   end program neqair_data
