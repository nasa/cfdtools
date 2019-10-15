!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program prepare_neqair_data

!     This adaptation of prepare_rapid_analysis prompts the user for specifics
!     about a flow solution and the body point(s) at which radiative heating
!     calculations are to be performed with NEQAIR, either in tangent-slab mode
!     or full angular integration mode.
!
!     The working directory should contain (or point to) the vertex-centered
!     volume grid and the cell-centered volume grid and flow data (three volume
!     files in PLOT3D multiblock form).  Whether they are 2D or 3D is automatic-
!     ally determined.  If full integration is specified for a 2D axisymmetric
!     case, 3D revolved grid and function files will be produced as part of the
!     automation procedures of which the present utility is the first step
!     following a flow calculation with the right choice of species.
!
!     The output from this utility is a shell script named prepare_LOS_data that
!     (when sourced) invokes the indicated utilities and control files (also
!     generated here) to set up line-of-sight data in NEQAIR's LOS.dat form as
!     LOS-1.dat, LOS-2.dat, ..., LOS-n.dat in the working directory.
!
!     A separate script (run_neqair) can then submit NEQAIR batch jobs for the
!     indicated range of line numbers (one line of sight per PBS job) via a
!     third script, neqair.pbs, which is expected to contain appropriate
!     walltime and nodes statements, load the intended NEQAIR module, and
!     launch NEQAIR on the appropriate number of processors.
!
!     Files neqair.inp and neqair.pbs are expected to be in the working
!     directory along with the body point file (indices or coordinates) and a
!     sample.LOS.dat file indicating the correct species, although the user
!     can choose to proceed if any of neqair.inp, neqair.pbs, or sample.LOS.dat
!     are missing.  The neqair_data step (los.* --> LOS.dat format) can then
!     be performed explicitly at a later time.
!
!     When NEQAIR is run on a range of lines, NEQAIR results will appear in
!     subdirectores /LINE-1, /LINE-2, ...
!
!     Full integration mode should work with a single body point per working
!     directory.  Since actually performing the integration cannot occur
!     until all the NEQAIR results are in hand, no attempt is made here to
!     deal with the NEQAIR_INTEGRATION step, which is prompt-driven.
!
!  Site-Specifics:
!
!     If this program is installed at another site, the following path should
!     be edited appropriately prior to compiling and linking:
!
!        /apps/pkgs/cfdtools/         ! See parameter constant "path" below
!
!     It is assumed that all the utilities invoked here reside in subdirectories
!     with uppercase names and this indicated path.
!
!     E.g., lines_of_sight, flow_interp, ... should be installed here:
!
!        /share/apps/cfdtools/LINES_OF_SIGHT/lines_of_sight
!        /share/apps/cfdtools/FLOW_INTERP/flow_interp   and so on.
!
!     The following utilities may be invoked by the script from this program:
!
!        lines_of_sight       ! One or more lines of sight in a 3D volume grid
!        lines_of_sight_2d    ! One or more lines of sight in a 2D volume grid
!        hemispheres_of_sight ! Many lines of sight for one body point
!        flow_interp          ! 3D flow interpolation at LOS coordinates
!        flow_interp_2D       ! 2D   "    "    "    "    "    "    "
!        neqair_data          ! PLOT3D *.g,*.f --> NEQAIR's LOS.dat form
!
!  Assumptions:
!
!     The line-of-sight calculations assume that the grid has a single layer of
!     blocks, with k = 1 at the geometry surface.
!
!  History:
!
!     07/01/14  D.A.Saunders  PREPARE_NEQAIR_DATA design, in the style of
!                             PREPARE_LOCAL_ANALYSIS.
!     07/07/14-  "     "      Initial implementation and ...
!     07/11/14                ... testing.
!     07/17/14   "     "      It wasn't clear to Dinesh that sample.LOS.dat
!                             needs the correct species in it.  Showing the
!                             path to all needed files should help.
!     07/31/14   "     "      Allowing for descriptive text on each body point
!                             input line means suppressing the check on whether
!                             the number of tokens = ndim.
!     11/26/14   "     "      Handle the wider choice of redistribution options
!                             now offered by NEQAIR_DATA (hybrid, relative).
!     12/01/14   "     "      A curvature-based redistribution option has been
!                             added to NEQAIR_DATA.
!     12/04/14   "     "      A simple 2x thinning option as also been added.
!     01/24/14   "     "      NEQAIR_DATA has an extra prompt to allow NEQAIR
!                             to integrate out from the body, as might be
!                             required for meteor studies.
!     02/26/18   "     "      Handled the HEMISPHERES_OF_SIGHT option to output
!                             less than a full hemisphere of lines (if a cone
!                             angle < 90 is entered on the body pt. x/y/z line).
!     04/04/18   "     "      Updated the path of the various utilities. Didn't
!                             eliminate it altogether because of the diagnostic
!                             that points to where certain ancillary files might
!                             be found.
!     05/24/18   "     "      "Same relative distribution" needs to write the
!                             requested number of points (unlike 2x thinning).
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, CA
!                  Now with AMA, Inc. at NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Constants:

   integer, parameter :: &
      lunin  = 1,        &  ! For checking if a file is formatted or unformatted
      lunout = 2,        &  ! Output control files
      lunsh  = 3,        &  ! For the output shell script
      lunkbd = 5,        &  ! For keyboard entries
      luncrt = 6,        &  ! For prompts and diagnostics
      nc     = 96           ! # characters in file names

   logical, parameter :: &
      false  = .false.,  &
      true   = .true.

!  The following path indicates where PREPARE_NEQAIR_DAT, LINES_OF_SIGHT, etc.
!  are expected to be installed.  Be sure to get its length right.

   character, parameter :: &
!!!   path*21          = '/share/apps/cfdtools/', &
      path*20          = '/apps/pkgs/cfdtools/', &
      blank*1          = ' ',   &
      n*1              = 'n',   &
      y*1              = 'y',   &
      file_rev_g*11    = 'revolved.gu', &
      file_rev_cc_g*14 = 'revolved_cc.gu', &
      file_rev_cc_f*14 = 'revolved_cc.fu'

!  Variables:

   integer :: &
      ier, ios, ismooth, ndim, ne, npts, ntokens

   real :: &
      cone_angle, power, s1, s2, xbp, ybp, zbp

   logical ::   &
      cr, eof, formatted_cc, formatted_vc, hemisphere, off_center, proceed, &
      relative, thin2x, towards_body, yes

   character (1) :: &
      answer, yesno

   character (nc) :: &
      filename_BP,   filename_cc_f,  filename_cc_g, &
      filename_grid, filename_los_f, filename_los_g, buffer

   character (nc-10) :: &
      identifier        ! May get '.1.lines.g' appended to make a filename <= nc

!  Execution:
!  !!!!!!!!!!

   write (luncrt, '(/, (a))') &
      ' Requirements summary:', ' ', &
      '    Grids and solutions may be 2D or 3D, formatted or not.', &
      '    If flow data are interpolated to grid vertices (interp = 2), the',  &
      '    associated grid can be used twice: once for the line-of-sight',     &
      '    coordinates and once for the flow interpolations.',                 &
      '    DPLR''s POSTFLOW should have been run once or twice already:',      &
      '      [(1) vertex-centered volume grid (only)    (interp = 2);]',       &
      '       (2) volume grid & flow data (interp = 1, but interp = 2 is OK);',&
      '           ivarp = 0 120 120 125 125 1200   (2 temperatures)   or',     &
      '           ivarp = 0 120 120 120 120 1200   (1 temperature)',           &
      '    The working directory should also contain or point to these files:',&
      '       (3) body pt. file (x/y/z or grid i/j/k [+ optional cone angle]', &
      '       (4) neqair.inp      (tangent-slab or not, regions, etc.)',       &
      '       (5) neqair.pbs      (# nodes, desired NEQAIR version, etc.)',    &
      '       (6) run_neqair      (first/last LOS numbers as arguments)',      &
      '       (7) sample.LOS.dat  (top portion indicating relevant species)'
   write (luncrt, '(/, 3a, //, 2a, /, a, /)') &
      ' Examples of files (4)-(7) should be found here: ', path, &
      'PREPARE_NEQAIR_DATA', &
      ' Be sure to enter the right species in sample.LOS.dat and to', &
      ' check neqair.inp and neqair.pbs.', &
      ' The *.chem file is the best place to get the right species from.', &
      ' It is OK to use vertex-centered flow data.  This will mean the', &
      ' same grid name is entered twice.'

   ier = 0
   call check_existence ('neqair.inp',     luncrt, ios);  if (ios /= 0) ier = 1
   call check_existence ('neqair.pbs',     luncrt, ios);  if (ios /= 0) ier = 1
   call check_existence ('run_neqair',     luncrt, ios);  if (ios /= 0) ier = 1
   call check_existence ('sample.LOS.dat', luncrt, ios);  if (ios /= 0) ier = 1

   if (ier == 1) write (luncrt, '(/, a)') &
      ' You may still go ahead with the set-up.'
   proceed = true
   call ready (luncrt, 'Proceed?  [y|n|^D=quit|<CR>=proceed]: ', &
               lunkbd, proceed, cr, eof)
   if (eof) go to 99
   if (.not. proceed) go to 99

!  (1) Vertex-centered volume grid (may be the same as for the flow data):

   filename_grid = 'volrad.gu';  ios = 1
   do while (ios /= 0)
      call reads (luncrt, &
          'Vertex-centered volume grid? [<CR> = ' // trim (filename_grid) // &
          ']: ', lunkbd, filename_grid, cr, eof)
      if (eof) go to 99

      call determine_grid_form (filename_grid, lunin, formatted_vc, ios)
   end do

   call determine_grid_dim (filename_grid, lunin, formatted_vc, ndim, ios)

!  (2) Cell- or vertex-centered volume grid and flow data:

   filename_cc_g = 'volrad.gu';  ios = 1
   do while (ios /= 0)
      call reads (luncrt, &
         'Flow volume grid?            [<CR> = ' // trim (filename_cc_g) // &
           ']: ', lunkbd, filename_cc_g, cr, eof)
      if (eof) go to 99

      call determine_grid_form (filename_cc_g, lunin, formatted_cc, ios)
   end do

   filename_cc_f = 'volrad.fu';  ios = 1
   do while (ios /= 0)
      call reads (luncrt, &
         'Flow volume data?            [<CR> = ' // trim (filename_cc_f) // &
           ']: ', lunkbd, filename_cc_f, cr, eof)
      if (eof) go to 99

      call determine_grid_form (filename_cc_f, lunin, formatted_cc, ios)
   end do

!  (3) Body point file:

   filename_BP = 'bp.inp';  ios = 1
   do while (ios /= 0)
      call reads (luncrt, &
            'Target body point file? [<CR> = ' // trim (filename_BP) // ']: ', &
                  lunkbd, filename_BP, cr, eof)
      if (eof) go to 99

      open (lunin, file=filename_BP, status='old', iostat=ios)

   end do

!  Catch a possible user error (originally, but now cone_angle may be present):

   read  (lunin, '(a)') buffer
   close (lunin)

   call token_count (buffer, blank, ntokens)

!  Suppress this because it prevents trailing descriptive text.
!  The LOS codes will read only ndim items per line of body point data, except
!  now the hemisphere option looks for an optional fourth (cone angle) token.

!! if (ntokens /= ndim) then
!!    write (luncrt, '(a, 2i4)') &
!!       '*** Likely error: ndim /= # body pt. coords./indices:', ndim, ntokens
!!    go to 99
!! end if

!  Hemisphere(s)-of-sight, tangent-slab towards body, or integrate away from
!  the body as for meteor studies?

   answer = 't'  ! Most likely tangent-slab towards body
   if (ndim == 2) then
      call readc (luncrt, &
                 'Tangent-slab/towards body (t), or out from body (o)? [t]: ', &
                  lunkbd, answer, cr, eof)
      if (eof) go to 99
   else
      call readc (luncrt, &
   'Tangent-slab (t), hemisphere lines (h), or out from body (o)? [<CR>=t]: ', &
                  lunkbd, answer, cr, eof)
      if (eof) go to 99
   end if

   towards_body = answer /= 'O'
   hemisphere   = answer == 'H'

   if (hemisphere) then
       off_center = false
       call ready (luncrt, 'Any off-center body points? [y|n; <CR> = n] ', &
                   lunkbd, off_center, cr, eof)
       if (eof) go to 99
   end if

!  It looks as though all set-up files are ready.  Go ahead with the scripting:

   open  (lunsh, file='prepare_LOS_data', status='unknown')
   write (lunsh, '(a)') '#! /bin/tcsh', &
      '#  Script for preparing line-of-sight data for NEQAIR analysis.'

   write (lunsh,  '(a)') &
      'echo '' ''',      &
      'echo '' Lines-of-sight calculations ...'''

!  LINES_OF_SIGHT[_2D]:
!  !!!!!!!!!!!!!!!!!!!!

   if (.not. hemisphere) then

      filename_los_g = 'los.g'  ! No need for them to be variable (?)
      filename_los_f = 'los.f'

      if (ndim == 2) then
         write (lunsh, '(a)') &
            path // 'LINES_OF_SIGHT_2D/lines_of_sight_2d << LINES_OF_SIGHT'
      else
         write (lunsh, '(a)') &
            path // 'LINES_OF_SIGHT/lines_of_sight << LINES_OF_SIGHT'
      end if

      answer = 'b'
      call reads (luncrt, &
       'Make lines body-normal (b), shock-normal (s), parallel to -Ox (x)?' // &
       ' [<CR> = b] ', lunkbd, answer, cr, eof)
      if (eof) go to 99

      if (formatted_vc) then
         yesno = y
      else
         yesno = n
      end if

      write (lunsh, '(a)') &
         answer,                &
         trim (filename_BP),    &
         trim (filename_grid),  &
         yesno,                 &
         trim (filename_los_g), &
         'LINES_OF_SIGHT',      &
         blank

   else  ! Full angular integration

      if (ndim == 2) then
         call set_up_revolve_grid ()
         if (ios /= 0) go to 99
      else
         if (off_center) then
            call set_up_reflect_blocks ()  ! For half a forebody at nonzero AoA
            if (ios /= 0) go to 99
         end if
      end if

      call set_up_hemispheres_of_sight ()
      if (ios /= 0) go to 99

   end if

!  FLOW_INTERP[_2D]:
!  !!!!!!!!!!!!!!!!!

   if (.not. hemisphere .and. ndim == 2) then  ! No control file

      write (lunsh, '(a)') &
         'echo '' ''',      &
         'echo '' 2-D flow interpolation at LOS points ...''', &
         path // 'FLOW_INTERP_2D/flow_interp_2d << FLOW_INTERP_2D', &
         trim (filename_cc_g),  &
         trim (filename_cc_f),  &
         trim (filename_los_g), &
         trim (filename_los_f), &
         'FLOW_INTERP_2D',      &
         blank

   else  ! 3-D and/or full angular integration

      call write_flow_interp_control_file ()

      write (lunsh, '(a)') &
         'echo '' ''',      &
         'echo '' 3-D flow interpolation at LOS points ...''', &
         path // 'FLOW_INTERP/flow_interp', &
         blank
   end if

!  NEQAIR_DATA:
!  !!!!!!!!!!!!

   call check_existence ('sample.LOS.dat', luncrt, ios)

   if (ios /= 0) then
      write (lunsh, '(a)') &
         'echo '' ''',     &
         'echo '' You still need to invoke neqair_data explicitly.'''
      write (luncrt, '(a)') &
       '*** sample.LOS.dat is missing. ***', &
       'Skip NEQAIR_DATA conversion to NEQAIR''s LOS.dat format for now.', &
       'You can enter "source prepare_LOS_data" now and "neqair_data" later,', &
       'before running NEQAIR.'
   else
      proceed = true
      write (luncrt, '(a)') &
         ' You are ready to convert all lines to NEQAIR''s LOS.dat format.'
      call ready (luncrt, &
                  'Do it now or run NEQAIR_DATA later? [y|n; <CR> = now] ', &
                  lunkbd, proceed, cr, eof)

      if (proceed .and. .not. eof) then

         write (lunsh, '(a)') &
            'echo '' ''',      &
            'echo '' Converting LOS data to NEQAIR format ...''', &
            path // 'NEQAIR_DATA/neqair_data << NEQAIR_DATA', &
            trim (filename_los_g), &
            trim (filename_los_f)
         write (luncrt, '(a)') &
          ' Redistribute the line of sight points?', &
          '    n = no', &
          '    y = yes via dT/ds (gradient-based, normally satisfactory)', &
          '    h = hybrid dT/ds + Vinokur scheme (forebody, extreme T, p)', &
          '    c = curvature-based along the (normalized) T profile', &
          '    r = same relative spacing'
         answer = y
         call reads (luncrt, &
           '   t = simple 2x thinning: [n|y|h|r|t; <CR> = y = dT/ds-based]: ', &
                     lunkbd, answer, cr, eof)
         if (eof) go to 99

         write (lunsh, '(a)') answer
         yes      = answer /= 'n'
         relative = answer == 'r'
         thin2x   = answer == 't'

         if (yes) then
            npts = 0
            if (.not.  thin2x) then
               do while (npts <= 0)
                  npts = 120
                  call readi (luncrt, &
                  'Desired # pts. along output lines of sight? [<CR> = 120] ', &
                              lunkbd, npts, cr, eof)
                  if (eof) go to 99
               end do
               write (lunsh, '(i3)') npts
            end if
            if (.not. (relative .or. thin2x)) then
               power = 1.
               call readr (luncrt, &
             'Initial shape function exponent to use [0. -> 2.; <CR> = 1.]: ', &
                           lunkbd, power, cr, eof)
               if (eof) go to 99
               ismooth = 1
               call readi (luncrt, &
                 'Smoothing control? [0|1|2|3; <CR> = 1 => shape fn. only]: ', &
                           lunkbd, ismooth, cr, eof)
               if (eof) go to 99
               write (lunsh, '(f5.2)') power
               write (lunsh, '(i1)')   ismooth
            end if
         end if

         if (hemisphere .or. towards_body) then
            write (lunsh, '(a)') 'y'
         else  ! Option to integrate out from the body
            write (lunsh, '(a)') 'n'
         end if

         write (lunsh, '(a)') &
            'NEQAIR_DATA',    &
            blank,            &
            'echo '' ''',     &
            'echo '' All done.''', &
            'echo '' You are ready to run NEQAIR on any range of lines.'''

         write (luncrt, '(2a)') &
            'Success.  Next, enter "source prepare_LOS_data" ', &
            'before running NEQAIR.'
      else
         write (lunsh, '(a)') &
            'echo '' ''',     &
            'echo '' You still need to invoke neqair_data explicitly.'''
         write (luncrt, '(2a)') &
            'OK.  Enter "source prepare_LOS_data" then "neqair_data" ', &
            'before running NEQAIR.'
      end if
   end if

   close (lunsh)

99 continue

!  Local procedure for program prepare_neqair_data:

   contains

!     .........................................................................
      subroutine set_up_hemispheres_of_sight ()
!     .........................................................................

      ios = 1
      ne = 25
      call readi (luncrt, &
                  '# hemisphere points, pole to equator? [<CR> = 25] ', &
                  lunkbd, ne, cr, eof)
      if (eof) go to 99

      identifier = 'stag_point'
      call reads (luncrt, &
   'Identifier for the hemisphere-related output files? [<CR> = stag-point] ', &
                  lunkbd, identifier, cr, eof)
      if (eof) go to 99

      write (luncrt, '(a)') &
        'Define the point distribution for the intersected hemisphere lines.', &
        'Remember that NEQAIR_DATA also has the option to redistribute the', &
        'points based on |dTv/ds| and change their number.', &
        'If you intend to do that, more is better here.'
      npts = 129
      call readi (luncrt, &
                'Desired # points along los.g lines of sight? [<CR> = 129] ', &
                  lunkbd, npts, cr, eof)
         if (eof) go to 99

      write (luncrt, '(a)') &
      'LOS discretization is controlled by multiples s1 & s2 of the first and',&
      'last relative grid spacing.  E.g., 1. & 1. (forebody), 2. & 2. (aft).'

      s1 = 1.
      call readr (luncrt, 's1 [<CR> = 1.]: ', lunkbd, s1, cr, eof)
      if (eof) go to 99
      s2 = 1.
      call readr (luncrt, 's2 [<CR> = 1.]: ', lunkbd, s2, cr, eof)
      if (eof) go to 99

      if (formatted_cc) then
         yesno = y
      else
         yesno = n
      end if

      write (lunsh, '(a)') &
         'echo '' ''',      &
         'echo '' Generating hemispherical lines of sight ...''', &
         path // 'HEMISPHERES_OF_SIGHT/hemispheres_of_sight << HEMISPHERE'
      write (lunsh, *) ne
      write (lunsh, '(a)') &
         trim (filename_BP), &
         trim (filename_grid), &
         yesno, &
         trim (identifier)
      write (lunsh, '(i3, /, f8.4, /, f8.4)') npts, s1, s2
      write (lunsh, '(a)') 'HEMISPHERE'

!     Radiometer option to process less than a full hemisphere of lines:

      if (ntokens == 3) then  ! Optional cone angle is not present
         filename_los_g = trim (identifier) // '.1.lines.g'
         filename_los_f = trim (identifier) // '.1.lines.f'
         ios = 0
      else  ! The cone angle appears in the reduced lines of sight file name
         read (buffer, *, iostat=ios) xbp, ybp, zbp, cone_angle
         if (ios /= 0) then
            write (luncrt, '(a, a)') &
               '*** Trouble reading cone_angle from body point line:', &
               buffer
            go to 99
         end if
         filename_los_g(1:11) = 'cone.xx.xx.'
         write (filename_los_g(6:10), '(f5.2)') cone_angle
         filename_los_g(12:)  = trim (identifier) // '.1.lines.g'
         filename_los_f(1:11) = filename_los_g(1:11)
         filename_los_f(12:)  = trim (identifier) // '.1.lines.f'
      end if

 99   return

      end subroutine set_up_hemispheres_of_sight

!     .........................................................................
!
      subroutine set_up_reflect_blocks ()
!
!     For the case of full angular integration at off-center body points on
!     half a forebody at AoA, set up reflection of both the half volume grid
!     and the half volume cell-centered data, ready for HEMISPHERES_OF_SIGHT.
!     Unsteady wakes are not symmetric, so this is inappropriate for full
!     bodies.
!     .........................................................................

      proceed = true
      call ready (luncrt, &
  'Are the volume data files whole bodies, not half-bodies? [y|n; <CR> = y] ', &
                  lunkbd, proceed, cr, eof)
      if (eof .or. .not. proceed) then
         write (luncrt, '(2a)') &
            'Reflecting of half-forebody data at angle of attack ', &
            'has not been automated.', &
      'Apply REFLECT_BLOCKS explicitly first and run PREPARE_NEQAIR_DATA again.'
         ios = 1
      else
         ios = 0
      end if

      end subroutine set_up_reflect_blocks

!     .........................................................................
      subroutine set_up_revolve_grid ()
!     .........................................................................

      integer :: npts_rotate
      real    :: total_angle

      ios = 1

      write (luncrt, '(a)') &
  'The 2-D datasets need to be rotated 180 or 360 deg. for angular integration.'

      total_angle = 180.
      call readr (luncrt, 'Total angle of revolution? [<CR> = 180 deg.] ', &
                  lunkbd, total_angle, cr, eof)
      if (eof) go to 99

      npts_rotate = 181
      call readr (luncrt, '# circumferential grid points? [<CR> = 181]  ', &
                  lunkbd, npts_rotate, cr, eof)
      if (eof) go to 99

      if (formatted_vc) then
         yesno = y
      else
         yesno = n
      end if

      write (lunsh, '(a)') &
         'echo '' ''',      &
         'echo '' Revolving the 2-D volume grid ...''', &
         path // 'REVOLVE_GRID/revolve_grid << REVOLVE_GRID'
      write (lunsh, '(a)') &
         trim (filename_grid), &
         yesno, &
         'n'
      write (lunsh, '(f4.0, /, i3)') total_angle, npts_rotate
      write (lunsh, '(a)') &
         file_rev_g, 'n', &
         'REVOLVE_GRID'

      filename_grid = file_rev_g
      formatted_vc  = false

!     Repeat for the cell-centered grid and flow data:

      if (formatted_cc) then
         yesno = y
      else
         yesno = n
      end if

      write (lunsh, '(a)') &
         'echo '' ''',      &
         'echo '' Revolving the cell-centered 2-D grid and flow data ...''', &
         path // 'REVOLVE_GRID/revolve_grid << REVOLVE_CC'
      write (lunsh, '(a)') &
         trim (filename_cc_g), &
         yesno, &
         'y', &
         trim (filename_cc_f)
      write (lunsh, '(f4.0, /, i3, /, i1)') total_angle, npts_rotate, 0
      write (lunsh, '(a)') &
         file_rev_cc_g, 'n', &
         file_rev_cc_f, &
         'REVOLVE_CC'

      filename_cc_g = file_rev_cc_g
      filename_cc_f = file_rev_cc_f
      formatted_cc  = false

      ios = 0

 99   return

      end  subroutine set_up_revolve_grid

!     .........................................................................
      subroutine write_flow_interp_control_file ()
!     .........................................................................

      open  (lunout, file='flow_interp.inp', status='unknown')
      write (lunout, '(a)') &
         'FLOW_INTERP controls', &
         '------------- INPUT SOLUTION GRID -------------', &
         trim (filename_cc_g)
      write (lunout, '(l1)') formatted_cc
      write (lunout, '(a)') &
         '------------- INPUT FLOW SOLUTION -------------', &
         trim (filename_cc_f)
      write (lunout, '(l1)') formatted_cc
      write (lunout, '(a)') &
         '------------- TARGET GRID ---------------------', &
         trim (filename_los_g)
      write (lunout, '(l1)') true  ! Formatted
      write (lunout, '(a)') &
         '------------- INTERPOLATED FLOW SOLUTION ------', &
         trim (filename_los_f)
      write (lunout, '(l1)') true  ! Formatted
      write (lunout, '(a)') &
         '------------- MISCELLANEOUS CONTROLS ----------', &
         '0.0001        Search in/out tolerance', &
         '------------- OPTIONAL CONTROLS ---------------', &
         'T             Suppress soln. blocks if possible', &
         '1             Plain ADT method'
      close (lunout)

      end subroutine write_flow_interp_control_file

   end program prepare_neqair_data
