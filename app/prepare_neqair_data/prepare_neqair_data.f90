!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program prepare_neqair_data

!     This adaptation of prepare_rapid_analysis prompts the user for specifics
!     about a flow solution and the body point(s) at which radiative heating
!     calculations are to be performed with NEQAIR, either in tangent-slab mode
!     or full angular integration mode.  STRUCTURED GRIDS IN PLOT3D FORMAT ARE
!     ASSUMED.  For the unstructured case, USLOS is available for constructing
!     individual or hemispherical collections of lines of sight, but at the time
!     of writing, interpolation of unstructured flow field data onto those lines
!     of sight remains an issue, because US3D's nearest-cell (KDTREE-based) form
!     of interpolation can choose wrong cell centroids in the presence of high-
!     aspect ratio cells near a flow-field shock.
!
!     The working directory should contain (or point to) the vertex-centered
!     volume grid and associated flow data (temperature(s) & number densities).
!     With the utilities that preceded SLOS, cell-centered volume grid and flow
!     data could also be used, but the vertex-centered grid (without the fn.
!     file) was still needed because the line-of-sight discretizations were
!     derived from its off-wall grid lines.  In SLOS, the discretizations follow
!     the two-stage scheme of DPLR (1-sided stretching with ds1, then 2-sided
!     restretching with ds2mult applied to the 1-sided outermost interval).
!     Use of VERTEX-CENTERED grid and function files is now recommended, but
!     cell-centered data with "halo" cells will produce much the same result,
!     as only the innner and outer layers at k = 1 and kmax are actually used.
!     All these files should be in PLOT3D multiblock form, formatted or not.
!     Whether they are 2D or 3D is automatically determined.
!
!     If full integration is specified for a 2D axisymmetric case, 3D revolved
!     grid and function files will be produced as part of the automation
!     procedures of which the present utility is the first step following a
!     flow calculation with the right choice of species.
!
!     The output from this utility is a shell script named prepare_LOS_data that
!     (when sourced) invokes the indicated utilities and control files (also
!     generated here) to set up line-of-sight data in NEQAIR's LOS.dat form as
!     LOS-1.dat, LOS-2.dat, ..., LOS-n.dat in the working directory.
!
!     After the line(s) of sight have been checked for possible flaws (such as
!     in the underlying line-surface intersection calculations), a separate
!     script, set_neqair should be invoked with two first and last line number
!     arguments.  This generates directories LINE-1, LINE-2 ... containing links
!     to LOS-1.dat, LOS-2.dat ... as the LOS.dat expected by NEQAIR, and links
!     to the appropriate neqair.inp control file, and to a neqair.pbs file.
!     Another shell script can then launch NEQAIR within each LINE-* directory,
!     and this script is system-dependent.
!
!     Note that on the NAS supercomputer facility at Ames Research Center, an
!     "array" job is recommended.  This automates looping over a range of
!     indices (which can be the n in the LINE-n directories), with a limit of
!     500 total at the time of writing.  E.g., neqair.1-325.pbs:
!
!        #PBS -W group_list=xxxxx
!        #PBS -q normal
!        #PBS -l walltime=0:10:00
!        #PBS -l select=2:ncpus=20:model=ivy
!        #PBS -N neqair
!        #PBS -j oe
!        #PBS -J 1-325
!
!        module use -a /home4/tsa/modulefiles
!        module load neqair/15.0
!
!        cd LINE-$PBS_ARRAY_INDEX
!        mpiexec neqair < neqair.inp  > neqair.out
!
!     Files neqair.inp and neqair.pbs are expected to be in the working
!     directory along with the body point file (indices or coordinates) and a
!     sample.LOS.dat file indicating the correct species, although the user
!     can choose to proceed if any of neqair.inp, neqair.pbs, or sample.LOS.dat
!     are missing.  The neqair_data step (los.* --> LOS.dat format) can then
!     be performed explicitly at a later time.
!
!     This neqair.pbs should not be confused with the array job (or alternative)
!     shown above.  It should be set to run a single line only, and might be
!     used in /LINE-1 as a check that all looks as it should before launching a
!     job that loops over many lines.
!
!     More about sample.LOS.dat (see NEQAIR_DATA):
!
!        NEQAIR_DATA can now write line-of-sight (LOS) data in both the original
!        rigid format read by NEQAIR versions up to and including v14 and the
!        simpler column-oriented format handled from NEQAIR v15 on.  Both data
!        formats involve a sample.LOS.dat file read here that indicates the
!        species names.  If a sample file is found with a NEQAIR version number
!        on line 1, then the later data format is indicated. The sample file is
!        transcribed as the header of all output LOS data files in either case.
!
!     When NEQAIR is run on a range of lines, its results will appear in
!     neqair.out in subdirectores /LINE-1, /LINE-2, ...  A command such as
!     grep Total LINE-*/neqair.out > qrad.dat followed by wc qrad.dat will
!     confirm that all lines have run or that NEQAIR may have failed on some.
!     The utility sort_rows is one way of ordering qrad.dat, using LINE- as
!     the sort string, although grep may have a switch that overrides the
!     default order of LINE-1, LINE-10, LINE-100, LINE-101, ...
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
!     Originally:
!        lines_of_sight       ! One or more lines of sight in a 3D volume grid
!        lines_of_sight_2d    ! One or more lines of sight in a 2D volume grid
!        hemispheres_of_sight ! Many lines of sight for one body point
!        flow_interp          ! 3D flow interpolation at LOS coordinates
!        flow_interp_2D       ! 2D   "    "    "    "    "    "    "
!        neqair_data          ! PLOT3D .g/.f --> NEQAIR's LOS.dat form
!
!     Currently:
!        slos                 ! 1 or more lines, tangent-slab or hemi., 3D vol.
!        lines_of_sight_2d    ! 1 or more tangent-slab lines in a 2D volume grid
!        flow_interp          ! 3D flow interpolation at LOS coordinates
!        flow_interp_2D       ! 2D   "    "    "    "    "    "    "
!        neqair_data          ! PLOT3D .g/.f --> NEQAIR's LOS.dat form
!
!        [The more recent SLOS utility was derived from USLOS (for unstructured
!        grids).  Both of these combine the functions of the LINES_OF_SIGHT and
!        HEMISPHERES_OF_SIGHT utilities, and make use of the diagonal of the
!        volume grid bounding box for the line-surface intersection search
!        intervals in place of earlier heuristics.  They also feature a more
!        fail-safe (?) way of recovering from an intersection calculation that
!        has failed because a line-surface distance isn't essentially zero.]
!
!  FLOW_INTERP Warning:
!
!        Past practice has been to use FLOW_INTERP's hybrid method 3, because
!        its original preferred method occasionally encounters a singular matrix
!        believed to be associated with a high aspect ratio cell in a boundary
!        layer.  The ADT search method now traps singularities and retries with
!        a perturbed system guaranteed to be full rank.  Method 3 (nearest cell
!        centroid search via KDTREE, followed by refinement within that cell)
!        has been found to give seriously wrong interpolations occasionally by
!        using cells with centroids outside the shock even though the target
!        points are actually inside the shock.  This explains occasional NEQAIR
!        failures:  low free-stream temperatures were being interpolated on to
!        the line of sight inside the shock.  USE FLOW_INTERP METHOD 1 WHENEVER
!        PRACTICAL; METHODS 2 OR 3 SHOULD BE CONFINED TO LAST RESORTS FOR REALLY
!        LARGE VOLUME GRIDS.
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
!     11/01/19   "     "      Invoke SLOS in place of the earlier LINES_OF_SIGHT
!                             and HEMISPHERES_OF_SIGHT.  See above for why.
!                             The 2D case (LINES_OF_SIGHT_2D) is still supported
!                             as is the option to revolve 2D data if hemi-
!                             spherical lines are required.  The option to
!                             reflect half-body data remains incomplete.  The
!                             documentation above has also been expanded to be
!                             as clear as possible.  These revisions proved
!                             much messier than expected, partly because of
!                             the dropping of certain options by SLOS.
!     11/06/19   "     "      The 2D/revolve case wasn't right.  The revolving
!                             has to happen before SLOS begins.
!     11/07/19   "     "      Warn about 2D/3D body point gotcha.
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
      mxchar = 128          ! File name and body pt. x/y/z[/cone] buffer limit

   logical, parameter :: &
      false  = .false.,  &
      true   = .true.

!  The following path indicates where PREPARE_NEQAIR_DAT, LINES_OF_SIGHT, etc.
!  are expected to be installed.  Be sure to get its length right.

   character, parameter :: &
!!!   path*21          = '/share/apps/cfdtools/', &
!!!   path*20          = '/apps/pkgs/cfdtools/', &
      path*16          = '/home5/dasaund1/', &
      blank*1          = ' ',   &
      n*1              = 'n',   &
      y*1              = 'y',   &
      file_los_f*5     = 'los.f', &
      file_los_g*5     = 'los.g', &
      file_rev_flow*11 = 'revolved.fu', &
      file_rev_grid*11 = 'revolved.gu'

!  Variables:

   integer :: &
      ier, ios, ismooth, ndim, nbps, ne, npts, ntokens

   real :: &
      cone_angle, ds1, ds2_mult, power, xbp, ybp, zbp

   logical ::   &
      any_angle, body_normal, cr, eof, formatted_in, formatted_out, &
      hemisphere, off_center, parallel_to_Ox, proceed, relative, shock_normal, &
      thin2x, towards_body, yes

   character (1) :: &
      answer, yesno

   character (mxchar) :: &
      filename_bp, filename_flow, filename_grid, filename_los_f, &
      filename_los_g, buffer

!  Execution:
!  !!!!!!!!!!

   write (luncrt, '(/, (a))') &
      ' Requirements summary:', ' ', &
      '    Grids and solutions may be 2D or 3D, formatted or not.', &
      '    3D volume grid and flow data may be vertex- or cell-centered, as',  &
      '    only the inner and outer layers are actually used to construct',    &
      '    the lines of sight.',                                               &
      '    2D volume data should be vertex-centered, as the radial grid line', &
      '    distributions determine the discretization of the lines of sight.', &
      ' ', &
      '    DPLR''s POSTFLOW should have been used for two volume files:',      &
      '     (1&2) vertex-centered volume grid + flow data (interp = 2) or' //  &
      ' (3D only)', &
      '           cell-centered grid + flow data with halo cells (interp = 1)',&
      '           ivarp = 0 120 125 1200  (1-4 temperatures, NEQAIR 15+) or',  &
      '           ivarp = 0 120 120 125 125 1200  (4 temperatures, NEQAIR 14)',&
      ' ', &
      '    The working directory should also contain or point to these files:',&
      '       (3) body pt. file (x/y/z [+ optional cone angle for hemi case]', &
      '       (4) neqair.inp      (tangent-slab or not, regions, etc.)',       &
      '       (5) neqair.pbs      (# nodes, desired NEQAIR version, etc.)',    &
      '       (6) set_neqair      (first/last LOS numbers as arguments)',      &
      '       (7) sample.LOS.dat  (indicates NEQAIR version & species names)'
   write (luncrt, '(/, 3a, //, 2a, /, a, /)') &
      ' Examples of files (4)-(7) should be found here: ', path, &
      'PREPARE_NEQAIR_DATA', &
      ' Be sure to enter the right species in sample.LOS.dat and to', &
      ' check neqair.inp and neqair.pbs.', &
      ' Plotting the los.g lines before running NEQAIR is advised.', &
      ' Use DPLR''s *.chem file to determine the right species names.'

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

!  Open the shell script to be generated:

   open  (lunsh, file='prepare_LOS_data', status='unknown')
   write (lunsh, '(a)') &
      '#! /bin/tcsh', &
      '#  Script for preparing line-of-sight data for NEQAIR analysis.', &
      'echo '' ''',      &
      'echo '' Line(s) of sight calculations ...'''

!  (1) Vertex- or cell-centered volume grid (3D); vertex-centered (2D):

   filename_grid = 'volrad.gu';  ios = 1
   do while (ios /= 0)
      call reads (luncrt, &
         'Volume grid? (Vertex-centered for 2D, vc or cc for 3D)? [<CR> = ' &
         // trim (filename_grid) // ']: ', lunkbd, filename_grid, cr, eof)
      if (eof) go to 99

      call determine_grid_form (filename_grid, lunin, formatted_in, ios)
   end do

!  2- or 3-D?

   call determine_grid_dim (filename_grid, lunin, formatted_in, ndim, ios)

   write (luncrt, '(a, i2)') ' Apparent volume grid dimension:', ndim

!  (2) Vertex- or cell-centered flow data (3D); vertex-centered (2D):

   filename_flow = 'volrad.fu';  ios = 1
   do while (ios /= 0)
      call reads (luncrt, &
         'Flow volume data file?                                  [<CR> = ' &
         // trim (filename_flow) // ']: ', lunkbd, filename_flow, cr, eof)
      if (eof) go to 99

      call determine_grid_form (filename_flow, lunin, formatted_in, ios)
   end do

!  (3) Body point file:

   filename_bp = 'bp.inp';  ios = 1
   do while (ios /= 0)
      call reads (luncrt, &
            'Target body point file? [<CR> = ' // trim (filename_bp) // ']: ', &
                  lunkbd, filename_bp, cr, eof)
      if (eof) go to 99

      open (lunin, file=filename_bp, status='old', iostat=ios)
   end do

!  Count the body points:

   nbps = 0
   do
      read (lunin, '(a)', iostat=ios) buffer
      if (ios /= 0) exit
      nbps = nbps + 1
   end do
   rewind (lunin)

   write (luncrt, '(/, a, i4)') ' # body points specified: ', nbps

!  More than one body point means this is not a hemisphere case.
!  What type of line(s) of sight?

   towards_body = true       ! But NEQAIR_DATA can reverse this if told to
   hemisphere   = nbps == 1  ! But allow for a single LOS

   if (hemisphere) then
      call ready (luncrt, &
     'One body point found. Do you want hemispherical lines? [y/n; <cr>=y]: ', &
                  lunkbd, hemisphere, cr, eof)
     if (eof) go to 99
   end if

   if (hemisphere) then
      if (ndim == 2) then
         proceed = true
         call ready (luncrt, &
                '*** WARNING: The body pt. file must have 3 coordinates. ' // &
                     'Continue?  [y|n=^D=quit; y=<cr>=continue]: ', &
                     lunkbd, proceed, cr, eof)
         if (eof .or. .not. proceed) go to 99
      end if
      body_normal    = false
      shock_normal   = false
      parallel_to_Ox = false
      any_angle      = false
      off_center     = false
      call ready (luncrt, 'Any off-center body points? [y|n; <CR> = n] ', &
                  lunkbd, off_center, cr, eof)
      if (eof) go to 99

      ne = 25
      call readi (luncrt, &
        '# hemisphere points, pole to equator? [<CR> = 25 = 650/1300 lines] ', &
                  lunkbd, ne, cr, eof)
      if (eof) go to 99
   else  ! Tangent-slab lines, or other possibilities from the past
      write (luncrt, '(/, (a))') &
          ' Make the line(s) body-normal                          (b)', &
          '                 shock-normal                          (s)', &
          '              parallel to -Ox                          (x)', &
          ' or at the angle(s) with Ox found in the body pt. file (a)?'
      answer = 'B'
      call readc (luncrt, 'Line type? [b|s|x|a; <cr>=b=body-normal]: ', &
                  lunkbd, answer, cr, eof)
      body_normal    = answer == 'B'
      shock_normal   = answer == 'S'
      parallel_to_Ox = answer == 'X'
      any_angle      = answer == 'A'
      call ready (luncrt, &
              'Should NEQAIR integrate towards the body? [y|n; <CR> = yes]: ', &
               lunkbd, towards_body, cr, eof)
      write (luncrt, '(a, l2)') 'Towards body:', towards_body
      if (eof) go to 99
   end if

!  Line discretization details (except for the 2D case), which still reuses the
!  volume grid's relative discretization of off-wall radial lines:

   if (hemisphere .or. ndim == 3) then

      write (luncrt, '(a)') &
        'Define the point distribution for the intersected hemisphere lines.', &
        'Remember that NEQAIR_DATA also has the option to redistribute the', &
        'points based on |dTv/ds| and/or change their number.', &
        'If you intend to do that, retain the flow resolution here.', &
        'Preserving full resolution is probably preferable to coarsening.'
      npts = 121
      call readi (luncrt, &
                'Desired # points along los.g lines of sight? [<CR> = 121] ', &
                  lunkbd, npts, cr, eof)
         if (eof) go to 99

      write (luncrt, '(a)') &
         'LOS discretization is controlled as in DPLR'' grid alignment:', &
      'wall spacing ds1 with 1-sided-stretching; ds2_mult applied to outer ds2.'

      ds1 = 1.e-5
      call readr (luncrt, 'ds1 [<CR> = 1.e-5]: ', lunkbd, ds1, cr, eof)
      if (eof) go to 99
      ds2_mult = 0.1
      call readr (luncrt, 'ds2_mult [<CR> = 0.1]: ', lunkbd, ds2_mult, cr, eof)
      if (eof) go to 99

   end if

   if (hemisphere) then

      if (off_center) then
         call set_up_reflect_blocks ()  ! Not complete
         if (ios /= 0) go to 99
      end if

      if (ndim == 2) then
         call set_up_revolve_grid ()
         if (ios /= 0) go to 99
      end if

   end if

   filename_los_g = 'los.g'  ! No need for them to be variable (?)
   filename_los_f = 'los.f'
   formatted_out  = true

!  Lines of sight (LINES_OF_SIGHT_2D or SLOS):
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (hemisphere) then

      write (lunsh, '(a)') &
         path // 'SLOS/slos << LOS', &
         trim (filename_bp), &
         y, &     ! "Yes" to the question about hemisphere|not hemi. (nbps = 1)
         trim (filename_grid)
      write (lunsh, '(i3, /, i3, /, es10.3, /, f6.3)') ne, npts, ds1, ds2_mult
      write (lunsh, '(a)') 'LOS'

   else ! Not hemispherical

      if (ndim == 2) then  ! LINES_OF_SIGHT_2D
         yesno = n
         if (formatted_in) yesno = y
         write (lunsh, '(a)')      &
            path // 'LINES_OF_SIGHT_2D/lines_of_sight_2d << LOS', &
            answer,                &  ! Type of line (B|S|X|A)
            trim (filename_bp),    &
            trim (filename_grid),  &
            yesno,                 &  ! Formatted or not
            trim (filename_los_g)
      else  ! SLOS and not hemispherical
         write (lunsh, '(a)') &
            path // 'SLOS/slos << LOS', &
            trim (filename_bp)
         if (nbps == 1) write (lunsh, '(a)') n  ! No the question hemi? question
         write (lunsh, '(a)') trim (filename_grid)
         write (lunsh, '(i3, /, es10.3, /, f6.3)') npts, ds1, ds2_mult
      end if
      write (lunsh, '(a)') 'LOS', blank

   end if

!  FLOW_INTERP[_2D]:
!  !!!!!!!!!!!!!!!!!

   if (.not. hemisphere .and. ndim == 2) then  ! No control file

      write (lunsh, '(a)')  &
         'echo '' ''',      &
         'echo '' 2-D flow interpolation at LOS points ...''', &
         path // 'FLOW_INTERP_2D/flow_interp_2d << FLOW_INTERP_2D', &
         trim (filename_grid),  &
         trim (filename_flow),  &
         file_los_g, &
         file_los_f, &
         'FLOW_INTERP_2D', &
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
            file_los_g, &
            file_los_f
         write (luncrt, '(a)') &
          ' Redistribute the line of sight points?', &
          '    n = no', &
          '    y = yes via dT/ds (gradient-based, unnormalized)', &
          '    h = hybrid dT/ds + Vinokur scheme (forebody, extreme T, p)', &
          '    c = curvature-based along the (normalized) T profile', &
          '    r = same relative spacing'
         answer = n
         call reads (luncrt, &
         '   t = simple 2x thinning: [n|y|h|r|t; <CR>=n=no redistribution]: ', &
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
                  npts = 121
                  call readi (luncrt, &
                  'Desired # pts. along output lines of sight? [<CR> = 121] ', &
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
            write (lunsh, '(a)') y
         else  ! Option to integrate out from the body
            write (lunsh, '(a)') n
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

!  Local procedures for program prepare_neqair_data:

   contains

!     .........................................................................
!
      subroutine set_up_reflect_blocks ()
!
!     For the case of full angular integration at off-center body points on
!     half a forebody at AoA, set up reflection of both the half volume grid
!     and the half volume flow data, ready for SLOS.  Unsteady wakes may not
!     be symmetric, so this may be inappropriate for full bodies.
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

      npts_rotate = 91
      call readr (luncrt, '# circumferential grid points? [<CR> = 91]  ', &
                  lunkbd, npts_rotate, cr, eof)
      if (eof) go to 99

      write (lunsh, '(a)') &
         'echo '' ''',      &
         'echo '' Revolving the 2-D volume data ...''', &
         path // 'REVOLVE_GRID/revolve_grid << REVOLVE'
      write (lunsh, '(a, /, l1, /, a, /, a)') &
         trim (filename_grid), &
         formatted_in, &
         y, &
         trim (filename_flow)

      write (lunsh, '(f4.0, /, i3, /, i1)') total_angle, npts_rotate, 0
      write (lunsh, '(a)') &
         file_rev_grid, &
         n, &
         file_rev_flow, &
         'REVOLVE'

      filename_grid = file_rev_grid  ! For the FLOW_INTERP control file
      filename_flow = file_rev_flow
      formatted_in  = false
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
         trim (filename_grid)
      write (lunout, '(l1)') formatted_in
      write (lunout, '(a)') &
         '------------- INPUT FLOW SOLUTION -------------', &
         trim (filename_flow)
      write (lunout, '(l1)') formatted_in
      write (lunout, '(a)') &
         '------------- TARGET GRID ---------------------', &
         trim (filename_los_g)
      write (lunout, '(l1)') formatted_out
      write (lunout, '(a)') &
         '------------- INTERPOLATED FLOW SOLUTION ------', &
         trim (filename_los_f)
      write (lunout, '(l1)') formatted_out
      write (lunout, '(a)') &
         '------------- MISCELLANEOUS CONTROLS ----------', &
         '0.0001        Search in/out tolerance', &
         '------------- OPTIONAL CONTROLS ---------------', &
         'T             Suppress soln. blocks if possible', &
         '1             Plain ADT method; use 2|3 judiciously'
      close (lunout)

      end subroutine write_flow_interp_control_file

   end program prepare_neqair_data
