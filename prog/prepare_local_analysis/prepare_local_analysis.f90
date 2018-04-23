!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program prepare_local_analysis

!     This is a variant of prepare_rapid_analysis with one extra prompt for the
!     number of species in the flow solution.  Renaming it avoids affecting
!     existing automation scripts.  This variant was prompted by turbulent flow
!     cases, where it is no longer feasible to handle 1 or 2 temperatures in air
!     without the prompt that affects existing automation.
!
!  Site-Specifics:
!
!     If this program is installed at another site, the following path should
!     be edited appropriately prior to compiling and linking:
!
!        /share/apps/cfdtools/         ! See parameter constant "path" below
!
!     It is assumed that all the utilities invoked here reside in subdirectories
!     with uppercase names and this indicated path.
!
!     E.g., thin_grid, adjust_grid, ... should be installed here:
!
!        /share/apps/cfdtools/THIN_GRID/thin_grid
!        /share/apps/cfdtools/ADJUST_GRID/adjust_grid   and so on.
!
!     The following utilities may be invoked by the script from this program:
!
!        thin_grid           ! Option to coarsen new input surface/cavity grids
!        adjust_grid         ! Option to convert input inches to meters
!        combine_blocks_turb ! Option to append/init. cavity blocks (lam|turb)
!        template            ! Generates dplr interface & control files
!        plug_interp         ! Specialized treatment of a wing leading edge plug
!
!  Description (mostly from the original prepare_rapid_analysis):
!
!     Program prepare_local_analysis generates a shell script, named
!     local_analysis_script, and control files to run the NASA Ames tools for
!     rapid CFD analysis of damage configurations (or other perturbations -
!     often a tile cavity) with minimal user interaction.  The intended flow
!     solver is DPLR.
!
!     The baseline grid is assumed to have a single layer of blocks, with
!     k = 1 at the geometry surface.
!
!     Another shell script may be used to run prepare_local_analysis, to make
!     its output script executable, and to invoke that script.  At NASA ARC, see
!
!        /share/apps/cfdtools/PREPARE_LOCAL_ANALYSIS/local_analysis.sh
!
!     which is also pointed to by a symbolic link that should be in your path:
!
!        /share/apps/cfdtools/bin/local_analysis
!
!     Then typing "local_analysis" performs everything between the Gridgen and
!     FCONVERT steps via answers to a few prompts.
!
!     Alternatively, run prepare_local_analysis directly, then use
!
!        source local_analysis_script
!
!     to execute its output script.
!
!  Usual sequence of steps:
!
!     Given a baseline grid (vertex-centered) and a flow solution with the
!     cell-centered form of the grid, including "halo" cells, plus a new
!     surface grid (probably confined to the locality of a surface perturb-
!     ation) and (commonly) an associated cavity grid block or plug/gap
!     filler collar blocks, the usual sequence is:
!
!     1:  Thin the new surface and volume grids (THIN_GRID with 2 2 1 inputs).
!
!     2:  Convert the thinned grids from inches to meters (ADJUST_GRID).
!
!     3:  Generate a volume grid from the new surface and baseline volume grids
!         (first application of RADIAL_INTERP, with vertex-centered data).
!
!     4:  Interpolate the flow solution to this volume grid (second application
!         of RADIAL_INTERP, with cell-centered data but same new surface as 3).
!
!     5:  Combine the main volume grid from 3 with the ancillary block(s) (first
!         application of COMBINE_BLOCKS, with vertex-centered data).
!
!     6:  Combine the main flow starting guess with estimates for the ancillary
!         block(s) calculated as a specialized option by a second application of
!         COMBINE_BLOCKS (cell-centered data except for the cavity grid).
!
!     7:  Generate DPLR control files (dplr.inputs and dplr.interfaces) using
!         TEMPLATE.  A sample.inputs file should be present so that the header
!         and trailer records can be transcribed to dplr.inputs.  Some editing
!         of dplr.inputs may still be needed:  the correct number of blocks, the
!         correct flight conditions, and the appropriate CFL schedule.
!
!  N.B.:  Note that TEMPLATE can now avoid the need to edit certain BCs for the
!         cavity blocks (and different BCs for blocks around a plug repair or a
!         protruding gap filler) via the 'template.inp.2' ancillary file.
!         Ideally, the present program would write such a file based on how many
!         blocks are in the (optional) second grid file.  Instead, if the file
!         is present, this program can use it to distinguish between a cavity
!         and a plug/gap filler case.  Otherwise, it assumes it is neither -
!         probably just a new surface grid.
!
!         Note that gap filler cases are treated as for plug cases, and should
!         contain 'plug' in template.inp.2 (case-insensitive).
!     
!     Two runs of FCONVERT could in principle be added, but how to split the
!     blocks may need to be decided at run time according to how many CPUs are
!     deemed appropriate.
!
!  History:
!
!     11/28/05  DAS  Initial step towards automating rapid cavity analysis.
!     11/29/05   "   Added the THIN_GRID, SCALE_GRID and TEMPLATE options.
!     01/09/06   "   Ryan McDaniel ran into a glitch in the no-cavity case,
!                    and suggested a prologue to summarize requirements.
!     03/31/06   "   Remind the user about editing the cavity wall BCs, and
!                    default the baseline volume grid and solution files.
!     04/14/06   "   Added handling of the leading edge plug repair case
!                    and exploited the new TEMPLATE option to make use of
!                    a 'template.inp.2' control file.
!     04/17/06   "   Fixed a couple of typos in the prologue.
!     08/04/06   "   RADIAL_INTERP's distance tolerance should be negative,
!                    meaning it will decay any surface discrepancies to
!                    zero at the outer boundary.
!     10/24/06   "   Changed the default Shuttle baseline file names from
!                    "right half" to "both halves".  Documented the fact
!                    that protruding gap filler cases are treated as for
!                    plug cases at this level.
!     07/11/08   "   The output script now begins with #! /bin/tcsh.
!     10/15/08   "   It now removes any RADIAL_INTERP log files that might
!                    be present from an earlier run.
!     10/19/08   "   A Todd White case broke the assumption that the new 
!                    surface and volume files are formatted the same way.
!                    That is no longer assumed.
!     09/09/09   "   PREPARE_LOCAL_ANALYSIS is just PREPARE_RAPID_ANALYSIS
!                    with one extra prompt for cavity-type cases (number of
!                    species in the flow solution).  Renaming it avoids
!                    impacting automation scripts.
!     09/11/09   "   The prologue shouldn't mention RADIAL_INTERP control 
!                    files because these are generated here once all the 
!                    prompts have been answered.
!     02/13/10   "   Decaying surface discrepancies to zero at the outer
!                    boundary is found to be the wrong choice for a WLE
!                    plug (working with an offset surface in a high- 
!                    curvature region).  Only shallow cavities really
!                    need the decay, so suppress it here via the +0.0002 
!                    inputs to RADIAL_INTERP that had been -0.0002.
!     01/18/12   "   Installation on Redwood-2 prompted use of a path to 
!                    the utilities.  But all of them have uppercase
!                    directories now, instead of the earlier mixture.
!     07/19/12   "   An extra prompt allows the decay/no-decay choice now.
!     01/13/14   "   ADJUST_GRID is preferable to SCALE_GRID now, so the
!                    option to convert inches to meters has been recoded.
!
!  Author:  David Saunders, ELORET Corp./NASA Ames Research Center, CA
!                           Now ERC, Inc./NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Constants:

   integer, parameter :: &
      lunin  = 1,        &  ! For checking if a file formatted o unformatted
      lunout = 2,        &  ! Output control files
      lunsh  = 3,        &  ! For the output shell script
      lunctl = 4,        &  ! For the (optional) template.inp.2 control file
      lunkbd = 5,        &  ! For keyboard entries
      luncrt = 6            ! For prompts

   logical, parameter :: &
      false  = .false.,  &
      true   = .true.

   character, parameter :: &
      path * 21 = '/share/apps/cfdtools/', &  ! Installation-dependent path to
      blank * 1 = ' ',   &                    ! PREPARE_LOCAL_ANALYSIS subdir.
      n * 1     = 'n',   &
      y * 1     = 'y'

!  Variables:

   integer :: &
      i, ios, nspecies

   logical ::   &
      cavity_case, cr, decay, eof, formatted, formatted_surface, &
      formatted_volume, proceed, plug_case, scale_grid, thin_grid

   character :: &
      answer * 1, filename * 64, filename_cavity * 64, filename_flow * 64,     &
      filename_surface * 64, gridname * 64, keyword * 4

!  Execution:
!  !!!!!!!!!!

   write (luncrt, '(/, (a))') &
      ' Requirements summary:', ' ', &
      '    Grids and solutions may be formatted or not.', &
      '    Thinning and scaling of input files are provided for.', &
      '    If a flow solution is being interpolated, it must be cell-centered',&
      '    with halo cells, and must be accompanied by a cell-centered grid.', &
      '    DPLR''s POSTFLOW utility should have been run already.', ' ',       &
      '    The appropriate 2-temperature inputs to POSTFLOW are:',             &
      '             ivarp = 0 1000 150 151 152 120 [125]',                     &
      '             Omit 125 for 1-temperature Shuttle applications.', ' ',    &
      '             The number of species need not be 5.',                     &
      '             Turbulent solutions should add 90 [+ 70 71 for SST].',' ', &
      '    The RADIAL_INTERP log files are scanned here for acceptable cell',  &
      '    search statistics.  Interrupt proceedings if large numbers of',     &
      '    points appear out of tolerance (probably mismatched units).', ' ',  &
      '    Specialized appending of cavity blocks is optional.',               &
      '    These may instead be collar blocks around plug/gap filler blocks.', &
      ' ',  &
      '    TEMPLATE is the last automated step, but its results will not be',  &
      '    valid if all grid blocks are not corner-to-corner matched.',        &
      '    Even if subfacing is present, TEMPLATE will warn of bad cells, so', &
      '    let it run anyway.', ' ',                                           &
      '    TEMPLATE also looks for two optional files:', ' ',                  &
      '      ''sample.inputs''  will have its header and trailer transcribed', &
      '                       to the body of ''dplr.inputs'' produced here;',  &
      '      ''template.inp.2'' allows cavity/plug/gap filler BCs to be set',  &
      '                       and also allows control of iseq, jseq, kseq.', ' '

   proceed = true
   call ready (luncrt, 'Proceed?  [y|n|^D=quit|<CR>=proceed]: ',               &
               lunkbd, proceed, cr, eof)
   if (eof) go to 99
   if (.not. proceed) go to 99

   open  (lunsh, file='local_analysis_script', status='unknown')

   write (lunsh, '(a)') '#! /bin/tcsh', &
    '#  Script for preparing a local CFD analysis of a perturbed configuration.'
   write (lunsh,  '(a)')
   write (lunsh,  '(a)') '\rm -f vertices.log centers.log'
   write (lunsh,  '(a)')
   write (luncrt, '(a)')

!  If TEMPLATE's optional control file is present, make possible use of it:

   cavity_case = false
   plug_case   = false

   open (lunctl, file='template.inp.2', status='old', iostat=ios)

   do while (ios == 0)

!     Look for 'CAVITY' or 'PLUG' on any line (case-insensitive)

      read (lunctl, *, iostat=ios) keyword

      if (ios /= 0) then ! EOF
         close (lunctl)
         exit
      end if

      call upcase (keyword)

      if (keyword == 'CAVI') cavity_case = true
      if (keyword == 'PLUG') plug_case   = true

   end do

   if (cavity_case) then
      filename_surface = 'cavity-surf.grd'
   else if (plug_case) then
      filename_surface = 'plug-surf.grd'
   else
      filename_surface = 'new-surf.grd'
   end if

   do while (ios /= 0)
      call reads (luncrt, &
         'New surface grid file name? [<CR> = ' // trim (filename_surface) //  &
         ']: ', lunkbd, filename_surface, cr, eof)
      if (eof) go to 99

      call determine_form (filename_surface, lunin, formatted_surface, ios)
   end do

   if (cavity_case .or. plug_case) then

      if (cavity_case) then
         filename_cavity = 'cavity-vol.grd'
         call reads (luncrt, &
                     'Cavity grid file name? [<CR> = cavity-vol.grd]: ',       &
                     lunkbd, filename_cavity, cr, eof)
      else ! Plug or gap filler case, but keep the same variable name
         if (index (trim (filename_surface), 'plug') > 0) then
            filename_cavity = 'plug-vol.grd'
            call reads (luncrt, &
                  'File name for blocks around plug? [<CR> = plug-vol.grd]: ', &
                        lunkbd, filename_cavity, cr, eof)
         else
            filename_cavity = 'gapfiller-vol.grd'
            call reads (luncrt, &
                  'File name for collar blocks? [<CR> = gapfiller-vol.grd]: ', &
                        lunkbd, filename_cavity, cr, eof)
         end if
      end if

      if (eof) go to 99

!     Check existence; allow for different formatting from the new surface.

      call determine_form (filename_cavity, lunin, formatted_volume, ios)

      if (ios /= 0) go to 99

   end if

!  Option to thin the new file(s):
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   thin_grid = true
   call ready (luncrt, 'Thin the new surface grid?       [y|n; <CR> = yes]: ', &
               lunkbd, thin_grid, cr, eof)
   if (eof) go to 99

   if (thin_grid) then

      write (lunsh, '(a)') &
         'echo '' ''',     &
         'echo '' Thinning the new surface grid ...''', &
         path // 'THIN_GRID/thin_grid << THIN_GRID_1'
      if (formatted_surface) then
         answer = y
      else
         answer = n
      end if
      write (lunsh, '(a)') &
         trim (filename_surface),   &
         answer
      write (lunsh, '(3i2)') 2, 2, 1
      filename_surface  = 'thinned_surface.gu'
      formatted_surface = false
      write (lunsh, '(a)') &
         trim (filename_surface),   &
         'n',                       &
         'THIN_GRID_1',             &
         blank

      if (cavity_case) then
         write (lunsh, '(a)')       &
            'echo '' ''',           &
            'echo '' Thinning the cavity grid ...''', &
            path // 'THIN_GRID/thin_grid << THIN_GRID_2'
         if (formatted_volume) then
            answer = y
         else
            answer = n
         end if
         write (lunsh, '(a)')       &
            trim (filename_cavity), &
            answer
         write (lunsh, '(3i2)') 2, 2, 1
         filename_cavity  = 'thinned_cavity.gu'
         formatted_volume = false
         write (lunsh, '(a)')       &
            trim (filename_cavity), &
            'n',                    &
            'THIN_GRID_2',          &
            blank
      end if

      if (plug_case) then
         write (lunsh, '(a)')       &
            'echo '' ''',           &
            'echo '' Thinning the blocks around the plug ...''', &
            path // 'THIN_GRID/thin_grid << THIN_GRID_2'
         if (formatted_volume) then
            answer = y
         else
            answer = n
         end if
         write (lunsh, '(a)')       &
            trim (filename_cavity), &
            answer
         write (lunsh, '(3i2)') 2, 2, 1
         filename_cavity  = 'thinned_plug.gu'
         formatted_volume = false
         write (lunsh, '(a)')       &
            trim (filename_cavity), &
            'n',                    &
            'THIN_GRID_2',          &
            blank
       end if

   end if

!  Option to change the units of the new file(s):
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   scale_grid = true
   call ready (luncrt, 'Convert new grid(s) to meters?   [y|n; <CR> = yes]: ', &
               lunkbd, scale_grid, cr, eof)
   if (eof) go to 99

   if (scale_grid) then

      write (lunsh, '(a)') &
         'echo '' ''',     &
         'echo '' Converting the surface grid to meters ...''', &
         path // 'ADJUST_GRID/adjust_grid << ADJUST_GRID_1',    &
         trim (filename_surface)
      write (lunsh, '(i2, /, i2)') 18, 99
      filename_surface  = 'new_surface_meters.gu'
      formatted_surface = false
      write (lunsh, '(a, /, a, /, a, //)') &
         trim (filename_surface), n,       &
         'ADJUST_GRID_1'

      if (cavity_case) then
         write (lunsh, '(a)') &
            'echo '' ''',     &
            'echo '' Converting the cavity grid to meters ...''', &
            path // 'ADJUST_GRID/adjust_grid << ADJUST_GRID_2',   &
            trim (filename_cavity)
         write (lunsh, '(i2, /, i2)') 18, 99
         filename_cavity  = 'cavity_meters.gu'
         formatted_volume = false
         write (lunsh, '(a, /, a, /, a, //)') &
            trim (filename_cavity), n,        &
            'ADJUST_GRID_2'
      end if

      if (plug_case) then
         write (lunsh, '(a)') &
            'echo '' ''',     &
            'echo '' Scaling the blocks around the plug  ...''', &
            path // 'ADJUST_GRID/adjust_grid << ADJUST_GRID_2',  &
            trim (filename_cavity)
         write (lunsh, '(i2, /, i2)') 18, 99
         filename_cavity  = 'plug_meters.gu'
         formatted_volume = false
         write (lunsh, '(a, /, a, /, a, //)') &
            trim (filename_cavity), n,        &
            'ADJUST_GRID_2'
      end if

   end if

!  RADIAL_INTERP (1):  Main volume grid.
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   decay = true
   call ready (luncrt, &
               'Decay any k=1 differences to 0. at k=nk? [y|n; <CR> = yes]: ', &
               lunkbd, decay, cr, eof)
   if (eof) go to 99

   open (lunout, file='radial_interp.inp.vertices', status='unknown')

   ios = 1
   do while (ios /= 0)
      filename = 'both-halves.cvertex.gu'  ! Shuttle-specific
      call reads (luncrt, &
        'Baseline volume grid (vertices)  [<CR> = both-halves.cvertex.gu]: ',  &
                  lunkbd, filename, cr, eof)
      if (eof) go to 99

      call determine_form (filename, lunin, formatted, ios)

   end do  ! Check existence

   write (lunout, '(a)') &
      'RADIAL_INTERP controls for vertex-centered volume grid',                &
      '---------------- INPUT VOLUME (or SURFACE) GRID ------------------',    &
      trim (filename)
   write (lunout, '(l1, T33, a)') formatted, 'Formatted? [T|F]'
   write (lunout, '(a1, T33, a)') 'F', 'Cell-centered? [T|F]'
   write (lunout, '(a)') &
      '--------------- INPUT FLOW SOLUTION (IF PRESENT) -----------------',    &
      'none                            Associated flow solution, or none.',    &
      'T                               Formatted? [T|F]',                      &
      'F                               Cell-centered? [T|F]',                  &
      '---------------- TARGET SURFACE or VOLUME GRID -------------------',    &
      trim (filename_surface)
   write (lunout, '(l1, T33, a)') formatted_surface, 'Formatted? [T|F]'
   write (lunout, '(a)') &
      '------------------ INTERPOLATED SURFACE GRID ---------------------',    &
      'cvertex.interp.surf.xyz  Output surface grid (check surface searches)', &
      'T                               Formatted? [T|F]',                      &
      '------------------ INTERPOLATED VOLUME GRID ----------------------',    &
      'cvertex.interp.volume.gu        Output volume grid, or none.',          &
      'F                               Formatted? [T|F]',                      &
      '----------------- INTERPOLATED FLOW SOLUTION ---------------------',    &
      'none                            Output flow field, or none.',           &
      'F                               Formatted? [T|F]',                      &
      '------------------- MISCELLANEOUS CONTROLS -----------------------',    &
      '1 Flow intp. method: 1 assumes consistent radial lines; 2 relaxes this',&
      'T T = allow srf. cell extrapln.; F = force surface (p,q)s into [0,1]'
   if (decay) then
      write (lunout, '(a)') &
      '-0.0002   Proj. dist. tol. for surf. cell searches (same units as grid)'
   else
      write (lunout, '(a)') &
      '+0.0002   Proj. dist. tol. for surf. cell searches (same units as grid)'
   end if

   close (lunout)

!  RADIAL_INTERP (2): cell-centered flow interpolation case:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   open (lunout, file='radial_interp.inp.centers', status='unknown')

   ios = 1
   do while (ios /= 0)
      filename = 'both-halves.ccenter.gu'
      call reads (luncrt, &
                 'Baseline volume grid (centers) [both-halves.ccenter.gu]: ',  &
                  lunkbd, filename, cr, eof)
      if (eof) go to 99

      call determine_form (filename, lunin, formatted, ios)

   end do

   write (lunout, '(a)') &
      'RADIAL_INTERP controls for cell-centered flow interpolation',           &
      '---------------- INPUT VOLUME (or SURFACE) GRID ------------------',    &
      trim (filename)
   write (lunout, '(l1, T33, a)') formatted, 'Formatted? [T|F]'
   write (lunout, '(a1, T33, a)') 'T', 'Cell-centered? [T|F]'

   ios = 1
   do while (ios /= 0) 
      filename = 'both-halves.ccenter.fu'
      call reads (luncrt, &
               'Associated cell-centered flow  [both-halves.ccenter.fu]: ',    &
                  lunkbd, filename, cr, eof)
      if (eof) go to 99

      call determine_form (filename, lunin, formatted, ios)

   end do

   write (lunout, '(a)') &
      '--------------- INPUT FLOW SOLUTION (IF PRESENT) -----------------',    &
      trim (filename)
   write (lunout, '(l1, T33, a)') formatted, 'Formatted? [T|F]'
   write (lunout, '(a)') &
      'T                               Cell-centered? [T|F]',                  &
      '---------------- TARGET SURFACE or VOLUME GRID -------------------',    &
      trim (filename_surface)
   write (lunout, '(l1, T33, a)') formatted_surface, 'Formatted? [T|F]'
   write (lunout, '(a)') &
      '------------------ INTERPOLATED SURFACE GRID ---------------------',    &
      'ccenter.interp.surf.xyz  Output surface grid (check surface searches)', &
      'T                               Formatted? [T|F]',                      &
      '------------------ INTERPOLATED VOLUME GRID ----------------------',    &
      'ccenter.interp.volume.gu        Output volume grid, or none.',          &
      'F                               Formatted? [T|F]',                      &
      '----------------- INTERPOLATED FLOW SOLUTION ---------------------',    &
      'ccenter.interp.volume.fu        Output flow field, or none.',           &
      'F                               Formatted? [T|F]',                      &
      '------------------- MISCELLANEOUS CONTROLS -----------------------',    &
      '1 Flow intp. method: 1 assumes consistent radial lines; 2 relaxes this',&
      'T T = allow srf. cell extrapln.; F = force surface (p,q)s into [0,1]'
   if (decay) then
      write (lunout, '(a)') &
      '-0.0002   Proj. dist. tol. for surf. cell searches (same units as grid)'
   else
      write (lunout, '(a)') &
      '+0.0002   Proj. dist. tol. for surf. cell searches (same units as grid)'
   end if

   close (lunout)

!  Script commands for running RADIAL_INTERP twice:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   write (lunsh, '(a)') &
     'rm -f radial_interp.inp',                                                &
     'ln -s radial_interp.inp.vertices radial_interp.inp',                     &
     'echo '' ''',                                                             &
     'echo '' Running RADIAL_INTERP to generate main volume grid blocks ...''',&
     path // 'RADIAL_INTERP/radial_interp > vertices.log',                     &
     'echo '' Main part of volume grid is done.  Check distance statistics:''',&
     'grep ''Fin'' vertices.log',                                              &
     'grep ''=   1   1   1'' vertices.log',                                    &
     blank

   write (lunsh, '(a)') &
     'rm -f radial_interp.inp',                                                &
     'ln -s radial_interp.inp.centers radial_interp.inp',                      &
     'echo '' ''',                                                             &
     'echo '' Running RADIAL_INTERP to interpolate main flow solution ...''',  &
     path // 'RADIAL_INTERP/radial_interp > centers.log',                      &
     'echo '' Main flow interpolation is done.   Check distance statistics:''',&
     'grep ''Fin'' centers.log',                                               &
     blank

!  Script commands for running COMBINE_BLOCKS twice?
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   gridname = 'cvertex.interp.volume.gu'  ! If no cavity grid

   if (cavity_case .or. plug_case) then

!     COMBINE_BLOCKS (1):  Vertex-centered grid
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (formatted_volume) then
         answer = y
      else
         answer = n
      end if

      write (lunsh, '(a)') &
         'echo '' Combining main volume & cavity or plug-related blocks ...''',&
        path // 'COMBINE_BLOCKS_TURB/combine_blocks_turb << COMBINE_BLOCKS_1', &
         'cvertex.interp.volume.gu', &
         n,                          &
         n,                          &
         trim (filename_cavity),     &
         answer

      gridname = 'combined.vertex.gu'
      call reads (luncrt, &
         'Name for the combined grid blocks   (<CR> = combined.vertex.gu): ',  &
                  lunkbd, gridname, cr, eof)
      if (eof) go to 99

      i = len_trim (gridname)
      formatted = gridname(i:i) /= 'u'

      if (formatted) then
         answer = y
      else
         answer = n
      end if

      write (lunsh, '(a)')   &
         gridname(1:i),      &
         answer,             &
         'COMBINE_BLOCKS_1', &
         blank,              &
         'echo '' Combined grid is complete.''', &
         blank

!     COMBINE_BLOCKS (2):  Cell centered main flow solution + cavity estimate
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      nspecies = 5
      call readi (luncrt, &
                  'Number of species in flow solution  (<CR> = 5): ',  &
                  lunkbd, gridname, cr, eof)
      if (eof) go to 99

      if (formatted_volume) then
         answer = y
      else
         answer = n
      end if

      write (lunsh, '(a)')   &
      'echo '' Combining main flow field with estimate for added blocks ...''',&
        path // 'COMBINE_BLOCKS_TURB/combine_blocks_turb << COMBINE_BLOCKS_2', &
         'ccenter.interp.volume.gu', &
         n,                          &
         y,                          &
         'ccenter.interp.volume.fu', &
         trim (filename_cavity),     &
         answer,                     &
         n
      write (lunsh, '(i2)') nspecies

      if (formatted) then     ! Same as for combined vertex grid
         answer   = y
         filename = 'combined.center.g'
      else
         answer   = n
         filename = 'combined.center.gu'
      end if

      write (lunsh, '(a)') &
         trim (filename),  &
         answer

      if (formatted) then
         filename_flow = 'combined.center.f'
      else
         filename_flow = 'combined.center.fu'
      end if

      write (lunsh, '(a)')     &
         trim (filename_flow), &
         answer,               &
         'COMBINE_BLOCKS_2',   &
         blank,                &
         'echo '' Combined flow starting guess is complete.''',                &
         blank

   end if

!  TEMPLATE:  Derive dplr.inputs and dplr.interfaces from the final volume grid:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   write (lunsh, '(a)')        &
      ' echo '' Generating DPLR control files from the new volume grid ...''', &
      path // 'TEMPLATE/template << TEMPLATE', &
      trim (gridname),         &
      '0.0001',                &
      'TEMPLATE',              &
      blank

   if (plug_case) then
      write (lunsh, '(a)')     &
      ' echo '' Adjusting the interpolated flow on the outer boundaries ...''',&
         path // 'PLUG_INTERP/plug_interp << PLUG_INTERP', &
         trim (filename),      &
         answer,               &
         trim (filename_flow)

      if (formatted) then
         filename_flow = 'combined.center.plug.adjusted.f'
         answer = y
      else
         filename_flow = 'combined.center.plug.adjusted.fu'
         answer = n
      end if

      write (lunsh, '(a)')     &
        trim (filename_flow),  &
        answer,                &
         'PLUG_INTERP',        &
         blank
   end if

   write (lunsh, '(a)')        &
      'echo '' All done.''',   &
      'echo '' Be sure to check dplr.inputs for the correct # blocks,''',      &
      'echo '' free-stream conditions, and CFL schedule,''',                   &
      'echo '' then run FCONVERT twice for the desired # CPUs.'''

   close (lunsh)

99 continue

!  Local procedure for program prepare_local_analysis

   contains

!     .........................................................................

      subroutine determine_form (multiblock_filename, lun, formatted, ios)

!     Open a multiblock grid or function file and attempt to read the number
!     of blocks, thus determining whether it appears to be unformatted or not.
!     Close the file upon return.  Using this in a while loop checking ios
!     allows retries in case of bad file names.
!     .........................................................................

      character, intent (in)  :: multiblock_filename * (*)
      integer,   intent (in)  :: lun
      logical,   intent (out) :: formatted
      integer,   intent (out) :: ios      ! ios /= 0 suggests a bad file name

      integer :: l, nb

      l = len_trim (multiblock_filename)

      open (lun, file=multiblock_filename(1:l), status='old', iostat=ios)

      if (ios /= 0) then
         write (*, '(2a)') ' Unable to open ', multiblock_filename(1:l)
         go to 99
      end if

      read (lun, *, iostat=ios) nb

      if (ios /= 0 .or. nb <= 0 .or. nb > 500) then
!!!      write (*,*) 'Formatted read gives nb & ios = ', nb, ios
!!!      write (*,*) 'Try again as unformatted.'
         close (lun)
         open (lun, file=multiblock_filename(1:l), status='old', &
               form='unformatted', iostat=ios)
         formatted = .false.
!!!      write (*, *) 'ios for unformatted open:', ios
         read  (lun, iostat=ios) nb
!!!      write (*, *) 'nb, ios:', nb, ios
      else
         formatted = .true.
!!!      write (*, *) 'Successful formatted open.'
!!!      write (*, *) 'nb, ios:', nb, ios
      end if

      close (lun, iostat=ios)

   99 return

      end subroutine determine_form

   end program prepare_local_analysis
