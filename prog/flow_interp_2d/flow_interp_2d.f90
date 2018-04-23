!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program flow_interp_2d

!  This program interpolates a 2-space flow field solution from one PLOT2D
!  multiblock grid to another.  The solution should be provided at all vertices
!  of the grid.  (Both files may in fact be cell-centered forms of the under-
!  lying computational grid, but they are always treated as though they are
!  vertex-centered.)  The same underlying geometry is assumed to be represented
!  by both the solution grid and the target grid, but the outer boundaries,
!  block counts, and/or point distributions may differ.
!
!  Lines of sight (perhaps from LINES_OF_SIGHT_2D) are likely to be the target
!  grid blocks (one 2-space line per block), as needed for radiative heating
!  calculations from axisymmetric solutions.  Tabulated forms of the output are
!  provided in this case (as determined by ni = 1 for all target blocks).
!  Otherwise, the interpolated flow is output as a PLOT2D function file.
!
!  Further Notes Inherited From The Earlier 3-Space FLOW_INTERP:
!
!     >  Ideally, all target grid points will be contained by the grid
!        associated with the given flow solution.
!     >  For target points outside the solution grid, the nearest solution cell
!        is used.  Thus, extrapolation is explicitly avoided.
!     >  The ADT (Alternating Digital Tree) technique used ensures efficiency
!        no matter how many of the target points are outside the searched grid.
!     >  Note that application to different geometries (even slightly perturbed
!        ones) is NOT appropriate:  boundary layer point distributions are all-
!        important, and the real-space interpolations performed here will not
!        preserve the boundary layer flow variables (or gradients, especially)
!        anywhere the two surface grids differ more than infinitesimally.
!     >  See SURFACE_INTERP for interpolating 3-space surface datasets in
!        Tecplot form.
!
!  Algorithm:
!
!     >  For each block in the target grid:
!        >  Read the target grid block.
!        >  If this is the first target block, perform some initialization:
!           >  Build the solution grid search tree from all cells of all blocks.
!              This requires insertion of z = 0 in the solution grid in order
!              to use the 3-space structured surface form of the ADT scheme.
!        >  For each point in the target block:
!           >  Locate the nearest solution cell by searching the ADT.
!           >  Interpolate the flow solution to the current target point using
!              the best solution cell found.  The interpolation is "bilinear"
!              within a single cell.  (This familiar formulation is not really
!              linear because it contains nonlinear terms, but the effect is
!              (bi)linear if the cell is perfectly rectangular.)
!        >  Output the interpolated solution.
!
!  Control:
!
!     A few prompts suffice.
!     This version allows for adding further options requiring further prompts.
!     First special option, prompted by shock-normal lines of sight:
!        Tabulate the component of velocity along each line of sight?
!
!  History:
!
!     03/24/04-08/22/08  D.A.Saunders  3-space FLOW_INTERP evolution.
!     02/24/12-02/27/12    "      "    2-space adaptation as FLOW_INTERP_2D.
!     07/24/13             "      "    Meeting a requirement for V.n along lines
!                                      of sight is best done here, but it calls
!                                      for further prompts.  Todd suggested that
!                                      we allow for other future possibilities
!                                      with the initial new prompt.
!     08/06/2013           "      "    All ADT variants have been merged into a
!                                      module with generic build & search calls.
!  Author:
!
!     David Saunders, ERC, Inc./NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:
!  !!!!!!!!

   use grid_block_structure  ! Data structure for one grid block

   implicit none

!  Constants:
!  !!!!!!!!!!

   integer, parameter :: &
      lunxyz  = 2,       &   ! Input solution grid
      luntarg = 3,       &   ! Input target grid
      lunflow = 4,       &   ! Input flow solution; output flow interpolations
      lunkbd  = 5,       &   ! Keyboard inputs
      luncrt  = 6            ! Prompts and diagnostics

   logical, parameter :: &
      false   = .false., &
      true    = .true.

!  Variables:
!  !!!!!!!!!!

   integer :: &
      i, ib, ios, j, n, nblocks_soln, nblocks_targ, ni, nj, numf

   logical :: &
      formatted_out, formatted_targ, initialize

   type (grid_type), pointer, dimension(:) :: &
      solution, target

!  Execution:
!  !!!!!!!!!!

   call control ()  ! Prompt for/read the solution to be searched, etc.

   if (ios /= 0) go to 99

!  Perform the flow field interpolation one target block at a time:
!  ----------------------------------------------------------------

   do ib = 1, nblocks_targ

      initialize = ib == 1

      ni = target(ib)%ni;  nj = target(ib)%nj

      allocate (target(ib)%x(ni,nj,1), target(ib)%y(ni,nj,1))

      if (formatted_targ) then
         read (luntarg, *, iostat=ios) target(ib)%x, target(ib)%y
      else
         read (luntarg,    iostat=ios) target(ib)%x, target(ib)%y
      end if

      if (ios /= 0) then
         write (luncrt, '(a, i4)') ' Trouble reading target block #', ib
         write (luncrt, '(a, 2i5)') ' Dimensions: ', ni, nj
         go to 99
      end if

      allocate (target(ib)%q(numf,ni,nj,1))

!     Interpolate the flow at the points of the current target block.
!     Some preprocessing of the input solution blocks is performed when ib = 1.

      call interp_q (nblocks_soln, numf, solution, ib, target(ib), initialize, &
                     luncrt)

!     Save results for this target block:

      if (formatted_out) then
         write (lunflow, '(1p, 6e18.10)', iostat=ios) &
            (((target(ib)%q(n,i,j,1), i = 1, ni), j = 1, nj), n = 1, numf)
      else
         write (lunflow, iostat=ios) &
            (((target(ib)%q(n,i,j,1), i = 1, ni), j = 1, nj), n = 1, numf)
      end if

      if (ios /= 0) then
         write (luncrt, '(/, a, i4)') &
            ' Trouble writing interpolated flow.  Block #:', ib
         write (luncrt, '(a, 3i5)') ' Dimensions & # functions: ', ni, nj, numf
         go to 99
      end if

      deallocate (target(ib)%x, target(ib)%y, target(ib)%q)

   end do ! Next target block

   close (luntarg)
   close (lunflow)

!  Done.

99 continue  ! Avoid system-dependent STOP behavior

!  Local procedures for program flow_interp_2d:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine control ()  ! Prompt for input/output files and other controls;
!                            ! read input files and set up output file.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local constants:

      real,      parameter :: zero = 0.
      character, parameter :: format*11 = 'unformatted'

!     Local variables:

      integer   :: i, i1, j, l, n
      logical   :: ascii, cr, eof, formatted_flow, formatted_xyz
      character :: filename*80

!     Execution:

!     Open and read the solution grid:
!     --------------------------------

      ios = 1
      do while (ios /= 0)  ! Check existence and format
         filename = 'ccenter.gu'
         call reads (luncrt, &
                     'Grid to be interpolated [<CR> = ccenter.gu]:   ',  &
                     lunkbd, filename, cr, eof)
         if (eof) go to 99
 
         call determine_grid_form (filename, lunxyz, formatted_xyz, ios)
      end do
 
      i1 = 1;  if (formatted_xyz) i1 = 3
 
      open (lunxyz, file=filename, form=format(i1:11), status='old')

      if (formatted_xyz) then
         read (lunxyz, *, iostat=ios) nblocks_soln
      else
         read (lunxyz,    iostat=ios) nblocks_soln
      end if

      if (ios /= 0) then
         write (luncrt, '(a)') 'Trouble reading # solution grid blocks.'
         go to 99
      end if

      allocate (solution(nblocks_soln))

      if (formatted_xyz) then
         read (lunxyz, *, iostat=ios) &
            (solution(ib)%ni, solution(ib)%nj, ib = 1, nblocks_soln)
      else
         read (lunxyz,    iostat=ios) &
            (solution(ib)%ni, solution(ib)%nj, ib = 1, nblocks_soln)
      end if

      if (ios /= 0) then
         write (luncrt, '(a)') 'Trouble reading solution grid dimensions.'
         go to 99
      end if

      do ib = 1, nblocks_soln
         ni = solution(ib)%ni;  nj = solution(ib)%nj;  solution(ib)%nk = 1
         allocate (solution(ib)%x(ni,nj,1), solution(ib)%y(ni,nj,1), &
                   solution(ib)%z(ni,nj,1))
         if (formatted_xyz) then
            read (lunxyz, *, iostat=ios) solution(ib)%x, solution(ib)%y
         else
            read (lunxyz,    iostat=ios) solution(ib)%x, solution(ib)%y
         end if

         if (ios /= 0) then
            write (luncrt, '(a, i5)') &
               'Trouble reading solution grid.  Block #:', ib
            go to 99
         end if

         solution(ib)%z = zero  ! The ADT search scheme is for 3-space data
      end do

      close (lunxyz)

!     Flow solution to be interpolated:
!     ---------------------------------

      ios = 1
      do while (ios /= 0)
         filename = 'ccenter.fu'
         call reads (luncrt, &
                     'Flow to be interpolated [<CR> = ccenter.fu]:   ',  &
                     lunkbd, filename, cr, eof)
         if (eof) go to 99

         call determine_grid_form (filename, lunflow, formatted_flow, ios)
      end do

      i1 = 1;  if (formatted_flow) i1 = 3

      open (lunflow, file=filename, form=format(i1:11), status='old')

      if (formatted_xyz) then
         read (lunflow, *, iostat=ios) nblocks_soln
         if (ios == 0) read (lunflow, *, iostat=ios) &
            (ni, nj, numf, ib = 1, nblocks_soln)
      else
         read (lunflow,    iostat=ios) nblocks_soln
         if (ios == 0) read (lunflow, iostat=ios) &
            (ni, nj, numf, ib = 1, nblocks_soln)
      end if

      if (ios /= 0) then
         write (luncrt, '(a)') 'Trouble reading solution file header.'
         go to 99
      end if

      do ib = 1, nblocks_soln
         ni = solution(ib)%ni;  nj = solution(ib)%nj
         allocate (solution(ib)%q(numf,ni,nj,1))

         if (formatted_flow) then
            read (lunflow, *, iostat=ios) &
               (((solution(ib)%q(n,i,j,1), i = 1, ni), j = 1, nj), n = 1, numf)
         else
            read (lunflow,    iostat=ios) &
               (((solution(ib)%q(n,i,j,1), i = 1, ni), j = 1, nj), n = 1, numf)
         end if

         if (ios /= 0) then
            write (luncrt, '(a, i5)') &
               'Trouble reading flow solution.  Block #:', ib
            go to 99
         end if
      end do

      close (lunflow)

!     Target grid:
!     ------------

      ios = 1
      do while (ios /= 0)  
         filename = 'los.g'
         call reads (luncrt, &
                     'Target grid to interpolate to [<CR> = los.g]:  ',  &
                     lunkbd, filename, cr, eof)
         if (eof) go to 99

         call determine_grid_form (filename, luntarg, formatted_targ, ios)
      end do

      i1 = 1;  if (formatted_targ) i1 = 3

      open (luntarg, file=filename, form=format(i1:11), status='old')

!     There is no need to store the entire target grid.  Read its header:

      if (formatted_targ) then
         read (luntarg, *, iostat=ios) nblocks_targ
      else
         read (luntarg,    iostat=ios) nblocks_targ
      end if

      if (ios /= 0) then
         write (luncrt, '(a)') 'Trouble reading # target grid blocks.'
         go to 99
      end if

      allocate (target(nblocks_targ))

      if (formatted_targ) then
         read (luntarg, *, iostat=ios) &
            (target(ib)%ni, target(ib)%nj, ib = 1, nblocks_targ)
      else
         read (luntarg,    iostat=ios) &
            (target(ib)%ni, target(ib)%nj, ib = 1, nblocks_targ)
      end if

      if (ios /= 0) then
         write (luncrt, '(a)') 'Trouble reading target grid dimensions.'
         go to 99
      end if

!     Set up the interpolated flow output file:
!     -----------------------------------------

      ios = 1
      do while (ios /= 0)
         filename = 'interpolated_flow.f'
         call reads (luncrt, &
                     'Output flow file [<CR> = interpolated_flow.f]: ',  &
                     lunkbd, filename, cr, eof)
         if (eof) go to 99

         formatted_out = true  ! Avoid the prompt unless it's really needed
         i1 = 1;  if (formatted_out) i1 = 3

         open (lunflow, file=filename, form=format(i1:11), status='unknown', &
               iostat=ios)
      end do

      if (formatted_out) then
         write (lunflow, '(i4)',  iostat=ios) nblocks_targ
         if (ios == 0) write (lunflow, '(3i5)', iostat=ios) &
            (target(ib)%ni, target(ib)%nj, numf, ib = 1, nblocks_targ)
      else
         write (lunflow, '(i4)',  iostat=ios) nblocks_targ
         if (ios == 0) write (lunflow, '(3i5)', iostat=ios) &
            (target(ib)%ni, target(ib)%nj, numf, ib = 1, nblocks_targ)
      end if

      if (ios /= 0) then
         write (luncrt, '(a)') 'Trouble writing interpolated flow header.'
         go to 99
      end if

 99   return

      end subroutine control

   end program flow_interp_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine interp_q (nblocks_soln, numf, solution, iblock_target, &
                        target_block, initialize, luncrt)
!
!  This routine interpolates the given 2-space flow solution onto the indicated
!  block (line?) of a target grid.  It is the argument-driven, reusable portion
!  of program flow_interp_2d, q.v.  Some of the calculations are in 3-space
!  because we don't have a 2-space analogue of the structured surface grid
!  search package.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:
!  !!!!!!!!

   use grid_block_structure  ! Data structure for one grid block
   use adt_utilities         ! All variants of the ADT routines

   implicit none

!  Arguments:
!  !!!!!!!!!!

   integer, intent (in)             :: &
      nblocks_soln                     ! # blocks in the solution being searched
   integer, intent (in)             :: &
      numf                             ! Number of flow variables
   type (grid_type), intent (inout) :: &
      solution(nblocks_soln)           ! Solution grid and flow blocks
   integer, intent (in)             :: &
      iblock_target                    ! Current target block #, for diagnosics
   type (grid_type), intent (inout) :: &
      target_block                     ! Current target block, input with its
                                       ! grid defined (only), and output with
                                       ! its flow field set via interpolation
                                       ! or from the nearest boundary point
   logical, intent (in)             :: &
      initialize                       ! T => build the search tree, etc.
   integer, intent (in)             :: &
      luncrt                           ! Logical unit for diagnostics

!  Local constants:
!  ----------------

   integer, parameter :: &
      lunkbd = 5                       ! For belated special option prompts

   real, parameter :: &
      big  = 1.e+32,  &
      one  = 1.0,     &
      zero = 0.0

   logical, parameter :: &
      true = .true.

!  Local variables:
!  !!!!!!!!!!!!!!!!

   integer :: i, ib, ic, icell, j, jc, n, ni, ninside, nj, nout
   logical :: cr, eof, profile
   real    :: dsqmin, dsqmax, dtol, p, pm1, q, qm1, xmax, xmin, Vdotn, &
              wall_distance
   real, dimension (2) :: un, V
   real, dimension (3) :: interp_coords, target_coords, wall_coords

!  Variables saved for further calls:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer,   save :: ispecial, iv, ncell
   integer,   dimension (:,:), allocatable, save :: conn
   real,      save :: dtolsq
   logical,   save :: special
   character, save :: pformat*39, pformat1*18

!  Execution:
!  !!!!!!!!!!

!  Pre-process the solution blocks?

   if (initialize) then

!     Build an Alternating Digital Tree for all cells of all solution blocks.
!     First, count the cells so the "connectivity" array (block number and
!     (i,j) for each cell) can be allocated:

      ncell = 0;  xmin = big;  xmax = -big

      do ib = 1, nblocks_soln
         ncell = (solution(ib)%ni - 1) * (solution(ib)%nj - 1) + ncell
         xmax  = max (xmax, maxval (solution(ib)%x))
         xmin  = min (xmin, maxval (solution(ib)%x))
      end do

      dtol     = (xmax - xmin) * max (10.0*epsilon (dtol), 1.e-8)
      dtolsq   = dtol * dtol

      allocate (conn(3,ncell))

      call build_adt (nblocks_soln, solution, ncell, conn)

!     Handle any number of functions in profile tabulations:

      pformat = '(i9, 4x, 2i4, 2f12.8, es15.7, nfes12.4)'
      write (pformat(31:32), '(i2)') numf

   end if        ! End of initialization


   ninside = 0   ! # target block points found inside a solution cell
   nout    = 0   ! # such points found more than the distance tolerance
                 ! from the nearest solution cell
   ni = target_block%ni
   nj = target_block%nj
   ib = iblock_target

   profile = ni == 1  ! A single radial line gets special treatment for
                      ! comparison with wind tunnel data perhaps
   if (profile) then
      if (initialize) then  ! Allow for adding special options
         special = .false.;  ispecial = 0
         call ready (luncrt, &
                     'Any special tabulation option(s)? [y|n; <CR> = n]: ', &
                     lunkbd, special, cr, eof)
         if (special) then
            write (luncrt, '(3x, a)') &
               'Option 1: Tabulate velocity component along current line.', &
               'Option 2: TBD'
            ispecial = 1
            call readi (luncrt, 'Choice? [<CR> = 1]: ', &
                        lunkbd, ispecial, cr, eof)
            select case (ispecial)
               case (1)  ! V.n
                  iv = 2
                  call readi (luncrt, 'Index of vx in the flow data: ', &
                              lunkbd, iv, cr, eof)
                  pformat1 = '(es15.7, nfes12.4)'
                  write (pformat1(10:11), '(i2)') numf + 1
               case default
                  write (luncrt, '(a, i4)') 'Invalid choice:', ispecial
                  special = .false.;  ispecial = 0
            end select
         end if
      end if
      select case (ispecial)
         case (0)  ! No special tabulation
            write (luncrt, '(/, a, i4, /, 2a)') '# Profile #:', ib,  &
               '# Soln. block  ic  jc           p           q  Wall distance', &
               '  Interpolated flow'
         case (1)  ! V.n
            write (luncrt, '(/, a, i4, /, a)') '# Profile #:', ib,  &
               '# Wall distance  Interpolated flow'
      end select
   else
      write (luncrt, '(/, a, i4, a, 2i5, a, i7)') &
         ' Beginning target grid block', ib, '.  Dimensions:', &
         ni, nj, '  Total:', ni * nj
   end if

   target_coords(3) = zero  ! Needed by 3-space surface search package
   dsqmax           = zero

!  For each point of the current target block:

   do j = 1, nj

      do i = 1, ni

         target_coords(1) = target_block%x(i,j,1)
         target_coords(2) = target_block%y(i,j,1)

         call search_adt (target_coords, icell, p, q, dsqmin, true, &
                          nblocks_soln, solution, ncell, conn, interp_coords)

         dsqmax = max (dsqmin, dsqmax)

         if (dsqmin < dtolsq) then
            ninside = ninside + 1
         else
            nout = nout + 1
         end if

!        Interpolate the flow solution within the closest cell found:

         n  = conn(1,icell) ! Block #
         ic = conn(2,icell) ! Lower left cell indices
         jc = conn(3,icell)

         pm1 = one - p
         qm1 = one - q

         target_block%q(:,i,j,1) = &
         qm1*(pm1*solution(n)%q(:,ic,jc,  1) + p*solution(n)%q(:,ic+1,jc,  1)) &
         + q*(pm1*solution(n)%q(:,ic,jc+1,1) + p*solution(n)%q(:,ic+1,jc+1,1))

         if (profile) then
           if (j == 1) then
               wall_coords(1) = &
             qm1*(pm1*solution(n)%x(ic,jc,  1) + p*solution(n)%x(ic+1,jc,  1)) &
             + q*(pm1*solution(n)%x(ic,jc+1,1) + p*solution(n)%x(ic+1,jc+1,1))
               wall_coords(2) = &
             qm1*(pm1*solution(n)%y(ic,jc,  1) + p*solution(n)%y(ic+1,jc,  1)) &
             + q*(pm1*solution(n)%y(ic,jc+1,1) + p*solution(n)%y(ic+1,jc+1,1))
             select case (ispecial)
               case (0)  ! No special tabulation
               case (1)  ! V.n needs a unit normal, assumed to be along the line
                 un(1) =   target_coords(1) - target_block%x(1,nj,1)
                 un(2) = -(target_coords(2) - target_block%y(1,nj,1)) ! Careful!
                 un(:) = un(:) / sqrt (un(1)**2 + un(2)**2)
             end select
           end if
           wall_distance = sqrt ((wall_coords(1) - target_coords(1))**2 + &
                                 (wall_coords(2) - target_coords(2))**2)
           select case (ispecial)
              case (0)  ! No special tabulation
                 write (luncrt, pformat) &
                    n, ic, jc, p, q, wall_distance, target_block%q(:,1,j,1)
              case (1)  ! V.n
                 V(:)  = target_block%q(iv:iv+1,1,j,1)
                 Vdotn = abs (dot_product (V, un))
                 write (luncrt, pformat1) &
                    wall_distance, target_block%q(:,1,j,1), Vdotn
           end select
         end if

      end do ! Next i for this target block

   end do ! Next j for this target block

   if (.not. profile) then
      write (luncrt, '(a, i4, a, i9, i8, a, 1p, e13.5)') &
         ' Finished  target grid block', ib, ';  in, out:', ninside, nout, &
         ';  largest distance from target:', sqrt (dsqmax)
   end if

   end subroutine interp_q
