!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program cavity_map
!
!  Description:
!
!     This program maps a CFD surface solution in and around a cavity or gouge
!     on to an idealized cavity grid for comparison with other calculations.
!
!     The earlier Ensemble program assumed the cavities were rectangular, but
!     no such assumption is made here.  Instead, an idealized rectilinear
!     form of the cavity is used as an ancillary input.  This has four cavity
!     walls and floor at arbitrary angles, plus surrounding surface patches.
!
!     The rectilinear outer envelope comes from the Boeing Cavity Heating Tool.
!     Mapping the CFD solution to the Cavity Heating Tool grid turns out to be
!     what was really required here (not the mapping to a unit cavity done in
!     Ensemble), so that CFD fidelity can indicate whether the engineering tool
!     is conservative or not.  Such a mapping has been added as a further
!     output, without eliminating the unit cavity mapping of both the CFD and
!     outer envelope solutions.  This version does more than interpolate at the
!     coarse mesh vertices: it averages the CFD solution over the areas of each
!     coarse cell.  It also generates slices of both the mapped/smoothed CFD
!     data and the input CHT data for line plot comparisons.
!
!  Strategy of original outer/inner envelope approach, and the alignment scheme:
!
!     >  Impose CFD-type dense grids on the idealized surface patch edges.
!        This is the outer envelope for the actual damage.
!     >  Construct a scaled-down, shallow form of the idealized walls and
!        floor, sitting on the smooth OML and centered sensibly.  This acts
!        as an inner envelope for the cavity portion.
!     >  Align the CFD surface grid with the working axes by calculating the
!        principal axes of the CFD cavity patches (masses <-> cell areas acting
!        at cell centroids lead to moments of inertia about the center of mass;
!        the eigenvectors of the inertia matrix form rows of the needed
!        rotation matrix).
!        N.B.:  This version no longer uses cavity L, W, D to distinguish the
!        principal axes, because as the dimensions become less distinct, and
!        with sloping walls in the picture, the area-based moments of inertia
!        are less likely to make the distinctions reliably.
!        Instead, all possible permutations of the tentative rotation matrix
!        are applied, and the one that best transforms a third (derived)
!        reference point and the center of mass is chosen.
!     >  For each line connecting corresponding points of the two sets of
!        idealized wall and floor patches, find the intersection with the
!        realigned CFD surface grid.  The CFD solution at that intersection is
!        the desired mapped solution on a unit cavity, which is topologically
!        the same as the dense, idealized cavity and can be derived from it by
!        normalizing the arc lengths as (u,v) in the unit square and using u
!        and v for the non-constant x, y, or z appropriately.
!     >  For the smooth OML patches downstream and adjacent to the unit cavity,
!        employ the ADT search and interpolation techniques used in Ensemble.
!     >  There appears to be no point in retaining the option to perform
!        ensemble averaging of multiple solutions.
!
!     N.B.: This version also permits use of the plain ADT search method for all
!     idealized patches.  It appears that the inner/outer envelope method for
!     the cavity matches may not be necessary, and it suffers in places where
!     the intersection lines are almost parallel to the CFD surface.  This
!     possibility forced suppression of possible line extrapolation during the
!     line/surface intersection calculations.  Method 2 is now recommended.
!
!  Assumptions:
!
!     >  Shuttle coordinate system:  X -> downstream, Y -> right wing, Z -> up
!     >  The idealized outer envelope comes from the Cavity Heating Tool and is
!        in a standard 10-patch form.  No assumptions are made about the (i,j)
!        convention, however.
!     >  The x,y,z coordinates have the same units.  A hidden scaling factor
!        applied to the CFD data is optional at the end of control file.
!
!  Control file format ('cavity_map.inp'):
!
!     Cavity_map control file
!     dplr_cavity_soln.dat  ! CFD surface solution (Tecplot ASCII file)
!     1                     ! Flow variable representing heating bump factor
!     1  10  33             ! CFD patch # and (i, j) of reference upstrm. pt.
!     8  10  33             ! CFD patch # and (i, j) of reference dnstrm. pt.
!     1.0                   ! Reference heating value; < 0 => diagnostics
!     ! List of CFD patches inside the cavity
!     1, 2, 10:25
!     heating_tool_grid.dat ! Idealized cavity grid (Tecplot ASCII)
!     2                     ! 1 = inner/outer envelope method; 2 = plain ADT
!    [-0.0254               ! Optional scale factor applied to the CFD xyzs;
!                           ! negative means divide rather than multiply.]
!
!  Notes:
!
!     (1) The reference CFD points should be at the upstream and downstream
!         cavity lips near the mid-points.  The line joining them after
!         realignment will be in the Z = 0 plane.
!     (2) The reference heating value allows flow values to be converted to
!         bump factors; use 1.0 if the flow values are already bump factors;
!         use < 0 to turn on output of intermediate grids and diagnostics
!     (3) The list of CFD cavity patches can use commas and colons.
!     (4) See the optional scale factor at the end of the control file above.
!
!  Input Surface Solution Format (Tecplot ASCII, structured, CFD & Idealized):
!
!        TITLE     = ""
!        VARIABLES = "x" "y" "z" "BF" ["..."]
!        ZONE T="G1"
!         I=66, J=66, K=1, ZONETYPE=Ordered
!         DATAPACKING=BLOCK
!         DT=(DOUBLE DOUBLE ... (variable # of them) ... DOUBLE DOUBLE )
!         6.647333145E+00 6.575638294E+00 6.489704609E+00 6.390774727E+00 6....
!          :               :               :               :               :
!     POINT data packing is also permitted.  One or more flow variables may be
!     present.
!
!  Output Results (Tecplot ASCII files):
!
!     (1) CFD solution mapped to a unit cavity with CFD-type grid density.
!     (2) CFD solution mapped to the plain outer envelope grid (CHT grid).
!         This was added belatedly as what was really wanted.
!         However, the realignment benefits from use of CFD-type resolution.
!         Results have been non-dimensionalized by the ideal cavity length.
!     (3) Two pairs of slices and dices of mapped/smoothed CFD data on the
!         CHT grid and of any functions input with the CHT grid.
!         The slice & dice abscissas are arc length, not X or Y, for better
!         display of results along vertical cavity walls.
!         Each cut is a separate zone.
!
!     The output formats are similar to the input formats (Tecplot ASCII).
!     The output file names are derived from the input CFD file name.
!
!  History:
!
!     05/13/05    D. Saunders  Final version of Ensemble (rectangular cavities).
!     6/14-6/21    "      "    Initial adaptation as Cavity_map.
!     06/22/05     "      "    Map the plain outer envelope (probably from the
!                              cavity heating tool) by making sure the indexing
!                              is in standard order then normalizing as for the
!                              form used with the interpolated CFD solution.
!     06/29/05     "      "    The transpose of the rotation matrix should be
!                              applied for realignment.  This was hidden during
!                              tests that mapped the outer envelope to itself.
!     07/01/05     "      "    Added the option to use the plain ADT mapping
!                              method everywhere.  Added mapping to the plain
!                              outer envelope grid.
!     07/03/05     "      "    Introduced a second CFD reference point in order
!                              to take out a misalignment that is inevitable if
!                              the ideal cavity walls are not symmetric.
!     09/21/05     "      "    Steve Alter urged a control-volume averaging
!                              approach for the interpolations to the coarse
!                              Cavity Heating Tool mesh. This is now done by
!                              dividing the coarse cells into many sub-cells.
!     10/19/05     "      "    Output line plots from the CHT grid (CFD + CHT).
!     10/29/05     "      "    Mapping a copy of the coarse, linear CHT data to
!                              the CHT grid highlighted the lack of continuity
!                              at patch edges caused by the area-based averaging
!                              to cell centers.  Therefore, force equality along
!                              common edges by averaging in map_to_envelope.
!     11/02/05     "      "    Make arc length the abscissa for slices & dices.
!                              This was a little awkward, requiring temporary
!                              copies of slice X & Ys, etc.
!     02/06/06     "      "    Internal procedure deallocate_blocks is now a
!                              part of Tecplot_io.f90.
!     02/16/06     "      "    Steve Alter's subscript check caught something in
!                              in subroutine slice_patch.  The slice arrays
!                              passed from the higher level are big enough but
!                              were wrongly dimensioned in the routine. Using
!                              (*) is simpler than passing the right length.
!     04/27/06     "      "    Using the same slice & dice arrays for the CFD
!                              and CHT solutions is bad inside the Tec I/O pkg.
!                              if the numbers of functions differ.  We need to
!                              deallocate %q arrays and reallocate them if a
!                              CHT solution is present.
!     05/02/06     "      "    Divide the dimensional outputs by cavity LENGTH.
!     06/01/06     "      "    Permute_R had a glitch affecting the case where
!                              cavity W > L > D.
!     08/14/06     "      "    Revised to use Tecplot 360 version of the
!                              Tecplot_io package.
!     08/24/06     "      "    The mapped.* file's variable names beyond 4
!                              were not being copied from the CFD file.
!     07/26/07     "      "    Use of L, D, W to distinguish the principal axes
!                              is not reliable if any two are close to each
!                              other, and sloping walls add to the fuzziness of
!                              the calculations.  Instead, all 48 permutations
!                              of the initial rotation matrix are now tried, and
!                              the one that best transforms a third, derived
!                              reference point is chosen.
!                              An optional scale factor may be included at the
!                              end of the control file.  It is applied to the
!                              CFD coordinates.  If it is negative, it is used
!                              to divide rather than multiply the coordinates.
!     07/29/07     "      "    Offset the derived reference point to help break
!                              ties, and include the transformed center of mass
!                              in the test for the best rotation.
!     07/30/07     "      "    Tabulate floor averages for 85/15% and 50/40/10%
!                              splits.
!     07/31/07     "      "    Tabulate floor averages for the CHT solution too.
!     08/05/07     "      "    CFD uncertainties are now added according to the
!                              Length/Depth ratio (floor and walls) or half the
!                              excess over 1 (outside, after averaging).  The
!                              line plot files now contain an extra variable.
!     08/10/07     "      "    Bill Wood revised the uncertainty definitions:
!                              walls and OML now use |BF - 1|.  Also, the short/
!                              long cavity cutoff changed from 10 to 15 L/D.
!     08/15/07     "      "    Before the principal axis calculations, align the
!                              reference points with the X axis.  This makes no
!                              difference in some cases (suggesting consistency)
!                              yet can help other cases (possibly with more of a
!                              rotation about OX).
!     08/16/07     "      "    Determine preliminary and final rotations about
!                              line joining reference points by 1D optimization.
!     08/20/07     "      "    Remove possible bad scaling of the CFD data by
!                              shifting it to the cavity patch center first.
!                              However, it doesn't help - moments of inertia
!                              work with X - Xcg, etc. anyway.
!                              Now that the third 2-D rotation can be found via
!                              1-D optimization, just go with that now.
!                              This does depend on the reference points' being
!                              on a line parallel to the length axis.
!     08/08/13     "      "    All ADT variants are now in one module with
!               Now ERC, Inc.  generic build_adt and search_adt interfaces.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_header_structure   ! Part of Tecplot_io.f90
   use grid_block_structure    !  "    "    "    "    "
   use tecplot_io_module       !  "    "    "    "    "
   use surface_patch_utilities ! Uses grid_block_structure
   use adt_utilities           !  "    "    "    "    "

   implicit none

!  Constants:

   integer, parameter :: &
      lunctl     =  1,   &  ! Control file
      lun_cfd    =  2,   &  ! Input CFD solution (Tecplot ASCII)
      lun_ideal  =  3,   &  ! Input idealized cavity grid (Tecplot ASCII)
      lun_norm   =  4,   &  ! Normalized form of idealized input cavity
      luncrt     =  6,   &  ! For diagnostics
      lun_interp =  7,   &  ! Interpolations of transformed CFD to CHT g.
      lun_interp2 = 8,   &  ! Area-averaged interpolations to CHT grid
      lun_out    =  9,   &  ! Unit cavity, CFD-like density + soln. (Tecplot)
      lun_out2   = 10,   &  ! CFD soln. mapped/avged. to ideal grid (Tecplot)
      lun_slice  = 11,   &  ! Line plots (CFD + CHT) along cavity   (Tecplot)
      lun_dice   = 12,   &  ! Line plots (CFD + CHT) across cavity  (Tecplot)
      namelength = 32,   &  ! Variable name limit hard-coded in Tecplot_io.f90
      nb_ideal   = 10,   &  ! # patches in idealized cavity grid
      ndim       = 3,    &  ! All Tecplot files are 3D
      nslices    = 3,    &  ! # streamwise cuts (Y stations)
      ndices     = 8,    &  ! # cross-cuts      (X stations)
      ndices1    = 4,    &  ! # dices through cavity
      ndices2    = 4,    &  ! "   "   aft of  "   "
      nedges     = 15       ! # common edges among CHT patches

!     N.B.: Keep both nslices and ndices odd in order to include cavity center.

   integer, parameter :: &
      ib_floor   = 5,    &  ! CHT grid patch on the floor
      nsubareas  = 5        ! Duplicated here after tabulation of floor averages
                            ! for CHT data was called for in addition to the
                            ! CFD floor averages implemented originally;
                            ! see subroutine lengthwise_floor_splits, where most
                            ! of this grubbiness is modularized
   real, parameter ::    &
      one  = 1.0,        &
      zero = 0.0,        &
      undefined = -999.     ! Should match internal variable of Tecplot_io.f90

!  Empirical uncertainty parameters used in procedures standard_cavity_grids,
!  lengthwise_splits and !  add_uncertainties:

   real, parameter ::    &
      OML_factor = 0.5,  &  ! Uncertainty = |BF - 1| * OML_factor outside cavity
      short_test = 15.,  &  ! "short" cavity <-> L/D < short_test
      uncertainty_long  = 0.40, & ! Added to BF inside "long"  cavities
      uncertainty_short = 0.15    ! Added to BF inside "short" cavities

   logical, parameter :: &
      false = .false.,   &
      true  = .true.

!  Variables:

   integer :: &
      datapacking, i, i_bump_var_in, i_ref_d, i_ref_u, ib, ib_ref_d, ib_ref_u, &
      ic, ios, iquad, j, j_ref_d, j_ref_u, jc, map_method, n, nb_cfd,          &
      nb_inner, nb_out, ni, ninside, nj, noutside, nquad,                      &
      num_cavity_patches, numf_cfd, numf_ideal, numf_slices

   integer :: &
      kf, kd, ku, ko, ki, kt  ! Patch numbers for the general case

   integer :: &
      medges(5,2,nedges) ! Set by save_common_edge_info; used in map_to_envelope

   integer, allocatable, dimension (:) :: &
      list_of_cavity_patches

   integer, allocatable, dimension (:,:) :: &
      conn               ! For patch & (i,j) of surface quads. being searched

   real :: &
      bump_factor_reference_value, bump_scale, cavity_length, cavity_depth,    &
      dmax, dmean, dsqmin, dtolsq, uncertainty_floor, xyz_scale,               &
      DEPTH, LENGTH, WIDTH

!     N.B.:  DEPTH, LENGTH, WIDTH are now AFTER dividing the idealized cavity
!            coordinates by the input ideal cavity length.  LENGTH is now 1.
!            These normalized values are no longer used.
!     See subroutine standard_cavity_grids for setting of uncertainty values.

   real, dimension (nsubareas) :: &
      floor_averages_CFD, floor_averages_CHT       ! See lengthwise_floor_splits

   real, allocatable, dimension (:,:) :: &
      dice_1_arcs, dice_1_ys, dice_2_arcs, dice_2_ys, slice_arcs, slice_xs

   logical :: &
      convert_to_BF, diagnostics, short_cavity

   character :: &
      buffer * 128

!  Derived data types:

   type (grid_header) :: &             ! These parallel the following grid_types
      header_cfd, header_ideal, header_norm, header_inner, header_out,         &
      header_out2, header_surf, header_surf2, header_slices, header_dices

   type (grid_type) :: &
      floor_centers_CHT              ! For cell-centered CHT floor function data

   type (grid_type), pointer, dimension (:) :: &
      xyzq_cfd,    & ! input CFD surface solution
      xyzq_ideal,  & ! idealized outer envelope surface grid (CHT grid + soln.)
      xyzq_norm,   & ! normalized form of input outer enveloped grid + solution
      xyzq_inner,  & ! inner envelope grid derived from regridded outer envelope
      xyzq_out,    & ! standard CFD-like grid + mapped soln.; unit cavity output
      xyzq_out2,   & ! CFD solution mapped to the plain outer envelope grid
      interp_surf, & ! real space interpolated surface points for plot checking
      interp_surf2,& ! ditto for interpolations to the plain outer envelope pts.
      slices,      & ! For slices along CHT grid (mapped CFD soln. & CHT soln.)
      dices          ! For slices across ......................................

!  !!!!!!!!!!
!  Execution:
!  !!!!!!!!!!

   open (lunctl, file='cavity_map.inp', status='old', iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') &
         ' Unable to open cavity_map.inp control file.'
      go to 99
   end if

!  Echo the control file to the output log:

   do ! Until EOF
      read (lunctl, '(a)', iostat=ios) buffer
      if (ios < 0) exit

      write (luncrt, '(1x, a)') trim (buffer)
   end do

   write (luncrt, '(a)');  rewind (lunctl)

   read (lunctl, *)                             ! Case description
   read (lunctl, *) header_cfd%filename         ! Input CFD solution file name
   read (lunctl, *) i_bump_var_in               ! Flow variable for bump factor
   read (lunctl, *) ib_ref_u, i_ref_u, j_ref_u  ! Block # & (i,j) on upstrm. lip
   read (lunctl, *) ib_ref_d, i_ref_d, j_ref_d  ! Block # & (i,j) on dnstrm. lip
   read (lunctl, *) bump_factor_reference_value ! Reference heating value

   diagnostics = bump_factor_reference_value < zero
   bump_factor_reference_value = abs (bump_factor_reference_value)

!  We need to know the number of patches to allocate the list of cavity patches,
!  so read the CFD data now:

   header_cfd%formatted = true
   header_cfd%ndim      = ndim

   ios = 0
   if (diagnostics) ios = 1 ! Turns on verbose mode

   call Tecplot_read (lun_cfd, header_cfd, xyzq_cfd, ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Trouble reading CFD Tecplot file ', &
         trim (header_cfd%filename)
      go to 99
   end if

   numf_cfd = header_cfd%numq
   nb_cfd   = header_cfd%nblocks

   read (lunctl, *)                    ! RDLIST can't handle a trailing comment

   allocate (list_of_cavity_patches(nb_cfd));  num_cavity_patches = nb_cfd ! Max

   call rdlist (luncrt, ' ', lunctl, num_cavity_patches, list_of_cavity_patches)

   read (lunctl, *) header_ideal%filename      ! Idealized cavity grid file name
   read (lunctl, *) map_method         ! 1 = inner/outer cavity envelope method;
                                       ! 2 = plain ADT method everywhere
   read (lunctl, *, iostat=ios) xyz_scale      ! Optional input

   if (ios /= 0) xyz_scale = one
   if (xyz_scale < zero) xyz_scale = -one / xyz_scale  ! [-]0.0254 is easier to
                                                       ! remember
   close (lunctl)

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!
!  End of control file inputs.
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!

   write (luncrt, '(/, a, /, (32i4))') &
      ' CFD cavity patches:', list_of_cavity_patches(1:num_cavity_patches)

   bump_scale  = one / bump_factor_reference_value
   convert_to_BF   =   bump_factor_reference_value /= one

   if (convert_to_BF) then
      if (i_bump_var_in < 0 .or. i_bump_var_in > numf_cfd) then
         write (luncrt, '(/, a, i3)') ' Bad i_bump_var_in:', i_bump_var_in
         go to 99
      end if
      do ib = 1, nb_cfd
         xyzq_cfd(ib)%q(i_bump_var_in,:,:,1) = bump_scale * &
            xyzq_cfd(ib)%q(i_bump_var_in,:,:,1)
      end do
      header_cfd%varname(ndim + i_bump_var_in) = 'BF'
   end if

   if (xyz_scale /= one) then
      do ib = 1, nb_cfd
         xyzq_cfd(ib)%x(:,:,1) = xyz_scale * xyzq_cfd(ib)%x(:,:,1)
         xyzq_cfd(ib)%y(:,:,1) = xyz_scale * xyzq_cfd(ib)%y(:,:,1)
         xyzq_cfd(ib)%z(:,:,1) = xyz_scale * xyzq_cfd(ib)%z(:,:,1)
      end do
   end if

!  Read the idealized cavity grid [probably with cavity heating tool soln.].
!  This is the outer envelope (full scale, but unrotated).

   header_ideal%formatted = true
   header_ideal%ndim      = ndim

   if (diagnostics) ios = 1 ! Turns on verbose mode

   call Tecplot_read (lun_ideal, header_ideal, xyzq_ideal, ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Trouble reading Tecplot file ', &
         trim (header_ideal%filename)
      go to 99
   end if

   numf_ideal = header_ideal%numq
   nb_out     = header_ideal%nblocks

   if (diagnostics) then
      write (luncrt, '(/, a, 2i4)') &
        ' # patches and variables in idealized cavity file:', nb_out, numf_ideal
      write (luncrt, '(a, 10(3x, a))') ' Names found:', &
         (trim (header_ideal%varname(i)), i = 1, ndim + numf_ideal)
   end if

   if (nb_out /= 10) then
      write (luncrt, '(/, 3a, i6)') &
       ' # patches found in ', trim (header_ideal%filename), &
       ' not the expected 10:', nb_out
      go to 99
   end if

!  Set up for mapping/averaging the CFD solution to the plain idealized grid.
!  First, the idealized grid coordinates are divided by the cavity length in
!  standard_cavity_grids, so adjust the variable names accordingly:

   header_ideal%varname(1) = 'X/L'
   header_ideal%varname(2) = 'Y/L'
   header_ideal%varname(3) = 'Z/L'

!  Clone_header is part of the Tecplot 360 update to Tecplot_io.f90.

   call clone_header (header_ideal, numf_cfd, 0, header_out2)  ! 0 = POINT

   n = ndim + numf_cfd;  header_out2%varname(4:n) = header_cfd%varname(4:n)

!  Clone_grid is part of the earlier surface_patch_utilities.f90.

   call clone_grid (nb_out, numf_cfd, xyzq_ideal, xyzq_out2)   ! All block dims.

!  Set up for a normalized form of the input outer envelope grid + solution.
!  The ni & nj values may have to be swapped around.  See STANDARD_CAVITY_GRIDS.

   call clone_header (header_ideal, numf_ideal, 0, header_norm)

   header_norm%varname(1) = 'X (unit cavity)'
   header_norm%varname(2) = 'Y (unit cavity)'
   header_norm%varname(3) = 'Z (unit cavity)'

   call clone_grid (nb_out, numf_ideal, xyzq_ideal, xyzq_norm) ! All block dims.

!  The CFD data mapped to a CFD-like standard grid imposed on the plain outer
!  envelope (CHT) grid has most dataset info. in common with the CFD data except
!  for the number of blocks.  First, the CFD coordinates are also divided by
!  the cavity length in standard_cavity_grids:

   header_cfd%varname(1) = 'X/L'
   header_cfd%varname(2) = 'Y/L'
   header_cfd%varname(3) = 'Z/L'

   call clone_header (header_cfd, numf_cfd, 0, header_out)

   header_out%nblocks = nb_out

!  Only the cavity faces are needed for the inner envelope derived from the
!  idealized outer envelope.  (This is still generated, but not used if the
!  recommended all-ADT mapping method is specified.)

   nb_inner = 5

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!  STANDARD_CAVITY_GRIDS is the first of numerous internal procedures appearing
!  below in order of use.

!  Impose a standardized grid with CFD-type density on the outer envelope.
!  Also, rectify the plain outer envelope grid and solution.
!  N.B.:  We now non-dimensionalize the ideal cavity first, by its length.
!         We therefore have to scale the CFD coordinates here too.

   call standard_cavity_grids () ! Generates xyzq_norm g + f and xyzq_out grid;
                                 ! xyzq_norm is rectified but not normalized yet


!  Awkward retrofit of lengthwise floor split averaging of CHT data:

   if (numf_ideal > 0) then

      ni = xyzq_norm(ib_floor)%ni
      nj = xyzq_norm(ib_floor)%nj

      allocate (floor_centers_CHT%q(numf_ideal,ni-1,nj-1,1))

      call vertices_to_centers (ni, nj, numf_ideal, xyzq_norm(ib_floor)%q,     &
                                floor_centers_CHT%q)

!     Assume the bump factors are the first functions:

      call lengthwise_floor_splits (ni, nj, numf_ideal, 1,                     &
                                    floor_centers_CHT%q, false,                &
                                    floor_averages_CHT, floor_averages_CFD)

      deallocate (floor_centers_CHT%q)

!     See subroutine map_to_envelope for the CFD floor averaging and tabulation

   else

      floor_averages_CHT(:) = 9.99

   end if


!  Set up coordinates for slicing & dicing the ideal cavity at constant X & Y:

   call set_up_cuts (nb_out, xyzq_norm)  ! In fields of slices & dices arrays;
                                         ! used below for CFD data & CHT data


!  Normalize the outer envelope grid in place and save it for possible slicing.
!  This is probably redundant now that slicing is done here, not via Tecplot.

   call normalize_patches (nb_out, xyzq_norm)

   header_norm%filename  = 'normalized.' // trim (header_ideal%filename)
   header_norm%formatted = true

   call Tecplot_write (lun_norm, header_norm, xyzq_norm, ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') &
         ' Trouble writing the normalized outer envelope solution.'
!!!   go to 99
   end if

   call deallocate_header (header_norm)
   call deallocate_blocks (1, nb_out, ndim, numf_ideal, xyzq_norm, ios)

   deallocate (xyzq_norm)

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!  Derive an inner envelope by scaling and shifting the outer envelope cavity:

   call inner_envelope_grid ()   ! Derives xyzq_inner grid from xyzq_out

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!  Unrotate/shift the CFD solution to match the idealized cavity as best we can:

   call align_cfd_grid ()        ! Moves xyzq_cfd in place

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!  Set up another variant of the standardized grid for storing and plotting
!  the forthcoming interpolations of the realigned CFD surface in real (but
!  non-dimensionalized by L) space:

   call clone_grid (nb_out, 0, xyzq_out, interp_surf)

   if (diagnostics) then

      call clone_header (header_out, numf_cfd, 0, header_surf)

      header_surf%filename  = 'interp_surf.dat'
      header_surf%formatted = true
      header_surf%title     = 'Interpolated CFD solution after realignment'

      call Tec_header_write (lun_interp, header_surf, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble writing header for interp_surf.dat.'
         go to 99
      end if

   end if

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!  Construct a search tree from all (realigned) CFD surface patches:

   nquad = 0
   do ib = 1, nb_cfd
      nquad = (xyzq_cfd(ib)%ni - 1) * (xyzq_cfd(ib)%nj - 1) + nquad
   end do

   allocate (conn(3,nquad)) ! For patch # and (i,j)

   call build_adt (nb_cfd, xyzq_cfd, nquad, conn)

   if (.not. diagnostics) write (luncrt, '(a)')

!  Tolerance for search diagnostics: use the length of the idealized cavity.

   dtolsq = (0.001 * (xyzq_ideal(3)%xmax - xyzq_ideal(1)%xmin)) ** 2

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!  Interpolate the realigned CFD data to the input idealized outer envelope.
!  This is really all that is needed for comparing with or correcting the
!  Cavity Heating Tool solution, but it has been retrofitted as an additional
!  output from the original unit cavity mapping scheme.

   call map_to_envelope ()


!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!  Perform the cavity surface portion of the mapping to the standardized grid.
!  Corresponding points of the outer and inner envelope surfaces (patches 1-5)
!  are intersected with the realigned CFD grid unless the plain ADT method
!  everywhere is indicated.

   if (map_method == 1) call map_cavity_patches ()

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!  Map the surrounding OML surface patches as for the earlier Ensemble.  If
!  map_method = 2, do the cavity patches this way as well.

   call map_OML_patches ()

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   deallocate (conn)

!! call deallocate_header (header_cfd)  ! No - needed for cloning slice headers

   call deallocate_blocks (1, nb_cfd, ndim, numf_cfd, xyzq_cfd, ios)

   deallocate (xyzq_cfd)

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!  Build another search tree from the idealized grid (CHT grid) for slicing:

   nquad = 0
   do ib = 1, nb_out
      nquad = (xyzq_ideal(ib)%ni - 1) * (xyzq_ideal(ib)%nj - 1) + nquad
   end do

   allocate (conn(3,nquad)) ! For patch # and (i,j)

   call build_adt (nb_out, xyzq_ideal, nquad, conn)

   write (luncrt, '(a)')

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!  Addition of CFD uncertainties to the line plots (only) is awkward!
!  We must insert an extra variable in the mapped/averaged xyzq_out2 structure,
!  before slicing, because the empirical uncertainties depend on being inside or
!  outside the cavity.  Figuring that out after slicing would be worse.

   call add_uncertainty ()

!  numf_cfd has been incremented by 1 now.


!  Complete the slicing of the mapped CFD data by interpolating xyzq_out2%q at
!  the slice and dice coordinates:

   call interpolate_at_cuts (nb_out, numf_cfd, xyzq_out2, nslices, slices)
   call interpolate_at_cuts (nb_out, numf_cfd, xyzq_out2,  ndices,  dices)

!  Setting up a dataset header for Tecplot 360 form of the I/O is a pain ...

   call clone_header (header_cfd, numf_cfd, 0, header_slices)

   header_slices%filename   = 'slices.' // trim (header_cfd%filename)
   header_slices%formatted  = true
   header_slices%nblocks    = nslices
   header_slices%varname(1) = 'arc/L'

!  Temporarily replace x with arc length:

   do ib = 1, nslices
      slice_xs(:,ib)      = slices(ib)%x(:,1,1)
      slices(ib)%x(:,1,1) = slice_arcs(:,ib)
   end do

   call Tecplot_write (lun_slice, header_slices, slices, ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Trouble writing the sliced CFD solution.'
      go to 99
   end if

   do ib = 1, nslices
      slices(ib)%x(:,1,1) = slice_xs(:,ib)
   end do

!  Similarly for the dices:

   call clone_header (header_cfd, numf_cfd, 0, header_dices)

   header_dices%filename   = 'dices.' // trim (header_cfd%filename)
   header_dices%formatted  = true
   header_dices%nblocks    = ndices
   header_dices%varname(2) = 'arc/L'

!  Temporarily replace y with arc length:

   do ib = 1, ndices1
      dice_1_ys(:,ib) = dices(ib)%y(:,1,1)
      dices(ib)%y(:,1,1) = dice_1_arcs(:,ib)
   end do

   do ib = 1, ndices2
      dice_2_ys(:,ib) = dices(ib+ndices1)%y(:,1,1)
      dices(ib+ndices1)%y(:,1,1) = dice_2_arcs(:,ib)
   end do

   call Tecplot_write (lun_dice, header_dices, dices, ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Trouble writing the diced CFD solution.'
      go to 99
   end if

   do ib = 1, ndices1
      dices(ib)%y(:,1,1) = dice_1_ys(:,ib)
   end do

   do ib = 1, ndices2
      dices(ib+ndices1)%y(:,1,1) = dice_2_ys(:,ib)
   end do


!  Likewise for slices of any data that came with the ideal envelope (CHT) grid:

   numf_slices = numf_cfd ! For deallocating below

   if (numf_ideal > 0) then

      if (numf_ideal /= numf_cfd) then  ! The Tecplot I/O pkg. assumes packing
         numf_slices = numf_ideal
         do ib = 1, nslices
            deallocate (slices(ib)%q)
            ni = slices(ib)%ni
            allocate   (slices(ib)%q(numf_ideal,ni,1,1))
         end do
         do ib = 1, ndices
            deallocate (dices(ib)%q)
            ni = dices(ib)%ni
            allocate   (dices(ib)%q(numf_ideal,ni,1,1))
         end do
      end if

      call interpolate_at_cuts (nb_out, numf_ideal, xyzq_ideal, nslices, slices)
      call interpolate_at_cuts (nb_out, numf_ideal, xyzq_ideal,  ndices,  dices)

      call deallocate_header (header_slices)
      call clone_header (header_ideal, numf_ideal, 0, header_slices)

      header_slices%filename   = 'slices.' // trim (header_ideal%filename)
      header_slices%formatted  = true
      header_slices%nblocks    = nslices
      header_slices%varname(1) = 'arc/L'

!     Replace x with s temporarily:

      do ib = 1, nslices
         slice_xs(:,ib)      = slices(ib)%x(:,1,1)
         slices(ib)%x(:,1,1) = slice_arcs(:,ib)
      end do

      call Tecplot_write (lun_slice, header_slices, slices, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble writing the sliced CHT solution.'
         go to 99
      end if

!     Repeat for CHT dices:

      call deallocate_header (header_dices)
      call clone_header (header_ideal, numf_ideal, 0, header_dices)

      header_dices%filename   = 'dices.' // trim (header_ideal%filename)
      header_dices%formatted  = true
      header_dices%nblocks    = ndices
      header_dices%varname(2) = 'arc/L'

!     Replace y with s temporarily:

      do ib = 1, ndices1
!!!      dice_1_ys(:,ib) = dices(ib)%y(:,1,1)
         dices(ib)%y(:,1,1) = dice_1_arcs(:,ib)
      end do

      do ib = 1, ndices2
!!!      dice_2_ys(:,ib) = dices(ib+ndices1)%y(:,1,1)
         dices(ib+ndices1)%y(:,1,1) = dice_2_arcs(:,ib)
      end do

      call Tecplot_write (lun_dice, header_dices, dices, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble writing the diced CHT solution.'
         go to 99
      end if

   end if

   call deallocate_header (header_cfd)
   call deallocate_header (header_slices)
   call deallocate_blocks (1, nslices, ndim, numf_slices, slices, ios)
   call deallocate_header (header_dices)
   call deallocate_blocks (1,  ndices, ndim, numf_slices,  dices, ios)

   deallocate (slices, dices, conn)
   deallocate (slice_arcs, slice_xs)
   deallocate (dice_1_arcs, dice_1_ys, dice_2_arcs, dice_2_ys)


!  Save the non-dimensionalized form of the ideal cavity [+ solution]:

   header_ideal%filename = 'nondimensional.' // trim (header_ideal%filename)

   call Tecplot_write (lun_ideal, header_ideal, xyzq_ideal, ios)

   if (ios /= 0) then
      write (luncrt, '(a)') ' Trouble writing nondimensionalized CHT data.'
   end if

   call deallocate_header (header_ideal)
   call deallocate_blocks (1, nb_out, ndim, numf_ideal, xyzq_ideal, ios)

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!  Finally, normalize the grid that goes with the mapped surface solution:

   call normalize_patches (nb_out, xyzq_out)

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   header_out%filename  = 'unit_cavity.' // trim (header_cfd%filename)
   header_out%formatted = true
   header_out%title     = 'Unit cavity form of mapped ' // &
                          trim (header_cfd%filename)
   header_out%varname(1) = 'X (unit cavity)'
   header_out%varname(2) = 'Y (unit cavity)'
   header_out%varname(3) = 'Z (unit cavity)'

   call Tecplot_write (lun_out, header_out, xyzq_out, ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Trouble writing the normalized mapped soln.'
      go to 99
   end if

   call deallocate_header (header_out)
   call deallocate_blocks (1, nb_out, ndim, numf_cfd, xyzq_out, ios)

   deallocate (xyzq_out)


99 continue

!  Internal procedures for program cavity_map in the order they're used:

   contains

!     --------------------------------------------------------------------------
      subroutine standard_cavity_grids ()
!     --------------------------------------------------------------------------

!     Construct a standard surface grid with CFD-type density on the outer
!     envelope patches, which are assumed to be full scale but are aligned
!     with the coordinate system.  N.B.:  We now non-dimensionalize the
!     outer envelope coordinates by its cavity length before proceeding so
!     that all related output coordinates are relative to the cavity length.
!     We therefore also have to nondimensionalize the CFD data likewise.
!
!               Y            Outer envelope starts with X = 0 upstream.
!                \           Z = 0 on the OML portion;  Z > 0 on the floor.
!                /|          Imposed grid dimensions:
!               X Z          l down side walls,  m streamwise,  n cross-stream
!
!     Patch    Location                         Dimensions
!
!       1      Front wall (upstream)               n x l
!       2      Side wall  (Y < 0)                  m x l
!       3      Back wall  (downstream)             n x l
!       4      Side wall  (Y > 0)                  m x l
!       5      Floor                               m x n
!       6      OML, adjacent to cavity (Y < 0)     m x n
!       7      OML, downstream (Y < 0)             m x n
!       8      OML, downstream of cavity           m x n
!       9      OML, downstream (Y > 0)             m x n
!      10      OML, adjacent to cavity (Y > 0)     m x n
!
!     Strategy:  Transfer appropriate corners to the denser patches.
!                All grid lines are straight, so constructing the edges
!                and interior points does not need TFI techniques.
!
!     A fully normalized form of this grid will become the "unit" cavity
!     inherited from the earlier ENSEMBLE (but probably not used now).
!
!     Retrofit:  Rectify the plain outer envelope too, ready for normalizing.

!     --------------------------------------------------------------------------

!     Local constants:

      integer, parameter :: l = 49, m = 65, n = 33

!     Local variables:

      integer :: i, ib, ier, j, j1, j2, ni, nj
      integer :: ilu, ild, iru, ird, jlu, jld, jru, jrd
      real    :: d1, d2, t(m)
      real, dimension (10) :: d1i, d2i, d1j, d2j

!     Patch:       1     2     3     4     5     6     7     8     9    10

      data d1i /.002, .002, .002, .002, .002, .002, .002, .002, .002, .002/
      data d2i /.002, .002, .002, .002, .002, .002, .050, .050, .050, .002/
      data d1j /.002, .002, .002, .002, .002, .050, .050, .002, .002, .002/
      data d2j /.002, .002, .002, .002, .002, .002, .002, .002, .050, .050/

!     Execution:

!     Calculate the data ranges of the outer envelope patches:

      do ib = 1, nb_out
         call patch_data_range (xyzq_ideal(ib))
         if (diagnostics) then
            if (ib == 1) write (luncrt, '(a)')
            write (luncrt, '(a, i3, 6f14.6)')                   &
               ' Data range of outer envelope (CHT) patch', ib, &
               xyzq_ideal(ib)%xmin, xyzq_ideal(ib)%xmax,        &
               xyzq_ideal(ib)%ymin, xyzq_ideal(ib)%ymax,        &
               xyzq_ideal(ib)%zmin, xyzq_ideal(ib)%zmax
         end if
      end do

!     Non-dimensionalize the ideal geometry by its cavity length:

      cavity_length = xyzq_ideal(6)%xmax - xyzq_ideal(6)%xmin
      cavity_depth  = xyzq_ideal(5)%zmax - xyzq_ideal(6)%zmin

!     Bump factor uncertainty parameters (see main program constants):

      short_cavity = (cavity_length / cavity_depth) < short_test

      if (short_cavity) then
         uncertainty_floor = uncertainty_short
      else
         uncertainty_floor = uncertainty_long
      end if

      do ib = 1, nb_out
         xyzq_ideal(ib)%x = xyzq_ideal(ib)%x / cavity_length
         xyzq_ideal(ib)%y = xyzq_ideal(ib)%y / cavity_length
         xyzq_ideal(ib)%z = xyzq_ideal(ib)%z / cavity_length

!        Play safe - recover the precise new data range:

         call patch_data_range (xyzq_ideal(ib))

      end do

!     Non-dimensionalize the CFD coordinates by the ideal cavity length too:

      do ib = 1, nb_cfd
         call patch_data_range (xyzq_cfd(ib))
         if (diagnostics) then
            if (ib == 1) write (luncrt, '(a)')
            write (luncrt, '(a, i3, 6f15.8)')                &
               ' Data range of CFD surface patch', ib,        &
               xyzq_cfd(ib)%xmin, xyzq_cfd(ib)%xmax,          &
               xyzq_cfd(ib)%ymin, xyzq_cfd(ib)%ymax,          &
               xyzq_cfd(ib)%zmin, xyzq_cfd(ib)%zmax
         end if
         xyzq_cfd(ib)%x = xyzq_cfd(ib)%x / cavity_length
         xyzq_cfd(ib)%y = xyzq_cfd(ib)%y / cavity_length
         xyzq_cfd(ib)%z = xyzq_cfd(ib)%z / cavity_length
      end do

      allocate (xyzq_out(nb_out)) ! Standardized dense patches

      xyzq_out( 1)%zone_title = 'Front wall'
      xyzq_out( 2)%zone_title = 'Side wall'
      xyzq_out( 3)%zone_title = 'Back wall'
      xyzq_out( 4)%zone_title = 'Side wall'
      xyzq_out( 5)%zone_title = 'Floor'
      xyzq_out( 6)%zone_title = 'OML alongside, left'
      xyzq_out( 7)%zone_title = 'OML downstream, left'
      xyzq_out( 8)%zone_title = 'OML downstream, mid'
      xyzq_out( 9)%zone_title = 'OML downstream, right'
      xyzq_out(10)%zone_title = 'OML alongside, right'

      do ib = 1, 4            ! Side walls
         xyzq_out(ib)%ni = m  ! Except for front and back
         xyzq_out(ib)%nj = l
      end do

      xyzq_out(1)%ni = n;  xyzq_out(3)%ni = n  ! Front and back

      do ib = 5, nb_out
         xyzq_out(ib)%ni = m
         xyzq_out(ib)%nj = n
      end do

      do ib = 1, nb_out
         xyzq_out(ib)%nk = 1
         xyzq_out(ib)%nzoneaux = 0
         xyzq_out(ib)%solutiontime = undefined

         call Tec_block_allocate (xyzq_out(ib), ndim, numf_cfd, ier)

         if (ier /= 0) then
            write (luncrt, '(a, i3)') &
               ' Trouble allocating standardized dense patch #', ib
            stop
         end if

         xyzq_norm(ib)%zone_title   = xyzq_out(ib)%zone_title
         xyzq_norm(ib)%nzoneaux     = 0
         xyzq_norm(ib)%solutiontime = undefined
      end do

!     All the desired standardized patches are defined by their corners
!     and by the edge distributions implied by d1i, d2i, d1j, and d2j and
!     the dimensions l, m, or n.
!     Therefore, the main complication is allowing for arbitrary (i,j)
!     indexing in the given outer envelope patches.
!
!     Saving of common edge info. has been retrofitted to enable averaging of
!     interpolated CFD along common edges (forced by the control volume scheme
!     for averaging CFD data over the coarse CHT grid cells then moving results
!     from cell centers to cell vertices).

      do ib = 1, nb_out

         call identify_input_corners (ib, xyzq_ideal(ib), &  ! l, r, u, d mean
                                      ilu, ild, iru, ird, &  ! left-upstream,
                                      jlu, jld, jru, jrd)    ! right downstrm ..

         call save_common_edge_info  (ib, nedges, medges, &  ! 15 edge pairs
                                      ilu, ild, iru, ird, &  ! defined in medges
                                      jlu, jld, jru, jrd)    ! (1:5,1:2,1:15)

         call rectify_idealized_grid (ib, numf_ideal, xyzq_ideal(ib), &
                                      ilu, ild, iru, ird, &
                                      jlu, jld, jru, jrd, &
                                      xyzq_norm(ib))

         call assign_output_corners  (ib, xyzq_ideal(ib), &  ! For densified g.
                                      ilu, ild, iru, ird, &
                                      jlu, jld, jru, jrd, &
                                      xyzq_out(ib))

         call trapezium_grid (xyzq_out(ib), d1i(ib), d2i(ib), d1j(ib), d2j(ib))

      end do

      if (diagnostics) then

         open (lun_out, file='regridded_outer_envelope.p3da', status='unknown')

         write (lun_out, '(i2)') nb_out
         write (lun_out, '(3i4)') &
            (xyzq_out(i)%ni, xyzq_out(i)%nj, 1, i = 1, nb_out)
         do i = 1, nb_out
            write (lun_out, '(1p, 6e19.11)') &
               xyzq_out(i)%x, xyzq_out(i)%y, xyzq_out(i)%z
         end do

         close (lun_out)

      end if

      end subroutine standard_cavity_grids

!     --------------------------------------------------------------------------

      subroutine identify_input_corners (ib, g,              &
                                         ilu, ild, iru, ird, &
                                         jlu, jld, jru, jrd)

!     Identify the left-upstream, left-downstream, and corresponding right
!     corner indices of the given surface patch g.
!     Surely there's an easier way than the awkward logic coded here!

!     --------------------------------------------------------------------------

!     Arguments:

      integer,          intent (in)  :: ib  ! Patch number
      type (grid_type), intent (in)  :: g   ! Patch grid
      integer,          intent (out) :: &
         ilu, ild, iru, ird,            &   ! Corner indices in g(i,j,1)
         jlu, jld, jru, jrd

!     Local constants:

      integer, parameter :: indef = 999
      real,    parameter :: big   = 1.e+10, half = 0.5, tol = 0.001

!     Local variables:

      integer :: i, j, ni, nj
      real    :: xmax, xmin, ymax, ymin, xmean, ymean, zmean

!     Execution:

      ilu = indef;  iru = indef;  ni = g%ni
      ild = indef;  ird = indef;  nj = g%nj

      xmean = half * (g%xmin + g%xmax)
      ymean = half * (g%ymin + g%ymax)
      zmean = half * (g%zmin + g%zmax)

      select case (ib)

         case (1) ! Upstream wall

            do j = 1, nj, nj - 1
               do i = 1, ni, ni - 1
                  if (g%z(i,j,1) < zmean) then  ! Corner on the OML
                     if (g%y(i,j,1) < ymean) then
                        ilu = i;  jlu = j
                     else
                        iru = i;  jru = j
                     end if
                  end if
               end do
            end do

            ymin = big;  ymax = -big  ! Guard against a floor way off-center
            do j = 1, nj, nj - 1
               do i = 1, ni, ni - 1
                  if (g%z(i,j,1) > g%zmin + tol) then
                     if (g%y(i,j,1) < ymin) then      ! Min. over floor corners
                        ymin = g%y(i,j,1);  ild = i;  jld = j
                     end if
                     if (g%y(i,j,1) > ymax) then
                        ymax = g%y(i,j,1);  ird = i;  jrd = j
                     end if
                  end if
               end do
            end do

         case (2) ! Left side wall

            do j = 1, nj, nj - 1
               do i = 1, ni, ni - 1
                  if (g%z(i,j,1) < zmean) then
                     if (g%x(i,j,1) < xmean) then
                        ilu = i;  jlu = j
                     else
                        ild = i;  jld = j
                     end if
                  end if
               end do
            end do

            xmin = big;  xmax = -big
            do j = 1, nj, nj - 1
               do i = 1, ni, ni - 1
                  if (g%z(i,j,1) > g%zmin + tol) then
                     if (g%x(i,j,1) < xmin) then
                        xmin = g%x(i,j,1);  iru = i;  jru = j
                     end if
                     if (g%x(i,j,1) > xmax) then
                        xmax = g%x(i,j,1);  ird = i;  jrd = j
                     end if
                  end if
               end do
            end do

         case (3) ! Downstream wall

            do j = 1, nj, nj - 1
               do i = 1, ni, ni - 1
                  if (g%z(i,j,1) < zmean) then
                     if (g%y(i,j,1) < ymean) then
                        ild = i;  jld = j
                     else
                        ird = i;  jrd = j
                     end if
                  end if
               end do
            end do

            ymin = big;  ymax = -big
            do j = 1, nj, nj - 1
               do i = 1, ni, ni - 1
                  if (g%z(i,j,1) > g%zmin + tol) then
                     if (g%y(i,j,1) < ymin) then
                        ymin = g%y(i,j,1);  ilu = i;  jlu = j
                     end if
                     if (g%y(i,j,1) > ymax) then
                        ymax = g%y(i,j,1);  iru = i;  jru = j
                     end if
                  end if
               end do
            end do

         case (4) ! Right side wall

            do j = 1, nj, nj - 1
               do i = 1, ni, ni - 1
                  if (g%z(i,j,1) < zmean) then
                     if (g%x(i,j,1) < xmean) then
                        iru = i;  jru = j
                     else
                        ird = i;  jrd = j
                     end if
                  end if
               end do
            end do

            xmin = big;  xmax = -big
            do j = 1, nj, nj - 1
               do i = 1, ni, ni - 1
                  if (g%z(i,j,1) > g%zmin + tol) then
                     if (g%x(i,j,1) < xmin) then
                        xmin = g%x(i,j,1);  ilu = i;  jlu = j
                     end if
                     if (g%x(i,j,1) > xmax) then
                        xmax = g%x(i,j,1);  ild = i;  jld = j
                     end if
                  end if
               end do
            end do

         case (5:10) ! Floor and all OML patches

            do j = 1, nj, nj - 1
               do i = 1, ni, ni - 1
                  if (g%x(i,j,1) < xmean) then
                     if (g%y(i,j,1) < ymean) then
                        ilu = i;  jlu = j
                     else
                        iru = i;  jru = j
                     end if
                  else
                     if (g%y(i,j,1) < ymean) then
                        ild = i;  jld = j
                     else
                        ird = i;  jrd = j
                     end if
                  end if
               end do
            end do

      end select

      if (ilu == indef .or. iru == indef .or. &
          ild == indef .or. ird == indef) then
         write (luncrt, '(/, a, i3)') &
            ' Unable to identify corners of outer envelope patch', ib
         stop
      end if

      end subroutine identify_input_corners

!     --------------------------------------------------------------------------

      subroutine save_common_edge_info (ib,  nedges, m,     &
                                        ilu, ild, iru, ird, &
                                        jlu, jld, jru, jrd)

!     Insert the common edge pointers for patch ib into the appropriate slots.
!     We assume we know which patch is which, but we're handling unknown index
!     ordering in the input CHT grid.  We're stuck with that grid (can't rectify
!     it in the output mapped CFD data), so it doesn't help that we are also
!     standardizing the grid indexing for other reasons.   :(
!
!     The edges are defined in the following particular order from the pointers
!     set up by identify_input_corners:
!
!        Edge #   Patch 1   i1  i2  j1  j2    Patch 2   i1  i2  j1  j2
!          1          1    ild ilu jld jlu        2    iru ilu jru jlu
!          2          2    ird ild jrd jld        3    ilu ild jlu jld
!          3          3    iru ird jru jrd        4    ild ird jld jrd
!          4          4    ilu iru jlu jru        1    ird iru krd jru
!          5          1    ild ird jld jrd        5    ilu iru jlu jru
!          6          2    iru ird jru jrd        5    ilu ild jlu jld
!          7          3    ilu iru jlu jru        5    ild ird jld jrd
!          8          4    ilu ild jlu jld        5    iru ird jru jrd
!          9          2    ilu ild jlu jld        6    iru ird jru jrd
!         10          4    iru ird jru jrd       10    ilu ild jlu jld
!         11          6    ild ird jld jrd        7    ilu iru jlu jru
!         12          3    ild ird jld jrd        8    ilu iru jlu jru
!         13         10    ild ird jld jrd        9    ilu iru jlu jru
!         14          7    iru ird jru jrd        8    ilu ild jlu jld
!         15          8    iru ird jru jrd        9    ilu ild jlu jld
!
!     The output array, medges (1:5, 1:2, 1:15) is declared in the main program
!     and used in map_to_envelope to ensure common edges have common f values.
!     Passing it as argument m(:,:,:) streamlines the code.

!     --------------------------------------------------------------------------

!     Arguments:

      integer, intent (in)    :: ib             ! CHT patch number
      integer, intent (in)    :: nedges         ! # common edges in CHT grid
      integer, intent (inout) :: m(5,2,nedges)  ! See above
      integer, intent (in)    :: &
         ilu, ild, iru, ird,     &              ! Corner indices in g(ib)(i,j,1)
         jlu, jld, jru, jrd

!     Execution:

      select case (ib)

      case (1)
        m(1,1, 1)=ib; m(2,1, 1)=ild; m(3,1, 1)=ilu; m(4,1, 1)=jld; m(5,1, 1)=jlu
        m(1,2, 4)=ib; m(2,2, 4)=ird; m(3,2, 4)=iru; m(4,2, 4)=jrd; m(5,2, 4)=jru
        m(1,1, 5)=ib; m(2,1, 5)=ild; m(3,1, 5)=ird; m(4,1, 5)=jld; m(5,1, 5)=jrd
      case (2)
        m(1,2, 1)=ib; m(2,2, 1)=iru; m(3,2, 1)=ilu; m(4,2, 1)=jru; m(5,2, 1)=jlu
        m(1,1, 2)=ib; m(2,1, 2)=ird; m(3,1, 2)=ild; m(4,1, 2)=jrd; m(5,1, 2)=jld
        m(1,1, 6)=ib; m(2,1, 6)=iru; m(3,1, 6)=ird; m(4,1, 6)=jru; m(5,1, 6)=jrd
        m(1,1, 9)=ib; m(2,1, 9)=ilu; m(3,1, 9)=ild; m(4,1, 9)=jlu; m(5,1, 9)=jld
      case (3)
        m(1,2, 2)=ib; m(2,2, 2)=ilu; m(3,2, 2)=ild; m(4,2, 2)=jlu; m(5,2, 2)=jld
        m(1,1, 3)=ib; m(2,1, 3)=iru; m(3,1, 3)=ird; m(4,1, 3)=jru; m(5,1, 3)=jrd
        m(1,1, 7)=ib; m(2,1, 7)=ilu; m(3,1, 7)=iru; m(4,1, 7)=jlu; m(5,1, 7)=jru
        m(1,1,12)=ib; m(2,1,12)=ild; m(3,1,12)=ird; m(4,1,12)=jld; m(5,1,12)=jrd
      case (4)
        m(1,2, 3)=ib; m(2,2, 3)=ild; m(3,2, 3)=ird; m(4,2, 3)=jld; m(5,2, 3)=jrd
        m(1,1, 4)=ib; m(2,1, 4)=ilu; m(3,1, 4)=iru; m(4,1, 4)=jlu; m(5,1, 4)=jru
        m(1,1, 8)=ib; m(2,1, 8)=ilu; m(3,1, 8)=ild; m(4,1, 8)=jlu; m(5,1, 8)=jld
        m(1,1,10)=ib; m(2,1,10)=iru; m(3,1,10)=ird; m(4,1,10)=jru; m(5,1,10)=jrd
      case (5)
        m(1,2, 5)=ib; m(2,2, 5)=ilu; m(3,2, 5)=iru; m(4,2, 5)=jlu; m(5,2, 5)=jru
        m(1,2, 6)=ib; m(2,2, 6)=ilu; m(3,2, 6)=ild; m(4,2, 6)=jlu; m(5,2, 6)=jld
        m(1,2, 7)=ib; m(2,2, 7)=ild; m(3,2, 7)=ird; m(4,2, 7)=jld; m(5,2, 7)=jrd
        m(1,2, 8)=ib; m(2,2, 8)=iru; m(3,2, 8)=ird; m(4,2, 8)=jru; m(5,2, 8)=jrd
      case (6)
        m(1,2, 9)=ib; m(2,2, 9)=iru; m(3,2, 9)=ird; m(4,2, 9)=jru; m(5,2, 9)=jrd
        m(1,1,11)=ib; m(2,1,11)=ild; m(3,1,11)=ird; m(4,1,11)=jld; m(5,1,11)=jrd
      case (7)
        m(1,2,11)=ib; m(2,2,11)=ilu; m(3,2,11)=iru; m(4,2,11)=jlu; m(5,2,11)=jru
        m(1,1,14)=ib; m(2,1,14)=iru; m(3,1,14)=ird; m(4,1,14)=jru; m(5,1,14)=jrd
      case (8)
        m(1,2,12)=ib; m(2,2,12)=ilu; m(3,2,12)=iru; m(4,2,12)=jlu; m(5,2,12)=jru
        m(1,2,14)=ib; m(2,2,14)=ilu; m(3,2,14)=ild; m(4,2,14)=jlu; m(5,2,14)=jld
        m(1,1,15)=ib; m(2,1,15)=iru; m(3,1,15)=ird; m(4,1,15)=jru; m(5,1,15)=jrd
      case (9)
        m(1,2,13)=ib; m(2,2,13)=ilu; m(3,2,13)=iru; m(4,2,13)=jlu; m(5,2,13)=jru
        m(1,2,15)=ib; m(2,2,15)=ilu; m(3,2,15)=ild; m(4,2,15)=jlu; m(5,2,15)=jld
      case (10)
        m(1,2,10)=ib; m(2,2,10)=ilu; m(3,2,10)=ild; m(4,2,10)=jlu; m(5,2,10)=jld
        m(1,1,13)=ib; m(2,1,13)=ild; m(3,1,13)=ird; m(4,1,13)=jld; m(5,1,13)=jrd
      end select

      end subroutine save_common_edge_info

!     --------------------------------------------------------------------------

      subroutine rectify_idealized_grid (ib, nf, g,          &
                                         ilu, ild, iru, ird, &
                                         jlu, jld, jru, jrd, &
                                         h)

!     Ensure that the indexing of the input outer envelope patch is in standard
!     order, ready for normalizing. This was adapted from ASSIGN_OUTPUT_CORNERS.
!     Passing the fourth corner indices is redundant.

!     --------------------------------------------------------------------------

!     Arguments:

      integer,          intent (in)    :: ib  ! Patch number
      integer,          intent (in)    :: nf  ! # functions with outer envelope
      type (grid_type), intent (in)    :: g   ! Outer envelope patch grid
      integer,          intent (in)    :: &
         ilu, ild, iru, ird,              &   ! Corner indices in g(i,j,1)
         jlu, jld, jru, jrd
      type (grid_type), intent (inout) :: h   ! Corresp. patch of rectified
                                              ! grid and solution
!     Execution:

      select case (ib)

         case (1) ! Upstream wall

            call rectify_patch (g, nf, ild, jld, ird, jrd, ilu, jlu, h)

         case (2) ! Left side wall

            call rectify_patch (g, nf, iru, jru, ird, jrd, ilu, jlu, h)

         case (3) ! Downstream wall

            call rectify_patch (g, nf, ilu, jlu, iru, jru, ild, jld, h)

         case (4:10) ! Right side wall, floor, and all OML patches

            call rectify_patch (g, nf, ilu, jlu, ild, jld, iru, jru, h)

         end select

      end subroutine rectify_idealized_grid

!     --------------------------------------------------------------------------

      subroutine rectify_patch (g, nf, i1, j1, i2, j2, i3, j3, h)

!     Rectify the given patch to a known (i,j) order.

!    --------------------------------------------------------------------------

!     Arguments:

      type (grid_type), intent (in)    :: g   ! Outer envelope patch grid
      integer,          intent (in)    :: nf  ! # functions with outer envelope
      integer,          intent (in)    :: &
         i1, j1, i2, j2, i3, j3               ! Corner indices in g(i,j,1)
                                              ! corresp. to (1,1), (mi,1),
                                              ! (1,mj) & (mi,mj) in h(i,j,1);
                                              ! the 4th corner is redundant
      type (grid_type), intent (inout) :: h   ! Corresp. standardized patch

!     Local variables:

      integer :: ig, ih, inc, jg, jh, jnc, m, ni, nj

!     Execution:

      ni = g%ni;  nj = g%nj

      if (j1 == j2) then

         h%ni = ni;  h%mi = ni
         h%nj = nj;  h%mj = nj
         h%nk =  1;  h%mk =  1

         allocate (h%x(ni,nj,1), h%y(ni,nj,1), h%z(ni,nj,1))
         if (nf > 0) allocate (h%q(nf,ni,nj,1))

         inc = sign (1, i2 - i1)
         jnc = sign (1, j3 - j1)
         jh  = 0
         do jg = j1, j3, jnc
            jh = jh + 1
            ih = 0
            do ig = i1, i2, inc
               ih = ih + 1
               h%x(ih,jh,1) = g%x(ig,jg,1)
               h%y(ih,jh,1) = g%y(ig,jg,1)
               h%z(ih,jh,1) = g%z(ig,jg,1)
               do m = 1, nf
                  h%q(m,ih,jh,1) = g%q(m,ig,jg,1)
               end do
            end do
         end do

      else ! i1 = i2

         h%ni = nj;  h%mi = nj
         h%nj = ni;  h%mj = ni
         h%nk =  1;  h%mk =  1

         allocate (h%x(nj,ni,1), h%y(nj,ni,1), h%z(nj,ni,1))
         if (nf > 0) allocate (h%q(nf,nj,ni,1))

         inc = sign (1, i3 - i1)
         jnc = sign (1, j2 - j1)
         jh  = 0
         do ig = i1, i3, inc
            jh = jh + 1
            ih = 0
            do jg = j1, j2, jnc
               ih = ih + 1
               h%x(ih,jh,1) = g%x(ig,jg,1)
               h%y(ih,jh,1) = g%y(ig,jg,1)
               h%z(ih,jh,1) = g%z(ig,jg,1)
               do m = 1, nf
                  h%q(m,ih,jh,1) = g%q(m,ig,jg,1)
               end do
            end do
         end do

      end if

      end subroutine rectify_patch

!     --------------------------------------------------------------------------

      subroutine assign_output_corners (ib, g,              &
                                        ilu, ild, iru, ird, &
                                        jlu, jld, jru, jrd, &
                                        h)

!     Transfer the appropriate corner coordinates from the outer envelope patch
!     to the corresponding patch to be regridded with CFD-type density.

!     --------------------------------------------------------------------------

!     Arguments:

      integer,          intent (in)    :: ib  ! Patch number
      type (grid_type), intent (in)    :: g   ! Outer envelope patch grid
      integer,          intent (in)    :: &
         ilu, ild, iru, ird,              &   ! Corner indices in g(i,j,1)
         jlu, jld, jru, jrd
      type (grid_type), intent (inout) :: h   ! Corresp. patch of standardized
                                              ! grid
!     Local variables:

      integer :: mi, mj, ni, nj

!     Execution:

      ni = g%ni;  nj = g%nj
      mi = h%ni;  mj = h%nj

      select case (ib)

         case (1) ! Upstream wall

            h%x( 1, 1,1) = g%x(ild,jld,1);  h%x(mi, 1,1) = g%x(ird,jrd,1)
            h%y( 1, 1,1) = g%y(ild,jld,1);  h%y(mi, 1,1) = g%y(ird,jrd,1)
            h%z( 1, 1,1) = g%z(ild,jld,1);  h%z(mi, 1,1) = g%z(ird,jrd,1)

            h%x( 1,mj,1) = g%x(ilu,jlu,1);  h%x(mi,mj,1) = g%x(iru,jru,1)
            h%y( 1,mj,1) = g%y(ilu,jlu,1);  h%y(mi,mj,1) = g%y(iru,jru,1)
            h%z( 1,mj,1) = g%z(ilu,jlu,1);  h%z(mi,mj,1) = g%z(iru,jru,1)

         case (2) ! Left side wall

            h%x( 1, 1,1) = g%x(iru,jru,1);  h%x(mi, 1,1) = g%x(ird,jrd,1)
            h%y( 1, 1,1) = g%y(iru,jru,1);  h%y(mi, 1,1) = g%y(ird,jrd,1)
            h%z( 1, 1,1) = g%z(iru,jru,1);  h%z(mi, 1,1) = g%z(ird,jrd,1)

            h%x( 1,mj,1) = g%x(ilu,jlu,1);  h%x(mi,mj,1) = g%x(ild,jld,1)
            h%y( 1,mj,1) = g%y(ilu,jlu,1);  h%y(mi,mj,1) = g%y(ild,jld,1)
            h%z( 1,mj,1) = g%z(ilu,jlu,1);  h%z(mi,mj,1) = g%z(ild,jld,1)

         case (3) ! Downstream wall

            h%x( 1, 1,1) = g%x(ilu,jlu,1);  h%x(mi, 1,1) = g%x(iru,jru,1)
            h%y( 1, 1,1) = g%y(ilu,jlu,1);  h%y(mi, 1,1) = g%y(iru,jru,1)
            h%z( 1, 1,1) = g%z(ilu,jlu,1);  h%z(mi, 1,1) = g%z(iru,jru,1)

            h%x( 1,mj,1) = g%x(ild,jld,1);  h%x(mi,mj,1) = g%x(ird,jrd,1)
            h%y( 1,mj,1) = g%y(ild,jld,1);  h%y(mi,mj,1) = g%y(ird,jrd,1)
            h%z( 1,mj,1) = g%z(ild,jld,1);  h%z(mi,mj,1) = g%z(ird,jrd,1)

         case (4:10) ! Right side wall, floor, and all OML patches

            h%x( 1, 1,1) = g%x(ilu,jlu,1);  h%x(mi, 1,1) = g%x(ild,jld,1)
            h%y( 1, 1,1) = g%y(ilu,jlu,1);  h%y(mi, 1,1) = g%y(ild,jld,1)
            h%z( 1, 1,1) = g%z(ilu,jlu,1);  h%z(mi, 1,1) = g%z(ild,jld,1)

            h%x( 1,mj,1) = g%x(iru,jru,1);  h%x(mi,mj,1) = g%x(ird,jrd,1)
            h%y( 1,mj,1) = g%y(iru,jru,1);  h%y(mi,mj,1) = g%y(ird,jrd,1)
            h%z( 1,mj,1) = g%z(iru,jru,1);  h%z(mi,mj,1) = g%z(ird,jrd,1)

         end select

      end subroutine assign_output_corners

!     --------------------------------------------------------------------------

      subroutine trapezium_grid (g, d1i, d2i, d1j, d2j)

!     Specialized 3-space planar grid face generator, where the corners are in
!     place and opposite edges have the same Vinokur-type distribution.  There
!     is no need to use the more general TFIQ3D, which recomputes edge arc
!     lengths unnecessarily.

!     --------------------------------------------------------------------------

!     Arguments:

      type (grid_type), intent (inout) :: g    ! Surface patch with dimensions
                                               ! and corner points in place
      real, intent (in) :: d1i, d2i, d1j, d2j  ! End-point increments on [0, 1]
                                               ! for i & j Vinokur distributions
!     Local constants:

      real, parameter :: one = 1., zero = 0.

!     Local variables:

      integer :: i, j, ier, ni, nj
      real    :: r, rm1
      real, allocatable, dimension (:) :: t

!     Execution:

      ni = g%ni;  nj = g%nj

      allocate (t(ni))

!     Vinokur expects the first and last point to be in place.

      t(1)  = zero
      t(ni) = one;   call vinokur (1, ni, d1i, d2i, t, luncrt, ier)

!     Edge distributions in the i direction:

      do i = 2, ni - 1
         r = t(i);  rm1 = one - r
         g%x(i, 1,1) = rm1 * g%x(1, 1,1) + r * g%x(ni, 1,1)
         g%y(i, 1,1) = rm1 * g%y(1, 1,1) + r * g%y(ni, 1,1)
         g%z(i, 1,1) = rm1 * g%z(1, 1,1) + r * g%z(ni, 1,1)

         g%x(i,nj,1) = rm1 * g%x(1,nj,1) + r * g%x(ni,nj,1)
         g%y(i,nj,1) = rm1 * g%y(1,nj,1) + r * g%y(ni,nj,1)
         g%z(i,nj,1) = rm1 * g%z(1,nj,1) + r * g%z(ni,nj,1)
      end do

      deallocate (t);  allocate (t(nj))

!     All interior lines in the j direction:

      t(1)  = zero
      t(nj) = one;   call vinokur (1, nj, d1j, d2j, t, luncrt, ier)

      do j = 2, nj - 1
         r = t(j);  rm1 = one - r
         do i = 1, ni
            g%x(i,j,1) = rm1 * g%x(i,1,1) + r * g%x(i,nj,1)
            g%y(i,j,1) = rm1 * g%y(i,1,1) + r * g%y(i,nj,1)
            g%z(i,j,1) = rm1 * g%z(i,1,1) + r * g%z(i,nj,1)
         end do
      end do

      deallocate (t)

      end subroutine trapezium_grid

!     --------------------------------------------------------------------------

      subroutine set_up_cuts (nb, g)

!     Simple-minded set-up of (x,y,z) coordinates for cuts at a handful of span-
!     wise and streamwise stations of the rectified cavity floor, including the
!     walls and surrounding patches.  Passing the rectified outer envelope grid
!     (xyzq_norm) as an argument cleans up the code.

!     --------------------------------------------------------------------------

!     Arguments:

      integer,          intent (in)    :: nb     ! # patches in outer envelope
      type (grid_type), intent (inout) :: g(nb)  ! Rectified outer envelope grid
                                                 ! (temporarily transposed here)
!     Local constants:

      real,      parameter :: half    = 0.50
      real,      parameter :: fourth  = 0.25
      character, parameter :: blank*1 = ' '

!     Local variables:

      integer :: i, ib, ii, ip, is, it, j, jinc, ni, nj, nlines, npts
      integer :: npts_per_dice_1, npts_per_dice_2, npts_per_slice
      integer :: ib_dices_1(5), ib_dices_2(3), ib_slices(4)

      real    :: average_length, average_width, dx, dy, x_centroid, y_centroid
      real    :: station, total, x, x0, y, y0
      real    :: xstations(ndices), ystations(nslices)

!     Storage:

      data ib_dices_1 /6, 2, 5, 4, 10/, &     ! Patch #s for cavity cross-cuts
           ib_dices_2 /7, 8, 9/,        &     ! ... and for dices aft of cavity
           ib_slices  /1, 5, 3, 8/            ! ... and streamwise slices

!     Execution:

!     Count the number of CHT grid points along the cuts, avoiding duplicates:

      npts_per_slice  = g(1)%nj + g(3)%nj + g(5)%ni + g(8)%ni - 3
      npts_per_dice_1 = g(2)%nj + g(4)%nj + g(5)%nj + g(6)%nj + g(10)%nj - 4
      npts_per_dice_2 = g(7)%nj + g(8)%nj + g(9)%nj - 2

      allocate (slices(nslices), dices(ndices))
      allocate (dice_1_arcs(npts_per_dice_1,ndices1))
      allocate (dice_1_ys  (npts_per_dice_1,ndices1))
      allocate (dice_2_arcs(npts_per_dice_2,ndices2))
      allocate (dice_2_ys  (npts_per_dice_2,ndices2))
      allocate (slice_arcs (npts_per_slice, nslices))
      allocate (slice_xs   (npts_per_slice, nslices))

!     %q is sized to the CFD solution here; it may have to be reallocated
!     for the CHT solution if there is one.
!     Allow for inserting a BF + uncertainty variable later.

      do ib = 1, nslices
         slices(ib)%ni = npts_per_slice
         slices(ib)%nj = 1
         slices(ib)%nk = 1
         slices(ib)%nzoneaux = 0
         call Tec_block_allocate (slices(ib), ndim, numf_cfd + 1, ios)
         if (ios /= 0) then
            write (luncrt, '(a, i3)') &
               ' Trouble allocating CFD slice #', ib
            stop
         end if
      end do

      do ib = 1, ndices
         if (ib <= ndices1) dices(ib)%ni = npts_per_dice_1
         if (ib >  ndices1) dices(ib)%ni = npts_per_dice_2
         dices(ib)%nj = 1
         dices(ib)%nk = 1
         dices(ib)%nzoneaux = 0
         call Tec_block_allocate (dices(ib), ndim, numf_cfd + 1, ios)
         if (ios /= 0) then
            write (luncrt, '(a, i3)') &
               ' Trouble allocating CFD dice #', ib
            stop
         end if
      end do

!     Average dimensions of the cavity floor:

      ni = g(5)%ni;  nj = g(5)%nj
      average_length = ((g(5)%x(ni,1,1) + g(5)%x(ni,nj,1))  -    &
                        (g(5)%x( 1,1,1) + g(5)%x( 1,nj,1))) * half
      x_centroid     =  (g(5)%x( 1,1,1) + g(5)%x(ni, 1,1)   +    &
                         g(5)%x(1,nj,1) + g(5)%x(ni,nj,1))  * fourth
      average_width  = ((g(5)%y(1,nj,1) + g(5)%y(ni,nj,1))  -    &
                        (g(5)%y(1, 1,1) + g(5)%y(ni, 1,1))) * half
      y_centroid     =  (g(5)%y( 1,1,1) + g(5)%y(ni, 1,1)   +    &
                         g(5)%y(1,nj,1) + g(5)%y(ni,nj,1))  * fourth

!     x and y stations of cuts.  Dicing downstream of the cavity needs care.
!     First, the cross-cuts centered on the cavity floor:

      dx = average_length / real (ndices1 + 1)
      x0 = x_centroid - dx * real ((ndices1 + 1)/2)
      do is = 1, ndices1
         x = x0 + dx * real (is)
         dices(is)%x(:,1,1) = x
         dices(is)%zone_title = blank
         station = real (is) / real (ndices1 + 1)
         write (dices(is)%zone_title(1:11), '(a, f4.1)') 'x/L =', station
      end do

!     Cross-cut stations aft of the cavity:

      x0 =  g(8)%x(1,1,1)
      ni =  g(8)%ni
      dx = (g(8)%x(ni,1,1) - x0) / real (ndices2 + 1)
      do is = 1, ndices2
         x = x0 + dx * real (is)
         it = is + ndices1
         dices(it)%x(:,1,1) = x
         dices(it)%zone_title = blank
         station = one + real (is) / real (ndices2 + 1)
         write (dices(it)%zone_title(1:11), '(a, f4.1)') 'x/L =', station
      end do

!     Streamwise cuts centered on the cavity floor:

      dy = average_width / real (nslices+1)
      y0 = y_centroid - dy * real ((nslices+1)/2)
      do is = 1, nslices
         y = y0 + dy * real (is)
         slices(is)%y(:,1,1) = y
         slices(is)%zone_title = blank
         station = real (is) / real (nslices + 1)
         write (slices(is)%zone_title(1:11), '(a, f5.2)') 'y/W =', station
      end do

!     Specialized location of remaining coordinates of slices.

      call transpose_patch (g(5), 0)  ! Temporarily
      call transpose_patch (g(8), 0)

      do is = 1, nslices

         ii = 1                       ! Slice point counter, downstream order
         station = slices(is)%y(1,1,1)
!!!      write (6, '(a, f8.4)') ' Slice station:', station

         do ip = 1, 4                 ! For slice patches 1, 5, 3, 8

            ib = ib_slices(ip)        ! Patch to slice
            npts   = g(ib)%ni
            nlines = g(ib)%nj
            jinc   = 1
            if (ib == 1) jinc = -1    ! Front wall; j always starts at the floor

!!!         write (6, '(a, 7i5)') ' is,ii,ip,ib,npts,nlines,jinc: ', &
!!!                                 is,ii,ip,ib,npts,nlines,jinc

            call slice_patch (ii, station, npts, nlines, jinc,                 &
                              g(ib)%y, g(ib)%z, g(ib)%x,                       &
                              slices(is)%z, slices(is)%x)

            ii = ii + nlines - 1      ! Avoid duplicating common edge points

         end do ! Next patch to slice

!        Replace x with arc length, except we can't overwrite yet.

         call chords2d (npts_per_slice, slices(is)%z, slices(is)%x, false,     &
                        total, slice_arcs(1,is))

      end do ! Next slice

      call transpose_patch (g(5), 0)  ! Restore the original indexing
      call transpose_patch (g(8), 0)

!     Likewise: specialized location of remaining coordinates of cross-cuts.
!     All rectified patches have i in the downstream direction.
!     First, the dices centered on the cavity floor:

      do is = 1, ndices1

         ii = 1                       ! Dice point counter, left to right
         station = dices(is)%x(1,1,1)
!!!      write (6, '(a, f8.4)') ' Dice station:', station

         do ip = 1, 5                 ! For dice patches 6, 2, 5, 4, 10

            ib = ib_dices_1(ip)       ! Patch to dice
            npts   = g(ib)%ni
            nlines = g(ib)%nj
            jinc   = 1
            if (ib == 2) jinc = -1    ! Left wall; j always starts at the floor

!!!         write (6, '(a, 7i5)') ' is,ii,ip,ib,npts,nlines,jinc: ', &
!!!                                 is,ii,ip,ib,npts,nlines,jinc

            call slice_patch (ii, station, npts, nlines, jinc,                 &
                              g(ib)%x, g(ib)%y, g(ib)%z,                       &
                              dices(is)%y, dices(is)%z)

            ii = ii + nlines - 1      ! Avoid duplicating common edge points

         end do ! Next patch to dice

!!!      write (6, '(3f12.7)') (dices(is)%x(i,1,1), dices(is)%y(i,1,1), &
!!!                             dices(is)%z(i,1,1), i = 1, npts_per_dice_1)

!        Replace y with arc length, except we can't overwrite yet.

         call chords2d (npts_per_dice_1, dices(is)%y, dices(is)%z, false,      &
                        total, dice_1_arcs(1,is))

      end do ! Next floor dice

!     Likewise for the rest of the cross-cuts aft of the cavity:

      do is = 1, ndices2

         ii = 1                       ! Dice point counter, left to right
         it = is + ndices1
         station = dices(it)%x(1,1,1)
!!!      write (6, '(a, f8.4)') ' Dice station:', station

         do ip = 1, 3                 ! For dice patches 7, 8, 9

            ib = ib_dices_2(ip)       ! Patch to dice
            npts   = g(ib)%ni
            nlines = g(ib)%nj
            jinc   = 1

!!!         write (6, '(a, 7i5)') ' it,ii,ip,ib,npts,nlines,jinc: ', &
!!!                                 it,ii,ip,ib,npts,nlines,jinc

            call slice_patch (ii, station, npts, nlines, jinc,                 &
                              g(ib)%x, g(ib)%y, g(ib)%z,                       &
                              dices(it)%y, dices(it)%z)

            ii = ii + nlines - 1      ! Avoid the duplicate

         end do ! Next patch to dice

!!!      write (6, '(3f12.7)') (dices(it)%x(i,1,1), dices(it)%y(i,1,1), &
!!!                             dices(it)%z(i,1,1), i = 1, npts_per_dice_2)

!        Replace y with arc length, except we can't overwrite yet.

         call chords2d (npts_per_dice_2, dices(it)%y, dices(it)%z, false,      &
                        total, dice_2_arcs(1,is))

      end do ! Next downstream dice

      end subroutine set_up_cuts

!     --------------------------------------------------------------------------

      subroutine slice_patch (i1, station, npts, nlines, jinc,                 &
                              abscissas, ordinates_1, ordinates_2,             &
                              slice_coords_1, slice_coords_2)

!     Slice all j lines of a patch at the indicated station, which could be
!     an x station or a y station.  Interpolate the missing x, y, or z coords.
!     More general methods exist, but we take advantage of the known cavity
!     grid topology at the higher level to do it this more straightforward way.

!     --------------------------------------------------------------------------

!     Arguments:

      integer, intent (in)  :: i1       ! Index in slice_coords_1/2 to start at
      real,    intent (in)  :: station  ! x or y station of a slice or dice
      integer, intent (in)  :: npts     ! # points i on each j line
      integer, intent (in)  :: nlines   ! # patch j lines to interpolate within
      integer, intent (in)  :: jinc     ! +/-1 to keep slice points in order
      real,    intent (in)  :: abscissas  (npts,nlines) ! Patch x or y coords.
      real,    intent (in)  :: ordinates_1(npts,nlines) ! Patch y &  z coords.
      real,    intent (in)  :: ordinates_2(npts,nlines) ! .. or z &  x coords.
      real,    intent (out) :: slice_coords_1(*)        ! Interpolated y & z or
      real,    intent (out) :: slice_coords_2(*)        ! interpolated z & x.

!     Local constants:

!     real, parameter :: one = 1.

!     Local variables:

      integer :: ic, ii, j1, j2
      real    :: r, rm1

!     Execution:

      ic = 1                       ! Interval pointer (in & out)
      ii = i1

      if (jinc == 1) then
         j1 = 1;  j2 = nlines
      else
         j1 = nlines;  j2 = 1
      end if

      do j = j1, j2, jinc

         call interval (npts, abscissas(1,j), station, one, ic)

         r   = (station           - abscissas(ic,j)) / &
               (abscissas(ic+1,j) - abscissas(ic,j))
         rm1 = one - r
         slice_coords_1(ii) = rm1 * ordinates_1(ic,j) + r * ordinates_1(ic+1,j)
         slice_coords_2(ii) = rm1 * ordinates_2(ic,j) + r * ordinates_2(ic+1,j)

!!!      write (6, '(a, 2i4, 3f12.7)') ' j,ic,xi,y,z:', j,ic, &
!!!         station, slice_coords_1(ii), slice_coords_2(ii)

         ii  = ii + 1

      end do

      end subroutine slice_patch

!     --------------------------------------------------------------------------

      subroutine inner_envelope_grid ()

!     Scale and shift the densified outer envelope grid to sit on the missing
!     OML, about half size.  Only the cavity walls and floor are transformed.
!     Half scale should suffice.

!     --------------------------------------------------------------------------

!     Local constants:

      real, parameter :: xyz_scale = 0.5,      & ! Applied to floor length, etc.
                         fourth = 0.25, half = 0.5
!     Local variables:

      integer :: ib, ier, ni, nj
      real    :: d, dx, dy, dz
      real    :: xid, xiu, xod, xou, yil, yir, yol, yor, zbi, zbo, zti, zto

!     Execution:

!     Make a copy of the densified outer envelope array details (patches 1-5):

      call clone_grid (nb_inner, numf_cfd, xyzq_out, xyzq_inner)

!     Calculate the shift that goes with scaling x.  Work with floor, not OML.

      ni  = xyzq_out(1)%ni
      xou = half * (xyzq_out(1)%x(1,1,1) + xyzq_out(1)%x(ni,1,1)) ! Outer upstr.
      xod = half * (xyzq_out(3)%x(1,1,1) + xyzq_out(3)%x(ni,1,1)) ! Outer dnstr.
      dx  = xod  - xou                                            ! Outer floor
      d   = half * (dx - xyz_scale * dx)
      xiu = xou  + d
      xid = xod  - d
      dx  = (xod * xiu - xou * xid) / dx
      LENGTH = xyzq_out(7)%x(1,1,1) - xyzq_out(6)%x(1,1,1)        ! At OML; now
                                                                  ! should be 1
!     Calculate the shift that goes with scaling y:

      nj  = xyzq_out(6)%nj
      yol = xyzq_out(6)%y(1,nj,1)                ! OML edges of cavity envelope
      yor = xyzq_out(9)%y(1, 1,1)
      dy  = yor  - yol;            WIDTH = dy    ! Cavity width at OML
      d   = half * (dy - xyz_scale * dy)
      yil = yol  + d
      yir = yor  - d
      dy  = (yor * yil - yol * yir) / dy

!     Calculate the shift that goes with scaling z:

      ni  = xyzq_out(5)%ni
      nj  = xyzq_out(5)%nj
      zbo = fourth * (xyzq_out(5)%z(1, 1,1) + xyzq_out(5)%z(ni, 1,1)  +        &
                      xyzq_out(5)%z(1,nj,1) + xyzq_out(5)%z(ni,nj,1)) ! Mid flr.
      zto = xyzq_out(6)%z(1,1,1)                 ! OML
      dz  = zbo  - zto;            DEPTH  = dz   ! Cavity depth
      zbi = zto
      zti = zbi  - xyz_scale * dz
      dz  = (zbo * zti - zto * zbi) / dz

!     Scale and shift the outer envelope to give the inner envelope:

      do ib = 1, nb_inner

         xyzq_inner(ib)%nzoneaux = 0

         call Tec_block_allocate (xyzq_inner(ib), ndim, 0, ier)

         if (ier /= 0) then
            write (luncrt, '(a, i3)') &
               ' Trouble allocating inner envelope patch #', ib
            stop
         end if

         call scale_and_shift (xyzq_out(ib), xyz_scale, dx, dy, dz,            &
                               xyzq_inner(ib))
      end do

      if (diagnostics) then

         open (lun_out, file='inner_envelope.p3da', status='unknown')

         write (lun_out, '(i1)') nb_inner
         write (lun_out, '(3i4)') &
            (xyzq_inner(i)%ni, xyzq_inner(i)%nj, 1, i = 1, nb_inner)
         do i = 1, nb_inner
            write (lun_out, '(1p, 6e19.11)') &
               xyzq_inner(i)%x, xyzq_inner(i)%y, xyzq_inner(i)%z
         end do

         close (lun_out)

      end if

      end subroutine inner_envelope_grid

!     --------------------------------------------------------------------------

      subroutine align_cfd_grid ()

!     Remove any rotations contained in the CFD surface grid, and shift the grid
!     to match the outer envelope as best we can.
!
!     Since the principal axis method becomes less reliable the closer any of
!     the cavity dimensions come to each other, and the axis order is lost when
!     the eigensystem calculation orders the eigenvalues (= moment of inertia)
!     into ascending order, this version tries all of the possible rotations
!     that can be derived from the tentative one by interchanging columns and
!     swapping signs.  The fact that the OML patch at the top of the cavity is
!     missing also contributes to the fuzziness when the dimensions are not very
!     distinct.  Original use of L, W, D at the OML to identify which principal
!     axis is which also ignores possible wall slopes that affect the moment of
!     inertia calculations.
!
!     The number of possible rotations is 3x2x1 column permutations times 2x2x2
!     sign permutations, giving 48 possibilites.
!
!     In order to avoid more than the two original input reference points, a
!     third reference point is derived from those.  The point at 0.6L from the
!     upstream ref. pt. is projected to the floor of the cavity of both the CFD
!     and the CHT grids.  (Actually, a wall point may be closer, but the idea
!     should still work, and offsetting the in-between point slightly should
!     break ties that bilateral symmetry can produced.)  The rotation that most
!     closely overlays the derived reference points and also the centers of mass
!     (including shifts and a 2-D rotation originally implemented to achieve
!     good alignment with the correct 3-D rotation) is chosen and applied to the
!     full CFD surface grid as originally.  This should be bulletproof if the
!     CFD cavity really matches the CHT cavity.  If the CFD cavity is more like
!     a gouge (as allowed for), it is probably still the best we can do.
!
!     Belated realization:  the area-based moments of inertia are affected by
!     misalignment with the axes.  Therefore, take out as much misalignment as
!     the two reference points allow before calculating the moments for the CFD
!     data.

!     --------------------------------------------------------------------------

      use trigd

      real, parameter :: half  = 0.5, one = 1.0, zero = 0.0, large = 1.e+10,   &
                         ratio = 0.6  ! For in-between reference pt.

      integer :: i, ib, icol_best, icol_perm, isign_best, isign_perm, j
      integer :: list_out(5), ni, nj, npts
      real    :: bias, bias_angle, dsq, dsq_min, dx, dy, dz
      real    :: phi, pratio, qratio, theta
      real    :: xmax_cav, xmin_cav, ymax_cav, ymin_cav, zmax_cav, zmin_cav
      real    :: CM_CFD(3), CM_CHT(3), CM_MAP(3)
      real    :: lambda(3), p(3), q(3), R(3,3), shift(3)
      real    :: S(3,3), T(3,3), R_transpose(3,3), R_transpose_best(3,3)
      real    :: ref_CHT_pts(3,3), derived_ref_CHT_pt(3), derived_cfd_best(3)

      type (grid_type) :: ref_cfd_pts(1), ref_pts(1)

      data list_out /1, 2, 3, 4, 5/ ! Patches 1-5 are the cavity walls & floor

!     Execution:
!     ----------

!     Derive a third reference point from the original two, for the CFD grid
!     and the idealized grid.  Avoid one that could be transformed correctly in
!     more than one way (i.e., avoid the mid-point).
!     First, for the idealized (CHT) grid:

      ref_CHT_pts(1,1)  =  xyzq_ideal(1)%xmin             ! Upstream lip, middle
      ref_CHT_pts(3,1)  =  xyzq_ideal(1)%zmin
      ref_CHT_pts(2,1)  = (xyzq_ideal(1)%ymin + xyzq_ideal(1)%ymax) * half

      ref_CHT_pts(1,2)  =  xyzq_ideal(3)%xmax           ! Downstream lip, middle
      ref_CHT_pts(3,2)  =  xyzq_ideal(3)%zmin
      ref_CHT_pts(2,2)  = (xyzq_ideal(3)%ymin + xyzq_ideal(3)%ymax) * half

      bias              = (xyzq_ideal(3)%ymax - xyzq_ideal(3)%ymin) * 0.01

      ref_CHT_pts(:,3)  = (one - ratio) * ref_CHT_pts(:,1) + &     ! Derived pt.
                                 ratio  * ref_CHT_pts(:,2)         ! on the OML
      ref_CHT_pts(2,3)  = ref_CHT_pts(2,3) + bias          ! Break possible ties

!     Build a search tree from the walls and floor of the idealized grid:

      nquad = 0
      do ib = 1, 5
         nquad = (xyzq_ideal(ib)%ni - 1) * (xyzq_ideal(ib)%nj - 1) + nquad
      end do

      allocate (conn(3,nquad)) ! For patch # and (i,j)

      call build_adt (5, xyzq_ideal, nquad, conn)

!     Project the derived ref. CHT pt. to the cavity floor (or possibly a wall):

      call search_adt (ref_CHT_pts(:,3), iquad, pratio, qratio, dsqmin, true,  &
                       5, xyzq_ideal, nquad, conn, derived_ref_CHT_pt)

      deallocate (conn)

      if (diagnostics) then
         write (luncrt, '(/, 2a, /)') &
            ' Reference ideal cavity points & derived ref. pt.',               &
            ' on floor or wall (/L):'
         write (luncrt, '(1p, 3e15.6)') &
            ((ref_CHT_pts(i,j), i = 1, 3), j = 1, 3), derived_ref_cht_pt(:)
      end if

!     Before repeating for the CFD reference points, improve the scaling by
!     shifting the coordinates, then align the CFD reference points with the
!     axes as best we can before doing the moments of inertia calculations.

      xmin_cav =  large;  ymin_cav =  large;  zmin_cav =  large
      xmax_cav = -large;  ymax_cav = -large;  zmax_cav = -large

      do i = 1, num_cavity_patches
         ib = list_of_cavity_patches(i)

         call patch_data_range (xyzq_cfd(ib))

         xmin_cav = min (xmin_cav, xyzq_cfd(ib)%xmin)
         xmax_cav = max (xmax_cav, xyzq_cfd(ib)%xmax)
         ymin_cav = min (ymin_cav, xyzq_cfd(ib)%ymin)
         ymax_cav = max (ymax_cav, xyzq_cfd(ib)%ymax)
         zmin_cav = min (zmin_cav, xyzq_cfd(ib)%zmin)
         zmax_cav = max (zmax_cav, xyzq_cfd(ib)%zmax)
      end do

      shift(1) = -(xmin_cav + xmax_cav) * half
      shift(2) = -(ymin_cav + ymax_cav) * half
      shift(3) = -(zmin_cav + zmax_cav) * half

      if (diagnostics) then
         write (luncrt, '(/, a, 1p, 3e15.6)') &
            ' CFD cavity (x,y,z) shifts to improve scaling:', shift
      end if

      do ib = 1, nb_cfd

        call scale_and_shift (xyzq_cfd(ib), one, shift(1), shift(2), shift(3), &
                              xyzq_cfd(ib))
      end do

!     Align the CFD reference points with the X axis:

      p(1) = xyzq_cfd(ib_ref_u)%x(i_ref_u,j_ref_u,1)
      p(2) = xyzq_cfd(ib_ref_u)%y(i_ref_u,j_ref_u,1)
      p(3) = xyzq_cfd(ib_ref_u)%z(i_ref_u,j_ref_u,1)

      dx = xyzq_cfd(ib_ref_d)%x(i_ref_d,j_ref_d,1) - p(1)
      dy = xyzq_cfd(ib_ref_d)%y(i_ref_d,j_ref_d,1) - p(2)
      theta = -atand (dy / dx)

      do ib = 1, nb_cfd

        npts = xyzq_cfd(ib)%ni * xyzq_cfd(ib)%nj

        call rotate2d (npts, xyzq_cfd(ib)%x, xyzq_cfd(ib)%y, theta, p(1), p(2))

      end do

      dx = xyzq_cfd(ib_ref_d)%x(i_ref_d,j_ref_d,1) - p(1)
      dz = xyzq_cfd(ib_ref_d)%z(i_ref_d,j_ref_d,1) - p(3)
      phi = -atand (dz / dx)

      do ib = 1, nb_cfd

        npts = xyzq_cfd(ib)%ni * xyzq_cfd(ib)%nj

        call rotate2d (npts, xyzq_cfd(ib)%x, xyzq_cfd(ib)%z, phi, p(1), p(3))

      end do

      if (diagnostics) then

         write (luncrt, '(/, a, 1p, 2e15.6)') &
            ' Rotations that align CFD ref. points with OX:', theta, phi
         write (luncrt, '(/, 2a, /)') &
            ' OX-aligned reference upstream & downstream CFD cavity pts. (/L):'
         write (luncrt, '(1p, 3e15.6)') &
            xyzq_cfd(ib_ref_u)%x(i_ref_u,j_ref_u,1), &
            xyzq_cfd(ib_ref_u)%y(i_ref_u,j_ref_u,1), &
            xyzq_cfd(ib_ref_u)%z(i_ref_u,j_ref_u,1), &
            xyzq_cfd(ib_ref_d)%x(i_ref_d,j_ref_d,1), &
            xyzq_cfd(ib_ref_d)%y(i_ref_d,j_ref_d,1), &
            xyzq_cfd(ib_ref_d)%z(i_ref_d,j_ref_d,1)

         open (lun_out, file='aligned_with_OX.g', status='unknown')

         write (lun_out, '(i3)') nb_cfd
         write (lun_out, '(3i4)') &
            (xyzq_cfd(i)%ni, xyzq_cfd(i)%nj, 1, i = 1, nb_cfd)

         do i = 1, nb_cfd
            write (lun_out, '(1p, 6e19.11)') &
               xyzq_cfd(i)%x, xyzq_cfd(i)%y, xyzq_cfd(i)%z
         end do

         close (lun_out)

      end if

!     Estimate a third rotation about the realigned reference pt. line via a
!     1-D optimization that minimizes |(zmax - zmin) - CHT depth|.

      call best_rotation_about_reference_axis (1)

      if (diagnostics) then

         open (lun_out, file='third_alignment.g', status='unknown')

         write (lun_out, '(i3)') nb_cfd
         write (lun_out, '(3i4)') &
            (xyzq_cfd(i)%ni, xyzq_cfd(i)%nj, 1, i = 1, nb_cfd)

         do i = 1, nb_cfd
            write (lun_out, '(1p, 6e19.11)') &
               xyzq_cfd(i)%x, xyzq_cfd(i)%y, xyzq_cfd(i)%z
         end do

         close (lun_out)

      end if

!     The CFD data should have all rotations removed now, but is still only
!     roughly shifted to improve the data scaling.

!     Derive a third ref. pt. for the largely-realigned CFD surface patches:
!     ----------------------------------------------------------------------

      ref_cfd_pts(1)%ni = 3;    allocate (ref_cfd_pts(1)%x(3,1,1))
      ref_cfd_pts(1)%nj = 1;    allocate (ref_cfd_pts(1)%y(3,1,1))
      ref_cfd_pts(1)%nk = 1;    allocate (ref_cfd_pts(1)%z(3,1,1))

      ref_cfd_pts(1)%x(1,1,1) = xyzq_cfd(ib_ref_u)%x(i_ref_u,j_ref_u,1)
      ref_cfd_pts(1)%y(1,1,1) = xyzq_cfd(ib_ref_u)%y(i_ref_u,j_ref_u,1)
      ref_cfd_pts(1)%z(1,1,1) = xyzq_cfd(ib_ref_u)%z(i_ref_u,j_ref_u,1)

      ref_cfd_pts(1)%x(2,1,1) = xyzq_cfd(ib_ref_d)%x(i_ref_d,j_ref_d,1)
      ref_cfd_pts(1)%y(2,1,1) = xyzq_cfd(ib_ref_d)%y(i_ref_d,j_ref_d,1)
      ref_cfd_pts(1)%z(2,1,1) = xyzq_cfd(ib_ref_d)%z(i_ref_d,j_ref_d,1)

      dy = ref_cfd_pts(1)%y(2,1,1) - ref_cfd_pts(1)%y(1,1,1)
      dx = ref_cfd_pts(1)%x(2,1,1) - ref_cfd_pts(1)%x(1,1,1)

!!!   if (diagnostics) &
!!!      write (6, '(/, a, 1p, 2e20.10)') ' Numerator & denominator:', dy, dx

      if (dy == zero ) then  ! Arc tan misbehavior observed
         bias_angle = zero
      else
         bias_angle = atan (dy / dx)
      end if

!     This angle should be zero now, unless the CFD pts. weren't truly aligned.

      dy =  bias * cos (bias_angle)
      dx = -bias * sin (bias_angle)

      if (diagnostics) then
         write (luncrt, '(/, a, 1p, 3e15.6)') &
            ' bias radians, dx, dy for CFD ref. points:', bias_angle, dx, dy
      end if

      p(1) = (one - ratio) * ref_cfd_pts(1)%x(1,1,1) + &  ! Derived pt. on OML
                    ratio  * ref_cfd_pts(1)%x(2,1,1) + dx ! Break ties
      p(2) = (one - ratio) * ref_cfd_pts(1)%y(1,1,1) + &
                    ratio  * ref_cfd_pts(1)%y(2,1,1) + dy
      p(3) = (one - ratio) * ref_cfd_pts(1)%z(1,1,1) + &
                    ratio  * ref_cfd_pts(1)%z(2,1,1)

!     We can't assume the floor and wall patches are contiguous, so use all.
!     The tree-building routine deallocates first if necessary.

      nquad = 0
      do ib = 1, nb_cfd
         nquad = (xyzq_cfd(ib)%ni - 1) * (xyzq_cfd(ib)%nj - 1) + nquad
      end do

      allocate (conn(3,nquad)) ! For patch # and (i,j)

      call build_adt (nb_cfd, xyzq_cfd, nquad, conn)

!     Project the derived ref. CFD pt. to the cavity floor (or maybe wall):

      call search_adt (p, iquad, pratio, qratio, dsqmin, true,  &
                       nb_cfd, xyzq_cfd, nquad, conn, q)

      ref_cfd_pts(1)%x(3,1,1) = q(1)  ! Derived reference pt. before realignment
      ref_cfd_pts(1)%y(3,1,1) = q(2)
      ref_cfd_pts(1)%z(3,1,1) = q(3)

      deallocate (conn)

      if (diagnostics) then
         write (luncrt, '(/, 2a, /)') &
            ' Partially-realigned reference upstream, downstream & derived',   &
            ' CFD cavity points (/L):'
         write (luncrt, '(1p, 3e15.6)') &
            (ref_cfd_pts(1)%x(i,1,1), ref_cfd_pts(1)%y(i,1,1), &
             ref_cfd_pts(1)%z(i,1,1), i = 1, 2), p(:), q(:)
      end if

!     Determine the center of mass of the DENSIFIED outer envelope cavity faces:
!     --------------------------------------------------------------------------

      call rotation_matrix (nb_out, xyzq_out, 5, list_out, CM_CHT, R, lambda)

      if (diagnostics) then
         write (luncrt, '(/, 2a, /, (1p, e15.6, 3x, 3e15.6, e18.6))')          &
            ' Idealized cavity patches CM, principal axis matrix, and',        &
            ' ordered eigenvalues:', (CM_CHT(i), R(i,:), lambda(i), i = 1, 3)
      end if

!     This routine also determines principal axes and hence a rotation matrix.
!     However, eigenvectors give relative directions only.  Also, the walls may
!     slope, and the top of the cavity at the OML is not in the surface data.
!     Thus, even if no two of LENGTH, WIDTH & DEPTH are equal, it is not always
!     feasible to make sure the ascending eigenvalues (which correspond to
!     moments of inertia) match the cavity dimensions in descending order.
!     Therefore, we calculate principal axes in some unknown order, and try all
!     possible permutations of the associated rotation matrix.

!     Determine a preliminary rotation matrix from the CFD cavity patch list
!     (walls and floor only), by calculating moments of inertia about the CM:

      call rotation_matrix (nb_cfd, xyzq_cfd, num_cavity_patches,              &
                            list_of_cavity_patches, CM_CFD, R, lambda)

      if (diagnostics) then
         write (luncrt, '(/, 2a, /, (1p, e15.6, 3x, 3e15.6, e18.6))')          &
            ' CFD cavity patches CM, principal axis matrix, and',              &
            ' ordered eigenvalues:', (CM_CFD(i), R(i,:), lambda(i), i = 1, 3)
!!!      write (luncrt, '(a)')
      end if

      shift = CM_CHT - CM_CFD  ! This is now the last shift needed; see go to ..

      do ib = 1, nb_cfd

        call scale_and_shift (xyzq_cfd(ib), one, shift(1), shift(2), shift(3), &
                              xyzq_cfd(ib))
      end do

!     We transform the reference CFD points the same way we will transform all
!     of them when the best rotation matrix is found.  We also track where the
!     transformed CM of the CFD points ends up.
!
!     First, shift the untransformed reference CFD points the same way shifting
!     all CFD points puts their cavity CM at the origin, because the rotations
!     below are about the CM.  This needs to be done only once.

      call scale_and_shift (ref_cfd_pts(1), one, shift(1), shift(2), shift(3), &
                            ref_cfd_pts(1))

      ref_pts(1)%ni = 3;    allocate (ref_pts(1)%x(3,1,1))
      ref_pts(1)%nj = 1;    allocate (ref_pts(1)%y(3,1,1))
      ref_pts(1)%nk = 1;    allocate (ref_pts(1)%z(3,1,1))

!     Since principal axis calculations can be thrown off by asymmetric cells
!     within symmetric cavities (e.g., 0.25 x 0.25 x 0.1" WLE chip modeled by
!     Peter Gnoffo), we abandon that approach in favor of the above simpler
!     realignment with the axes.  No more rotation matrix calculations ...

      if (nb_cfd > 0) go to 500   ! Kludge preserves the original for posterity.
                                  ! --------------------------------------------


!     Loop over all column permutations & sign permutations for each col. perm.:
!     --------------------------------------------------------------------------

      dsqmin = 1.e+6

      do icol_perm = 1, 6

         select case (icol_perm)

            case (1)  ! 1 2 3
              S = R
            case (2)  ! 1 3 2
              S(:,1) = R(:,1);  S(:,2) = R(:,3);  S(:,3) = R(:,2)
            case (3)  ! 2 3 1
              S(:,1) = R(:,2);  S(:,2) = R(:,3);  S(:,3) = R(:,1)
            case (4)  ! 2 1 3
              S(:,1) = R(:,2);  S(:,2) = R(:,1);  S(:,3) = R(:,3)
            case (5)  ! 3 1 2
              S(:,1) = R(:,3);  S(:,2) = R(:,1);  S(:,3) = R(:,2)
            case (6)  ! 3 2 1
              S(:,1) = R(:,3);  S(:,2) = R(:,2);  S(:,3) = R(:,1)

         end select

         do isign_perm = 1, 8

            T = S

            select case (isign_perm)

               case (1)  ! + + +
                  ! T = S
               case (2)  ! - + +
                  T(:,1) = -T(:,1)
               case (3)  ! + - +
                  T(:,2) = -T(:,2)
               case (4)  ! - - +
                  T(:,1:2) = -T(:,1:2)
               case (5)  ! + + -
                  T(:,3) = -T(:,3)
               case (6)  ! - + -
                  T(:,1) = -T(:,1);  T(:,3) = -T(:,3)
               case (7)  ! + - -
                  T(:,2:3) = -T(:,2:3)
               case (8)  ! - - -
                  T(:,:) = -T(:,:)

            end select

!           Apply to the reference points what will be applied to the
!           shifted CFD points when the best rotation is found:

            R_transpose = transpose (T)

            do i = 1, 3
               p(1) = ref_cfd_pts(1)%x(i,1,1)
               p(2) = ref_cfd_pts(1)%y(i,1,1)
               p(3) = ref_cfd_pts(1)%z(i,1,1)
               q    = matmul (R_transpose, p)
               ref_pts(1)%x(i,1,1) = q(1)   ! Because these are overwritten
               ref_pts(1)%y(i,1,1) = q(2)
               ref_pts(1)%z(i,1,1) = q(3)
            end do

!           Apply the shift that makes the CFD and dense/ideal grid CMs match:

            call scale_and_shift (ref_pts(1), one, CM_CHT(1), CM_CHT(2),       &
                                  CM_CHT(3), ref_pts(1))

!           Take out any misalignment of the length-wise principal axis with the
!           Z = 0 plane that can result from asymmetric cavity walls.  Rotate
!           about the CM of the densified idealized cavity, using the (current
!           form of the) reference CFD cavity lip pts.

            dx = ref_pts(1)%x(2,1,1) - ref_pts(1)%x(1,1,1)
            dz = ref_pts(1)%z(2,1,1) - ref_pts(1)%z(1,1,1)
            theta = -atand (dz / dx)

            call rotate2d (3, ref_pts(1)%x, ref_pts(1)%z, theta, CM_CHT(1),    &
                           CM_CHT(3))

            CM_MAP = CM_CHT

            call rotate2d (1, CM_MAP(1), CM_MAP(3), theta, CM_CHT(1), CM_CHT(3))

!           Finally, shift so the realigned cavity lips are at Z = 0:

            shift(3) = xyzq_out(6)%z(1,1,1) - ref_pts(1)%z(2,1,1)
            ref_pts(1)%z(:,1,1) = ref_pts(1)%z(:,1,1) + shift(3)
            CM_MAP(3) = CM_MAP(3) + shift(3)

!           How close are the transformed reference points and CM?

            dsq = zero

            do i = 1, 2
               q(1) = ref_pts(1)%x(i,1,1) - ref_CHT_pts(1,i)
               q(2) = ref_pts(1)%y(i,1,1) - ref_CHT_pts(2,i)
               q(3) = ref_pts(1)%z(i,1,1) - ref_CHT_pts(3,i)

               dsq  = dot_product (q, q) + dsq
            end do

            q(1) = ref_pts(1)%x(3,1,1) - derived_ref_CHT_pt(1)
            q(2) = ref_pts(1)%y(3,1,1) - derived_ref_CHT_pt(2)
            q(3) = ref_pts(1)%z(3,1,1) - derived_ref_CHT_pt(3)

            dsq  = dot_product (q, q) + dsq

            q(:) = CM_MAP(:) - CM_CHT(:)

            dsq  = dot_product (q, q) + dsq

            if (dsq < dsqmin) then
               dsqmin = dsq
               R_transpose_best = R_transpose
               derived_cfd_best(1) = ref_pts(1)%x(3,1,1)
               derived_cfd_best(2) = ref_pts(1)%y(3,1,1)
               derived_cfd_best(3) = ref_pts(1)%z(3,1,1)
               icol_best  = icol_perm
               isign_best = isign_perm
            end if

            if (diagnostics) then
               write (luncrt, '(a, 2i3, 1p, e16.8)') &
                  '   icolperm, isignperm, dsq: ', icol_perm, isign_perm, dsq
            end if

         end do  ! Next sign permutation for this column permutation

      end do  ! Next column permutation

!     Below is the original code following the attempted permutation of R
!     based on the relative sizes of L, W, D.

!     Align the CFD grid with the idealized cavity by applying the rotation
!     matrix, which really contains the direction cosines of the principal axes.
!     Actually, we multiply the CFD solution by the transpose of R.

      R_transpose = R_transpose_best

      write (luncrt, '(/, a, 2i3)') ' Best column & sign permutations: ',      &
         icol_best, isign_best

!     First, shift the CM to the origin, because the rotation is about the CM:

      shift = -CM_CFD

      do ib = 1, nb_cfd

        call scale_and_shift (xyzq_cfd(ib), one, shift(1), shift(2), shift(3), &
                              xyzq_cfd(ib))

         do j = 1, xyzq_cfd(ib)%nj
            do i = 1, xyzq_cfd(ib)%ni
               p(1) = xyzq_cfd(ib)%x(i,j,1)
               p(2) = xyzq_cfd(ib)%y(i,j,1)
               p(3) = xyzq_cfd(ib)%z(i,j,1)
               q    = matmul (R_transpose, p)
               xyzq_cfd(ib)%x(i,j,1) = q(1)
               xyzq_cfd(ib)%y(i,j,1) = q(2)
               xyzq_cfd(ib)%z(i,j,1) = q(3)
            end do
         end do
      end do

!     Shift the CFD grid so the CM of the cavity matches the outer envelope CM:

      shift = CM_CHT

      do ib = 1, nb_cfd

        call scale_and_shift (xyzq_cfd(ib), one, shift(1), shift(2), shift(3), &
                              xyzq_cfd(ib))
      end do

!     Take out any misalignment of the length-wise principal axis with the Z = 0
!     plane that can result from asymmetric cavity walls.  Rotate about the CM
!     of the densified idealized cavity, using the reference CFD cavity lip pts.

      dx = xyzq_cfd(ib_ref_d)%x(i_ref_d,j_ref_d,1) - &
           xyzq_cfd(ib_ref_u)%x(i_ref_u,j_ref_u,1)
      dz = xyzq_cfd(ib_ref_d)%z(i_ref_d,j_ref_d,1) - &
           xyzq_cfd(ib_ref_u)%z(i_ref_u,j_ref_u,1)
      theta = -atand (dz / dx)

!!!   write (6, '(/, a)') ' Suppressing 2-D rotation.'

      do ib = 1, nb_cfd

        npts = xyzq_cfd(ib)%ni * xyzq_cfd(ib)%nj

        call rotate2d (npts, xyzq_cfd(ib)%x, xyzq_cfd(ib)%z, theta, &
                       CM_CHT(1), CM_CHT(3))
      end do

!     (Almost) finally, shift so the realigned cavity lips are at Z = 0:

      shift(1) = zero
      shift(2) = zero
      shift(3) = xyzq_out(6)%z(1,1,1) - xyzq_cfd(ib_ref_d)%z(i_ref_d,j_ref_d,1)

!!!   write (6, '(/, a)') ' Suppressing dZ shift.'

      do ib = 1, nb_cfd

        call scale_and_shift (xyzq_cfd(ib), one, shift(1), shift(2), shift(3), &
                              xyzq_cfd(ib))
      end do

      if (diagnostics) then

         write (luncrt, '(/, a, 1p, 2e15.6)') &
            ' 2-D rotation and vertical shift to align OML:', theta, shift(3)
         write (luncrt, '(/, 2a, /)') &
            ' Penultimate realigned upstream, downstream & derived CFD cavity',&
            ' grid points:'
         write (luncrt, '(1p, 3e15.6)') &
            xyzq_cfd(ib_ref_u)%x(i_ref_u,j_ref_u,1), &
            xyzq_cfd(ib_ref_u)%y(i_ref_u,j_ref_u,1), &
            xyzq_cfd(ib_ref_u)%z(i_ref_u,j_ref_u,1), &
            xyzq_cfd(ib_ref_d)%x(i_ref_d,j_ref_d,1), &
            xyzq_cfd(ib_ref_d)%y(i_ref_d,j_ref_d,1), &
            xyzq_cfd(ib_ref_d)%z(i_ref_d,j_ref_d,1), &
            derived_cfd_best(:)

      end if

!     Some cases can still be off.  Therefore, find the rotation about the line
!     joining the CHT reference points that best aligns the CFD cavity patches
!     from where they've been transformed so far.  This is a one-variable
!     minimization problem.  Since we haven't actually mapped to the ideal (CHT)
!     grid yet, we just use the mins. and maxs. of the transformed cavity as a
!     measure of goodness.

      call best_rotation_about_reference_axis (2)


 500  continue  ! The original rotation matrix scheme is no longer used
                ! Nor is the derived reference point

      if (diagnostics) then

         write (luncrt, '(/, a, /)') &
            ' Final realigned upstream, downstream & derived CFD cavity points:'
         write (luncrt, '(1p, 3e15.6)') &
            xyzq_cfd(ib_ref_u)%x(i_ref_u,j_ref_u,1), &
            xyzq_cfd(ib_ref_u)%y(i_ref_u,j_ref_u,1), &
            xyzq_cfd(ib_ref_u)%z(i_ref_u,j_ref_u,1), &
            xyzq_cfd(ib_ref_d)%x(i_ref_d,j_ref_d,1), &
            xyzq_cfd(ib_ref_d)%y(i_ref_d,j_ref_d,1), &
            xyzq_cfd(ib_ref_d)%z(i_ref_d,j_ref_d,1), &
            ref_cfd_pts(1)%x(3,1,1), ref_cfd_pts(1)%y(3,1,1), &
            ref_cfd_pts(1)%z(3,1,1)
!!!         derived_cfd_best(:)

         open (lun_out, file='realigned_cfd_grid.p3da', status='unknown')

         write (lun_out, '(i3)') nb_cfd
         write (lun_out, '(3i4)') &
            (xyzq_cfd(i)%ni, xyzq_cfd(i)%nj, 1, i = 1, nb_cfd)

         do i = 1, nb_cfd
            write (lun_out, '(1p, 6e19.11)') &
               xyzq_cfd(i)%x, xyzq_cfd(i)%y, xyzq_cfd(i)%z
         end do

         close (lun_out)

      end if

      end subroutine align_cfd_grid

!     --------------------------------------------------------------------------

      subroutine best_rotation_about_reference_axis (mode)

      integer, intent (in) :: mode

!     Original usage (mode 2):
!     Determine via 1-variable optimization the angle of rotation about the line
!     joining the ideal cavity reference points that produces the best match
!     between the CFD cavity patches (as realigned so far) with the ideal cavity
!     patches, as measured by patch (x,y,z) data ranges.
!
!     Similar preliminary usage (mode 1):
!     Estimate a third 2-D rotation about the line joining the CFD reference
!     points (now parallel to OX) that "levels" the CFD cavity, by minimizing
!     the difference between the CFD cavity-patch Z data range and the CHT
!     cavity depth.

!     --------------------------------------------------------------------------

!     Local constants:

      integer,   parameter :: nfmax  = 60    ! Limit on # function evaluations
      real,      parameter :: half = 0.5, large  = 1.e+10
      character, parameter :: caller * 6 = 'BEST_R'

!     Local variables:

      integer :: i, ib, istat, lunerr, ni, nj, npts, numfun
      real    :: depth_over_L, dsq, r, ra, rb, tol, dlip(3), ulip(3)
      real    :: xmax_CHT, xmin_CHT, ymax_CHT, ymin_CHT, zmax_CHT, zmin_CHT
      real    :: xmax_rot, xmin_rot, ymax_rot, ymin_rot, zmax_rot, zmin_rot

      type (grid_type), allocatable, dimension (:) :: xyz_rotate

!     Execution:

      allocate (xyz_rotate(num_cavity_patches))

      do i = 1, num_cavity_patches
         ib = list_of_cavity_patches(i)
         ni = xyzq_cfd(ib)%ni
         nj = xyzq_cfd(ib)%nj
         xyz_rotate(i)%ni = ni
         xyz_rotate(i)%nj = nj
         xyz_rotate(i)%nk = 1

         allocate (xyz_rotate(i)%x(ni,nj,1), xyz_rotate(i)%y(ni,nj,1),         &
                   xyz_rotate(i)%z(ni,nj,1))
      end do

      if (mode == 1) then ! Estimate third preliminary rotation

         ulip(1) = xyzq_cfd(ib_ref_u)%x(i_ref_u,j_ref_u,1)
         ulip(2) = xyzq_cfd(ib_ref_u)%y(i_ref_u,j_ref_u,1)
         ulip(3) = xyzq_cfd(ib_ref_u)%z(i_ref_u,j_ref_u,1)

!!!      dlip(2) = xyzq_cfd(ib_ref_d)%y(i_ref_d,j_ref_d,1)
!!!      dlip(3) = xyzq_cfd(ib_ref_d)%z(i_ref_d,j_ref_d,1)

         depth_over_L = cavity_depth / cavity_length

      else ! mode == 2: Final correction about line joining CHT reference pts.

         ulip(1) =  xyzq_ideal(1)%xmin
         ulip(3) =  xyzq_ideal(1)%zmin
         ulip(2) = (xyzq_ideal(1)%ymin + xyzq_ideal(1)%ymax) * half

         dlip(1) =  xyzq_ideal(3)%xmax
         dlip(3) =  xyzq_ideal(3)%zmin
         dlip(2) = (xyzq_ideal(3)%ymin + xyzq_ideal(3)%ymax) * half

!        Idealized cavity data range:

         xmin_CHT = ulip(1)
         xmax_CHT = dlip(1)
         ymin_CHT = xyzq_ideal(2)%ymin
         ymax_CHT = xyzq_ideal(4)%ymax
         zmin_CHT = xyzq_ideal(1)%zmin
         zmax_CHT = xyzq_ideal(5)%zmax

         if (diagnostics) then
            write (6, '(/, a, 1p, 3e20.10)') ' ulip: ', ulip
            write (6, '(   a, 1p, 3e20.10)') ' dlip: ', dlip
            write (6, '(/, a, 1p, 3e20.10)') ' CHT min (x,y,z): ', &
               xmin_CHT, ymin_CHT, zmin_CHT
            write (6, '(   a, 1p, 3e20.10)') ' CHT max (x,y,z): ', &
               xmax_CHT, ymax_CHT, zmax_CHT
         end if

      end if

!     A minimum can be found to within sqrt (machine eps.), but avoid the sqrt:

      if (epsilon (tol) < 1.e-10) then
         tol = 1.e-8
      else
         tol = 1.e-4
      end if

      ra = -100.          ! Degrees of rotation for search interval
      rb =  100.          ! [not obvious what to choose]
      numfun = nfmax      ! Limit; starting guess may be poor
      lunerr = luncrt     ! < 0 suppresses iteration printout
      istat = 2           ! Initialize the minimization

10    continue

         call fminrc (ra, rb, r, dsq, tol, numfun, caller, lunerr, istat)

         if (istat < -1) then

            write (luncrt, '(/, 2a)') caller, ': FMINRC fatal error'
            stop

         else if (istat < 0) then ! Iteration limit; may be usable

            write (luncrt, '(/, 2a)') caller, ': Iteration limit.'

         else if (istat > 0) then ! Evaluate the objective function at r

!           Internal procedures can't use internal procedures, so do it in-line.

            xmin_rot =  large;  ymin_rot =  large;  zmin_rot =  large
            xmax_rot = -large;  ymax_rot = -large;  zmax_rot = -large

            do i = 1, num_cavity_patches
               ib = list_of_cavity_patches(i)
               xyz_rotate(i)%x = xyzq_cfd(ib)%x
               xyz_rotate(i)%y = xyzq_cfd(ib)%y
               xyz_rotate(i)%z = xyzq_cfd(ib)%z

               npts = xyz_rotate(i)%ni * xyz_rotate(i)%nj

               call rotate2d (npts, xyz_rotate(i)%y, xyz_rotate(i)%z, r,       &
                              ulip(2), ulip(3))

               call patch_data_range (xyz_rotate(i))

               xmin_rot = min (xmin_rot, xyz_rotate(i)%xmin)
               xmax_rot = max (xmax_rot, xyz_rotate(i)%xmax)
               ymin_rot = min (ymin_rot, xyz_rotate(i)%ymin)
               ymax_rot = max (ymax_rot, xyz_rotate(i)%ymax)
               zmin_rot = min (zmin_rot, xyz_rotate(i)%zmin)
               zmax_rot = max (zmax_rot, xyz_rotate(i)%zmax)
            end do

            if (mode == 1) then ! Preliminary 3rd rotation

               dsq = ((zmax_rot - zmin_rot) - depth_over_L)**2

!!!            write (6, '(10x, a, 1p, 3e15.6)') &
!!!               'zmax_rot, zmin_rot, dsq:', zmax_rot, zmin_rot, dsq

            else ! mode ==2 (final correction)

               dsq = (xmin_rot - xmin_CHT)**2 + (xmax_rot - xmax_CHT)**2 +     &
                     (ymin_rot - ymin_CHT)**2 + (ymax_rot - ymax_CHT)**2 +     &
                     (zmin_rot - zmin_CHT)**2 + (zmax_rot - zmax_CHT)**2
            end if

            go to 10

         else ! istat = 0 (success)

         end if

!     Apply the apparently best 2-D rotation to all the transformed CFD patches:

      do ib = 1, nb_cfd

         npts = xyzq_cfd(ib)%ni * xyzq_cfd(ib)%nj

         call rotate2d (npts, xyzq_cfd(ib)%y, xyzq_cfd(ib)%z, r,               &
                        ulip(2), ulip(3))
      end do

      do i = 1, num_cavity_patches
         deallocate (xyz_rotate(i)%x, xyz_rotate(i)%y, xyz_rotate(i)%z)
      end do

      deallocate (xyz_rotate)

      end subroutine best_rotation_about_reference_axis

!     --------------------------------------------------------------------------

      subroutine map_to_envelope ()

!     Map the realigned CFD solution to the plain outer envelope grid.
!     This has been added as an afterthought.  Initially, only the plain ADT
!     search method was implemented for each coarse Cavity Heating Tool grid pt.
!     This version breaks each coarse cell into many subcells, interpolates at
!     the center of each subcell, and averages these to give coarse cell center
!     flow values.  These in turn are converted to coarse cell vertex values.
!     Common edge values must then be averaged to force continuity.

!     --------------------------------------------------------------------------

!     Local constants:

      integer, parameter :: nsubi = 8, nsubj = 10   ! # subcell vertices
      real,    parameter :: half = 0.5

!     Local variables:

      integer :: i, i1, i2, ib, ic, ie, ier, ii, inc, ip, iquad, is,           &
                 j, j1, j2, jc, jj, jnc, js, l, l1, l2, lnc, m, m1, m2, mnc,   &
                 n, n1, n2, ni, nj, nsearches_for_this_patch

      logical :: lvarying

      real    :: p, pm1, q, qm1, interp_xyz(3), target_xyz(3)

      real, allocatable :: sum_fa(:),             & ! For accumulating fs * area
                           f_interp(:)              ! For subcell centroid check

      type (grid_type)  :: centers_and_areas,     & ! For one coarse grid patch
                           coarse_centers,        & ! Coarse cell center fs only
                           subgrid,               & ! For subcell vertex x/y/zs
                           subcenters_and_areas     ! For subcell cntrs & areas
!     Execution:

      subgrid%ni = nsubi;  subgrid%nj = nsubj;  subgrid%nk = 1
      subgrid%nzoneaux = 0

      call Tec_block_allocate (subgrid, ndim, 0, ier)

      if (ier /= 0) then
         write (luncrt, '(a)') &
            ' Trouble allocating subgrid for mapping to idealized cavity.'
         stop
      end if

      allocate (sum_fa(numf_cfd), f_interp(numf_cfd))

!     Save the interpolations to subcell centroids?

      if (diagnostics) allocate (interp_surf2(nb_out))

!     Surface patch utility "centroids_areas" does further allocations.

      write (luncrt, '(/, 2a, /)') ' Interpolating transformed CFD solution',  &
         ' at subcell centers of ideal coarse grid for area-based averaging.'

!     For each coarse patch of the outer envelope ...

      do ib = 1, nb_out

         xyzq_out2(ib)%zone_title = xyzq_out(ib)%zone_title
         xyzq_out2(ib)%nzoneaux   = 0

         call Tec_block_allocate (xyzq_out2(ib), ndim, numf_cfd, ier) ! Vertices

         if (ier /= 0) then
            write (luncrt, '(a, i3)') &
               ' Trouble allocating xyzq_out2 patch #', ib
            stop
         end if

         ni = xyzq_ideal(ib)%ni
         nj = xyzq_ideal(ib)%nj

         do j = 1, nj
            do i = 1, ni
               xyzq_out2(ib)%x(i,j,1) = xyzq_ideal(ib)%x(i,j,1)
               xyzq_out2(ib)%y(i,j,1) = xyzq_ideal(ib)%y(i,j,1)
               xyzq_out2(ib)%z(i,j,1) = xyzq_ideal(ib)%z(i,j,1)
            end do
         end do

         nsearches_for_this_patch = (ni - 1)*(nj - 1) * (nsubi - 1)*(nsubj - 1)

         if (diagnostics) then

            interp_surf2(ib)%ni = (ni - 1) * (nsubi - 1)
            interp_surf2(ib)%nj = (nj - 1) * (nsubj - 1)
            interp_surf2(ib)%nk = 1
            interp_surf2(ib)%zone_title   = xyzq_out(ib)%zone_title
            interp_surf2(ib)%nzoneaux     = 0
            interp_surf2(ib)%solutiontime = undefined

            call Tec_block_allocate (interp_surf2(ib), ndim, numf_cfd, ier)

            if (ier /= 0) then
               write (luncrt, '(a, i3)') &
                  ' Trouble allocating interp_surf2 patch #', ib
               stop
            end if

         end if

         ninside  = 0;  dmax  = zero
         noutside = 0;  dmean = zero

         call centroids_areas (xyzq_ideal(ib), centers_and_areas) ! Surface util

         allocate (coarse_centers%q(numf_cfd,ni-1,nj-1,1)) ! Already have x/y/zs

!        For each coarse cell ...

         do j = 1, nj - 1

            do i = 1, ni - 1

               call densify_cell (i, j, xyzq_ideal(ib), subgrid)    ! Local proc

               call centroids_areas (subgrid, subcenters_and_areas) ! Surf. util

               sum_fa(:) = zero

!              For each subcell center ...

               do js = 1, nsubj - 1
                  do is = 1, nsubi - 1
                     target_xyz(1) = subcenters_and_areas%x(is,js,1)
                     target_xyz(2) = subcenters_and_areas%y(is,js,1)
                     target_xyz(3) = subcenters_and_areas%z(is,js,1)

                     call search_adt (target_xyz, iquad, p, q, dsqmin, true,   &
                                      nb_cfd, xyzq_cfd, nquad, conn, interp_xyz)

                     if (dsqmin < dtolsq) then ! Nearest quad was within tol.
                        ninside  = ninside + 1
                     else
                        noutside = noutside + 1
                     end if

                     dmax  = max (dmax, dsqmin)
                     dmean = dmean + sqrt (dsqmin)

!                    Interpolate the aligned surface flow at this target point:

                     n   = conn(1,iquad) ! Block #
                     ic  = conn(2,iquad) ! Lower left quad. indices
                     jc  = conn(3,iquad)
                     pm1 = one - p
                     qm1 = one - q

                     f_interp(:) = qm1*(pm1*xyzq_cfd(n)%q(:,ic,  jc,  1)  +    &
                                          p*xyzq_cfd(n)%q(:,ic+1,jc,  1)) +    &
                                     q*(pm1*xyzq_cfd(n)%q(:,ic,  jc+1,1)  +    &
                                          p*xyzq_cfd(n)%q(:,ic+1,jc+1,1))
                     sum_fa(:) = f_interp(:)*subcenters_and_areas%q(1,is,js,1) &
                                 + sum_fa(:)

                     if (diagnostics) then
!!!                     if (ib == 3) then
!!!                        if (is == 3 .and. js == 3) then
!!!                           write (6,'(a, 3i5)') ' ib,i,j: ', ib, i, j
!!!                           write (6,'(a, 2i5)') ' is, js: ', is, js
!!!                           write (6,'(a, 3i5)') ' n, ic, jc: ', n, ic, jc
!!!                           write (6,'(4f16.8)') &
!!!                              xyzq_cfd(n)%q(1,ic:ic+1,jc:jc+1,1)
!!!                           write (6,'(a, (3f16.8))') ' p, q, finterp: ',    &
!!!                              p, q, f_interp(:)
!!!                        end if
!!!                     end if

                        ii = (i-1)*(nsubi-1) + is
                        jj = (j-1)*(nsubj-1) + js
                        interp_surf2(ib)%x(ii,jj,1) = interp_xyz(1)
                        interp_surf2(ib)%y(ii,jj,1) = interp_xyz(2)
                        interp_surf2(ib)%z(ii,jj,1) = interp_xyz(3)
                        interp_surf2(ib)%q(:,ii,jj,1) = f_interp(:)
                     end if

                  end do ! Next subcell is

               end do ! Next subcell js

               coarse_centers%q(:,i,j,1) = sum_fa(:) / &
                                           centers_and_areas%q(1,i,j,1)

            end do ! Next cell i for this target block

         end do ! Next cell j for this target block

         dmax  = sqrt (dmax)
         dmean = dmean / real (nsearches_for_this_patch)

         write (luncrt, '(a, i3, a, 2i7, a, 1p, 2e12.5)')                      &
            '   Outer envelope patch', ib,                                     &
            ':  # points in/outside tolerance:', ninside, noutside,            &
            ';  max/mean devn.:', dmax, dmean


!        Awkward determination of averages over portions of the floor:
!        -------------------------------------------------------------

         if (ib == ib_floor) then  ! Calculate but don't print the averages

            call lengthwise_floor_splits (ni, nj, numf_cfd, i_bump_var_in,     &
                                         coarse_centers%q, false,              &
                                         floor_averages_CFD, floor_averages_CHT)
         end if

         if (ib == nb_out) then  ! Tabulate the averages

            call lengthwise_floor_splits (ni, nj, numf_cfd, i_bump_var_in,     &
                                         coarse_centers%q, true,               &
                                         floor_averages_CHT, floor_averages_CFD)
         end if


!        Move the interpolated flow from coarse cell centers to cell vertices:

         call centers_to_vertices (ni, nj, numf_cfd, &
                                   coarse_centers%q, xyzq_out2(ib)%q)

         deallocate (centers_and_areas%x, centers_and_areas%y,                 &
                     centers_and_areas%z, centers_and_areas%q, coarse_centers%q)

      end do ! Next coarse outer envelope patch

      deallocate (subgrid%x, subgrid%y, subgrid%z, sum_fa)
      deallocate (subcenters_and_areas%x, subcenters_and_areas%y,              &
                  subcenters_and_areas%z, subcenters_and_areas%q)

      if (diagnostics) then

         header_surf2%filename    = 'interp_surf2.dat'
         header_surf2%formatted   = true
         header_surf2%ndim        = ndim
         header_surf2%numq        = numf_cfd
         header_surf2%nblocks     = nb_out
         header_surf2%datapacking = 0
         header_surf2%title     = 'Subcell interpolations of realigned CFD data'
         allocate (header_surf2%varname(ndim + numf_cfd))
         header_surf2%varname(:)  = header_cfd%varname(:)
         header_surf2%ndatasetaux = 0

         call Tecplot_write (lun_interp2, header_surf2, interp_surf2, ios)

         if (ios /= 0) then
            write (luncrt, '(/, a)') &
               ' Trouble writing subcell interpolations to interp_surf2.dat.'
            stop
         end if

         call deallocate_header (header_surf2)
         call deallocate_blocks (1, nb_out, ndim, numf_cfd, interp_surf2, ios)

      end if

!     Force continuity along common patch edges.  See save_common_edge_info.

      do ie = 1, nedges
         n1 = medges(1,1,ie);  n2 = medges(1,2,ie)  ! Patch numbers
         i1 = medges(2,1,ie);  l1 = medges(2,2,ie)
         i2 = medges(3,1,ie);  l2 = medges(3,2,ie)
         inc= i2 - i1;         lnc= l2 - l1
         j1 = medges(4,1,ie);  m1 = medges(4,2,ie)
         j2 = medges(5,1,ie);  m2 = medges(5,2,ie)
         jnc= j2 - j1;         mnc= m2 - m1

!        One of the indices takes only one value for each patch.
!        Check for programmer error:

         ier = 0
         if (inc == 0 .and. jnc == 0) ier = 1
         if (inc /= 0 .and. jnc /= 0) ier = 1
         if (lnc == 0 .and. mnc == 0) ier = 1
         if (lnc /= 0 .and. mnc /= 0) ier = 1

         if (ier /= 0) then
            write (luncrt, '(/, (a, 7i5))') &
               ' ?? n1, i1, i2, inc, j1, j2, jnc: ', &
                    n1, i1, i2, inc, j1, j2, jnc,    &
               '    n2, l1, l2, lnc, m1, m2, mnc: ', &
                    n2, l1, l2, lnc, m1, m2, mnc,    &
               '    Map_to_envelope: Bad edge info.  Edge number, ie: ', ie
            write (luncrt, '(/, a, /, (i3, 2(i6, 2x, 4i4)))') &
               ' medges(1:5,1:2,1:nedges):', &
               (ii, medges(1:5,1:2,ii), ii = 1, nedges)
            stop
         end if

         if (inc == 0) then
            inc = 1  ! Can't have a zero increment
         else
            inc = sign (1, inc)
         end if

         if (jnc == 0) then
            jnc = 1
         else
            jnc = sign (1, jnc)
         end if

         lvarying = lnc /= 0

         if (lnc == 0) then
            lnc = 1
         else
            lnc = sign (1, lnc)
         end if

         if (mnc == 0) then
            mnc = 1
         else
            mnc = sign (1, mnc)
         end if

         l = l1;  m = m1

         do j = j1, j2, jnc
            do i = i1, i2, inc
               f_interp(:) = (xyzq_out2(n1)%q(:,i,j,1)  + &
                              xyzq_out2(n2)%q(:,l,m,1)) * half
               xyzq_out2(n1)%q(:,i,j,1) = f_interp(:)
               xyzq_out2(n2)%q(:,l,m,1) = f_interp(:)

               if (lvarying) then
                  l = l + lnc
               else
                  m = m + mnc
               end if

            end do
         end do

      end do ! Next common edge

      deallocate (f_interp)

!     Save the CFD data mapped to the CHT grid:

      header_out2%filename  = 'mapped.' // trim (header_cfd%filename)
      header_out2%formatted = true

      call Tecplot_write (lun_out2, header_out2, xyzq_out2, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') &
            ' Trouble writing solution mapped to outer envelope grid.'
         stop
      end if

!     Don't deallocate xyzq_out2 yet - we need to slice it still.

      end subroutine map_to_envelope

!     --------------------------------------------------------------------------

      subroutine densify_cell (i, j, coarse, subgrid)

!     Subdivide the indicated coarse grid cell, using the numbers of new grid
!     points already set in the subgrid data structure.  We need cell centers
!     for the new subcells, but it is convenient to generate cell vertices first
!     as is done here via TFI.
!
!     --------------------------------------------------------------------------

!     Arguments:

      integer, intent (in)             :: i, j      ! Coarse cell to densify
      type (grid_type), intent (in)    :: coarse    ! One coarse surface patch
      type (grid_type), intent (inout) :: subgrid   ! Finer dimensions input;
                                                    ! finer vertices output

!     Local variables:

      integer :: ic, is, jc, js, nsubi, nsubj
      real    :: dx, dy, dz, fraction
      real, allocatable :: arcs (:)

!     Execution:

      nsubi = subgrid%ni;  nsubj = subgrid%nj

!     For each coarse cell edge in the i direction:

      fraction = 1./ real (nsubi - 1);  js = 1

      do jc = j, j + 1
         dx = (coarse%x(i+1,jc,1) - coarse%x(i,jc,1)) * fraction
         dy = (coarse%y(i+1,jc,1) - coarse%y(i,jc,1)) * fraction
         dz = (coarse%z(i+1,jc,1) - coarse%z(i,jc,1)) * fraction
         do is = 1, nsubi
            subgrid%x(is,js,1) = coarse%x(i,jc,1) + real (is - 1) * dx
            subgrid%y(is,js,1) = coarse%y(i,jc,1) + real (is - 1) * dy
            subgrid%z(is,js,1) = coarse%z(i,jc,1) + real (is - 1) * dz
         end do
         js = nsubj
      end do

!     For each coarse cell edge in the j direction:

      fraction = 1./ real (nsubj - 1);  is = 1

      do ic = i, i + 1
         dx = (coarse%x(ic,j+1,1) - coarse%x(ic,j,1)) * fraction
         dy = (coarse%y(ic,j+1,1) - coarse%y(ic,j,1)) * fraction
         dz = (coarse%z(ic,j+1,1) - coarse%z(ic,j,1)) * fraction
         do js = 2, nsubj - 1  ! No need to do the corners again
            subgrid%x(is,js,1) = coarse%x(ic,j,1) + real (js - 1) * dx
            subgrid%y(is,js,1) = coarse%y(ic,j,1) + real (js - 1) * dy
            subgrid%z(is,js,1) = coarse%z(ic,j,1) + real (js - 1) * dz
         end do
         is = nsubi
      end do

!     Fill the refined interior via transfinite interpolation:

      allocate (arcs(2*(nsubi + nsubj)))

      call tfiq3d (nsubi, 1, nsubi, 1, nsubj, subgrid%x, subgrid%y, subgrid%z, &
                   arcs)
      deallocate  (arcs)

      end subroutine densify_cell

!     --------------------------------------------------------------------------

      subroutine lengthwise_floor_splits (ni, nj, nf, ibf, bf, tabulate, av,   &
                                          av2)

!     Calculate (and later print) certain averages over cavity floor portions.
!
!     Who came up with the indicated arbitrary splits in the face of the 29 x 4
!     cells on the CHT floor grid ???!
!
!     The underlying CHT grid cells are assumed to be uniform, so their areas
!     aren't actually needed.
!
!     Arranging to do the same for the CHT data and printing them side-by-side
!     with the CFD averages spoiled the original modularization.

!     --------------------------------------------------------------------------

!     Arguments:

      integer, intent (in) :: ni, nj, nf, ibf   ! ibf = fn. # for bump factors
      real,    intent (in) :: bf(nf,ni-1,nj-1)  ! Averaged CFD flow, CHT centers
      logical, intent (in) :: tabulate          ! Allows for printing later
      real,    intent (inout) :: av(*)          ! Dimension should be nsubareas
      real,    intent (inout) :: av2(*)         ! Likewise; print if tabulating

!     Local constants:

      integer, parameter :: nsubareas = 5
      real,    parameter :: small_fraction = 0.0001

!     Local variables:

      integer       :: i, j, n
      integer, save :: mi, mj
      integer       :: i1(nsubareas), i2(nsubareas)  ! Indices nearest to splits
      real          :: asum, bf_plus_uncertainty, f1, f2, ri, rj, split
      real,    save :: s1(nsubareas), s2(nsubareas)  ! Fractional splits
      character     :: string * 4

      data    s1 / 0.00, 0.50, 0.90, 0.00, 0.85/
      data    s2 / 0.50, 0.90, 1.00, 0.85, 1.00/

!     Execution:

      if (.not. tabulate) then

         mi = ni - 1;  ri = real (mi)  ! Floor cell counts
         mj = nj - 1;  rj = real (mj)

         do n = 1, nsubareas

            split = ri * s1(n)
            i1(n) = int (split)  ! Truncates
            if (s1(n) == zero) then
               i1(n) = 1
               f1 = zero       ! Fraction of cell in front of i1(n)
            else
               f1 = split - real (i1(n))
               if (f1 < small_fraction) then
                  f1 = zero
                  i1(n) = i1(n) + 1
               else
                  f1 = one - f1
                  i1(n) = i1(n) + 2
               end if
            end if

!!!         write (6, '(/, a, 2i5)') '   *** n, i1(n): ', n, i1(n)
!!!         write (6, '(a, 2f10.5)') '   split, f1: ', split, f1

            split = ri * s2(n)
            i2(n) = int (split)
            if (s2(n) == one) then
               f2 = zero
            else
               f2 = split - real (i2(n))
               if (f2 < small_fraction) f2 = zero
            end if

            av(n) = zero;  asum = zero

            if (f1 > zero) then
               i = i1(n) - 1
               do j = 1, mj
                  av(n) = av(n) + f1 * bf(ibf,i,j)
               end do
!!!            write (6, '(/, a, 4f10.5, /)') &
!!!               '   bf(ibf,i,:): ', bf(ibf,i,:)
               asum = asum + f1
            end if

            do i = i1(n), i2(n)
               do j = 1, mj
                  av(n) = av(n) + bf(ibf,i,j)
               end do
!!!            write (6, '(a, 4f10.5)') &
!!!               '   bf(ibf,i,:): ', bf(ibf,i,:)
               asum  = asum + one
            end do

            if (f2 > zero) then
               i = i2(n) + 1
               do j = 1, mj
                   av(n) = av(n) + f2 * bf(ibf,i,j)
               end do
!!!            write (6, '(/, a, 4f10.5, /)') &
!!!               '   bf(ibf,i,:): ', bf(ibf,i,:)
               asum = asum + f2
            end if

!!!         write (6, '(a, 3i5)') '   n, i1(n), i2(n): ', n, i1(n), i2(n)
!!!         write (6, '(a, 5f12.5)') '   split, f1, f2, asum, sum a*f: ', &
!!!                                      split, f1, f2, asum, av(n)

            av(n) = av(n) / (asum * rj)

!!!         write (6, '(a, f10.5)') '   average: ', av(n)

         end do  ! Next sub-area

      else  ! Tabulate results later to keep printout orderly

         rj = cavity_length / cavity_depth

         write (luncrt, '(/, a, f6.2, a)') &
            ' Bump Factor Averages on CHT Grid Floor (Length / Depth =', rj, ')'

         do n = 1, nsubareas
            j = nint (100. * (s2(n) - s1(n)))

            if (s1(n) == zero) then
               write (luncrt, '(/, 15x, a, /)') 'CHT   CFD   CFD + Uncertainty'
               string = 'Fore'
            else if (s2(n) == one) then
               string = 'Last'
            else
               string = 'Next'
            end if

            bf_plus_uncertainty = av2(n) + uncertainty_floor

            write (luncrt, '(4x, a, i3, a, 2f6.2, f13.2)') &
               string, j, '%', av(n), av2(n), bf_plus_uncertainty
         end do

      end if

      end subroutine lengthwise_floor_splits

!     --------------------------------------------------------------------------

      subroutine centers_to_vertices (ni, nj, nf, centers, vertices)

!     Move cell center function values to cell vertices for a surface patch.

!     --------------------------------------------------------------------------

      integer, intent (in)  :: ni, nj                ! # vertices
      integer, intent (in)  :: nf                    ! # functions per cell
      real,    intent (in)  :: centers(nf,ni-1,nj-1) ! Cell-centered fun. values
      real,    intent (out) :: vertices(nf,ni,nj)    ! Vertex-centered funs.

      real, parameter :: fourth = 0.25

      integer :: i, ia, ib, j, ja, jb

      do j = 1, nj
         ja = max (1, j - 1);  jb = min (j, nj - 1)
         do i = 1, ni
            ia = max (1, i - 1);  ib = min (i, ni - 1)
            vertices(:,i,j) = (centers(:,ia,ja) + centers(:,ib,ja)  + &
                               centers(:,ia,jb) + centers(:,ib,jb)) * fourth
         end do
      end do

      end subroutine centers_to_vertices

!     --------------------------------------------------------------------------

      subroutine vertices_to_centers (ni, nj, nf, vertices, centers)

!     Move function values from vertices to cell centers for a surface patch.

!     --------------------------------------------------------------------------

      integer, intent (in)  :: ni, nj                ! # vertices
      integer, intent (in)  :: nf                    ! # functions per cell
      real,    intent (in)  :: vertices(nf,ni,nj)    ! Vertex-centered funs.
      real,    intent (out) :: centers(nf,ni-1,nj-1) ! Cell-centered fun. values

      real, parameter :: fourth = 0.25

      integer :: i, j

      do j = 1, nj - 1
         do i = 1, ni - 1
            centers(:,i,j) = (vertices(:,i,j)   + vertices(:,i+1,j)    + &
                              vertices(:,i,j+1) + vertices(:,i+1,j+1)) * fourth
         end do
      end do

      end subroutine vertices_to_centers

!     --------------------------------------------------------------------------

      subroutine map_cavity_patches ()

!     Specialized mapping of CFD surface data to an idealized cavity. The cavity
!     modeled by the CFD has an arbitrary geometry and grid topology,  while the
!     idealized cavity has 4 planar walls and a planar floor at arbitrary angles
!     with fixed patch order and standardized point distributions, from which an
!     inner envelope has been derived.  Straight lines joining the corresponding
!     points of the inner and outer envelopes are intersected with the CFD sur-
!     face, and the flow solution is interpolated at the intersection points.
!     The 2-point lines may be extrapolated within reason.

!     --------------------------------------------------------------------------

!     Local constants:

      integer, parameter :: nl = 2                    ! 2-point lines always
      real,    parameter :: tmin = -0.01, tmax = 1.01 ! Allow tiny extrapolation
      character, parameter :: lcs_method * 1 = 'L'    ! LCSFIT not needed though

!     Local variables:

      integer :: i, ib, ic, iquad, j, jc, m, n, ni, nj
      real    :: p, pm1, q, qm1, t, interp_xyz(3)

      real, dimension (2) :: tl, xl, yl, zl

!     Execution:

      tl(1) = tmin;  tl(nl) = tmax

      do ib = 1, 5 ! The outer and inner envelope cavities have 5 patches

         ni = xyzq_out(ib)%ni
         nj = xyzq_out(ib)%nj

         allocate (interp_surf(ib)%x(ni,nj,1), interp_surf(ib)%y(ni,nj,1),     &
                   interp_surf(ib)%z(ni,nj,1))

         if (diagnostics) then
            interp_surf(ib)%zone_title   = xyzq_out(ib)%zone_title
            interp_surf(ib)%nzoneaux     = 0
            interp_surf(ib)%solutiontime = undefined
            allocate (interp_surf(ib)%q(numf_cfd,ni,nj,1))
         end if

         ninside  = 0;  dmax  = zero
         noutside = 0;  dmean = zero

         do j = 1, nj

            do i = 1, ni

               xl(1) = xyzq_out(ib)%x(i,j,1);  xl(nl) = xyzq_inner(ib)%x(i,j,1)
               yl(1) = xyzq_out(ib)%y(i,j,1);  yl(nl) = xyzq_inner(ib)%y(i,j,1)
               zl(1) = xyzq_out(ib)%z(i,j,1);  zl(nl) = xyzq_inner(ib)%z(i,j,1)

               call intsec6 (nb_cfd, xyzq_cfd, nquad, conn, nl, xl, yl, zl, tl,&
                             lcs_method, iquad, p, q, t, interp_xyz, dsqmin)

               if (dsqmin < dtolsq) then ! The nearest quad was within tolerance
                  ninside  = ninside + 1
               else
                  noutside = noutside + 1
               end if

               dmax  = max (dmax, dsqmin)
               dmean = dmean + sqrt (dsqmin)

               interp_surf(ib)%x(i,j,1) = interp_xyz(1)
               interp_surf(ib)%y(i,j,1) = interp_xyz(2)
               interp_surf(ib)%z(i,j,1) = interp_xyz(3)

!              Interpolate the surface flow at this target point:

               n   = conn(1,iquad) ! Block #
               ic  = conn(2,iquad) ! Lower left quad. indices
               jc  = conn(3,iquad)
               pm1 = one - p
               qm1 = one - q

               xyzq_out(ib)%q(:,i,j,1) = &
                  qm1 * (pm1 * xyzq_cfd(n)%q(:,ic,  jc,  1)  + &
                           p * xyzq_cfd(n)%q(:,ic+1,jc,  1)) + &
                    q * (pm1 * xyzq_cfd(n)%q(:,ic,  jc+1,1)  + &
                           p * xyzq_cfd(n)%q(:,ic+1,jc+1,1))

            end do

         end do

         dmax  = sqrt (dmax)
         dmean = dmean / real (ni * nj)
         write (luncrt, '(a, i3, a, 2i7, a, 1p, 2e12.5)')                      &
            '   Standardized patch #', ib,                                     &
            ':  # points in/outside tolerance:', ninside, noutside,            &
            ';  max/mean devn.:', dmax, dmean

         if (diagnostics) then

            interp_surf(ib)%q = xyzq_out(ib)%q

            call Tec_block_write (lun_interp, header_surf, interp_surf(ib), ios)

            if (ios /= 0) then
               write (luncrt, '(/, a, i3)') &
                  ' Trouble writing interpolated surface block', ib
               stop
            end if

         end if

         deallocate (interp_surf(ib)%x, interp_surf(ib)%y, interp_surf(ib)%z)

         if (diagnostics) deallocate (interp_surf(ib)%q)

      end do ! Next cavity block

      end subroutine map_cavity_patches

!     --------------------------------------------------------------------------

      subroutine map_OML_patches ()

!     Interpolate the CFD surface solution to the OML patches of the idealized
!     cavity grid.  The method here is the same as that for the earlier Ensemble
!     program.  Would it be better to include "OML" patches with the inner
!     envelope employed for mapping cavity walls and floor?
!
!     (Later:)  Do the cavity patches the straightforward way if map_method = 2.

!     --------------------------------------------------------------------------

!     Local variables:

      integer :: i, ib, ib1, ic, iquad, j, jc, m, n, ni, nj
      real    :: p, pm1, q, qm1, interp_xyz(3), target_xyz(3)

!     Execution:

      if (map_method == 1) then ! OML patches only
         ib1 = 6
      else                      ! All patches
         ib1 = 1
      end if

      write (luncrt, '(/, 2a, /)') ' Interpolating transformed CFD solution',  &
         ' at vertices of densified idealized cavity grid.'

      do ib = ib1, nb_out

         ni = xyzq_out(ib)%ni
         nj = xyzq_out(ib)%nj

         allocate (interp_surf(ib)%x(ni,nj,1), interp_surf(ib)%y(ni,nj,1),  &
                   interp_surf(ib)%z(ni,nj,1))

         if (diagnostics) then
            interp_surf(ib)%zone_title   = xyzq_out(ib)%zone_title
            interp_surf(ib)%nzoneaux     = 0
            interp_surf(ib)%solutiontime = undefined
            allocate (interp_surf(ib)%q(numf_cfd,ni,nj,1))
         end if


         ninside  = 0;  dmax  = zero
         noutside = 0;  dmean = zero

         do j = 1, nj

            do i = 1, ni

               target_xyz(1) = xyzq_out(ib)%x(i,j,1)
               target_xyz(2) = xyzq_out(ib)%y(i,j,1)
               target_xyz(3) = xyzq_out(ib)%z(i,j,1)

               call search_adt (target_xyz, iquad, p, q, dsqmin, true,         &
                                nb_cfd, xyzq_cfd, nquad, conn, interp_xyz)

               if (dsqmin < dtolsq) then ! The nearest quad was within tolerance
                  ninside  = ninside + 1
               else
                  noutside = noutside + 1
               end if

               dmax  = max (dmax, dsqmin)
               dmean = dmean + sqrt (dsqmin)

               interp_surf(ib)%x(i,j,1) = interp_xyz(1)
               interp_surf(ib)%y(i,j,1) = interp_xyz(2)
               interp_surf(ib)%z(i,j,1) = interp_xyz(3)

!              Interpolate the surface flow at this target point:

               n   = conn(1,iquad) ! Block #
               ic  = conn(2,iquad) ! Lower left quad. indices
               jc  = conn(3,iquad)
               pm1 = one - p
               qm1 = one - q

               xyzq_out(ib)%q(:,i,j,1) = &
                  qm1 * (pm1 * xyzq_cfd(n)%q(:,ic,  jc,  1)  + &
                           p * xyzq_cfd(n)%q(:,ic+1,jc,  1)) + &
                    q * (pm1 * xyzq_cfd(n)%q(:,ic,  jc+1,1)  + &
                           p * xyzq_cfd(n)%q(:,ic+1,jc+1,1))

            end do ! Next i for this target block

         end do ! Next j for this target block

         dmax  = sqrt (dmax)
         dmean = dmean / real (ni * nj)
         write (luncrt, '(a, i3, a, 2i7, a, 1p, 2e12.5)')                      &
            '   Standardized patch #', ib,                                     &
            ':  # points in/outside tolerance:', ninside, noutside,            &
            ';  max/mean devn.:', dmax, dmean

         if (diagnostics) then

            interp_surf(ib)%q = xyzq_out(ib)%q

            call Tec_block_write (lun_interp, header_surf, interp_surf(ib), ios)

            if (ios /= 0) then
               write (luncrt, '(/, a, i3)') &
                  ' Trouble writing interpolated surface block', ib
               stop
            end if

         end if

         deallocate (interp_surf(ib)%x, interp_surf(ib)%y, interp_surf(ib)%z)

         if (diagnostics) deallocate (interp_surf(ib)%q)

      end do ! Next target block

      if (diagnostics) close (lun_interp)

      end subroutine map_OML_patches

!     --------------------------------------------------------------------------

      subroutine add_uncertainty ()

!     Make room in the xyzq_out2 data structure for an extra CFD variable,
!     BF + uncertainty, and insert it in each CHT-related block.
!     The input numf_cfd is incremented by 1 here prior to slicing.

!     --------------------------------------------------------------------------

!     Local variables:

      integer :: ib, ni, nj, numf_new
      real    :: bf
      real,      allocatable :: variables (:,:,:)
      character, allocatable :: names (:) * (namelength)

!     Execution:

!     The xyzq_out2 data structure uses header_cfd even though its grid is that
!     of the CHT (idealized cavity).

      numf_new = numf_cfd + 1

      allocate (names(ndim + numf_cfd))
      names(:) =  header_cfd%varname(:)
      deallocate (header_cfd%varname)
      allocate   (header_cfd%varname(ndim + numf_new))
      header_cfd%varname(1:ndim + numf_cfd) = names(:)
      header_cfd%varname(ndim + numf_new) = 'BF_plus_uncertainty'
      header_cfd%numq =  numf_new
      deallocate (names)

!     For each CHT grid patch ...

      do ib = 1, nb_out

         ni = xyzq_out2(ib)%ni
         nj = xyzq_out2(ib)%nj

         allocate (variables(numf_cfd,ni,nj))
         variables(:,:,:) = xyzq_out2(ib)%q(:,:,:,1)
         deallocate (xyzq_out2(ib)%q)
         allocate   (xyzq_out2(ib)%q(numf_new,ni,nj,1))

         do j = 1, nj
            do i = 1, ni
               xyzq_out2(ib)%q(1:numf_cfd,i,j,1) = variables(:,i,j)
            end do
         end do

         deallocate (variables)

!        Bump factor uncertainty parameters are defined in the main program.
!        uncertainty_floor is set in standard_cavity_grids.

         if (ib == 5) then  ! Floor patch

            xyzq_out2(ib)%q(numf_new,:,:,1) = uncertainty_floor + &
               xyzq_out2(ib)%q(i_bump_var_in,:,:,1)

         else  ! Wall patch or outside the cavity

            do j = 1, nj
               do i = 1, ni
                  bf = xyzq_out2(ib)%q(i_bump_var_in,i,j,1)
                  xyzq_out2(ib)%q(numf_new,i,j,1) = bf + &
                     abs (bf - one) * OML_factor
               end do
            end do

         end if

      end do  ! Next idealized grid patch

      numf_cfd = numf_new

      end subroutine add_uncertainty

!     --------------------------------------------------------------------------

      subroutine interpolate_at_cuts (nb, nf, xyzq, nslices, slices)

!     Interpolate the given function data (on the outer envelope/cavity heating
!     tool grid) at the given slices (streamwise or crosswise).  The search tree
!     has already been constructed for use on both the CFD data & any CHT data.

!     --------------------------------------------------------------------------

!     Arguments:

      integer, intent (in)             :: nb, &      ! # blocks in surface data
                                          nf         ! # funcs. in surface data
      type (grid_type), intent (in)    :: xyzq(nb)   ! surface data structure
      integer, intent (in)             :: nslices    ! # cuts being interpolated
      type (grid_type), intent (inout) :: slices(nb) ! Input with target slice
                                                     ! (x,y,z)s; output with
                                                     ! interpolated functions
!     Local variables:

      integer :: i, ib, ic, iquad, j, jc, n, ni
      real    :: p, pm1, q, qm1, interp_xyz(3), target_xyz(3)

!     Execution:

      ninside  = 0;  dmax  = zero
      noutside = 0;  dmean = zero

      do ib = 1, nslices

         ni = slices(ib)%ni  ! nj = nk = 1

         do i = 1, ni

            target_xyz(1) = slices(ib)%x(i,1,1)
            target_xyz(2) = slices(ib)%y(i,1,1)
            target_xyz(3) = slices(ib)%z(i,1,1)

            call search_adt (target_xyz, iquad, p, q, dsqmin, true,            &
                             nb, xyzq, nquad, conn, interp_xyz)

            if (dsqmin < dtolsq) then ! The nearest quad was within tolerance
               ninside  = ninside + 1
            else
               noutside = noutside + 1
            end if

            dmax  = max (dmax, dsqmin)
            dmean = dmean + sqrt (dsqmin)

            n   = conn(1,iquad) ! Block #
            ic  = conn(2,iquad) ! Lower left quad. indices
            jc  = conn(3,iquad)
            pm1 = one - p
            qm1 = one - q

!!!         write (6,'(a, 2i5)') ' ib,i: ', ib, i
!!!         write (6,'(a, 3i5)') '    n, ic, jc: ', n, ic, jc
!!!         write (6,'(3x,4f16.8)') xyzq(n)%q(1,ic:ic+1,jc:jc+1,1)

            slices(ib)%q(1:nf,i,1,1) = &
               qm1 * (pm1 * xyzq(n)%q(1:nf,ic,  jc,  1)  + &
                        p * xyzq(n)%q(1:nf,ic+1,jc,  1)) + &
                 q * (pm1 * xyzq(n)%q(1:nf,ic,  jc+1,1)  + &
                        p * xyzq(n)%q(1:nf,ic+1,jc+1,1))

!!!         write (6,'(a, (3f16.8))') ' p, q, finterp: ', p, q, &
!!!            slices(ib)%q(1:nf,i,1,1)

         end do ! Next i for this slice

      end do ! Next target block

      dmax  = sqrt (dmax)
      dmean = dmean / real (ni * nslices)
      write (luncrt, '(a, a, i6, i5, a, 1p, 2e12.5)')                          &
         '   Slice/dice search statistics:',                                   &
         '  # pts. in/outside tolerance:', ninside, noutside,                  &
         ';  max/mean devn.:', dmax, dmean

      end subroutine interpolate_at_cuts

!     --------------------------------------------------------------------------

      subroutine normalize_patches (nb, g)

!     Convert the given standardized surface patches to unit cavity form by
!     normalizing the arc lengths of each patch and setting the third coordinate
!     constant.

!     --------------------------------------------------------------------------

!     Arguments:

      integer,          intent (in)    :: nb      ! # patches
      type (grid_type), intent (inout) :: g(nb)   ! Outer envelope or densified
                                                  ! form, both standardized
!     Local constants:

      real, parameter :: half = 0.5, onept5 = 1.5

!     Local variables:

      integer :: i, ib, j, ni, nj

      real, allocatable, dimension (:,:,:) :: u, v

!     Execution:

      do ib = 1, nb

         ni = g(ib)%ni;  nj =  g(ib)%nj

         allocate (u(ni,nj,1), v(ni,nj,1))

         u(1,1,1) = zero  ! -999. would suppress normalization

         call param2d (ni, nj, 1, ni, 1, nj, g(ib)%x, g(ib)%y, g(ib)%z, u, v)

         select case (ib)

         case (1) ! Upstream wall

            g(ib)%x = zero
            g(ib)%y = u - half
            g(ib)%z = v

         case (2) ! Side wall, left

            g(ib)%x = u
            g(ib)%y = -half
            g(ib)%z = v

         case (3) ! Downstream wall

            g(ib)%x = one
            g(ib)%y = u - half
            g(ib)%z = v

         case (4) ! Side wall, right

            g(ib)%x = u
            g(ib)%y = half
            g(ib)%z = v

         case (5) ! Floor

            g(ib)%x = u
            g(ib)%y = v - half
            g(ib)%z = zero

         case (6) ! Alongside, left

            g(ib)%x = u
            g(ib)%y = v - onept5
            g(ib)%z = one

         case (7) ! Downstream, left

            g(ib)%x = u + one
            g(ib)%y = v - onept5
            g(ib)%z = one

         case (8) ! Downstream, mid

            g(ib)%x = u + one
            g(ib)%y = v - half
            g(ib)%z = one

         case (9) ! Downstream, right

            g(ib)%x = u + one
            g(ib)%y = v + half
            g(ib)%z = one

         case (10) ! Alongside, right

            g(ib)%x = u
            g(ib)%y = v + half
            g(ib)%z = one

         end select

         deallocate (u, v)

      end do ! Next patch

      end subroutine normalize_patches

   end program cavity_map
