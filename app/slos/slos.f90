!
   program slos
!
!  Description:
!
!     SLOS has been adapted from USLOS (unstructured [grid] lines of sight),
!     itself an adaptation of LINES_OF_SIGHT and HEMISPHERES_OF_SIGHT, which
!     are now superseded by this combination of the two purposes (body-normal
!     lines for tangent-slab radiation calculations, or hemispherical lines
!     for just one body point at a time) for the 3D structured grid case.
!     (See also LINES_OF_SIGHT_2D.)
!
!     The volume grid is expected to be right-handed, and for the hemisphere
!     case, the y >= 0 half is expected (with x pointing downstream).  An off-
!     center body point requires both halves for the hemisphere case.
!
!     The inner surface is expected to represent a mostly-convex surface.  Any
!     concavities may cause difficulties with hemisphere lines that encounter
!     the body.  Ancillary utility ADJUST_HEMISPHERE_LOS is the initial answer
!     for such cases.
!
!  Body Point Inputs (one entry per line; trailing comments ignored):
!
!     Body-normal (tangent-slab radiation calculation) case:
!
!        0.  0.  0.
!        1.628659741950E-01  5.741545696080E-01  0.  ! Nose/cone tangency
!        8.985316209400E-01  1.837811469800E+00  0.  ! Cone/shoulder tangency
!         :                   :                   :
!
!        If only one body point is detected, a prompt still allows for the
!        tangent-slab case.
!
!     Hemispherical lines (full angular integration) case:
!
!        1.0742  -1.81  1.083  [10.]  ! x, y, z[, cone angle]
!
!        Only one body point is handled for this case.  The optional semi-cone
!        angle allows for writing fewer than the full hemisphere's worth of
!        lines to avoid unnecessary radiation cone angles outside the viewing
!        angle of a radiometer instrument.  Note that NEQAIR_INTEGRATION can
!        still integrate over a narrower angle than the indicated data angle.
!        If the optional angle is omitted, it is taken to be 90 degrees.
!
!  Input Structured Volume Grid (PLOT3D Format):
!
!     Formatted or unformatted multiblock files are handled, and a single layer
!     of grid blocks is assumed (normally the case for hypersonic vehicles) so
!     that the outer surface of all blocks serves as the shock-aligned boundary
!     with which intersections of lines from the body point(s) are calculated.
!
!  Outputs:
!
!     While the precursors of SLOS discretize the lines of sight with the same
!     relative point distributions as found in the nearest off-body grid line,
!     this variant follows the USLOS approach of using DPLR-type controls nk,
!     d1, and ds2_fraction that are prompted for.
!
!     The lines of sight are written to a single multiblock PLOT3D-type file,
!     one line per block.  These should be checked visually before proceeding
!     with radiation calculations.  Ensuing procedures will interpolate flow
!     data onto the discretized lines of sight and write results compatible with
!     the radiation solver.  The file name is hard-coded as 'los.g' as part
!     of simplifying the original nomenclature.
!
!  Method:
!
!     >  A handful of prompts serve to drive the program -- no control file.
!     >  Read the indicated body point (x,y,z) file -- one or more points.
!     >  If a single body point is found, prompt in case hemispherical lines are
!        NOT intended.
!     >  If all body points are on the centerline (y = 0.), the input volume
!        grid is assumed to be the starboard (right) half.  That half (y > 0.)
!        also serves for off-center body-normal lines, but for a hemisphere at
!        an off-center body point, a full volume grid is expected, with the
!        outer surface closed and convex.
!     >  Build a search tree from the inner multiblock surface, and search
!        for each body point.
!     >  Generate and store the corresponding unit body-normal vector(s).
!     >  Build a new search tree from the outer boundary, as needed for line-
!        surface intersection calculations.
!
!     >  Body-normal line case:
!           >  For each body point B:
!              >  The body normal unit vector is already in hand.
!              >  A search interval along the body-normal line for minimizing
!                 the distance to the outer surface is provided by s = 0 at the
!                 low end and s = the diagonal of the bounding box for the
!                 outer surface mesh.
!              >  The intersection point is then conveniently calculated via
!                 a 1D minimization w.r.t. distance s from B along the body-
!                 normal line (INTSEC6).
!              >  Impose a CFD-type point distribution on the intersected
!                 interval BP, as needed for interpolated flow field data.
!                 DPLR-type use of inputs nk, ds1, ds2fraction will suffice.
!           >  Write the discretized line(s) of sight to the indicated PLOT3D-
!              type output file.
!
!     >  Hemisphere lines of sight case:
!           >  The vertices of a triangulated unit hemispherical quadrant define
!              the lines of sight between the body point and the outer grid
!              boundary.  The number of uniform points specified along an edge
!              determines the resolution.  Typically, 25 such points are used,
!              giving 650 lines of sight for a centerline body point, or 1300
!              lines if the body point is off-center.  Centerline points require
!              just half the volume data, and the eventual integration over
!              solid angle of radiances from the lines of sight can be doubled.
!           >  The initial unit quadrant is transformed so that line 1/vertex 1
!              corresponds to the body point normal. After The lines for the
!              first quadrant are generated, the quadrant is rotated 90 degrees
!              about the normal and the next set of lines is generated.  Then if
!              the body point is off-center, this is repeated twice more.
!           >  Each intersection involves a 2-point line from the body point
!              through a quadrant vertex, of length derived from the bounding
!              box of the outer shock surface.  This length is overkill for many
!              intersection search intervals, but that doesn't matter.  Each
!              hemisphere line is intersected and discretized as for the body
!              normal (tangent-slab) case described above.
!  History:
!
!     11/02/2018  D.A.Saunders  Started adapting USLOS for the structured case.
!                               The option for entering body point indices in
!                               lieu of x/y/z coordinates has been dispensed
!                               with.  Use GU to coarsen a surface grid if
!                               necessary.  Generation of shock-normal lines
!                               and lines parallel to Ox have also been omitted.
!                               Use the original LINES_OF_SIGHT for these.
!     11/07/2018    "    "      Implementation and testing completed.
!     11/09/2018    "    "      Use 'los.g' for the output file name as part of
!                               adjusting NEQAIR_Integration to handle the
!                               original nomenclature and this simplified form.
!     10/04/2019    "    "      Resetting the lower intersection search interval
!                               to zero was not being done for each new line.
!                               Earlier testing must not have encountered any
!                               intersection failures that are caused by outer
!                               grid boundaries that are not strictly convex!
!     10/05/2019    "    "      Failed intersections prompted another way of
!                               retrying: discretize the part of the line beyond
!                               the failed solution and evaluate line-surface
!                               distances from those points; pick the interval
!                               containing the smallest distance, and redo the
!                               line-surface intersection.  REMEMBER THAT THE
!                               ARC LENGTH INTERVAL PASSED TO INTSEC6 SHOULD BE
!                               NORMALIZED BY THE LINE LENGTH IF THAT IS THE
!                               INTERVAL EXPECTED TO CONTAIN THE INTERSECTION,
!                               AS IT IS HERE, USING THE BOUNDING BOX DIAGONAL.
!     07/21/2020    "    "      A "small" tolerance in transform_quadrant was
!                               too big.  Body points near the nose (normally
!                               not the case for full angular integration)
!                               produced inaccurate rotation axes for the second
!                               quadrant of the transformed unit hemisphere.
!                               This led to irregular heat flux results near the
!                               heatshield apex on the centerline, now fixed.
!                               Avoiding a cross product of parallel unit
!                               vectors (zero length vector) was the root cause.
!     07/26/2020    "    "      Reluctant special-casing of the (0,0,0) body pt.
!                               that proves pathological in the sense that the
!                               associated surface grid cell is almost certain
!                               to be a corner cell for which the existing
!                               adjustment that forces the body normal to be in
!                               the symmetry place still doesn't force it to be
!                               along -Ox as it should be for a body of revolu-
!                               tion centered on Ox.  Remember that we normally
!                               don't apply the hemisphere LOS method to the
!                               nose region anyway, but it has been needed as
!                               part of comparing radiative heat flux results
!                               with those from a vasdtly more efficient method.
!     10/08/2020    "    "      Replace intsec6 with intsec9 as introduced for
!                               line_surface.f90 (interp3d library) used by the
!                               Stardust_Lines utility.  This packages a retry
!                               method that is virtually bulletproof.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure        ! Employed by the ADT module structured cases
   use xyzq_io_module              ! PLOT3D-type file I/O
   use tri_header_structure        ! Data structures used by triangulation_io
   use tri_zone_structure
   use triangulation_io            ! Tecplot-type triangulation I/O package
   use adt_utilities               ! All variants of ADT build & search routines

   implicit none

!  Local constants:

   integer, parameter :: &
      lunbps     = 1,    &   ! Body point (x,y,z) coordinates
      lunvolg    = 2,    &   ! Input volume grid (PLOT3D)
      lunlos     = 3,    &   ! Output line(s) of sight (PLOT3D)
      lunkbd     = 5,    &   ! Keyboard inputs
      luncrt     = 6,    &   ! Screen
      lunhemi1   = 7,    &   ! Output triangulated unit hemisphere quadrant
      lunhemi2   = 8,    &   ! Output body-tangent [half-]hemisphere triangulns.
      lunhemi3   = 9,    &   ! Output triangulation of LOS/boundary intersectns.
      luncone    = 10,   &   ! Output lines within the indicated cone angle
      ngeometric = 2,    &   ! Meaning none in the boundary layer
      nl         = 2         ! 2-point lines are passed to INTSEC6

   real, parameter ::    &
      half      = 0.5,   &
      one       = 1.0,   &
      ninety    = 90.0,  &
      rblayer   = 1.05,  &   ! Unused geometric growth rate in boundary layer
      zero      = 0.0

   logical, parameter ::   &
      false     = .false., &
      true      = .true.

   character (11), parameter :: &
      format = 'unformatted'

   character (1), parameter :: &
      lcs_method = 'L'         ! But LCSFIT isn't needed (hook for general case)

!  Local variables:

   integer :: &
      i, i1, ib, ic, in, ios, iquad, j, jc, jn, k, l, lbp, lid, m, n, nblocks, &
      nbps, ne, nelements, nf, ni, ninside, nj, nk, nlines, nnodes, noutside, &
      nquad, nquadrants

   real :: &
      area, bbox_diagonal, davg, dmax, dmean, d1norm, d2norm, ds1, &
      ds2_fraction, dsq, dsqmin, dtolsq, length, p, q, s, xmax, xmin, &
      xyz_interp(3), xyz_intersect(3), ymax, ymin, zmax, zmin

   logical :: &
      at_nose, cell_centered, cr, eof, formatted, hemisphere, oncenter

   integer, allocatable :: &
      conn(:,:)            ! For (patch,i,j) of surface quads. to search

   real, dimension (nl) :: &
      sl, xl, yl, zl       ! 2-pt. intersection line between body/outer boundary

   real, dimension (3) ::  &
      un, v1, v2           ! For a unit normal and rotation axis end points

   real, allocatable, dimension (:) :: &
      cone_angle, distsq,              & ! 90 or missing => full hemisphere
      sline, xline, yline, zline         ! For one radial volume grid line

   real, allocatable, dimension (:,:) :: &
      unit_norm, xyz_bp                  ! Unit normals, surface (x,y,z)s of bps

   character (128) :: filename

   type (grid_type), dimension(:), pointer :: &
      inner_surface, outer_surface, radial_lines, volume_grid
   type (tri_header_type) :: &
      tri_header                         ! One header for all transformations
   type (tri_type), pointer, dimension (:) :: &
      intersected_quadrant, transformed_quadrant, unit_hemi_quadrant

!  Execution:
!  ----------

!  Read the target body point coordinates [and optional cone angles?]:
!  -------------------------------------------------------------------

   call get_bp_data ()  ! Sets hemisphere = T if hemisphere lines are indicated
   if (ios /= 0) go to 99

!  Read the volume grid.  Only the k = 1 and k = nk planes are used.
!  -----------------------------------------------------------------

   ios = 1
   do while (ios /= 0)
      call reads (luncrt, 'Input volume grid file name (PLOT3D /mgrid): ', &
                  lunkbd, filename, cr, eof)
      if (eof) go to 99  ! Quit
      if (cr) cycle

      call determine_grid_form (filename, lunvolg, formatted, ios)
   end do

   i1 = 1;  if (formatted) i1 = 3
   open (lunvolg, file=filename, form=format(i1:11), status='OLD')

   call xyzq_read (lunvolg, -lunvolg, formatted, nblocks, nf, cell_centered,   &
                   volume_grid, ios)
   if (ios /= 0) go to 99

!  Extract the inner and outer boundaries as multiblock surface grids.

   allocate (inner_surface(nblocks), outer_surface(nblocks))

   do ib = 1, nblocks
      inner_surface(ib)%ni = volume_grid(ib)%ni
      outer_surface(ib)%ni = volume_grid(ib)%ni
      inner_surface(ib)%nj = volume_grid(ib)%nj
      outer_surface(ib)%nj = volume_grid(ib)%nj
      inner_surface(ib)%nk = 1
      outer_surface(ib)%nk = 1
   end do

   do ib = 1, nblocks
      call xyz_allocate (inner_surface(ib), ios)

      inner_surface(ib)%x(:,:,1) = volume_grid(ib)%x(:,:,1)
      inner_surface(ib)%y(:,:,1) = volume_grid(ib)%y(:,:,1)
      inner_surface(ib)%z(:,:,1) = volume_grid(ib)%z(:,:,1)

      call xyz_allocate (outer_surface(ib), ios)

      nk = volume_grid(ib)%nk
      outer_surface(ib)%x(:,:,1) = volume_grid(ib)%x(:,:,nk)
      outer_surface(ib)%y(:,:,1) = volume_grid(ib)%y(:,:,nk)
      outer_surface(ib)%z(:,:,1) = volume_grid(ib)%z(:,:,nk)
   end do

   call get_data_range ()  ! Bounding box diagonal, for intersection purposes
   if (ios /= 0) go to 99  ! May be wrong half flow soln. for a hemisphere case

   if (hemisphere) then
      call get_hemisphere_details ()  ! # edge points, etc.
      if (ios /= 0) go to 99
   else
      nlines = nbps
   end if

!  Set up the output line(s) of sight file:
!  ----------------------------------------

   open (lunlos, file='los.g', status='unknown') ! Assumed by NEQAIR_Integration

   nk = 121  ! Suited to 120 NEQAIR processors
   call readi (luncrt, '# pts. on each discretized line of sight? <cr>=121: ', &
               lunkbd, nk, cr, eof)
   if (eof) go to 99

   ds1 = 1.e-5
   call readr (luncrt, 'Wall spacing, ds1? <cr>=1.e-5: ', &
               lunkbd, ds1, cr, eof)
   d1norm = ds1/bbox_diagonal  ! We work with normalized distributions

   ds2_fraction = 0.1
   call readr (luncrt, &
               'Fraction of 1-sided ds2 for outer spacing? <cr>=0.1: ', &
               lunkbd, ds2_fraction, cr, eof)
   if (eof) go to 99

!  Count the surface quads (same for inner and outer surfaces):
!  ------------------------------------------------------------

   nquad = 0;  ymin = 1.e+30;  ymax = -ymin

   do ib = 1, nblocks
      ni = inner_surface(ib)%ni
      nj = inner_surface(ib)%nj
      nquad = (ni - 1)*(nj - 1) + nquad
   end do

   allocate (conn(3,nquad)) ! For (patch,i,j) of each surface quad.

!  Build a search tree for the inner surface and do the body point searching:

   call build_adt (nblocks, inner_surface, nquad, conn)

   ninside  = 0;  dmax  = zero
   noutside = 0;  dmean = zero
   dtolsq = (bbox_diagonal*1.e-6)**2

   do l = 1, nbps

      write (luncrt, '(a, 3es16.8)') 'BP:', xyz_bp(:,l)
      call search_adt (xyz_bp(:,l), iquad, p, q, dsqmin, true, nblocks, &
                       inner_surface, nquad, conn, xyz_interp)

      if (dsqmin < dtolsq) then ! The nearest triangle was within tolerance
         ninside  = ninside + 1
      else
         noutside = noutside + 1
      end if

      write (luncrt, '(a, es16.8)') 'dsqmin:', dsqmin
      dmax  = max (dmax, dsqmin)
      dmean = dmean + sqrt (dsqmin)
      oncenter = xyz_bp(2,l) == zero ! We ensure a normal in the symmetry plane
      at_nose  = oncenter .and. xyz_bp(1,l) == zero .and. xyz_bp(3,l) == zero
      if (.not. at_nose) xyz_bp(:,l) = xyz_interp(:) ! Else ensure that the body
                                                     ! pt. is on the surface
      write (luncrt, '(a, l2)') 'oncenter:', oncenter
      write (luncrt, '(a, l2)') 'at_nose: ', at_nose

!     Construct a carefully-calculated unit normal at the current body point.

      n = conn(1,iquad) ! Patch # at inner surface for this target surface point
      i = conn(2,iquad) ! Corresponding cell indices (lower left)
      j = conn(3,iquad)

      write (luncrt, '(a, 3i4)') 'n,i,j:', n,i,j
      write (luncrt, '(a, 2es16.8)') 'p,q:', p,q

      call surface_normal (inner_surface(n)%ni, inner_surface(n)%nj, &
                           inner_surface(n)%x,  inner_surface(n)%y,  &
                           inner_surface(n)%z,  i, j, p, q, unit_norm(:,l))

      write (luncrt, '(a, 3es19.11)') 'Body normal:', unit_norm(:,l)
      if (oncenter) then
         unit_norm(2,l) = sqrt (unit_norm(1,l)**2 + unit_norm(3,l)**2)
         unit_norm(1,l) = unit_norm(1,l) / unit_norm(2,l)
         unit_norm(3,l) = unit_norm(3,l) / unit_norm(2,l)
         unit_norm(2,l) = zero
         xyz_bp(2,l)    = zero
         write (luncrt, '(a, 3f12.8)') ' Adjusted body point:  ', xyz_bp(:,l)
         write (luncrt, '(a, 3f12.8)') ' Adjusted body normal: ', unit_norm(:,l)

!        Special handling of the pathological nose-point case (that assumes the
!        geometry is a body of revolution centered about Ox):

         if (at_nose) then
            unit_norm(1,l) = -one
            unit_norm(2,l) = zero  ! (Already so)
            unit_norm(3,l) = zero
            write (luncrt, '(a, 3f12.8)') ' Corrected nose normal:', &
               unit_norm(:,l)
         end if

      end if

   end do ! Next target surface point

   dmax = sqrt (dmax);  dmean = dmean / real (nbps)

   write (luncrt, '(/, a, 2i5, /, a, 1p, 2e12.5)') &
      ' # surface points inside/outside tolerance:', ninside, noutside, &
      ' max & mean distance:', dmax, dmean

!  Set up storage for the output lines of sight.
!  They are saved in (1,1,nk) multiblock PLOT3D form.
!  --------------------------------------------------

   allocate (radial_lines(nlines), stat=ios)
   if (ios /= 0) then
      write (luncrt, '(/, a, i10)') &
         ' Allocation trouble with radial_lines. nbps:', nbps
      go to 99
   end if

   do l = 1, nlines
      radial_lines(l)%ni = 1
      radial_lines(l)%nj = 1
      radial_lines(l)%nk = nk

      call xyz_allocate (radial_lines(l), ios)
      if (ios /= 0) then
         write (luncrt, '(/, a, i10)') &
            ' Allocation trouble with radial_lines(l). l:', l
         write (luncrt, '(/, a, 3i10)') ' ni, nj, nk: ', 1, 1, nk
         go to 99
      end if

   end do

!  Work space for point distributions:

   allocate (xline(nk), yline(nk), zline(nk), sline(nk), distsq(nk), stat=ios)
   if (ios /= 0) then
      write (luncrt, '(/, a, i10)') &
         ' Allocation trouble with x/y/z/sline. nk:', nk
      go to 99
   end if

!  Build a new search tree from the outer shock boundary:
!  ------------------------------------------------------

   call release_adt ()    ! Deallocate work-space from first ADT

   call build_adt (nblocks, outer_surface, nquad, conn)

!  Calculate each line of sight of the indicated type:
!  ---------------------------------------------------

   if (.not. hemisphere) then  ! Body-normal line(s) only -- do it in-line

      ninside  = 0;  dmax  = zero
      noutside = 0;  dmean = zero

      do l = 1, nbps

         xl(1)  = xyz_bp(1,l);  xl(nl) = xl(1) + bbox_diagonal*unit_norm(1,l)
         yl(1)  = xyz_bp(2,l);  yl(nl) = yl(1) + bbox_diagonal*unit_norm(2,l)
         zl(1)  = xyz_bp(3,l);  zl(nl) = zl(1) + bbox_diagonal*unit_norm(3,l)
         sl(1)  = zero          ! Lower search interval limit
!!!      sl(nl) = bbox_diagonal ! Upper limit on intersection search along line;
         sl(nl) = one           ! overkill for forebody; safe for aft body
                                ! These are NORMALIZED arc lengths

         call intsec6 (nblocks, outer_surface, nquad, conn, nl, xl, yl, zl,    &
                       sl, lcs_method, iquad, p, q, s, xyz_intersect, dsqmin)

!        The output s is normalized.
         write (luncrt, '(a,4i10)') 'iquad,conn(:,iquad):', iquad,conn(:,iquad)
         write (luncrt, '(a,4es19.11)') 'p, q, s, bbox:', p, q, s, bbox_diagonal
         write (luncrt, '(a,4es19.11)') 'xyzint, dsqmin:', xyz_intersect, dsqmin

         if (dsqmin < dtolsq) then ! The nearest quad was within tolerance
            ninside  = ninside + 1
         else
            noutside = noutside + 1
         end if

         dmax  = max (dmax, dsqmin)
         dmean = dmean + sqrt (dsqmin)

!        Discretize the intersected line.
!        s is the normalized arc length along the line at the intersection.
!        EXPDIS5 is equivalent to Vinokur's 1-sided tanh-based method.
!        MODE 1 clusters towards the low end.

         call expdis5 (1, zero, s, d1norm, nk, sline, -luncrt)

         d2norm = (s - sline(nk-1)) * ds2_fraction

!        Pure 2-sided Vinokur if ngeometric = 2 else geometric in b.layer:

         call blgrid (nk, d1norm, d2norm, ngeometric, rblayer, sline, luncrt, &
                      ios)
         do k = 1, nk
            radial_lines(l)%x(1,1,k) = xl(1) + ((xl(nl) - xl(1)))*sline(k)
            radial_lines(l)%y(1,1,k) = yl(1) + ((yl(nl) - yl(1)))*sline(k)
            radial_lines(l)%z(1,1,k) = zl(1) + ((zl(nl) - zl(1)))*sline(k)
         end do

      end do  ! Next line of sight

      dmax = sqrt (dmax);  dmean = dmean / real (nbps)
      write (luncrt, '(/, a, 2i5, /, a, 1p, 2e12.5)') &
         ' # outer intersection points inside/outside tolerance:', &
         ninside, noutside, ' max & mean distance:', dmax, dmean

   else  ! One body point/many lines of sight

      call hemisphere_lines ();  if (ios /= 0) go to 99
!     ------------------------

   end if

!  Save the lines of sight as a PLOT3D formatted file:
!  ---------------------------------------------------

   call xyz_write (lunlos, true, nlines, radial_lines, ios)

!  Let Fortran do all the deallocations.

99 continue

!  Internal procedures for program SLOS:
!  -------------------------------------

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_bp_data ()  ! The body pt. file may include cone angles
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: lbuf, n, nitems
      logical :: ascii
      character (3)  :: seps
      character (64) :: buffer

!     Execution:

      ascii = true
      call file_prompt (lunbps, 0, 'body point x/y/z[/cone angle]', 'old', &
                        false, ascii, ios)
      if (ios /= 0) go to 99

      nbps = 0
      do ! Count the target points till EOF
         read (lunbps, '(a)', iostat=ios) buffer
         if (ios /= 0) exit
         nbps = nbps + 1
      end do
      rewind (lunbps)

      write (luncrt, '(/, a, i4)') ' # target points specified: ', nbps

      allocate (xyz_bp(3,nbps), unit_norm(3,nbps), cone_angle(nbps))

      lbuf = index (buffer, '!') - 1  ! Avoid counting comment tokens
      if (lbuf < 0) lbuf = len_trim (buffer)

      seps = ' ,' // char (9)  ! Space, comma, or tab

      call token_count (buffer(1:lbuf), seps, nitems)  ! Should be 3 or 4

      if (nitems == 3) then
         write (luncrt, '(a)') ' Not expecting cone angle with BP inputs.'
      else
         write (luncrt, '(a)') ' Expecting cone angle with BP inputs.'
      end if

      cone_angle(:) = 90.
      do n = 1, nbps
         if (nitems == 3) then
            read (lunbps, *, iostat=ios) xyz_bp(:,n)
         else
            read (lunbps, *, iostat=ios) xyz_bp(:,n), cone_angle(n)
         end if
         if (ios /= 0) then
            write (luncrt, '(a, i5)') 'Trouble reading body point', n
            go to 99
         end if
      end do
      close (lunbps)

      hemisphere = nbps == 1
      if (hemisphere) then  ! Allow for a single body-normal line of sight
         call ready (luncrt, &
             '1 body pt. => hemisphere lines? ' // &
             'y|n; y=<cr>=yes; n=body-normal line: ', &
             lunkbd, hemisphere, cr, eof)
         if (eof) then
           ios = 1;  go to 99
         end if
      end if

 99   return

      end subroutine get_bp_data

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_data_range ()  ! .. of the shock boundary, for intersectns.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      xmin = 1.e+10
      xmax = -xmin
      ymin =  xmin
      ymax =  xmax
      zmin =  xmin
      zmax =  xmax

      do ib = 1, nblocks
         do j = 1, outer_surface(ib)%nj
            do i = 1, outer_surface(ib)%ni
               xmax = max (xmax, outer_surface(ib)%x(i,j,1))
               xmin = min (xmin, outer_surface(ib)%x(i,j,1))
               ymax = max (ymax, outer_surface(ib)%y(i,j,1))
               ymin = min (ymin, outer_surface(ib)%y(i,j,1))
               zmax = max (zmax, outer_surface(ib)%z(i,j,1))
               zmin = min (zmin, outer_surface(ib)%z(i,j,1))
            end do
         end do
      end do

      bbox_diagonal = sqrt ((xmax-xmin)**2 + (ymax-ymin)**2 + (zmax-zmin)**2)
      length = bbox_diagonal  ! Shorter name

      write (luncrt, '(a, 6es14.6)') &
         'xmin, xmax, ...', xmin, xmax, ymin, ymax, zmin, zmax, &
         'Bounding box diagonal length:', length

!     For the hemisphere case, it is too awkward to handle either half of the
!     geometry when the body point is on the centerline.  Insist on the right
!     half, or both halves for an off-center body point:

      ios = 0
      if (hemisphere) then
         oncenter = xyz_bp(2,1) == zero
         if (ymax < (ymax - ymin)*0.1) then
            if (oncenter) then
               write (luncrt, '(/, (a))') &
                  ' The geometry appears to be the left half (y < 0).', &
                  ' A centerline body point requires the right half (y >= 0).'
               ios = -1
            end if
          end if
          if (.not. oncenter) then
            if (min (abs (ymax), abs (ymin)) < (ymax - ymin)*0.1) then
               write (luncrt, '(/, (a))') &
                  ' The geometry appears to be a half body, but', &
                  ' an off-center body point requires the whole body.'
               ios = -1
            end if
         end if
      end if

      end subroutine get_data_range

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_hemisphere_details ()
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ne = 0
      do while (ne <= 1 .or. ne > 1000)
         ne = 25
         call readi (luncrt, &
               '# hemisphere pts., pole to equator? <cr>=25=650/1300 lines: ', &
                     lunkbd, ne, cr, eof)
         ios = -1;  if (eof) exit

         ios = 0
         nelements  = (ne - 1)**2    ! # triangles in a hemisphere quadrant
         nnodes     = (ne + 1)*ne/2  ! # vertices   "    "    "    "    "
         nquadrants = 4;  if (oncenter) nquadrants = 2
         nlines     = nnodes*nquadrants  ! Too hard to suppress edge duplicates
         write (luncrt, '(/, a, i5)') '# lines of sight indicated:', nlines
      end do

 99   return

      end subroutine get_hemisphere_details

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine hemisphere_lines ()
!
!     See the main program header for the method description.  A single body
!     point is assumed.  A triangulated quadrant of a hemisphere defines the
!     lines of sight generated.  Some variables are necessarily global.  Where
!     feasible, they are local.  File names for hemisphere-related outputs are
!     hard-coded rather than using a run-time identifier.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local constants:

      character (4), parameter :: identifier = 'hemi'

!     Local variables:

      integer :: &
         i, m, nn
      real, dimension (3) ::  &
         un, v1, v2              ! For unit normals and rotation axis end pts.

      type (tri_header_type) :: &
         tri_header              ! One header for various hemi-related files

!     Execution:
!     ----------

!     Generate the indicated triangulation of a quadrant of a unit hemisphere:

      allocate (unit_hemi_quadrant(1))  ! Not a scalar because it's a pointer

      call spherical_triangulation (ne, unit_hemi_quadrant(1))

!     This unit quadrant has all of x, y, z in [0, 1], and the rotational trans-
!     formations that complete the hemisphere should be about its Oz axis, not
!     Ox, because the triangulation is produced via constant-z slices.  At the
!     nose it needs to start outside the body, not inside as initially, then be
!     shifted and rotated to align its original Oz axis with the body normal.
!     It is less confusing if we rotate the initial quadrant 90 degrees so its
!     original Oz axis aligns with -Ox.

      v1(:) = zero  ! End points of rotation axis
      v2(:) = zero
      v2(2) = -one

      call rotate_xyz (nnodes, unit_hemi_quadrant(1)%xyz, ninety, v1, v2)

      unit_hemi_quadrant(1)%zone_title   = 'Unit hemisphere quadrant'
      unit_hemi_quadrant(1)%element_type = 'TRIANGLE'

!     Save it for visualization as a Tecplot surface triangulation:

      tri_header%filename    = 'unit_hemi_quadrant_vertices.dat'
      tri_header%fileform    = 1        ! Vertex-centered
      tri_header%formatted   = true
      tri_header%nvertices   = 3        ! Triangles, not tets
      tri_header%numf        = 0
      tri_header%nzones      = 1
      tri_header%datapacking = 0        ! Point order
      tri_header%title       = '# edge points:'
      write (tri_header%title(15:18), '(i4)') ne
      allocate (tri_header%varname(3))
      tri_header%varname(1)  = 'x'
      tri_header%varname(2)  = 'y'
      tri_header%varname(3)  = 'z'

      call tri_write (lunhemi1, tri_header, unit_hemi_quadrant, ios)
      if (ios /= 0) go to 99

!     Two further triangulations are produced, each with either 2 or 4 quadrant
!     zones: (a) a unit [half-]hemisphere, tangent to the body at the body pt.;
!     (b) the surface triangulation produced by intersecting lines through each
!     vertex of (a) with the outer boundary mesh.

      allocate (transformed_quadrant(1))
      allocate (intersected_quadrant(1))

      transformed_quadrant(1)%nelements = nelements
      intersected_quadrant(1)%nelements = nelements
      transformed_quadrant(1)%nnodes    = nnodes
      intersected_quadrant(1)%nnodes    = nnodes

      call tri_zone_allocate (tri_header, transformed_quadrant(1), ios)
      if (ios /= 0) go to 99

      call tri_zone_allocate (tri_header, intersected_quadrant(1), ios)
      if (ios /= 0) go to 99

      transformed_quadrant(1)%conn = unit_hemi_quadrant(1)%conn
      intersected_quadrant(1)%conn = unit_hemi_quadrant(1)%conn

      deallocate (unit_hemi_quadrant(1)%conn)

      tri_header%filename = 'transformed_unit_hemi.dat'
      call tri_header_write (lunhemi2, tri_header, transformed_quadrant, ios)
      if (ios /= 0) go to 99

      tri_header%filename = 'intersected_boundary.dat'
      call tri_header_write (lunhemi3, tri_header, transformed_quadrant, ios)
      if (ios /= 0) go to 99

!     The output lines of sight have been allocated as for the body-normal case.
!     The (one) unit body normal has also been constructed as for that case.

      un(:) = unit_norm(:,1)  ! Simpler notation
      xl(1) = xyz_bp(1,1)
      yl(1) = xyz_bp(2,1)
      zl(1) = xyz_bp(3,1)

      transformed_quadrant(1)%zone_title   = 'Quadrant x'
      transformed_quadrant(1)%element_type = 'TRIANGLE'
      intersected_quadrant(1)%element_type = 'TRIANGLE'

      ninside  = 0;  dmax  = zero
      noutside = 0;  dmean = zero

      l = 0
      do m = 1, nquadrants

         write (transformed_quadrant(1)%zone_title(10:10), '(i1)') m
         intersected_quadrant(1)%zone_title = transformed_quadrant(1)%zone_title

         call transform_quadrant (un, m)

         call tri_zone_write (lunhemi2, tri_header, transformed_quadrant(1), &
                              ios)
         if (ios /= 0) go to 99

         do nn = 1, nnodes
            l = l + 1  ! LOS count over all quadrants

!           2-pt. line through this vertex/node more than long enough:

            xl(nl) = xl(1) + (transformed_quadrant(1)%xyz(1,nn) - xl(1))*length
            yl(nl) = yl(1) + (transformed_quadrant(1)%xyz(2,nn) - yl(1))*length
            zl(nl) = zl(1) + (transformed_quadrant(1)%xyz(3,nn) - zl(1))*length
            sl(1)  = zero    ! Intersection search interval limits
            sl(nl) = one     ! (normalized)

!!!         call xyz_allocate (radial_lines(l), ios) ! Done as for non-hemi case
!!!         if (ios /= 0) go to 99

            call intersect_line (l, m, nn)

            intersected_quadrant(1)%xyz(1,nn) = radial_lines(l)%x(1,1,nk)
            intersected_quadrant(1)%xyz(2,nn) = radial_lines(l)%y(1,1,nk)
            intersected_quadrant(1)%xyz(3,nn) = radial_lines(l)%z(1,1,nk)

         end do  ! Next line of this quadrant

         call tri_zone_write (lunhemi3, tri_header, intersected_quadrant(1), &
                              ios)
         if (ios /= 0) go to 99

      end do  ! Next quadrant

      dmax = sqrt (dmax);  dmean = dmean / real (nlines)
      write (luncrt, '(/, a, 2i5, /, a, 2es12.5)') &
         '# outer intersection points inside/outside tolerance:', &
         ninside, noutside, ' max & mean distance:', dmax, dmean

!!!   call xyz_write (lunlos, true, nlines, radial_lines, ios) ! Common to body-
                                                               ! normal case
!     Allow for running NEQAIR on only the lines within the cone angle:

      if (cone_angle(1) < ninety) call save_cone_only ()

 99   return

      end subroutine hemisphere_lines

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine transform_quadrant (un, m)
!
!     A quadrant of a triangulated unit hemisphere is transformed 2 or 4 times
!     to define all the lines of sight for one body point.  Processing for each
!     quadrant is completed before the next.  For the first quadrant (m = 1),
!     the unit quadrant is transformed so that its primary axis along Ox is
!     coincident with the unit normal at the body point. Each next quadrant is
!     obtained by rotating the previous quadrant 90 degrees about the body
!     normal.  The un(:) argument is needed because subroutine hemisphere_lines
!     is a local procedure here but it is the essence of the main program in
!     HEMISPHERES_OF_SIGHT, much of which has been transcribed.  The muddle of
!     global and local variables is the result of sharing code with the portion
!     that has been transcribed from LINES_OF_SIGHT.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use trigd

!     Argument:

      real,    intent (in) :: un(3)  ! Unit normal at the body point
      integer, intent (in) :: m      ! Quadrant number (1:nquadrants)

!     Local constants:

      real, parameter :: small_value = 0.000001

!     Local variables:

      integer       :: i
      real, save    :: theta
      logical, save :: parallel

!     Execution:

      if (m == 1) then  ! Align the unit hemisphere quadrant w/ the body normal

         v1(1) = -one;  v1(2) = zero;  v1(3) = zero  ! Unit vector along -Ox

         write (luncrt, '(a, 3f12.8)') ' Body pt. x,y,z:', xyz_bp(:,1)
         write (luncrt, '(a, 3f12.8)') ' Unit normal:   ', un(:)

         call cross (v1, un, v2)     ! v2 = v1 x un provides a second point on
                                     ! a rotation axis through the body point,
                                     ! and also the initial angle of rotation

         write (luncrt, '(a, 3f12.8)') ' -Ox x un:      ', v2(:)

!        Avoid theta = 90, where tan (theta) is infinite:

         if (abs (un(1)) < small_value) then
            theta = sign (ninety, un(3))
         else                                                      ! |v1 x un| /
            theta = atand (sqrt (un(2)**2 + un(3)**2) / (-un(1)))  !  v1 . un
            if (un(1) > zero) theta = theta + ninety + ninety      ! Aft body
         end if

         write (luncrt, '(a, 3f14.8)') ' theta  :', theta

!        Avoid parallel vectors, with zero-magnitude cross-product:

         parallel = dot_product (v2, v2) < small_value  ! Normal is along Ox

         v1(:) = xyz_bp(:,1)         ! First point on initial rotation axis

         if (parallel) then          ! Second  "   "   "   "   "   "   "
            v2(:) = v1(:)
            v2(2) = v2(2) + one      ! Axis is parallel to Oy
            write (luncrt, '(a)') 'Evidently the body normal is along Ox.'
         else
            v2(:) = v1(:) + v2(:)
         end if


         write (luncrt, '(a, 3f12.8)') ' Axis p2:', v2(:)

!        The first transformation requires an initial shift:

         do i = 1, nnodes
            transformed_quadrant(1)%xyz(:,i) = &
               unit_hemi_quadrant(1)%xyz(:,i) + v1(:)
         end do

      else if (m == 2) then     ! The extra quadrants rotate the previous
                                ! quadrant about the body normal
         v2(:) = v1(:) + un(:)  ! Second point along normal
         theta = ninety

      end if

!     Rotate in-place:

      call rotate_xyz (nnodes, transformed_quadrant(1)%xyz, theta, v1, v2)

      end subroutine transform_quadrant

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine intersect_line (l, m, nn)
!
!     Intersect the line of sight defined by the current vertex of the
!     current transformed quadrant, and impose the specified distribution.
!     The arguments are necessary for diagnostics because hemisphere_lines is
!     also a local procedure.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in) :: l  ! Line number over all lines
      integer, intent (in) :: m  ! Quadrant number
      integer, intent (in) :: nn ! Node number of quadrant m

!     Local constants:

!     Local variables:

      real :: fraction

!     Execution:

      call intsec9 (l, nblocks, outer_surface, nquad, conn, nl, xl, yl, zl, &
                    sl, lcs_method, iquad, p, q, s, xyz_intersect, dsqmin)

!!!   write (luncrt, '(a, 2i4, 4i10)') &
!!!      ' vertex, iquad, conn(:,iquad):', nn, iquad, conn(:,iquad)
!!!   write (luncrt, '(a, 3es19.11)') ' p, q, s:   ', p, q, s, &
!!!                                   ' xyz_intersect:', xyz_intersect

      if (dsqmin < dtolsq) then ! The nearest cell was within tolerance
         ninside  = ninside + 1
      else
         noutside = noutside + 1
         write (luncrt, '(a, i6, 4es13.6)') &
         ' *** Intersection trouble.  Line #, sa, sb, dsq, s:', l, sl, dsqmin, s

         ! Keep going
      end if

      dmax  = max (dmax, dsqmin)
      dmean = dmean + sqrt (dsqmin)

!     Discretize the intersected line.
!     EXPDIS5 is equivalent to Vinokur's 1-sided tanh-based method.
!     s is the normalized arc length along the line at the intersection.
!     MODE 1 clusters towards the low end.

      call expdis5 (1, zero, s, d1norm, nk, sline, -luncrt)

      d2norm = (s - sline(nk-1)) * ds2_fraction

!     Pure 2-sided Vinokur if ngeometric = 2 else geometric in b.layer:

      call blgrid (nk, d1norm, d2norm, ngeometric, rblayer, sline, luncrt, ios)

      do k = 1, nk
         radial_lines(l)%x(1,1,k) = xl(1) + ((xl(nl) - xl(1)))*sline(k)
         radial_lines(l)%y(1,1,k) = yl(1) + ((yl(nl) - yl(1)))*sline(k)
         radial_lines(l)%z(1,1,k) = zl(1) + ((zl(nl) - zl(1)))*sline(k)
      end do

      end subroutine intersect_line

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine save_cone_only ()  ! Suppress lines outside the cone angle
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i, ncone
      logical, allocatable :: keep(:)
      real :: angle, vnormal(3), v2(3)
      type (grid_type), pointer, dimension (:) :: cone_lines
      character (12) :: filename

!     Execution:

      filename = 'cone_xx.xx.g'
      write (filename(6:10), '(f5.2)') cone_angle(1)
      open (luncone, file=filename, status='unknown')

      vnormal(1) = radial_lines(1)%x(1,1,nk) - radial_lines(1)%x(1,1,1)
      vnormal(2) = radial_lines(1)%y(1,1,nk) - radial_lines(1)%y(1,1,1)
      vnormal(3) = radial_lines(1)%z(1,1,nk) - radial_lines(1)%z(1,1,1)
      allocate (keep(nlines))
      keep(1)  = true
      keep(2:) = false

      ncone = 1
      do i = 2, nlines
         v2(1) = radial_lines(i)%x(1,1,nk) - radial_lines(1)%x(1,1,1)
         v2(2) = radial_lines(i)%y(1,1,nk) - radial_lines(1)%y(1,1,1)
         v2(3) = radial_lines(i)%z(1,1,nk) - radial_lines(1)%z(1,1,1)
         call angle_between_vectors (vnormal, v2, angle)
         if (angle > cone_angle(1)) cycle
            ncone = ncone + 1
            keep(i) = true
      end do

      allocate (cone_lines(ncone))
      cone_lines(:)%ni = 1
      cone_lines(:)%nj = 1
      cone_lines(:)%nk = nk

      ncone = 0
      do i = 1, nlines
         if (.not. keep(i)) cycle
            ncone = ncone + 1
            allocate (cone_lines(ncone)%x(1,1,nk), &
                      cone_lines(ncone)%y(1,1,nk), &
                      cone_lines(ncone)%z(1,1,nk))
            cone_lines(ncone)%x(:,:,:) = radial_lines(i)%x(:,:,:)
            cone_lines(ncone)%y(:,:,:) = radial_lines(i)%y(:,:,:)
            cone_lines(ncone)%z(:,:,:) = radial_lines(i)%z(:,:,:)
      end do

      call xyz_write (luncone, true, ncone, cone_lines, ios)

      do i = 1, ncone
         deallocate (cone_lines(i)%x, cone_lines(i)%y, cone_lines(i)%z)
      end do
      deallocate (keep, cone_lines)

      end subroutine save_cone_only

   end program slos
