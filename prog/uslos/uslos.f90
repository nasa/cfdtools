!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program uslos
!
!  Description:
!
!     USLOS (unstructured grid lines of sight) combines the functions of the
!     LINES_OF_SIGHT and HEMISPHERES_OF_SIGHT utilities for the case of 3-space
!     volume grids that are NOT multiblock volume grids.  Since these are purely
!     geometric functions, it is assumed that the needed inner and outer (shock)
!     boundaries can be provided as triangulated surfaces, and that these are
!     compatible with Tecplot (a de facto standard at NASA).  The lines of sight
!     are expected to be used for radiative heat flux calculations.
!
!     Anticipated usage is for either a list of forebody points or a single body
!     point on the aft body or perhaps near the shoulder of a typical atmos-
!     pheric entry vehicle.  More than one set of hemisphere lines of sight is
!     NOT recommended, so if the body point list contains multiple entries, it
!     is assumed that only the body-normal lines are required.
!
!     The inner surface is expected to be right-handed.  This means that for
!     vertices 1, 2, 3 of a triangular element, the body normal vector is de-
!     fined by the cross product vector w = (v2 - v1) x (v3 - v1) and this
!     should point into the flow field (i.e., outwards from what is presumed to
!     be a mostly-convex surface).  Any concavities may cause difficulties such
!     as those handled by ancillary utility ADJUST_HEMISPHERE_LOS for the
!     structured surface case, which will need an analogue for the unstructured
!     case.  If the outer surface is derived from the same (right-handed) volume
!     grid, this means that a vector defined similarly to be normal to the outer
!     surface would point away from the body. But no such normal is ever needed,
!     so the outer surface need not be right-handed.  It is used only for line-
!     surface intersection calculations.
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
!  Input Surface Triangulation Format (inner and outer, in separate files):
!
!     TITLE = "Dataset title"
!     VARIABLES = "x", "y", "z"[, "f1"[, "f2"] ...]
!     ZONE T="ZONE 1", NODES=450816, ELEMENTS=150272,
!     DATAPACKING=POINT, ZONETYPE=FETRIANGLE
!      3.5443134E-01  1.0169792E+00  5.8715301E-01 [optional functions, ignored]
!      3.5443952E-01  1.0356432E+00  5.5356261E-01     "         "         "
!           :              :              :
!           :              :              :
!      2.8716209E+00 -1.1242014E-01  6.4854415E-02  0.0000000E+00
!           1       2       3
!           4       5       6
!           :       :       :
!           :       :       :
!
!      450814  450815  450816
!
!     Multizone ASCII files are permitted.  Unformatted file input is not sup-
!     ported.  Any function data are ignored.  One way of triangulating a multi-
!     block PLOT3D-type structured PLOT3D-type surface is to run HYPER_AERO and
!     use the Tecplotable output of Newtonian Cp predictions.  The Cps will be
!     ignored.  They are vertex-centered, so tri_header%fileform = 1 here.
!
!  Outputs:
!
!     PLOT3D-type output is retained as for the structured case (one or more
!     blocks, one line of sight per block).  Ensuing procedures will interpolate
!     flow data onto the discretized lines of sight and write results compatible
!     with the radiation solver.  The file name is hard-coded as 'los.g' as part
!     of simplifying the original nomenclature.
!
!  Method:
!
!     >  A handful of prompts serve to drive the program -- no control file.
!     >  Read the indicated body point (x,y,z) file -- one or more points.
!     >  If a single body point is found, prompt in case hemispherical lines are
!        NOT intended.
!     >  If all body points are on the centerline (y = 0.), the input surfaces
!        are assumed to be the starboard (right) half.  That half (y > 0.)
!        also serves for off-center body-normal lines, but for a hemisphere at
!        an off-center body point, a full volume grid is expected, with the
!        outer surface closed and convex.
!     >  Read an inner surface triangulation then an outer surface triangulation
!        from what is assumed to be a shock-aligned volume grid.
!     >  Build a search tree from the inner surface triangulation, and search
!        for each body point.
!     >  Generate and store the corresponding unit body-normal vector(s).
!     >  Read the indicated outer boundary triangulation, and build a new search
!        tree as needed for line-surface intersection calculations.
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
!                 normal line (INTSEC8).
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
!     09/27/2018  D.A.Saunders  Initial design.
!     10/20/2018     "    "     After a hiatus, finished the body-normal case.
!     10/22/2018     "    "     Allowed for multizone surfaces.
!     10/31/2018     "    "     Finished adding & testing the hemisphere option.
!     11/01/2018     "    "     Filled in the missing hemisphere method outline.
!     11/06/2018     "    "     Corrected some of the description above.  Note
!                               SLOS, the structured-grid analogue of USLOS, is
!                               now an alternative to the original utilities.
!     11/09/2018     "    "     Use 'los.g' for the output file name as part of
!                               adjusting NEQAIR_Integration to handle the
!                               original nomenclature and this simplified form.
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
   use adt_utilities               ! ADT build & search routines

   implicit none

!  Local constants:

   integer, parameter :: &
      lunbps     = 1,    &   ! Body point (x,y,z) coordinates
      lunbody    = 2,    &   ! Input body surface triangulation (Tecplot)
      lunshock   = 3,    &   ! Input shock boundary triangulation (Tecplot)
      lunlos     = 4,    &   ! Output line(s) of sight (PLOT3D)
      lunkbd     = 5,    &   ! Keyboard inputs
      luncrt     = 6,    &   ! Screen
      lunhemi1   = 7,    &   ! Output triangulated unit hemisphere quadrant
      lunhemi2   = 8,    &   ! Output body-tangent [half-]hemisphere triangulns.
      lunhemi3   = 9,    &   ! Output triangulation of LOS/boundary intersectns.
      luncone    = 10,   &   ! Output lines within the indicated cone angle
      ngeometric = 2,    &   ! Meaning none in the boundary layer
      nl         = 2         ! 2-point lines are passed to INTSEC8

   real, parameter ::    &
      half      = 0.5,   &
      one       = 1.0,   &
      ninety    = 90.0,  &
      rblayer   = 1.05,  &   ! Unused geometric growth rate in boundary layer
      zero      = 0.0

   logical, parameter ::   &
      false     = .false., &
      true      = .true.

   character, parameter :: &
      lcs_method * 1 = 'L'    ! But LCSFIT isn't needed (hook for general case)

!  Local variables:

   integer :: &
      i1, i2, i3, in, ios, itri, k, l, nbps, ne, nelements, ninside, nk, &
      nlines, nnode, nnodes, noutside, nquadrants, ntri

   real :: &
      area, bbox_diagonal, davg, dmax, dmean, d1norm, d2norm, ds1, &
      ds2_fraction, dsq, dsqmin, dtolsq, p, q, r, s, xmax, xmin, &
      xyz_interp(3), xyz_intersect(3), ymax, ymin, zmax, zmin

   logical :: &
      cr, eof, formatted, hemisphere, oncenter

   real, dimension (nl) :: &
      sl, xl, yl, zl       ! 2-pt. intersection line between body/outer boundary

   real, dimension (3) ::  &
      un, v1, v2           ! For a unit normal and rotation axis end points

   real, allocatable, dimension (:) :: &
      cone_angle,                      & ! 90 or missing => full hemisphere
      sline, xline, yline, zline         ! For one radial volume grid line

   real, allocatable, dimension (:,:) :: &
      unit_norm, xyz_bp                  ! Unit normals, surface (x,y,z)s of bps

   character (128) :: filename

   type (grid_type), dimension(:), pointer :: &
      radial_lines
   type (tri_header_type)                  :: &
      tri_header_body, tri_header_shock
   type (tri_type), dimension (:), pointer :: &
      tri_body, tri_shock
   type (tri_type), pointer, dimension (:) :: &
      intersected_quadrant, transformed_quadrant, unit_hemi_quadrant

!  Execution:
!  ----------

!  Read the target body point coordinates [and optional cone angles?]:
!  -------------------------------------------------------------------

   call get_bp_data ()  ! Sets hemisphere = T if hemisphere lines are indicated
   if (ios /= 0) go to 99

!  Prompt for, check existence of, and read the inner and outer volume surfaces:
!  -----------------------------------------------------------------------------

   call reads (luncrt, 'Body surface triangulation?  ', lunkbd, &
               tri_header_body%filename, cr, eof)
   if (cr .or. eof) go to 99

   open (lunbody, file=tri_header_body%filename, status='old', iostat=ios)
   if (ios /= 0) then
      write (luncrt, '(2a)') 'Cannot open ', trim (tri_header_body%filename)
      go to 99
   end if

   close (lunbody)  ! tri_read opens it
   ios = 1  ! Verbose mode
   tri_header_body%fileform  = 1  ! Vertex-centered is assumed for any fns.
   tri_header_body%formatted = true
   tri_header_body%nvertices = 3
   tri_header_body%combine_zones = true
   call tri_read (lunbody, tri_header_body, tri_body, ios)
   if (ios /= 0) go to 99

   call reads (luncrt, 'Shock surface triangulation? ', lunkbd, &
               tri_header_shock%filename, cr, eof)
   if (cr .or. eof) go to 99

   open (lunshock, file=tri_header_shock%filename, status='old', iostat=ios)
   if (ios /= 0) then
      write (luncrt, '(2a)') 'Cannot open ', trim (tri_header_shock%filename)
      go to 99
   end if

   close (lunshock)  ! tri_read opens it
   ios = 1  ! Verbose mode
   tri_header_shock%fileform  = 1  ! Vertex-centered is assumed for any fns.
   tri_header_shock%formatted = true
   tri_header_shock%nvertices = 3
   tri_header_shock%combine_zones = true
   call tri_read (lunshock, tri_header_shock, tri_shock, ios)
   if (ios /= 0) go to 99

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

!  Build a search tree for the inner surface and do the body point searching:

   nnode = tri_header_body%nnodes
   ntri  = tri_header_body%nelements

   write (luncrt, '(a, i10)') 'Inner surface node count:    ', nnode
   write (luncrt, '(a, i10)') 'Inner surface triangle count:', ntri

   call build_adt (nnode, ntri, tri_header_body%conn, tri_header_body%xyz)

   ninside  = 0;  dmax  = zero
   noutside = 0;  dmean = zero
   dtolsq = (bbox_diagonal*1.e-6)**2

   do l = 1, nbps

      call search_adt (xyz_bp(:,l), itri, p, q, r, dsqmin, true, nnode, ntri, &
                       tri_header_body%conn, tri_header_body%xyz, xyz_interp)

      if (dsqmin < dtolsq) then ! The nearest triangle was within tolerance
         ninside  = ninside + 1
      else
         noutside = noutside + 1
      end if

      dmax  = max (dmax, dsqmin)
      dmean = dmean + sqrt (dsqmin)
      oncenter = xyz_bp(2,l) == zero ! We ensure a normal in the symmetry plane
      xyz_bp(:,l) = xyz_interp(:)    ! Ensure that the body pt is on the surface

      call tri_normal_and_area (nnode, ntri, tri_header_body%xyz, &
                                tri_header_body%conn, itri, area, &
                                unit_norm(:,l))

      write (luncrt, '(a, 3es19.11)') 'Body normal:', unit_norm(:,l)
      if (oncenter) then
         unit_norm(2,l) = sqrt (unit_norm(1,l)**2 + unit_norm(3,l)**2)
         unit_norm(1,l) = unit_norm(1,l) / unit_norm(2,l)
         unit_norm(3,l) = unit_norm(3,l) / unit_norm(2,l)
         unit_norm(2,l) = zero
         xyz_bp(2,l)    = zero
         write (luncrt, '(a, 3f12.8)') ' Adjusted body point: ', xyz_bp(:,l)
         write (luncrt, '(a, 3f12.8)') ' Adjusted body normal:', unit_norm(:,l)
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

   allocate (xline(nk), yline(nk), zline(nk), sline(nk), stat=ios)
   if (ios /= 0) then
      write (luncrt, '(/, a, i10)') &
         ' Allocation trouble with x/y/z/sline. nk:', nk
      go to 99
   end if

   sl(1)  = zero          ! For all 2-point lines to be intersected
   sl(nl) = bbox_diagonal

!  Build a new search tree from the outer shock boundary:
!  ------------------------------------------------------

   call release_adt ()              ! Deallocate work-space from first ADT
   nnode = tri_header_shock%nnodes  ! One or more zones are handled now
   ntri  = tri_header_shock%nelements

   write (luncrt, '(a, i10)') 'Outer boundary node count:    ', nnode
   write (luncrt, '(a, i10)') 'Outer boundary triangle count:', ntri

   call build_adt (nnode, ntri, tri_header_shock%conn, tri_header_shock%xyz)

!  Calculate each line of sight of the indicated type:
!  ---------------------------------------------------

   if (.not. hemisphere) then  ! Body-normal line(s) only -- do it in-line

      ninside  = 0;  dmax  = zero
      noutside = 0;  dmean = zero

      do l = 1, nbps

         xl(1)  = xyz_bp(1,l);  xl(2) = xl(1) + bbox_diagonal*unit_norm(1,l)
         yl(1)  = xyz_bp(2,l);  yl(2) = yl(1) + bbox_diagonal*unit_norm(2,l)
         zl(1)  = xyz_bp(3,l);  zl(2) = zl(1) + bbox_diagonal*unit_norm(3,l)
         sl(nl) = bbox_diagonal ! Upper limit on intersection search along line;
                                ! overkill for forebody; safe for aft body

         call intsec8 (nnode, ntri, tri_header_shock%conn, &
                       tri_header_shock%xyz, nl, xl, yl, zl, sl, lcs_method, &
                       itri, p, q, r, s, xyz_intersect, dsqmin)

!!!      write (luncrt, '(a,4i10)') ' itri,conn(:,itri):', itri, &
!!!                                 tri_header_shock%conn(:,itri)
!!!      write (luncrt, '(a,4es19.11)') 'p, q, r, s:', p, q, r, s
         write (luncrt, '(a,4es19.11)') 'xyzint, dsqmin:', xyz_intersect, dsqmin

         if (dsqmin < dtolsq) then ! The nearest triangle was within tolerance
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
            radial_lines(l)%x(1,1,k) = xl(1) + ((xl(2) - xl(1)))*sline(k)
            radial_lines(l)%y(1,1,k) = yl(1) + ((yl(2) - yl(1)))*sline(k)
            radial_lines(l)%z(1,1,k) = zl(1) + ((zl(2) - zl(1)))*sline(k)
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

!  Internal procedures for program USLOS:
!  --------------------------------------

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

      do in = 1, tri_header_shock%nnodes
         xmax = max (xmax, tri_header_shock%xyz(1,in))
         xmin = min (xmin, tri_header_shock%xyz(1,in))
         ymax = max (ymax, tri_header_shock%xyz(2,in))
         ymin = min (ymin, tri_header_shock%xyz(2,in))
         zmax = max (zmax, tri_header_shock%xyz(3,in))
         zmin = min (zmin, tri_header_shock%xyz(3,in))
      end do

      bbox_diagonal = sqrt ((xmax-xmin)**2 + (ymax-ymin)**2 + (zmax-zmin)**2)

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
      real :: &
         length
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
!     (b) the surface triangulation produced by intersecting lines throught each
!     vertex of (a) with outer boundary mesh.

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
      length   = bbox_diagonal  ! Shorter name

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

            xl(2) = xl(1) + (transformed_quadrant(1)%xyz(1,nn) - xl(1))*length
            yl(2) = yl(1) + (transformed_quadrant(1)%xyz(2,nn) - yl(1))*length
            zl(2) = zl(1) + (transformed_quadrant(1)%xyz(3,nn) - zl(1))*length

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

      real, parameter :: small_value = 0.001

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

!     Execution:

      call intsec8 (nnode, ntri, tri_header_shock%conn, &
                    tri_header_shock%xyz, nl, xl, yl, zl, sl, lcs_method, &
                    itri, p, q, r, s, xyz_intersect, dsqmin)

!!!   write (luncrt, '(a,3i5)') ' Quadrant, node, line:', m, nn, l
!!!   write (luncrt, '(a,4i10)') ' itri,conn(:,itri):', itri, &
!!!                              tri_header_shock%conn(:,itri)
!!!   write (luncrt, '(a,4es19.11)') 'p, q, r, s:', p, q, r, s
      write (luncrt, '(a,4es19.11)') 'xyzint, dsqmin:', xyz_intersect, dsqmin

      if (dsqmin < dtolsq) then ! The nearest triangle was within tolerance
         ninside  = ninside + 1
      else
         noutside = noutside + 1
         write (luncrt, '(a, i6, 4es13.6)') &
         ' *** Intersection trouble.  Line #, sa, sb, dsq, s:', l, sl, dsqmin, s

         write (luncrt, '(a)') &
         '     Assume local min. due to boundary concavity.  Do one retry.'
         sl(1) = s*1.000001

         call intsec8 (nnode, ntri, tri_header_shock%conn, &
                       tri_header_shock%xyz, nl, xl, yl, zl, sl, lcs_method, &
                       itri, p, q, r, s, xyz_intersect, dsqmin)

!!!      write (luncrt, '(a,4i10)') ' itri,conn(:,itri):', itri, &
!!!                                 tri_header_shock%conn(:,itri)
!!!      write (luncrt, '(a,4es19.11)') 'p, q, r, s:', p, q, r, s

         if (dsqmin < dtolsq) then ! The nearest triangle was within tolerance
            ninside  = ninside + 1;  noutside = noutside - 1
            write (luncrt, '(a, i6, 4es13.6)') &
               '     Recovered.       Line #, sa, sb, dsq, s:', l, sl, dsqmin, s
         else
            write (luncrt, '(a, i6, 4es13.6)') &
         ' *** Intersection failure.  Line #, sa, sb, dsq, s:', l, sl, dsqmin, s

            ! Keep going
         end if
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
         radial_lines(l)%x(1,1,k) = xl(1) + ((xl(2) - xl(1)))*sline(k)
         radial_lines(l)%y(1,1,k) = yl(1) + ((yl(2) - yl(1)))*sline(k)
         radial_lines(l)%z(1,1,k) = zl(1) + ((zl(2) - zl(1)))*sline(k)
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

   end program uslos
