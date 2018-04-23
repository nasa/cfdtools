!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program surface_project
!
!  Description:
!
!  This is an adaptation of SURFACE_INTERP for the case where the target surface
!  points need to be projected onto a somewhat different surface grid.  Only the
!  (x,y,z) coordinates are treated - no function data are involved.  Further, it
!  has the option to retain a target point rather than the projected point if it
!  finds that the interpolation coefficients imply extrapolation - i.e., if the
!  grid being projected to does not encompass some of the target points.
!
!  This functionality was prompted by the need to perturb a smooth surface grid
!  to match the (noisy) as-built coordinates of a mesh of heating sensors on a
!  wind tunnel model.  The as-built shape can be blended into the grid for the
!  ideal shape this way, and the flow solution can be reconverged in the hope of
!  matching the wind tunnel data better.  (This was on the wing leading edge,
!  where the flow is expected to be extremely sensitive to curvature anomalies.)
!
!  The program reads a target Tecplot dataset of one or more surface patches
!  (ni x 1 x 1 if unstructured), followed by a structured surface dataset (also
!  Tecplot ASCII) onto which the target points will be projected.  Any function
!  data present are ignored.
!
!  Actually, if the first line of the target point file is purely numeric, it is
!  read as a list of (x,y,z) coordinates, one point per line.
!
!  Optionally, points outside that surface can be left unchanged.  Otherwise,
!  they will be output on the boundary of that surface.
!
!  Results are output as a Tecplot surface dataset in BLOCK order (for ease of
!  adapting as a PLOT3D file).  A single output function variable contains the
!  projection distance for each target point, except this is set to zero (for
!  plotting convenience) if extrapolation is implied and the option to retain
!  input coordinates of target points outside the desired surface is specified.
!
!  Control file format (standard input):
!
!     SURFACE_PROJECT control file
!     target_points.dat     ! Target pts.: x,y,z list | Tecplot ASCII; 1+ zones
!     WLE-sensor-mesh.dat   ! Structured surface to be projected to; Tec. ASCII
!     0        ! For target pts. outside mesh:  0 = untouched; 1 = best edge pt.
!     projected_points.dat  ! Plottable results: Tecplot ASCII, BLOCK order
!     2 1 1                 ! Target block, i, j for debug output; 0 suppresses
!
!  History:
!
!     02/15/09  David Saunders  Initial adaptation from SURFACE_INTERP, for the
!                               application outlines above.
!     02/22/09    "      "      If the target points file starts with a purely
!                               numeric line, it is interpreted as an (x,y,z)
!                               list, one point per line.
!     08/08/13    "      "      All ADT variants are in one module now with
!               Now ERC, Inc.   generic build_adt and search_adt interfaces.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_header_structure  ! See Tecplot_io.f90 for these three modules
   use grid_block_structure
   use tecplot_io_module
   use adt_utilities          ! All ADT variants

   implicit none

!  Constants:

   integer, parameter :: &
      luntarget  = 1,    &  ! Input target surface points (x,y,z list | Tecplot)
      lunsurface = 2,    &  ! Input surface grid to project to   (Tecplot ASCII)
      lunresults = 3,    &  ! Output plottable results           (Tecplot ASCII)
      lunctl     = 5,    &  ! Control file (standard input)
      luncrt     = 6,    &  ! For diagnostics
      ndim       = 3

   real, parameter ::    &
      one  = 1.0,        &
      zero = 0.0

   logical, parameter :: &
      false = .false.,   &
      true  = .true.

!  Variables:

   integer :: &
      extrapolation_type, ib, idebug, ios, jdebug, m, nb_surface, nb_target,   &
      ndebug, nf_out, nf_surface, npts, nquad

   integer, allocatable, dimension (:,:) :: &
      conn               ! For patch & (i,j) of surface quads. being searched

   logical :: &
      list_of_targets

   character :: &
      buffer * 132       ! For checking the first line of the target point file

!  Composite data types:

   type (grid_header) :: &
      header_target, header_surface, header_results

   type (grid_type), pointer, dimension (:) :: &
      target, surface, results      ! (x,y,z,f) for the input surfaces & results

!  Procedures:

   external :: alpha
   logical  :: alpha

!  Execution:
!  !!!!!!!!!!

   read (lunctl, *)
   read (lunctl, *) header_target%filename

   header_target%formatted = true
   header_target%ndim      = ndim

!  Check for just an (x,y,z) list instead of target surface patch(es):

   open (luntarget, file=header_target%filename, status='old', iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') &
         ' Trouble opening ', trim (header_target%filename)
      go to 99
   end if

   read (luntarget, '(a)') buffer

   write (luncrt, '(a)') trim (buffer)

   list_of_targets = .not. alpha (trim (buffer))

   write (luncrt, '(a, l1)') 'list_of_targets: ', list_of_targets

   if (list_of_targets) then

      rewind (luntarget)

      call read_list ()  ! Internal procedure below counts/assigns npts points

   else  ! Tecplot ASCII assumed

      close (luntarget)

      call Tecplot_read (luntarget, header_target, target, ios)

   end if

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Trouble reading target points file ',        &
         trim (header_target%filename)
      go to 99
   end if

   nb_target = header_target%nblocks

   read (lunctl, *) header_surface%filename ! Tecplot surface grid to project to

   header_surface%formatted = true
   header_surface%ndim      = ndim
   ios                      = 1    ! Verbose input mode

   call Tecplot_read (lunsurface, header_surface, surface, ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Trouble reading surface grid Tecplot file ', &
         trim (header_surface%filename)
      go to 99
   end if

   nb_surface = header_surface%nblocks
   nf_surface = header_surface%numq

   if (nf_surface > 0) then
      write (luncrt, '(a, i3)') ' # functions in surface solution:',           &
      nf_surface, ' These will NOT be interpolated.  SURFACE_INTERP does that.'
   end if

!  Read the choice for what to do with target points outside this surface:

   read (lunctl, *) extrapolation_type

   if (extrapolation_type < 0 .or. extrapolation_type > 1) then
      write (luncrt, '(/, a, i5)') ' Bad extrapolation control:', &
         extrapolation_type, ' 0 = preserve input (x,y,z); 1 = move to edge.'
      go to 99
   end if

!  Output file specs:

   read (lunctl, *) header_results%filename
   header_results%datapacking = 1  ! BLOCK order
   header_results%numq        = 1  ! Projection distance
   nf_out                     = 1

!  Debug controls:

   read (lunctl, *) ndebug, idebug, jdebug ! ndebug <= 0 suppresses debug output

   close (lunctl)

!  Allocate interpolated results as though they are structured (may be ni x 1):

   allocate (results(nb_target))

   do ib = 1, nb_target

      results(ib)%ni = target(ib)%ni
      results(ib)%nj = target(ib)%nj
      results(ib)%nk = target(ib)%nk
      results(ib)%mi = target(ib)%ni
      results(ib)%mj = target(ib)%nj
      results(ib)%mk = target(ib)%nk
      results(ib)%zone_title = target(ib)%zone_title
      results(ib)%nzoneaux   = target(ib)%nzoneaux

      call Tec_block_allocate (results(ib), ndim, nf_out, ios)

      if (results(ib)%nzoneaux > 0) then
          results(ib)%zoneauxname(:)  = target(ib)%zoneauxname(:)
          results(ib)%zoneauxvalue(:) = target(ib)%zoneauxvalue(:)
      end if

   end do

!  Construct a search tree from all patches of the surface data:

   nquad = 0
   do ib = 1, nb_surface
      nquad = (surface(ib)%ni - 1) * (surface(ib)%nj - 1) + nquad
   end do

   allocate (conn(3,nquad)) ! For patch # and (i,j)

   call build_adt (nb_surface, surface, nquad, conn)


!  Perform the interpolations:

   call interpolate_surface ()


!  Save results:

   allocate (header_results%varname(ndim + nf_out))

   header_results%varname(1:ndim) = header_surface%varname(1:ndim)
   header_results%varname(ndim+1:ndim+1) = 'Distance projected'

!  Save results in BLOCK order, formatted:

   header_results%formatted   = true
   header_results%ndim        = ndim
   header_results%numq        = nf_out
   header_results%nblocks     = header_target%nblocks
   header_results%title       = 'Projected surface values'
   m                          = header_target%ndatasetaux
   header_results%ndatasetaux = m

   if (m > 0) then
      allocate (header_results%datasetauxname(m))
      header_results%datasetauxname(:)  = header_target%datasetauxname(:)
      allocate (header_results%datasetauxvalue(m))
      header_results%datasetauxvalue(:) = header_target%datasetauxvalue(:)
   end if

   call Tecplot_write (lunresults, header_results, results, ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') &
         ' Trouble writing the projected results in Tecplot form.'
   end if

   close (lunresults)

99 continue  ! Avoid system-dependent STOP issues; let F90 deallocate everything


!  Local procedures for program SURFACE_PROJECT
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   contains

!     --------------------------------------------------------------------------
      subroutine interpolate_surface ()

!     Locate the target surface points in the surface grid being projected to.
!     If requested, trap points that appear to be outside that grid, and just
!     output the target coordinates untouched.
!     --------------------------------------------------------------------------

!     Local variables:

      integer :: i, ib, ib1, ic, ier, iquad, j, jc, m, n, ni, ninside, nj,     &
                 noutside
      real    :: dmax, dmean, dmin, dsqmin, dtolsq, p, pm1, q, qm1,            &
                 interp_xyz(3), target_xyz(3)
      logical :: trap, use_projection

!     Execution:

!     Tolerance for search diagnostics (refine later):

      dtolsq = (0.001) ** 2
      trap   = extrapolation_type == 0

      do ib = 1, nb_target

         ninside  = 0;  dmax  = zero
         noutside = 0;  dmean = zero
         npts     = 0

         ni = target(ib)%ni
         nj = target(ib)%nj

         do j = 1, nj

            do i = 1, ni

               target_xyz(1) = target(ib)%x(i,j,1)
               target_xyz(2) = target(ib)%y(i,j,1)
               target_xyz(3) = target(ib)%z(i,j,1)

               call search_adt (target_xyz, iquad, p, q, dsqmin, true,         &
                                nb_surface, surface, nquad, conn, interp_xyz)

               use_projection = dsqmin < dtolsq  ! Always, presumably

               if (trap) then  ! Above tol may still be inside searched surface
                  if (.not. use_projection) then
                     if (p /= zero .and. p /= one .and. &  ! Exactness unlikely
                         q /= zero .and. q /= one) then
                        use_projection = true
                     else  ! Most likely outside (NOT FAIL-SAFE THOUGH!)
                        noutside = noutside + 1  ! Outside searched surface
                     end if
                  end if
               else
                  if (use_projection) then
                     ninside = ninside + 1
                  else
                     noutside = noutside + 1  ! Outside tolerance
                  end if
                  use_projection = true  ! Always
               end if

               if (use_projection) then
                  results(ib)%x(i,j,1)   = interp_xyz(1)
                  results(ib)%y(i,j,1)   = interp_xyz(2)
                  results(ib)%z(i,j,1)   = interp_xyz(3)
                  results(ib)%q(1,i,j,1) = sqrt (dsqmin)
                  dmax  = max (dmax, dsqmin)
                  dmean = dmean + dsqmin
               else
                  noutside = noutside + 1  ! Don't keep their statistics
                  results(ib)%x(i,j,1)   = target_xyz(1)
                  results(ib)%y(i,j,1)   = target_xyz(2)
                  results(ib)%z(i,j,1)   = target_xyz(3)
                  results(ib)%q(1,i,j,1) = zero  ! For plotting purposes
               end if

               npts = npts + 1

!              Debug output:

               if (ib == ndebug) then
                  if (i == idebug .and. j == jdebug) then
                     write (6, '(a, 3f20.8)') ' target xyz: ', target_xyz
                     write (6, '(a, 3f20.8)') ' interp xyz: ', interp_xyz
                     n   = conn(1,iquad) ! Block #
                     ic  = conn(2,iquad) ! Lower left quad. indices
                     jc  = conn(3,iquad)
                     write (6, '(a, 3i5)') &
                        ' Search result: block, ic, jc: ', n, ic, jc
                     write (6, '(a, 2f20.10)') ' p, q: ', p, q,     &
                                               ' dsqmin, dtolsq: ', &
                                                 dsqmin, dtolsq
                     dmin = sqrt (dsqmin)
                     write (6, '(a, f20.10)') ' shortest distance: ', dmin
                     write (6, '(a)') ' Cell vertex coordinates:'
                     write (6, '(4f20.8)') &
                       surface(n)%x(ic:ic+1,jc:jc+1,1), &
                       surface(n)%y(ic:ic+1,jc:jc+1,1), &
                       surface(n)%z(ic:ic+1,jc:jc+1,1)
                  end if
               end if

            end do ! Next i for this target block

         end do ! Next j for this target block

         write (luncrt, '(a, i5, a, 2i7, a, 1p, 2e12.5)')                      &
            '  Target patch', ib,                                              &
            ':  # points in/outside tolerance:', ninside, noutside,            &
            ';  max/mean devn.:', sqrt (dmax), sqrt (dmean / real (npts))

      end do ! Next target block

      end subroutine interpolate_surface

!     --------------------------------------------------------------------------
      subroutine read_list ()

!     Count the lines in an opened file, allocate (x,y,z) storage, read the list
!     of coordinates (one per line) and close the file.
!     --------------------------------------------------------------------------

      integer :: i

      npts = -1
      do while (ios == 0) ! Until EOF (ios < 0)
         npts = npts + 1
         read (luntarget, *, iostat=ios)
      end do

      write  (luncrt, '(/, a, i7)') '# target points found:', npts

      rewind (luntarget)

      header_target%numq        = 0
      header_target%nblocks     = 1
      header_target%datapacking = 0  ! POINT order
      header_target%title       = 'Target points'
      header_target%ndatasetaux = 0
      allocate (header_target%varname(3))
      header_target%varname(1) = 'x'
      header_target%varname(2) = 'y'
      header_target%varname(3) = 'z'

      allocate (target(1))
      allocate (target(1)%x(npts,1,1), &
                target(1)%y(npts,1,1), &
                target(1)%z(npts,1,1))
      target(1)%ni = npts
      target(1)%nj = 1
      target(1)%nk = 1
      target(1)%mi = npts
      target(1)%mj = 1
      target(1)%mk = 1
      target(1)%nzoneaux   = 0
      target(1)%zone_title = header_target%title

      do i = 1, npts
         read (luntarget, *, iostat=ios) &
            target(1)%x(i,1,1), target(1)%y(i,1,1), target(1)%z(i,1,1)
         if (ios /= 0) exit
      end do

      close (luntarget)

      end subroutine read_list

   end program surface_project
