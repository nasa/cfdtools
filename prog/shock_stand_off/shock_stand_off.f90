!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program shock_stand_off

!  Description:
!
!     For an entire CFD grid and function file containing just the function used
!     to detect shock location (normally Temperature), calculate shock stand-off
!     distance at every inner surface point.  A single function is stipulated
!     because the normals off the surface may cross block boundaries, meaning
!     all the volume data must be stored.  This precludes use of the volume data
!     processed by the boundary layer edge detection utility, results from which
!     the stand-off distances are to be merged with in a subsequent step that
!     requires the stand-off results to be saved in Tecplot surface file form.
!
!     This version handles 2-D cases as well.  It is most expedient to turn a
!     single plane of data into two parallel planes then apply the existing
!     3-D technique.
!
!  Initial assumptions (probably generalizable, but they may never need to be):
!
!     o  The structured volume grid contains one layer of blocks, with k = 1
!        at the wall.  This simplifies determination of the inner and outer
!        boundary patches.  (To overcome these restrictions, one could use the
!        boundary condition data employed by the relevant flow solver.)
!
!     o  For the 2-D case, assume j = 1 at the wall initially.  The program
!        converts 2-D inputs to 3-D form with k = 1 at the wall.
!
!  Strategy:
!
!     o  Read the entire volume data and extract the inner and outer boundaries
!        as multiblock surface grids.
!
!        o  In order to handle 2-D cases transparently, we need to parse the
!           grid file header, so 2-D files must be formatted.
!
!        o  For a 2-D case, the 3rd coordinates are assumed to be missing, even
!           though this is incompatible with the I/O package used here.
!
!        o  For a 2-D case, replicate the given plane of data with a parallel
!           plane close apart in the third dimension, then proceed as for 3-D.
!
!     o  Build a search tree from the outer boundary.
!
!     o  For each inner surface point:
!
!        > Construct a two-point line normal to the wall with length that of
!          a local radial line.  This should be at least as long as the straight
!          line distance to the outer boundary.
!
!        > Intersect the line with the outer boundary and store the new end pt.
!
!     o  Release the outer surface search tree & build one for the volume grid.
!
!     o  For each two-point normal, if the surface x exceeds a specified limit,
!        set the stand-off distance to zero (assumed to be on the aft body),
!        else impose a radial-grid-line-like distribution on it and for k = nk,
!        1, -1 interpolate the flow quantity until a spike is apparent, making
!        further (relatively costly) interpolations redundant.  Linearly
!        interpolate to find the arc length where the flow quantity is (say)
!        1.1 x free-stream to define the stand-off distance.
!
!     o  Save results in Tecplot form for possible merging with BLAYER results
!        via the MERGE_FILES utility.
!
!  History:
!
!     04/20/07  D.A.Saunders  Initial adaptation of LINES_OF_SIGHT and the other
!                             single-function tools originally strung together
!                             to provide data for QRAD.  They are better suited
!                             to a thinned set of results.
!     05/08/07    "     "     Another XYZ convention forced a new control input.
!     08/28/09    "     "     Option to handle 2-D cases added transparently.
!     08/08/13    "     "     All ADT variants have been merged into a module
!               Now ERC, Inc. with generic build_adt and search_adt interfaces.
!     08/28/17  Now AMA, Inc. Cases from John Theisinger at NASA LaRC were so
!                             tightly aligned that uniform spacing along the
!                             discretized body-normal profiles was a bad choice.
!                             Worse: at low Mach, with little or no chemistry,
!                             there may not be a temperature spike--just a
!                             rise/no obvious drop, so seek the shock onset
!                             (rise) rather than a drop that may not be there.
!                             This precludes using Mach as the function.
!                             Mach, in retrospect, with a different search,
!                             would have been preferable.
!     08/30/17    "     "     Followed John's suggestion of using 1.1xTinf
!                             to interpolate a stand-off distance smoothly.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure  ! Derived data type for one grid block
   use adt_utilities         ! All ADT variants
   use xyzq_io_module        ! PLOT3D-type I/O package

   implicit none

!  Local constants:

   integer, parameter :: &
      lunvolg   = 1,     &   ! Input volume grid ...
      lunvolf   = 2,     &   ! ... and function file
      lunout    = 3,     &   ! Output file (Tecplot surface with standoffs)
      lunkbd    = 5,     &   ! Keyboard inputs
      luncrt    = 6,     &   ! Screen
      lunshockg = 7,     &   ! For plotting shock envelope intersection data
      lunshockf = 8,     &
      nl        = 2          ! 2-point lines are passed to INTSEC7

   real, parameter ::    &
      fs_scale  = 1.1,   &   ! For Tinf; cf. 0.9 x Minf scheme in DPLR
      half      = 0.5,   &
      one       = 1.0,   &
      zero      = 0.0

   logical, parameter :: &
      false = .false.,   &
      true  = .true.

   character, parameter :: &
      lcs_method * 1 = 'L'    ! But LCSFIT isn't needed (hook for general case)

!  Local variables:

   integer :: &
      i, ib, ic, ier, ihex, ios, iquad, j, jc, k, kc, n, nblocks, nf, nhex,    &
      ni, nj, nk, ninside, noutside, npts, nquad

   integer, allocatable :: &
      conn(:,:)               ! For (patch,i,j) of boundary quads. to search

   real :: &
      cutoff, dmax, dmean, dsqmin, dt, dt1, dt2, dtolsq, dx, dy, dz,           &
      freestream, p, pm1, q, qm1, r, rm1, rnk, streamwise_coord, t, total,     &
      ymax, ymin

   real, dimension (nl) :: &
      tl, xl, yl, zl          ! For 2-point lines passed to INTSEC7

   real, dimension (3) ::  &
      un, vec1, vec2,      &  ! For unit normals
      xyz_interp, xyz_target

   real, allocatable, dimension (:) :: &
      fline, tline, xline, yline, zline  ! For one discretized line of sight

   logical :: &
      cell_centered, formatted, formatted_output, not_2D, show, unused

   character :: &
      streamwise_axis * 1

   type (grid_type), pointer, dimension (:) :: &
      inner_surface, intersections, outer_surface, volume_grid

!  Execution:

!  Deal with the volume data (and initializing the output results file):
!  ---------------------------------------------------------------------

   call file_prompt (lunvolg, 0, 'volume grid', 'old', true, formatted, ios)
   if (ios /= 0) go to 99

   call file_prompt (lunvolf, 0, 'volume function', 'old', false, formatted,   &
                     ios)
   if (ios /= 0) go to 99

   formatted_output = true
   call file_prompt (lunout, 0, 'output stand-off results', 'unknown', false,  &
                     formatted_output, ios)

   write (luncrt, '(a)', advance='no') ' Streamwise axis [x|y|z]: '
   read  (lunkbd, *) streamwise_axis

   write (luncrt, '(3a)', advance='no') &
      ' Cut-off ', streamwise_axis, ' defining aft body: '
   read  (lunkbd, *) cutoff

!  Check for handling a 2-D case:
!  ------------------------------

   call make_3D_if_needed ()  ! Internal procedure below

   if (ios /= 0) go to 99

   if (not_2D) then

!     Read the entire volume data:

      call xyzq_read (lunvolg, lunvolf, formatted, nblocks, nf, cell_centered, &
                      volume_grid, ios)
      if (ios /= 0) go to 99

   else
      ! The datasets are in 3-D form now
   end if

!  Extract the inner and outer boundaries as multiblock surface grids:
!  -------------------------------------------------------------------

   allocate (inner_surface(nblocks), outer_surface(nblocks),                   &
             intersections(nblocks))

   npts = 0;  nquad = 0;  nhex = 0
   ymin = 1.e+30;  ymax = -ymin
   dt2  = 1.e+30               ! Need smallest outermost radial spacing

   do ib = 1, nblocks

      ni = volume_grid(ib)%ni
      nj = volume_grid(ib)%nj
      nk = volume_grid(ib)%nk
      npts  = ni * nj + npts
      nquad = (ni - 1) * (nj - 1) + nquad
      nhex  = (ni - 1) * (nj - 1) * (nk - 1) + nhex
      inner_surface(ib)%ni = ni
      outer_surface(ib)%ni = ni
      intersections(ib)%ni = ni
      inner_surface(ib)%nj = nj
      outer_surface(ib)%nj = nj
      intersections(ib)%nj = nj
      inner_surface(ib)%nk = 1
      outer_surface(ib)%nk = 1
      intersections(ib)%nk = 1

      call xyz_allocate (inner_surface(ib), ios)

      inner_surface(ib)%x(:,:,1) = volume_grid(ib)%x(:,:,1)
      inner_surface(ib)%y(:,:,1) = volume_grid(ib)%y(:,:,1)
      inner_surface(ib)%z(:,:,1) = volume_grid(ib)%z(:,:,1)

      call xyz_allocate (outer_surface(ib), ios)

      nk = volume_grid(ib)%nk
      outer_surface(ib)%x(:,:,1) = volume_grid(ib)%x(:,:,nk)
      outer_surface(ib)%y(:,:,1) = volume_grid(ib)%y(:,:,nk)
      outer_surface(ib)%z(:,:,1) = volume_grid(ib)%z(:,:,nk)

      do j = 1, nj
         do i = 1, ni
            ymin = min (ymin, outer_surface(ib)%y(i,j,1))
            ymax = max (ymax, outer_surface(ib)%y(i,j,1))
            dx   = volume_grid(ib)%x(i,j,nk) - volume_grid(ib)%x(i,j,nk-1)
            dy   = volume_grid(ib)%y(i,j,nk) - volume_grid(ib)%y(i,j,nk-1)
            dz   = volume_grid(ib)%z(i,j,nk) - volume_grid(ib)%z(i,j,nk-1)
            dt   = sqrt (dx**2 + dy**2 + dz**2)
            dt2  = min (dt2, dt)  ! For discretizing all 2-pt. body-normal lines
         end do
      end do

   end do  ! Next block

   dtolsq = ((ymax - ymin) * 0.00001) ** 2  ! Tolerance for search diagnostics
   dx   = volume_grid(1)%x(i,j,2) - volume_grid(1)%x(i,j,1)  ! Wall spacing is
   dy   = volume_grid(1)%y(i,j,2) - volume_grid(1)%y(i,j,1)  ! not critical --
   dz   = volume_grid(1)%z(i,j,2) - volume_grid(1)%z(i,j,1)  ! block 1 is OK
   dt1  = sqrt (dx**2 + dy**2 + dz**2)
   dt2  = dt2 + dt2  ! 2x is safe; avoid excessive smallness for many lines

   allocate (conn(3,nquad)) ! For (patch,i,j) of each boundary quad.

!  Build a search tree from the outer boundary patches:
!  ----------------------------------------------------

   call build_adt (nblocks, outer_surface, nquad, conn)

   ninside  = 0;  dmax  = zero
   noutside = 0;  dmean = zero

!  For each inner surface point ...
!  --------------------------------

   do ib = 1, nblocks

      ni = volume_grid(ib)%ni;  intersections(ib)%mi = ni; ! For stand-off dist.
      nj = volume_grid(ib)%nj;  intersections(ib)%mj = nj
      nk = volume_grid(ib)%nk;  intersections(ib)%mk = 1;  rnk = real (nk-1)

      allocate (xline(nk), yline(nk), zline(nk), tline(nk))

      call xyz_allocate (intersections(ib), ios)
      call   q_allocate (intersections(ib), nf, ios)  ! For intersection arcs

      do j = 1, nj

         do i = 1, ni

!!          show = i == 60 .and. j == nj

!           Derive a two-point line length from the volume grid radial line:

            do k = 1, nk
               xline(k) = volume_grid(ib)%x(i,j,k)
               yline(k) = volume_grid(ib)%y(i,j,k)
               zline(k) = volume_grid(ib)%z(i,j,k)
            end do

!           Unnormalized arc lengths and total length:

            call chords3d (nk, xline, yline, zline, false, total, tline)

!           Carefully calculated unit normal:

!!!         if (i < ni) then  ! (ic,jc) must point to "lower left"
!!!            ic = i         ! No need:  this is done inside the
!!!            p = zero       ! surface_normal utility now.
!!!         else
!!!            ic = ni - 1
!!!            p = one
!!!         end if

!!!         if (j < nj) then
!!!            jc = j
!!!            q = zero
!!!         else
!!!            jc = nj - 1
!!!            q = one
!!!         end if

            call surface_normal (ni, nj, inner_surface(ib)%x,                  &
                                 inner_surface(ib)%y, inner_surface(ib)%z,     &
                                 i, j, zero, zero, un)

!           Two-point line off the wall:

            xl(1) = xline(1);  xl(2) = xl(1) + total * un(1)
            yl(1) = yline(1);  yl(2) = yl(1) + total * un(2)
            zl(1) = zline(1);  zl(2) = zl(1) + total * un(3)

!!          if (show) then
!!             write (luncrt, '(a, 3f12.6, a, 3f11.7)') &
!!                ' xl,yl,zl:', xl(1), yl(1), zl(1), '  un(:): ', un(:), &
!!                '          ', xl(2), yl(2), zl(2), '  total: ', total
!!          end if

!           Intersect the two-point line with the outer boundary:

            tl(1)  = 0.001       ! Tell INTSEC6 the normalized t search range
            tl(nl) = 10.         ! Allow for extrapolation

            call intsec6 (nblocks, outer_surface, nquad, conn, nl, xl, yl, zl, &
                          tl, lcs_method, iquad, p, q, t, xyz_interp, dsqmin)

!!          if (show) then
!!    write (luncrt, '(a, 4i10)') ' iquad, conn(:,iquad):', iquad, conn(:,iquad)
!!    write (luncrt, '(a, 1p, 4e19.11)') ' p, q, t, total:', p, q, t, total
!!    write (luncrt, '(a, 1p, 4e19.11)') ' xyzint, dsqmin:', xyz_interp, dsqmin
!!          end if

            if (t > 0.99 * tl(nl)) write (luncrt, '(a, 3i5, 3f11.7)') &
               ' Intersection warning. ib, i, j:', ib, i, j, t, tl(nl), total

            if (dsqmin < dtolsq) then ! The nearest quad was within tolerance
               ninside  = ninside + 1
            else
               noutside = noutside + 1
            end if

            dmax  = max (dmax, dsqmin)
            dmean = dmean + sqrt (dsqmin)

            intersections(ib)%x(i,j,1)   = xyz_interp(1)
            intersections(ib)%y(i,j,1)   = xyz_interp(2)
            intersections(ib)%z(i,j,1)   = xyz_interp(3)
            intersections(ib)%q(1,i,j,1) = t  ! For plotting purposes

         end do  ! Next i

      end do  ! Next j

      deallocate (xline, yline, zline, tline)

   end do  ! Next block


   dmax = sqrt (dmax);  dmean = dmean / real (npts)

   write (luncrt, '(/, a, 2i6, /, a, 1p, 2e12.5)') &
      ' # outer intersection points inside/outside tolerance:',                &
      ninside, noutside, ' max & mean distance:', dmax, dmean

   open (lunshockg, file='shock_envelope.g', status='unknown')
   open (lunshockf, file='shock_envelope.f', status='unknown')

   call xyz_write (lunshockg, true, nblocks,     intersections, ios)
   call   q_write (lunshockf, true, nblocks, nf, intersections, ios)

!  Release the surface search tree data in readiness for volume interpolations:

   call release_adt ()

   deallocate (conn)


!  Build the volume grid search tree:
!  ----------------------------------

   allocate (conn(4,nhex))

   call build_adt (nblocks, volume_grid, nhex, conn, unused)

!  For each surface normal, discretize it like a radial grid line,
!  interpolate the function data, and look for the shock jump.
!  -----------------------------------------------------------

   do ib = 1, nblocks

      ni = volume_grid(ib)%ni
      nj = volume_grid(ib)%nj
      nk = volume_grid(ib)%nk

      allocate (xline(nk), yline(nk), zline(nk), tline(nk), fline(nk))

      rnk = real (nk-1)

      do j = 1, nj

         do i = 1, ni

!!          show = i == 60 .and. j == nj

            xline(1)  = inner_surface(ib)%x(i,j,1)
            yline(1)  = inner_surface(ib)%y(i,j,1)
            zline(1)  = inner_surface(ib)%z(i,j,1)

!           Avoid trying to find a shock in the wake region:

            select case (streamwise_axis)

               case ('x', 'X')

                  streamwise_coord = xline(1)

               case ('y', 'Y')

                  streamwise_coord = yline(1)

               case default

                  streamwise_coord = zline(1)

            end select

            if (streamwise_coord > cutoff) then
               intersections(ib)%q(1,i,j,1) = zero
               cycle
            end if

            freestream = volume_grid(ib)%q(1,i,j,nk) ! Outermost value
            xline(nk) = intersections(ib)%x(i,j,1)
            yline(nk) = intersections(ib)%y(i,j,1)
            zline(nk) = intersections(ib)%z(i,j,1)

            dx = xline(nk) - xline(1)
            dy = yline(nk) - yline(1)
            dz = zline(nk) - zline(1)
            dt = sqrt (dx**2 + dy**2 + dz**2)  ! Total arc length

!           Two-sided stretching that won't miss the shock details:

            call htdis4 (true, zero, dt, dt1, dt2, nk, tline, -luncrt, ier)
            if (ier /= 0) then
               write (luncrt, '(a, 3i5)') 'HTDIS4 trouble at ib,i,j:', ib, i, j
               stop
            end if

            do k = 1, nk
               r = tline(k) / dt
               xline(k) = xline(1) + r*dx
               yline(k) = yline(1) + r*dy
               zline(k) = zline(1) + r*dz
            end do

!!          if (show) &
!!             write (luncrt, '(4f12.6)') &
!!                (xline(k), yline(k), zline(k), tline(k), k = 1, nk)

!           Search from the outer boundary for the jump at the shock:

            fline(nk) = zero  ! To pass the shock test

            do k = nk, 1, -1

               xyz_target(1) = xline(k)
               xyz_target(2) = yline(k)
               xyz_target(3) = zline(k)

               call search_adt (xyz_target, ihex, p, q, r, dsqmin, true, &
                                nblocks, volume_grid, nhex, conn, xyz_interp)

               n  = conn(1,ihex) ! Block #
               ic = conn(2,ihex) ! Lower left cell indices
               jc = conn(3,ihex)
               kc = conn(4,ihex)

               pm1 = one - p
               qm1 = one - q
               rm1 = one - r

               fline(k) = &
                  rm1 * (qm1 * (pm1 * volume_grid(n)%q(1,ic,  jc,  kc  )   +   &
                                  p * volume_grid(n)%q(1,ic+1,jc,  kc  ))  +   &
                           q * (pm1 * volume_grid(n)%q(1,ic,  jc+1,kc  )   +   &
                                  p * volume_grid(n)%q(1,ic+1,jc+1,kc  ))) +   &
                    r * (qm1 * (pm1 * volume_grid(n)%q(1,ic,  jc,  kc+1)   +   &
                                  p * volume_grid(n)%q(1,ic+1,jc,  kc+1))  +   &
                           q * (pm1 * volume_grid(n)%q(1,ic,  jc+1,kc+1)   +   &
                                  p * volume_grid(n)%q(1,ic+1,jc+1,kc+1)))

!!             if (show) then
!!                write (luncrt, '(a, 3f12.6)') ' xyz_target: ', xyz_target(:)
!!                write (luncrt, '(a, 4i4, a, 4i4, f12.3)') &
!!               ' ib,i,j,k:', ib,i,j,k, '  n,ic,jc,kc,f:', n,ic,jc,kc, fline(k)
!!             end if

               if (k == nk) cycle  ! Just need fline(nk)
!!             if (show) then
!!                write (luncrt, '(a, i4, 3f12.6)')  &
!!                   'k, fline(k), fline(k+1), fline(k)*fraction:', &
!!                    k, fline(k), fline(k+1), fline(k)*fraction
!!             end if

!!!            if (fline(k) - fline(k+1)  >  fline(k)*fraction) then ! Prior fix
!!!               intersections(ib)%q(1,i,j,1) = tline(k+1)  ! Stand-off dist.
!!!               exit
!!!            end if

               q = freestream*fs_scale
               if (fline(k) >= q) then  ! Linearly interpolate
                  r = (q - fline(k+1))/(fline(k) - fline(k+1))
                  intersections(ib)%q(1,i,j,1) = (one-r)*tline(k+1) + r*tline(k)
                  exit
               end if

            end do  ! Next k

         end do  ! Next i

      end do  ! Next j

      deallocate (xline, yline, zline, tline, fline)

   end do  ! Next block

   call release_adt ()

   deallocate (conn)


!  Save the stand-off data as a Tecplot file:
!  ------------------------------------------

   call save_results ()  ! Local procedure below


!! call xyz_write (lunout, true, nblocks, inner_surface, ios)

!! open (lunout, file='shock_stand_off.f', status='unknown')

!! call   q_write (lunout, true, nblocks, nf, intersections, ios)


!  Let F90 do the other deallocates.

99 continue

!  Local procedure for program shock_stand_off:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine make_3D_if_needed ()

!     Read the (already opened) grid file header to see if it is 2-D.  Any
!     read error is interpreted as being due to an unformatted 3-D grid file
!     that will be processed normally upon return.
!     If the grid has only two dimensions, reorganize the data to have a
!     dimension of 2 with k = 1 at the wall, which is assumed to be at j = 1
!     in the input grid.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

      character :: numi*8, numj*8, numk*8

!     Execution:

      read (lunvolg, *, iostat=ios) nblocks

      not_2D = ios /= 0  ! Probably an unformatted file;  assume it's 3-D

      if (not_2D) go to 98

!!!   read (lunvolg, *, iostat=ios) (ni, nj, nk, ib = 1, nblocks)

!!!   Reading 1.e-30 as an integer has been observed to give nk = 0 & ios = 0.

      read (lunvolg, *, iostat=ios) (numi, numj, numk, ib = 1, nblocks)

      if (ios /= 0) then
         write (luncrt, '(/, 2a)') &
            'Trouble reading an apparently formatted grid file header ', &
            'as character data.'
         go to 99
      end if

      not_2D = index (numk, '.') == 0  ! Dimensions appear to be 3-D

!!!   write (luncrt, '(a, i5, l2)') 'ios, not_2D:', ios, not_2D

      if (not_2D) go to 98

!     It appears to be a 2-D formatted file.  Read it the tedious way:

      rewind (lunvolg);  read (lunvolg, *)  ! nblocks

      allocate (volume_grid(nblocks), stat=ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble allocating grid (2-D case).'
         go to 99
      end if

      read (lunvolg, *, iostat=ios) &
         (volume_grid(ib)%ni, volume_grid(ib)%nk, ib = 1, nblocks)

      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading block dimensions (2-D case).'
         go to 99
      end if

      nj = 2;  volume_grid(:)%nj = nj  ! By the time we're done

      do ib = 1, nblocks  ! Read the grid and turn it into 2 planes

         ni = volume_grid(ib)%ni;  nk = volume_grid(ib)%nk
         write (luncrt, '(a, 4i5)') 'nb,ni,nj,nk:', nblocks, ni,nj,nk

         allocate (volume_grid(ib)%x(ni,nj,nk), volume_grid(ib)%y(ni,nj,nk), &
                   volume_grid(ib)%z(ni,nj,nk), stat=ios)

         if (ios /= 0) then
            write (luncrt, '(/, a)') 'Trouble allocating grid x/y/z (2-D case).'
            go to 99
         end if

         read (lunvolg, *, iostat=ios) &
            volume_grid(ib)%x(:,1,:), volume_grid(ib)%z(:,1,:)

         if (ios /= 0) then
            write (luncrt, '(/, a, i5)') &
               'Trouble reading grid x/z (2-D case).  Block #:', ib
            go to 99
         end if

         volume_grid(ib)%y(:,1,:) = zero

         volume_grid(ib)%x(:,2,:) = volume_grid(ib)%x(:,1,:)
         volume_grid(ib)%y(:,2,:) = 0.1
         volume_grid(ib)%z(:,2,:) = volume_grid(ib)%z(:,1,:)

      end do

      close (lunvolg)

      read (lunvolf, *, iostat=ios) nf  ! # blocks

      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading # blocks in function file.'
         go to 99
      else if (nf /= nblocks) then
         write (luncrt, '(/, a, 2i5)') &
            'Mismatch g and f block counts:', nblocks, nf
         ios = 1
         go to 99
      end if

      read (lunvolf, *, iostat=ios) &
         (volume_grid(ib)%mi, volume_grid(ib)%mk, nf, ib = 1, nblocks)

      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading 2-D function file dimens.'
         go to 99
      end if

      volume_grid(:)%mj = nj  ! When we're done

      do ib = 1, nblocks

         ni = volume_grid(ib)%mi;  nk = volume_grid(ib)%mk
         allocate (volume_grid(ib)%q(nf,ni,nj,nk), stat=ios)

         if (ios /= 0) then
            write (luncrt, '(/, a)') 'Trouble allocating grid x/y/z (2-D case).'
            go to 99
         end if

         read (lunvolf, *, iostat=ios) &
            (((volume_grid(ib)%q(n,i,1,k), i = 1, ni), k = 1, nk), n = 1, nf)

         if (ios /= 0) then
            write (luncrt, '(/, a, i5)') &
               'Trouble reading function file (2-D case).  Block #:', ib
            go to 99
         end if

         volume_grid(ib)%q(:,:,2,:) = volume_grid(ib)%q(:,:,1,:)

      end do

      close (lunvolf)
      go to 99

98    rewind (lunvolg);  ios = 0

99    return

      end subroutine make_3D_if_needed

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine save_results ()

!     Write the CFD surface grid with shock stand-off distance as a single
!     function, in formatted Tecplot format for possible merging with surface
!     and boundary layer edge data from BLAYER.  The available Tecplot_io pkg.
!     could be used, but avoid that clutter for now.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write (lunout, '(a)') &
         'TITLE = "Shock stand-off distance"', &
         'VARIABLES = "x, m", "y, m", "z, m", "Stand-off, m"'

      do ib = 1, nblocks

         write (lunout, '(a, i2, a)') 'ZONE T="Zone', ib, '"'
         write (lunout, '(a,i3,a,i3,a)') &
            'I=', intersections(ib)%mi, ', J=', intersections(ib)%mj, &
            ', K=1, ZONETYPE=Ordered'

         if (not_2D) then

            write (lunout, '(a)') 'DATAPACKING=BLOCK'

            write (lunout, '(1p, 8e14.6)') inner_surface(ib)%x
            write (lunout, '(1p, 8e14.6)') inner_surface(ib)%y
            write (lunout, '(1p, 8e14.6)') inner_surface(ib)%z
            write (lunout, '(1p, 8e14.6)') intersections(ib)%q

         else  ! Point order is more convenient for 2D; allow contouring still

            write (lunout, '(a)') 'DATAPACKING=POINT'

            do j = 1, intersections(ib)%mj
               do i = 1, intersections(ib)%mi
                  write (lunout, '(1p, 5e14.6)') &
                     inner_surface(ib)%x(i,j,1), inner_surface(ib)%y(i,j,1), &
                     inner_surface(ib)%z(i,j,1), intersections(ib)%q(1,i,j,1)
               end do
            end do

         end if

      end do

      end subroutine save_results

   end program shock_stand_off
