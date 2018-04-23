!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program project_to_sphere
!
!  Description:
!
!     Project a structured multiblock surface grid to the intended surface of a
!  sphere of specified center and radius.
!
!  History:
!
!     11/19/05  D. A. Saunders  Adaptation of SURFACE_CURVATURE, prompted partly
!                               by the test hemisphere & partly by the baseline
!                               CEV capsule.
!     08/06/10     "     "      The corrections were in the wrong direction!
!
!  Author:  David Saunders, ELORET/NASA Ames Research Ctr, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure
   use xyzq_io_module

   implicit none

!  Constants:

   integer, parameter :: &
      lunin  = 1,        &
      lunout = 2,        &
      lunkbd = 5,        &
      luncrt = 6

   real, parameter :: &
      zero = 0.

   character, parameter :: &
      blank * 1   = ' ',   &
      format * 11 = 'unformatted'

!  Variables:

   integer :: &
      i, i1, ib, ios, j, nblocks, ni, nj, npts

   real :: &
      dr, dx, dy, dz, drmax, drmin, drmean, phi, r, radius, theta,       &
      x, y, z, xc, yc, zc, xnew, ynew, znew

   logical :: &
      formatted_in, formatted_out

   character :: &
      answer * 1, filename * 80

!  Composite data type:

   type (grid_type), pointer, dimension (:) :: &
      grid

!  Execution:

   write (luncrt, '(/, a)', advance='no') ' Input grid name:  '
   read  (lunkbd, *) filename
   write (luncrt, '(a)', advance='no') ' Formatted?  [y|n]: '
   read  (lunkbd, *) answer
   formatted_in = answer == 'y' .or. answer == 'Y'
   i1 = 1;  if (formatted_in) i1 = 3

   open (lunin, file=filename, form=format(i1:11), status='old', iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') &
         ' Unable to open input grid: ', filename(1:len_trim(filename))
      go to 99
   end if

!  Allocate grid coordinate arrays and read the header records:

   call xyz_header_io (1, lunin, formatted_in, nblocks, grid, ios)
   if (ios /= 0) go to 99

   write (luncrt, '(a)', advance='no') ' Output grid name: '
   read  (lunkbd, *) filename
   write (luncrt, '(a)', advance='no') ' Formatted?  [y|n]: '
   read  (lunkbd, *) answer
   formatted_out = answer == 'y' .or. answer == 'Y'
   i1 = 1;  if (formatted_out) i1 = 3

   open (lunout, file=filename, form=format(i1:11), status='unknown',iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') &
         ' Unable to open output grid: ', filename(1:len_trim(filename))
      go to 99
   end if

   call xyz_header_io (2, lunout, formatted_out, nblocks, grid, ios)
   if (ios /= 0) go to 99

   write (luncrt, '(a)', advance='no') ' Intended sphere radius: '
   read  (lunkbd, *) radius
   write (luncrt, '(a)', advance='no') ' Intended center coords: '
   read  (lunkbd, *) xc, yc, zc

!  Process one block at a time:

   do ib = 1, nblocks

      call xyz_allocate (grid(ib), ios)
      if (ios /= 0) go to 99

      ni = grid(ib)%ni;  nj = grid(ib)%nj;  npts = ni * nj

      call xyz_block_io (1, lunin, formatted_in, npts,          &
                         grid(ib)%x, grid(ib)%y, grid(ib)%z, ios)
      if (ios /= 0) go to 99

      drmax = zero;  drmin = 999.;  drmean = zero

      do j = 1, nj
         do i = 1, ni
            x = grid(ib)%x(i,j,1) - xc
            y = grid(ib)%y(i,j,1) - yc
            z = grid(ib)%z(i,j,1) - zc
            r = sqrt (x*x + y*y + z*z)

            dr     = r - radius
            drmin  = min (dr, drmin);  drmax = max (dr, drmax)
            drmean = abs (dr) + drmean
            phi    = atan2 (y, x)
            theta  = acos  (z / r)

            xnew   = x - dr * sin (theta) * cos (phi)
            ynew   = y - dr * sin (theta) * sin (phi)
            znew   = z - dr * cos (theta)

            grid(ib)%x(i,j,1) = xnew + xc
            grid(ib)%y(i,j,1) = ynew + yc
            grid(ib)%z(i,j,1) = znew + zc
         end do
      end do

      drmean = drmean / real (npts)

      write (luncrt, '(a, i4, a, 1p, 3e15.7)') &
         ' Block:', ib, '  min, max, average dr:', drmin, drmax, drmean
 
      call xyz_block_io (2, lunout, formatted_out, npts,        &
                         grid(ib)%x, grid(ib)%y, grid(ib)%z, ios)
      if (ios /= 0) go to 99

      deallocate (grid(ib)%x, grid(ib)%y, grid(ib)%z)

   end do ! Next block

99 continue

! *** stop ! Avoid system dependencies.

  end program project_to_sphere
