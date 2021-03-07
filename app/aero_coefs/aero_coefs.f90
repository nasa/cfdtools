!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program aero_coefs
!
!  Description:
!
!     For a multiblock surface dataset in PLOT3D form, with vertex-centered
!     pressures or Cps in the function file (nf >= 1), calculate aerodynamic
!     force and moment coefficients.  The assumed right-handed coordinate
!     system is as follows, although it probably doesn't matter what the
!     convention is.
!
!        Ox points downstream; Oy points horizontally; Oz points up.
!
!     SI units are assumed.
!
!     Prompts suffice for the input files.  Providing pressure coefficients
!     rather than surface pressures suppresses prompts for free stream density
!     and velocity and associated divides by half*rhofree*Vfree^2.
!
!     Results go to standard output.
!
!     See also the accompanying LoverD utility for angle attack effects.
!
!  History:
!
!     06/01/2018  D.A.Saunders  Initial implementation, to help Dinesh Prahbu
!                               with RCS thruster calculations.
!     06/04/2018    "     "     Dinesh recommended handling of either pressures
!                               or pressure coefficients, and possibly shear
!                               stresses for viscous forces some day (not now).
!     06/05/2018    "     "     Dinesh found a couple of glitches.
!     02/28/2019     "    "     Surface normals point outward from a convex
!                               right-handed surface, so they are now negated
!                               in order for the x component of the forces to
!                               be positive as one would normally expect.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure   ! Data structure for one grid bloc
   use xyzq_io_module         ! PLOT3D I/O package

   implicit none

!  Constants:

   integer, parameter :: &
      luning = 1,        &  ! Input structured surface grid  (PLOT3D)
      luninf = 2,        &  ! Corresponding pressures
      lunkbd = 5,        &  ! Keyboard inputs
      luncrt = 6            ! Screen outputs

   real, parameter :: &
      half = 0.5, zero = 0.

   logical, parameter :: &
      false = .false.,   &
      true  = .true.

!  Variables:

   integer :: &
      i, ib, ios, j, ip, nblocks, ni, nj, nf, npts

   real :: &
      Lref, rhofree, Sref, xyzm(3), Vfree

   logical :: &
      cr, eof, formatted, pressures

   character (1) :: &
      answer

   character (256) :: &
      filename

   real, allocatable, dimension (:) :: &
      total
   real, allocatable, dimension (:,:) :: &
      coefs

!  Derived data types:

   type (grid_type), pointer, dimension (:) :: &
      surface               ! (x,y,z,f) for the input multiblock surface

!  Execution:

   call reads (luncrt, 'Input surface grid (PLOT3D multiblock): ', &
               lunkbd, filename, cr, eof)
   if (eof) go to 99

   call determine_grid_form (filename, luning, formatted, ios)
   if (ios /= 0) go to 99

   if (formatted) then
      open  (luning, file=trim (filename), status='old', iostat=ios)
   else
      open  (luning, file=trim (filename), status='old', &
             form='unformatted', iostat=ios)
   end if
   if (ios /= 0) then
      write (luncrt, '(2a)') '*** Unable to open ', trim (filename)
      go to 99
   end if

   call reads (luncrt, 'Vertex-centered function file: ', &
               lunkbd, filename, cr, eof)
   if (eof) go to 99

   if (formatted) then
      open (luninf, file=trim (filename), status='old', iostat=ios)
   else
      open (luninf, file=trim (filename), status='old', &
            form='unformatted', iostat=ios)
   end if
   if (ios /= 0) then
      write (luncrt, '(2a)') '*** Unable to open ', trim (filename)
      go to 99
   end if

!  Read the input file header(s).

   call xyz_header_io (1, luning, formatted, nblocks, surface, ios)
   if (ios /= 0) go to 99

   call q_header_io (1, luninf, formatted, nblocks, nf, surface, ios)
   if (ios /= 0) go to 99

   write (luncrt, '(a, i3)') ' Number of functions found:', nf, &
      ' Entering pressures requires free stream rho & V; Cps do not.'
   answer='Y'
   call readc (luncrt, 'Assume pressures? [y|n; <cr>=y]: ', &
               lunkbd, answer, cr, eof)
   if (eof) go to 99
   pressures = answer == 'Y'

   ip = 1
   if (pressures) then
      call readi (luncrt, 'Function index for pressure [<cr>=1]: ', &
                  lunkbd, ip, cr, eof)
      if (eof) go to 99
      call readr (luncrt, 'Free stream density:  ', lunkbd, rhofree, cr, eof)
      if (eof) go to 99
      call readr (luncrt, 'Free stream velocity: ', lunkbd, Vfree, cr, eof)
      if (eof) go to 99
   else
      call readi (luncrt, 'Function index for Cp [<cr>=1]: ', &
                  lunkbd, ip, cr, eof)
      if (eof) go to 99
   end if

   write (luncrt, '(a)', advance='no') ' Lref and Sref (2 values):  '
   read  (lunkbd, *) Lref, Sref
   write (luncrt, '(a)', advance='no') ' Moment center (3 coords.): '
   read  (lunkbd, *) xyzm(:)

!  Process one block at a time:

   allocate (coefs(6,nblocks), total(6));  total(:) = zero
 
   write (luncrt, '(/, 2a)') &
      '  Block       cx            cy            cz', &
      '              mx            my            mz'
   do ib = 1, nblocks

      call xyz_allocate (surface(ib), ios)
      if (ios /= 0) go to 99

      ni = surface(ib)%ni;  nj = surface(ib)%nj
      npts = ni * nj

      call xyz_block_io (1, luning, formatted, npts, &
                         surface(ib)%x, surface(ib)%y, surface(ib)%z, ios)
      if (ios /= 0) go to 99

      call q_allocate (surface(ib), nf, ios)
      if (ios /= 0) go to 99

      call q_block_io (1, luninf, formatted, nf, ni, nj, 1, surface(ib)%q, ios)
      if (ios /= 0) go to 99

      call patch_coefs (ni, nj, nf, ip, &
                        surface(ib)%x, surface(ib)%y, &
                        surface(ib)%z, surface(ib)%q, &
                        xyzm, Lref, Sref, coefs(:,ib))

      if (pressures) coefs(:,ib) = coefs(:,ib)/(half*rhofree*Vfree**2)
      coefs(1:3,ib) = coefs(1:3,ib)
      coefs(4:6,ib) = coefs(4:6,ib)
      total(:) = coefs(:,ib) + total(:)

      write (luncrt, '(i7, 3es14.6, 3x, 3es14.6)') ib, coefs(:,ib)

      deallocate (surface(ib)%x, surface(ib)%y, surface(ib)%z, surface(ib)%q)

   end do  ! Next surface block

   write (luncrt, '(/, a, /, 7x, 3es14.6, 3x, 3es14.6)') ' Totals:', total(:)

99 continue

   end program aero_coefs
