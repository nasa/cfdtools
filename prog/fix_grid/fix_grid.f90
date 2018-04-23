!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program fix_grid
!
!  Description:
!
!     This is an attempt at repairing a defective surface grid by 1-D spline
!     interpolation at specified indices.  Initially, each defect is assumed
!     to be only 1 point wide in either the i or j direction, but this could
!     be generalized.
!
!  Control File Format:
!
!     One or more pairs of i and j index ranges should be listed with the patch
!     number preceding each pair.  For a given pair, the shorter range of the
!     two determines the direction of interpolation.  E.g.,:
!
!        1
!        198
!        2:3
!        1
!        195:231
!        1
!        3
!        183
!        66
!        :
!
!     Here, the first triplet means patch 1 will be interpolated in the i
!     direction to fix points (198,2) and (198,3).
!
!     If the interpolation direction implies crossing a patch boundary,
!     periodicity is assumed.
!
!  History:
!
!     11/13/12  D.A.Saunders  Simple implementation prompted by a SIAD grid
!                             that started as a flawed IGES file.
!
!  Author:  David Saunders, ERC, Inc., at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure   ! Derived data type for one grid block
   use xyzq_io_module         ! PLOT3D multiblock file I/O package (3-space)

!  Constants:

   implicit none

   integer, parameter :: &
      lunin  = 1,        &
      lunout = 2,        &
      lunctl = 3,        &
      lunkbd = 5,        &
      luncrt = 6

   character, parameter :: &
      format * 11 = 'unformatted'

!  Variables:

   integer :: i1, ios, maxint, n, nblock, ni, nj, numf

   logical :: cell_centered, cr, eof, formatted

   character :: answer*1, filename*128, filename2*128

   integer, allocatable :: ilist(:), jlist(:)

   real,    allocatable :: s(:), xline(:), yline(:), zline(:)

!  Derived data types:

   type (grid_type), pointer, dimension (:) :: &
      block                                    ! Grid block array

!  Execution:

   open (lunctl, file='fix_grid.inp', status='old', iostat=ios)
   if (ios /= 0) then
      write (luncrt, '(a)') &
         'You need (patch,irange,jrange) triples in a fix_grid.inp file.'
      go to 99
   end if

!  Prompt for the input file and check for its existence and apparent form:

   ios = 1
   do while (ios /= 0)

      call reads (luncrt, 'Input surface grid (PLOT3D /mgrid): ', &
                  lunkbd, filename, cr, eof)
      if (eof) go to 99   ! Quit
      if (cr) cycle

      call determine_grid_form (filename, lunin, formatted, ios)

   end do

   i1 = 1;  if (formatted) i1 = 3

   open (lunin, file=filename, form=format(i1:11), status='OLD')

!  Read all blocks.

   call xyzq_read (lunin, -1, formatted, nblock, numf, cell_centered, &
                   block, ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Trouble reading the grid.'
      go to 99
   end if

!  Display the block dimensions to reassure the user:

   write (luncrt, '(/, (i6, 2x, 3i5))') &
      (n, block(n)%ni, block(n)%nj, block(n)%nk, n = 1, nblock)

!  The largest patch dimension bounds the possible number of interpolations.

   maxint = 1
   do n = 1, nblock
      maxint = max (maxint, block(n)%ni, block(n)%nj)
   end do

   allocate (ilist(maxint), jlist(maxint), s(maxint), &
             xline(maxint), yline(maxint), zline(maxint))

!  Loop over possibly several sets of interpolations:
!  --------------------------------------------------

   do  ! Until EOF

      read (lunctl, *, iostat=ios) n  ! Patch number
      if (ios < 0) exit

      ni = maxint
      call rdlist (luncrt, 'Reading i range: ', lunctl, ni, ilist)
      nj = maxint
      call rdlist (luncrt, 'Reading j range: ', lunctl, nj, jlist)
      write (luncrt, '(20i4)') ilist(1:ni)
      write (luncrt, '(20i4)') jlist(1:nj)

      call fix_this_range ()

   end do

   deallocate (ilist, jlist, s, xline, yline, zline)

80 continue

!  Save repaired grid?

   filename2 = ' '
   call reads (luncrt, 'Output file name [^D = quit without saving]: ', &
               lunkbd, filename2, cr, eof)

   if (.not. eof) then

      if (cr) filename2 = filename
      if (filename2 == filename) then

         call readc (luncrt, 'Overwrite the input file?  [Y|N; <cr>=Yes]: ', &
                     lunkbd, answer, cr, eof)
         if (eof) go to 80

         cr = cr .or. answer == 'Y'

         if (.not. cr) go to 80

      end if

      formatted = .false.
      call ready (luncrt, 'Formatted (Y) or unformatted (N); <cr>=No: ', &
                  lunkbd, formatted, cr, eof)
      if (eof) go to 80

      i1 = 1;  if (formatted) i1 = 3

      open (lunout, file=filename2, form=format(i1:11), status='unknown')

      call xyz_write (lunout, formatted, nblock, block, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble saving the grid.'
      end if

   end if

   do n = 1, nblock
      deallocate (block(n)%x, block(n)%y, block(n)%z)
   end do
 
99 continue

! *** STOP ! Avoid system dependencies.


!  FIX_GRID internal procedure:

   contains

!     ----------------------------
      subroutine fix_this_range ()
!     ----------------------------

!     Fix a set of bad surface grid points defined by ilist(1:ni) & jlist(1:nj).
!     The shorter list defines the direction of each interpolation.
!     For now, assume the shorter list has length 1.

      real, parameter :: half = 0.5, silly = 1.e+10
      integer :: i, ii, ieval, imx, j, jj, jeval, jmx, nline
      real    :: derivs(3), sint, total, xint, yint, zint

      if (ni <= nj) then  ! Interpolate in the i direction
         imx   = block(n)%ni
         ieval = 1
         i     = ilist(1)
         nline = imx - 1
         do j = 1, nj
            jj = jlist(j)
            xline(1:i-1)   = block(n)%x(1:i-1,jj,1)
            xline(i:nline) = block(n)%x(i+1:nline+1,jj,1)
            yline(1:i-1)   = block(n)%y(1:i-1,jj,1)
            yline(i:nline) = block(n)%y(i+1:nline+1,jj,1)
            zline(1:i-1)   = block(n)%z(1:i-1,jj,1)
            zline(i:nline) = block(n)%z(i+1:nline+1,jj,1)
            call chords3d (nline, xline, yline, zline, .false., total, s)
            sint = (s(i-1) + s(i))*half
            call plscrv3d (nline, xline, yline, zline, s, 'B', .true., &
                           .false., sint, ieval, xint, yint, zint, derivs)
            block(n)%x(i,jj, 1) = xint
            block(n)%y(i,jj, 1) = yint
            block(n)%z(i,jj, 1) = zint
         end do
      else  ! Interpolate in the j direction
         jmx   = block(n)%nj
         jeval = 1
         j     = jlist(1)

         if (j == 1) then  ! Treat SIAD case, where j = 1 matches j = jmx
            nline = 4  ! Points j = jmx-2, jmx-1, 2, 3 for 4-point local spline
            do i = 1, ni
               ii = ilist(i)
               xline(1:2) = block(n)%x(ii,jmx-2:jmx-1,1)
               xline(3:4) = block(n)%x(ii,2:3,1)
               yline(1:2) = block(n)%y(ii,jmx-2:jmx-1,1)
               yline(3:4) = block(n)%y(ii,2:3,1)
               zline(1:2) = block(n)%z(ii,jmx-2:jmx-1,1)
               zline(3:4) = block(n)%z(ii,2:3,1)
               call chords3d (nline, xline, yline, zline, .false., total, s)
               sint = (s(2) + s(3))*half
               call plscrv3d (nline, xline, yline, zline, s, 'B', .true., &
                              .false., sint, jeval, xint, yint, zint, derivs)
               block(n)%x(ii,1,1) = xint;  block(n)%x(ii,jmx,1) = xint
               block(n)%y(ii,1,1) = yint;  block(n)%y(ii,jmx,1) = yint
               block(n)%z(ii,1,1) = zint;  block(n)%z(ii,jmx,1) = zint
            end do

         else  ! No awkward boundary condition
            nline = jmx - 1
            do i = 1, ni
               ii = ilist(i)
               xline(1:j-1)   = block(n)%x(ii,1:j-1,1)
               xline(j:nline) = block(n)%x(ii,j+1:nline+1,1)
               yline(1:j-1)   = block(n)%y(ii,1:j-1,1)
               yline(j:nline) = block(n)%y(ii,j+1:nline+1,1)
               zline(1:j-1)   = block(n)%z(ii,1:j-1,1)
               zline(j:nline) = block(n)%z(ii,j+1:nline+1,1)
               call chords3d (nline, xline, yline, zline, .false., total, s)
               sint = (s(j-1) + s(j))*half
               call plscrv3d (nline, xline, yline, zline, s, 'B', .true., &
                              .false., sint, jeval, xint, yint, zint, derivs)
               if (abs (xint) > silly .or. abs (yint) > silly .or. &
                   abs (zint) > silly) write (luncrt, '(a, 3i4)') &
                   '   Trouble at (n,i,j) =', n, i, j
               block(n)%x(ii,j,1) = xint
               block(n)%y(ii,j,1) = yint
               block(n)%z(ii,j,1) = zint
            end do
         end if
      end if

      end subroutine fix_this_range

   end program fix_grid
