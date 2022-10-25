!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine unique_xyz_list (npts, xyz, nf, f, nvertices, ncells, conn)
!
!  Description:
!
!     For a list of (x,y,z) points with possible repeats, and a list of pointers
!  in to this list, search for non-unique (x,y,z)s, suppress repeats, and adjust
!  the pointer list appropriately.  The input arguments are updated in place.  A
!  brute force method is used initially.
!
!  History:
!
!     09/26/2022  D.A.Saunders  Initial implementation as a work-around for the
!                               straightforward conversion of a structured grid
!                               to unstructured form.  Evidently, non-uniqueness
!                               is not handled by the US3D flow solver.
!     10/14/2022     "     "    The bugs seem to be out after long delays and a
!                               lot of puzzling results.
!     10/19/2022     "     "    The 1-D x, y, z arguments became 2-D xyz(3,npts)
!                               as part of handling function data. This seems to
!                               be twice as slow.
!     10/20/2022     "     "    Handle function data now.
!
!  Author:  David A. Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (inout) :: npts    ! Length of input and output (x,y,z) list
   real,    intent (inout), dimension (3,npts) :: xyz  ! List of 3-space coords.
   integer, intent (in) :: nf         ! Number of function variables, >= 0
   real,    intent (inout), dimension (nf,npts) :: f   ! Optional function data
   integer, intent (in) :: nvertices  ! Number of vertices in each cell
   integer, intent (in) :: ncells     ! Number of cells in the unstructured grid
   integer, intent (inout), dimension (nvertices,ncells) :: conn  ! Connectivity
                                                        ! info, updated in place
!  Local constants:

   real, parameter :: tol = 1.e-7  ! Fraction of average (x,y,z) magnitude for a
                                   ! test of an x, y, or z match
!  Local variables:

   integer :: i, icell, icon, ipt, npts_new
   real    :: avgmag, dsq, eps, epssq, sumsq

   real, allocatable :: xyzsq(:)

!  Execution:

!  Calculate a tolerance for x, y, z comparisons (squared magnitudes being
!  necessary but not sufficient for uniqueness tests).

   allocate (xyzsq(npts))

   sumsq = 0.
   do i = 1, npts
      xyzsq(i) = xyz(1,i)**2 + xyz(2,i)**2 + xyz(3,i)**2
      sumsq    = xyzsq(i) + sumsq
   end do

!!!avgmag = sqrt (sum (xyzsq(:))/real (npts))
   avgmag = sqrt (sumsq/real (npts))
   write (*, '(/, a, i9)') 'npts (not unique):', npts
   write (*, '(a, es15.7)') 'avg. magnitude:', avgmag

   eps = avgmag*tol
   epssq = eps*eps

!  Search from the top for a repeated (x,y,z), and adjust conn(:,:) as needed:

   npts_new = npts
   ipt = 0  ! Index of point being compared to those below it

   do  ! Until ipt = the end of a list that may be shortening

      ipt = ipt + 1
      if (ipt == npts_new) exit

      if (mod (ipt, 100000) == 0) write (*, '(a, i9)') 'Processing point ', ipt

      dsq = xyzsq(ipt)
      i = ipt + 1  ! Index of point to compare with
!!!   write (70, '(a, 2i9)') '#i, ipt:', i, ipt

      do  ! To the end of a list that may be shortening
         if (i > npts_new) exit

!!!      write (70, '(a, 2i9)') '   i, ipt:', i, ipt

         if (abs (xyzsq(i) - dsq) < epssq) then
            if (abs (xyz(1,i) - xyz(1,ipt)) < eps .and. &
                abs (xyz(2,i) - xyz(2,ipt)) < eps .and. &
                abs (xyz(3,i) - xyz(3,ipt)) < eps) then ! Duplicate; suppress it

!!!             write (50, '(a, i9, 3es16.8)') 'ipt, x/y/z/ipt:', &
!!!                                             ipt, xyz(:,ipt)

                call move_up () ! Slide x/y/z/dsq up & adjust affected pointers

            end if
         end if

         i = i + 1

      end do
   end do

   deallocate (xyzsq)
   npts = npts_new
   write (*, '(/, a, i9)') 'npts (unique):', npts

!  Internal procedure for subroutine unique_xyz_list

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine move_up ()  ! Suppress a duplicate by moving up in place.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: ie, ie2, iv, iv2, j

!!!   write (52, '(a, 3i9)') 'i, ipt, npts_new:', i, ipt, npts_new
!!!   write (52, '(a, 3es16.8)') 'x/y/z(i):', xyz(:,i)

!!    x(i:npts_new-1) = x(i+1:npts_new)  ! Doing all in the same loop
!!    y(i:npts_new-1) = y(i+1:npts_new)  ! should be more efficient (?)
!!    z(i:npts_new-1) = z(i+1:npts_new)

      do j = i, npts_new - 1
         xyz(:,j) = xyz(:,j+1)
         xyzsq(j) = xyzsq(j+1)
      end do

      if (nf > 0) then
         do j = i, npts_new - 1
            f(:,j) = f(:,j+1)
         end do
      end if

      npts_new = npts_new - 1

      do ie = 1, ncells
         do iv = 1, nvertices
            if (conn(iv,ie) == i) then
!!!            write (99, '(a, 2i8, 8i9)') 'Start ipt, ie:', ipt, ie, conn(:,ie)
               conn(iv,ie) = ipt
!!!            write (99, '(a, 2i8, 8i9)') 'End   ipt, ie:', ipt, ie, conn(:,ie)
            end if
         end do
      end do

!     Further pointers may need to move up 1 too:

      do ie2 = 1, ncells
         do iv2 = 1, nvertices
            if (conn(iv2,ie2) > i) conn(iv2,ie2) = conn(iv2,ie2) - 1
         end do
!!!      write (91, '(a, i5, 4i9)') '  Further ptrs', ie2, conn(:,ie2)
      end do

      end subroutine move_up

   end subroutine unique_xyz_list
