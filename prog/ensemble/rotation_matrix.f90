!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine rotation_matrix (nblocks, grid, nlist, list, R)

!  This routine is prompted by the need to transform the grid blocks defining an
!  arbitrary cavity or gouge in a surface to a normalized cavity for comparisons
!  with other cases.  For the gouge-related subset of grid blocks, the principal
!  axes needed to align the cavity with the working axes are taken to be defined
!  by the eigenvectors of the inertia matrix established here.  The eigenvectors
!  in turn are stacked as three unit column vectors of a rotation matrix R which
!  describes the principal axes relative to the working axes.
!
!  Treating the nonuniform volume grid points as though they represent particles
!  uniformly distributed with unit mass is not very precise. Nor is inclusion of
!  common face points more than once.   However, mapping an arbitrary shape to a
!  cube is inherently approximate, so these considerations are ignored, at least
!  initially.
!
!  References:  http://kwon3d.com/ (Theories and Practices of Motion Analysis);
!               Introduction to Robotics Mechanics & Control, John J. Craig
!
!  06/02/05  David Saunders  Initial implementation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure ! Module from the XYZQ_IO package defining grid_type

   implicit none

!  Arguments:

   integer, intent (in)          :: nblocks        ! Dimension of grid(:)
   type (grid_type), intent (in) :: grid(nblocks)  ! Grid blocks
   integer, intent (in)          :: nlist          ! Active length of list(:)
   integer, intent (in)          :: list(nlist)    ! Blocks to be processed
   real,    intent (out)         :: R(3,3)         ! Rotation matrix describing
                                                   ! the principal axes of the
                                                   ! cavity w.r.t. Ox, Oy, Oz;
!  Normal usage should produce columns of R in the order that corresponds to
!  principal axes along the cavity length, width, and depth respectively.  The
!  application should swap columns 2 & 3 if cavity width is less than depth.

!  Local constants:

   real, parameter :: zero = 0.

!  Local variables:

   integer             :: i, j, k, l, ib, ier
   real                :: xx, yy, zz, xy, yz, zx
   real, dimension (3) :: w,                    &  ! Eigenvalues, descending
                          fv1, fv2                 ! Work-space
   real                :: A(3,3)                   ! Inertia matrix; columns
                                                   ! correspond to eigenvectors
!  Execution:

!  Generate the inertia matrix from unit masses at the volume grid points:

   xx = zero;  yy = zero;  zz = zero
   xy = zero;  yz = zero;  zx = zero

   do l = 1, nlist
      ib = list(l)
      do k = 1, grid(ib)%nk
         do j = 1, grid(ib)%nj
            do i = 1, grid(ib)%ni
               x = grid(ib)%x(i,j,k)
               y = grid(ib)%y(i,j,k)
               z = grid(ib)%z(i,j,k)
               xx = y*y + z*z + xx
               yy = z*z + x*x + yy
               zz = x*x + y*y + zz
               xy = x*y + xy
               yz = y*z + yz
               zx = z*x + zx
            end do
         end do
      end do
   end do

   A(1,1) = xx;  A(1,2) = xy;  A(1,3) = zx
   A(2,1) = xy;  A(2,2) = yy;  A(2,3) = yz
   A(3,1) = zx;  A(3,2) = yz;  A(3,3) = zz

!  Calculate the eigenvectors of the inertia matrix:

   call rs (3, 3, A, w, 1, R, fv1, fv2, ier)

   if (ier /= 0) then
      write (*, *) ' Eigenvector routine "rs" reports "tql2" error: ', ier
      stop
   end if

   end subroutine rotation_matrix
