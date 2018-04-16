!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine gaussian_curvature (ni, nj, x, y, z, curv_gauss, curv_mean)
!
!  Calculate the Gaussian curvature and mean curvature distributions for one
!  regular surface patch.  The bulk of the calculations are the same for both
!  types of curvature, so returning both is more efficient than returning one
!  or the other from separate calls.
!
!  Assumptions:
!
!  o  Differencing along grid lines is not the same as fixing v and differencing
!     in the u direction, and vice versa, unless the patch is rectangular.
!     Therefore, the closer to rectangular the patch is, the less approximate
!     the derivatives and hence curvature values will be.
!  o  A degenerate patch with one edge collapsed would seem to be a problem for
!     the arc-length-based derivatives used here.  However, PARAM2D handles such
!     edges by inserting a uniform distribution for the normalized arc lengths,
!     so a result is still produced.  How meaningful it is on the collapsed edge
!     is another matter.
!  o  In these days of dynamic allocation, assuming packed arrays should not be
!     a drawback, and it keeps the argument list down.
!
!  Definitions:
!
!     At any point (i,j), the two types of curvature may be written as
!
!                  L N  - M**2
!        K  =  ------------------              (Gaussian curvature)
!              g11 g22  -  g12**2
!
!              1  g11 N  -  2 g12 M  +  g22 L
!        H  =  -  ---------------------------    (mean curvature)
!              2      g11 g22  -  g12**2
!
!     where (using ' to mean transpose):
!
!        g11  =  t1't1  =   Xu**2  +  Yu**2   +  Zu**2;    t1  =  (Xu  Yu  Zu)'
!        g12  =  t1't2  =   Xu Xv  +  Yu Yv   +  Zu Zv;    t2  =  (Xv  Yv  Zv)'
!        g22  =  t2't2  =   Xv**2  +  Yv**2   +  Zv**2;     n  =   unit normal
!
!         L   =  n't11  =  n1 Xuu  +  n2 Yuu  +  n3 Zuu;  t11  = (Xuu Yuu Vuu)'
!         M   =  n't12  =  n1 Xuv  +  n2 Yuv  +  n3 Zuv;  t12  = (Xuv Yuv Zuv)'
!         N   =  n't22  =  n1 Xvv  +  n2 Yvv  +  n3 Zvv;  t22  = (Xvv Yvv Zvv)'
!
!     Also, the relationships to the two principal curvatures k1, k2 are:
!
!         K   =   k1 k2                 and              H  = (k1 + k2) / 2
!
!     where k1 and k2 are the roots of  k**2  -  2Hk  +  K  =  0:
!
!         k1  =  H  +  sqrt (H**2 - K)  and  k2  =  H  -  sqrt (H**2 - K)
!
!     Note that right-handed patches on convex solid surfaces (normals pointing
!     out) have negative curvature.  Thus, for right-handed patches on a sphere
!     of radius r, the principal and mean curvatures should all be close to -1/r
!     if the grid lines are close to orthogonal.  The Gaussian curvatures should
!     all be close to 1/r**2.
!
!  References:
!
!     Mathematica and various web sites, none of which appear to offer software
!     for structured grids.
!
!  History:
!
!     09/07/05  D.A.Saunders  Initial implementation.
!     10/15/05    "     "     Safeguard extreme cases such as at the corner of
!                             a spherical cap patch (virtually 180 deg angle).
!                             Simply recalculating the unit normal from the cell
!                             diagonals doesn't help much, so kludge it: use
!                             values from the opposite corner of the cell, after
!                             processing all grid points first.  (There must be
!                             a better way.)
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: ni, nj                        ! Structured surface
                                                          ! patch dimensions ...
   real, intent (in),  dimension (ni,nj) :: x, y, z       ! ... and coordinates

   real, intent (out), dimension (ni,nj) :: curv_gauss, & ! Two forms of surface
                                            curv_mean     ! curvature distribn.
!  Local constants:

   real, parameter :: eps  = 1.e-10, &                    ! For g11 g22 - g12**2
                      half = 0.5, zero = 0.0
!  Local variables:

   integer :: i, ii, j, jj, maxpts

   real                                 :: denom, g11, g12, g22, L, M, N, vsize
   real, dimension (3)                  :: un, vec1, vec2
   real, allocatable, dimension (:)     :: abscissa, ordinate, d1, d2, unused
   real, allocatable, dimension (:,:)   :: u, v
   real, allocatable, dimension (:,:,:) :: t1, t2, t11, t12, t22, tuv, tvu

   logical, allocatable :: bad(:,:)

!  Execution:

!  Normalized arc lengths:

   allocate (u(ni,nj), v(ni,nj))

   u(1,1) = zero  ! -999. would suppress normalization

   call param2d (ni, nj, 1, ni, 1, nj, x, y, z, u, v)

!  To make use of 1-D differencing utilities, we process columns then rows.
!  The cross-derivatives t12 are estimated as the averages of two possibilities.

   maxpts = max (ni, nj)
   allocate (unused(maxpts))  ! For unneeded 1-D curvatures

   allocate (t1(ni,nj,3), t2(ni,nj,3), t11(ni,nj,3), t12(ni,nj,3),             &
             t22(ni,nj,3), tuv(ni,nj,3), tvu(ni,nj,3))

   do j = 1, nj

!     1st & 2nd partial derivs. w.r.t. u for column j:

      call fd12k (ni, u(1,j), x(1,j), t1(1,j,1), t11(1,j,1), unused)
      call fd12k (ni, u(1,j), y(1,j), t1(1,j,2), t11(1,j,2), unused)
      call fd12k (ni, u(1,j), z(1,j), t1(1,j,3), t11(1,j,3), unused)

   end do ! Next column

   allocate (abscissa(maxpts), ordinate(maxpts), d1(maxpts), d2(maxpts))

   do i = 1, ni

!     1st & 2nd partial derivs. w.r.t. v for row i:

      call getrow (ni, i, 1, nj, v, abscissa)

      call getrow (ni, i, 1, nj, x, ordinate)
      call fd12k  (nj, abscissa, ordinate, d1, d2, unused)
      call putrow (ni, i, 1, nj,  t2(1,1,1), d1)
      call putrow (ni, i, 1, nj, t22(1,1,1), d2)

      call getrow (ni, i, 1, nj, y, ordinate)
      call fd12k  (nj, abscissa, ordinate, d1, d2, unused)
      call putrow (ni, i, 1, nj,  t2(1,1,2), d1)
      call putrow (ni, i, 1, nj, t22(1,1,2), d2)

      call getrow (ni, i, 1, nj, z, ordinate)
      call fd12k  (nj, abscissa, ordinate, d1, d2, unused)
      call putrow (ni, i, 1, nj,  t2(1,1,3), d1)
      call putrow (ni, i, 1, nj, t22(1,1,3), d2)

!     d(Xu)/dv for row i:

      call getrow (ni, i, 1, nj,  t1(1,1,1), ordinate)
      call fd12k  (nj, abscissa, ordinate, d1, unused, unused)
      call putrow (ni, i, 1, nj, tvu(1,1,1), d1)

      call getrow (ni, i, 1, nj,  t1(1,1,2), ordinate)
      call fd12k  (nj, abscissa, ordinate, d1, unused, unused)
      call putrow (ni, i, 1, nj, tvu(1,1,2), d1)

      call getrow (ni, i, 1, nj,  t1(1,1,3), ordinate)
      call fd12k  (nj, abscissa, ordinate, d1, unused, unused)
      call putrow (ni, i, 1, nj, tvu(1,1,3), d1)

   end do ! Next row

   do j = 1, nj

!     d(Xv)/du for column j:

      call fd12k (ni, u(1,j), t2(1,j,1), tuv(1,j,1), unused, unused)
      call fd12k (ni, u(1,j), t2(1,j,2), tuv(1,j,2), unused, unused)
      call fd12k (ni, u(1,j), t2(1,j,3), tuv(1,j,3), unused, unused)

   end do

   deallocate (u, v, abscissa, ordinate, d1, d2, unused)

!  Complete the calculations at every point.

   allocate (bad(ni,nj));  bad(:,:) = .false.

   do j = 1, nj

      do i = 1, ni

         vec1 = t1(i,j,:);   vec2 = t2(i,j,:)
         g11  = dot_product (vec1, vec1)
         g12  = dot_product (vec1, vec2)
         g22  = dot_product (vec2, vec2)

         call cross (vec1, vec2, un)         ! Normal vector; |un| could be ~ 0
         vsize = sqrt (dot_product (un, un)) ! |un| = |vec1| |vec2| sin (theta)

         if (vsize / max (sqrt (g11*g22), eps) < 0.01) then
            vsize = eps       ! sin (1 degree) ~ 0.017
            bad(i,j) = .true.
         end if

         un = un / vsize                            ! Unit normal

         vec2 =  t11(i,j,:);                        L = dot_product (un, vec2)
         vec2 = (tuv(i,j,:) + tvu(i,j,:)) * half;   M = dot_product (un, vec2)
         vec2 =  t22(i,j,:);                        N = dot_product (un, vec2)

         denom = g11 * g22 - g12 * g12
         denom = sign (max (abs (denom), eps), denom)
         curv_gauss(i,j) = (L * N - M * M) / denom
         curv_mean (i,j) = ((g11*N + g22*L)*half - g12*M) / denom

!!!      if (abs (curv_gauss(i,j)) > 100.) then
!!!         write (10, '(a, 2i5, a, 1p, 3e12.4, a, 3e12.4, a, 3e12.4)') &
!!!         ' i,j:', i,j, '  gs:', g11,g12,g22, '  LMN:', L,M,N, &
!!!         '  denom,K,H:', denom, curv_gauss(i,j), curv_mean(i,j)
!!!      end if

      end do ! Next i

   end do ! Next j

   deallocate (t1, t2, t11, t12, t22, tuv, tvu)

!  Check for fudging bad values seen in hemisphere cap "corners":

   do j = 1, nj
      do i = 1, ni
         if (bad(i,j)) then
            if (i < ni) then
               ii = i + 1
            else
               ii = i - 1
            end if
            if (j < nj) then
               jj = j + 1
            else
               jj = j - 1
            end if
            curv_gauss(i,j) = curv_gauss(ii,jj)
            curv_mean (i,j) = curv_mean (ii,jj)
         end if
      end do
   end do

   deallocate (bad)

   end subroutine gaussian_curvature
