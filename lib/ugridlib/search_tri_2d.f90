!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine search_tri_2d (xytarg, nvert, xyvert, ntri, ivert, &
                             itri, p, q, r, xyinterp, dsqmin)
!  Description:
!
!     Rapid searching of a 3-space triangulated surface is well handled by an
!  existing ADT (Alternating Digital Tree) method which uses bounding box
!  techniques and a search tree to avoid more than data-range testing for the
!  bulk of cells in a given mesh.  See the earlier /adt/triangles and/or
!  /adt/all_grids directories accompanying this /ugridlib utility.
!
!     For the case of a 2-space triangulated mesh (possibly from a Delaunay
!  triangulation of a list of (x,y) points), we could either carry around a
!  third coordinate (z = 0 for all points) or we could implement a 2-space
!  variant of the ADT scheme that doesn't need to do so.  If mesh size becomes
!  an issue, these are still viable options, but for modest numbers of points,
!  they are overkill.  Therefore, this utility simply tests all triangles for
!  containing a target (x,y) point, and returns a triangle and associated
!  interpolation coefficients representing the best possible result whether the
!  target point is contained by the mesh or not.  More precisely, it returns a
!  best-possible point of interpolation that is not outside the mesh, and the
!  corresponding squared distance from the target. If dsq > 0, then the target
!  is outside the mesh and a best point on an outer edge of the mesh is returned
!  as xyinterp that is as close as possible to the target.
!
!  History:
!
!     07/25/2017  D.A.Saunders  Initial implementation as part of identifying
!                               target points too far for meaningful interpo-
!                               lation within scattered 2-space data points
!                               via the Delaunay triangulation.
!
!  Author:  David Saunders, AMA, Inc./NASA Ames Research Ctr., Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real,    intent (in)  :: xytarg(2)        ! Target point in 2-space
   integer, intent (in)  :: nvert            ! Number of vertices in the mesh
   real,    intent (in)  :: xyvert(2,nvert)  ! x,y coords. of (unique) mesh pts.
   integer, intent (in)  :: ntri             ! Number of triangular elements
   integer, intent (in)  :: ivert(3,ntri)    ! Pointers into xyvert(:,:) for
                                             ! each triangular element
   integer, intent (out) :: itri             ! Best triangle found
   real,    intent (out) :: p, q, r          ! Interp. coefs. for that triangle;
                                             ! p + q + r = 1;
                                             ! p x1 + q x2 + r x3 = xyinterp,
                                             ! where x1(:) = xyvert(i1,itri),
                                             ! and i1 = ivert(1,itri), etc.
   real,    intent (out) :: xyinterp(2)      ! Point closest to xytarg(:) that
                                             ! is not outside the triangulation
   real,    intent (out) :: dsqmin           ! The corresponding squared dist-
                                             ! ance from the target, >= 0.
!  Local variables:

   integer :: i, i1, i2, i3 
   real    :: dsq, pi, qi, ri, xnear(2)

!  Execution:

   dsqmin = 1.e+30

   do i = 1, ntri                                             
      i1 = ivert(1,i)
      i2 = ivert(2,i)
      i3 = ivert(3,i)

      call nearest_tri_pt_2d (xyvert(:,i1), xyvert(:,i2), xyvert(:,i3), &
                              xytarg(:), xnear(:), pi, qi, ri, dsq)
      if (dsq < dsqmin) then
          dsqmin = dsq
          itri = i
          p = pi;  q = qi;  r = ri
          xyinterp(:) = xnear(:)
          if (dsq == 0.) exit
      end if
   end do  ! Next triangle

   end subroutine search_tri_2d
