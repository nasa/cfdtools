!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine triangulate_patches (npatches, nf, patches, ref_length,          &
                                   nnode, ntri, tri_xyz, tri_f, conn)

!  Convert a structured surface grid to triangulated form.  Function data may
!  accompany the grid (or not).  This variation was needed for compatibility
!  with the "grid_type" derived data type now favored by the author.  Earlier
!  implementations treated packed lists of (x,y,z) triples that meant no need
!  to copy the data, but function data were not handled and compatibility with
!  current I/O utilities is desired here.
!
!  The handedness of each given patch is preserved in each pair of triangles
!  (normally) derived from each quad cell.  Quad. cells are split along the
!  shorter diagonals.  Degenerate quad. cells that are really single triangles
!  are checked for, so that all output triangles have nonzero areas.
!
!  Rather than expecting the calling program to provide enough storage for the
!  triangulation, pointer arguments allow this routine to allocate that storage.
!  However, to avoid two passes through the structured data, this allocation
!  assumes there are no degenerate cells.  Thus, upon return, the output arrays
!  are at least big enough, but will not be fully packed if any singular points
!  appear in the data (such as at the apex of a sharp nose cone).
!
!  Note that the output triangles with vertices 1:3 are such that cross product
!  (v2-v1) x (v3-v1) is in the same direction as cross product (v2-v1) x (v4-v1)
!  from the quad cell with vertices 1=(i,j), 2=(i+1,j), 3=(i+1,j+1), 4=(i,j+1),
!  or (v4-v3) x (v2-v3) if (i,j) is at a collapsed edge.  This clarifies what
!  convention to use for triangle normals, namely (v2-v1) x (v3-v1).
!
!  07/19/06  D.A.Saunders  Initial implementation as part of determining the
!                          projected area of a capsule from a given view point.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Derived data type definition:

   use grid_block_structure            ! From xyzq_io package or equivalent

   implicit none

!  Arguments:

   integer, intent (in)  :: &
      npatches,             &          ! # input surface patches
      nf                               ! # functions accompanying (x,y,z);
                                       ! nf >= 0
   type (grid_type), pointer :: &
      patches(:)                       ! Array of surface patches

   real,    intent (in)  :: &
      ref_length                       ! Reference length used to help
                                       ! identify collapse edges
   integer, intent (out) :: &
      nnode                            ! # points in all surface patches, with
                                       ! no attempt to avoid duplicates
   integer, intent (out) :: &
      ntri                             ! # triangle elements derived from the
                                       ! quad. cells, including only one if a
                                       ! quad. cell is collapsed on one edge
   real, pointer ::         &
      tri_xyz(:,:)                     ! The (x,y,z) coordinates of the patches
                                       ! are repacked as tri_xyz(1:3,1:nnode)
   real, pointer ::         &
      tri_f(:,:)                       ! Any associated functions are returned
                                       ! in tri_f(1:nf,1:nnode)
   integer, pointer ::      &
      conn(:,:)                        ! Vertices 1:3 of triangle n are pointed
                                       ! to in tri_xyz(:,*) by conn(1:3,n);
                                       ! vertices 1, 2, 3 are always such that
                                       ! cross product (v2-v1) x (v3-v1) points
                                       ! "out" if the normal to the associated
                                       ! quad. cell pointed "out", & vice versa
!  Local constants:

   real, parameter :: tol = 1.e-10     ! Applied to reference length to identify
                                       ! collapsed edges

!  Local variables:

   integer :: &
      i, i1, i2, i3, i4, ib, ip, it, j, n, ni, nj, npoints, ntrimax

   real :: &
      diag_sq_13, diag_sq_24, dsq_neglect, &
      edge_sq_12, edge_sq_23, edge_sq_34, edge_sq_41

   real, dimension (4) :: &
      x, y, z

!  Execution:

   dsq_neglect = (tol * ref_length)**2

!  Count the number of input quad. cells so we can allocate storage for the
!  triangulation that is at least large enough.  (We avoid the two passes that
!  precise allocation would require.)

   ntrimax = 0;  npoints = 0
   do ib = 1, npatches
      ni = patches(ib)%ni;  nj = patches(ib)%nj
      npoints = npoints + ni * nj
      ntrimax = ntrimax + (ni - 1) * (nj - 1)
   end do
   ntrimax = 2 * ntrimax;  nnode = npoints

   allocate (tri_xyz(3,npoints), conn(3,ntrimax))

   if (nf > 0) allocate (tri_f(nf,ntrimax))

!  Process one patch (block) at a time:

   ip = 0  ! Point number
   it = 0  ! Triangle number

   do ib = 1, npatches

      ni = patches(ib)%ni;  nj = patches(ib)%nj

      do j = 1, nj - 1

         do i = 1, ni - 1

            ip = ip + 1;  i1 = ip  ! (i,j) for patch ib
            i2 = i1 + 1            ! (i+1,j)
            i3 = i2 + ni           ! (i+1,j+1)
            i4 = i3 - 1            ! (i+1,j)

            x(1) = patches(ib)%x(i,j,1);    tri_xyz(1,ip) = x(1)
            y(1) = patches(ib)%y(i,j,1);    tri_xyz(2,ip) = y(1)
            z(1) = patches(ib)%z(i,j,1);    tri_xyz(3,ip) = z(1)
            x(2) = patches(ib)%x(i+1,j,1) 
            y(2) = patches(ib)%y(i+1,j,1)
            z(2) = patches(ib)%z(i+1,j,1)
            x(3) = patches(ib)%x(i+1,j+1,1) 
            y(3) = patches(ib)%y(i+1,j+1,1)
            z(3) = patches(ib)%z(i+1,j+1,1)
            x(4) = patches(ib)%x(i,j+1,1) 
            y(4) = patches(ib)%y(i,j+1,1)
            z(4) = patches(ib)%z(i,j+1,1)

            edge_sq_12 = (x(1) - x(2))**2 + (y(1) - y(2))**2 + (z(1) - z(2))**2
            edge_sq_23 = (x(2) - x(3))**2 + (y(2) - y(3))**2 + (z(2) - z(3))**2
            edge_sq_34 = (x(3) - x(4))**2 + (y(3) - y(4))**2 + (z(3) - z(4))**2
            edge_sq_41 = (x(4) - x(1))**2 + (y(4) - y(1))**2 + (z(4) - z(1))**2

            diag_sq_13 = (x(1) - x(3))**2 + (y(1) - y(3))**2 + (z(1) - z(3))**2
            diag_sq_24 = (x(2) - x(4))**2 + (y(2) - y(4))**2 + (z(2) - z(4))**2

            if (diag_sq_13 < diag_sq_24) then

               if (edge_sq_12 > dsq_neglect .and. edge_sq_23 > dsq_neglect) then
                  it = it + 1
                  conn(1,it) = i2
                  conn(2,it) = i3
                  conn(3,it) = i1
               end if

               if (edge_sq_34 > dsq_neglect .and. edge_sq_41 > dsq_neglect) then
                  it = it + 1
                  conn(1,it) = i4
                  conn(2,it) = i1
                  conn(3,it) = i3
               end if

            else ! Diagonal 2-4 is shorter

               if (edge_sq_12 > dsq_neglect .and. edge_sq_41 > dsq_neglect) then
                  it = it + 1
                  conn(1,it) = i1
                  conn(2,it) = i2
                  conn(3,it) = i4
               end if

               if (edge_sq_23 > dsq_neglect .and. edge_sq_34 > dsq_neglect) then
                  it = it + 1
                  conn(1,it) = i3
                  conn(2,it) = i4
                  conn(3,it) = i2
               end if

            end if

            if (nf > 0) tri_f(:,ip) = patches(ib)%q(:,i,j,1)

         end do ! Next i

         ip = ip + 1
         tri_xyz(1,ip) = patches(ib)%x(ni,j,1)
         tri_xyz(2,ip) = patches(ib)%y(ni,j,1)
         tri_xyz(3,ip) = patches(ib)%z(ni,j,1)

         if (nf > 0) tri_f(:,ip) = patches(ib)%q(:,ni,j,1)

      end do ! Next j

      j = ip
      do i = 1, ni
         ip = ip + 1
         tri_xyz(1,ip) = patches(ib)%x(i,nj,1)
         tri_xyz(2,ip) = patches(ib)%y(i,nj,1)
         tri_xyz(3,ip) = patches(ib)%z(i,nj,1)
      end do

      if (nf > 0) then
         ip = j
         do i = 1, ni
            ip = ip + 1
            tri_f(:,ip) = patches(ib)%q(:,i,nj,1)
         end do
      end if

   end do ! Next patch

!! nnode = ip  ! Should match npoints above
   ntri  = it

   end subroutine triangulate_patches
