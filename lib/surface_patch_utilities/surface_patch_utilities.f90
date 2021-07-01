!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module surface_patch_utilities

!  This module packages a few utilities for manipulating surface grid patches
!  and optional surface data.  It was prompted by a need to reorganize patches
!  into standard form for comparison of corresponding results.  Most of the
!  utilities treat one patch at a time, but a recent surface_connectivity addi-
!  tion to the collection opens further possibilities involving use of neighbor-
!  ing patch data.
!
!     centroids_areas   (patchin, patchout)    patchout%xyzf = centroids & areas
!     clone_grid        (np, nf, grid1, grid2) allocates and copies header info.
!     match_patch_edges (np, grid, tol)        averages common edge coordinates
!     patch_data_range  (patch)                assigns patch%xmin/xmax, etc.
!     reflect_x_y_or_z  (patch, n)             reflects about any coord. plane
!     reverse_patch_i   (patch, nf)            reverses the i indices
!     reverse_patch_j   (patch, nf)            reverses the j indices
!     scale_and_shift   (patch, scale, dx, dy, dz, pout)  scales/shifts x, y, z
!     surface_connectivity (np, grid, tol, icon)  sets up patch abutment info
!     surface_vectors   (np,grid,ip,i,j,tani,tanj,unorm) -> unit norm @ ip(i,j)
!     swap_coordinates  (patch, m, n)          swaps coordinates m and n
!     transpose_patch   (patch, nf)            swaps i and j in place
!     update_patch      (patch1, patch2, nf)   transfers patch1 to patch2
!     ????                                     what else?
!
!  04/09/05  D. A. Saunders   Initial implementation for Shuttle cavity grids.
!  06/16/05     "      "      Added scale and/or shift option.
!  06/20/05     "      "      Added calculation of cell centroids and areas.
!  02/22/06     "      "      Premature STOPs on Steve Alter's Opteron were
!                             traced to an ier /= 0 test in centroids_areas that
!                             should not have been there.
!  03/31/09     "      "      Added update_patch option.
!  07/26/10     "      "      %mi, %mj, %mk are now set in centroids_areas
!                             in case of I/O via xyzq_io.
!  09/18/15     "      "      Started adding a surface_vectors utility that
!                             produces unique results along patch edges for unit
!                             normal vectors (at grid points only, for now).
!  09/22/15     "      "      Started adding a surface_connectivity utility that
!                             is needed by surface_vectors.
!  09/28/15     "      "      Started adding a match_patch_edges utility, which
!                             is really needed first for surface_vectors to do
!                             what it promises (unique results at edge points).
!  10/06/15     "      "      The surface_vectors utility appears to be working
!                             in a version of SURFACE_CURVATURE, including
!                             averaged unit normals at corner points.  We still
!                             don't have more than a stub for match_patch_edges.
!  11/30/15     "      "      Completed match_patch_edges with match_corners and
!                             y0_edges.
!  05/19/15     "      "      Two glitches were detected when the CONVEX_GRID
!                             application made use of surface_vectors:
!                             patches reflected about y = 0 DO need to be
!                             rectified, and "if (jev == 2)" had iev, not jev.
!  06/21/21     "      "      The cloning was missing index (ib).
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!           Now with AMA, Inc., at NASA Ames Research Center.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure ! From Tecplot_io module, or equivalent

   implicit none  ! Applies to all procedures

!  Internal working variables:

   integer, allocatable, dimension (:,:,:), save, private :: &
      icon        ! Connectivity data for wgrid in match_patch_edges and
                  ! surface_vectors

   type (grid_type), allocatable, dimension (:), save, private :: &
      wgrid       ! Internal working grid employed by match_patch_edges and
                  ! surface_vectors

!  Public variables:

   logical, save, public :: first_unit_normal  ! Allows serial access to
                                               ! 2 or more surfaces
!  Public procedures:

   public :: centroids_areas
   public :: clone_grid
   public :: match_patch_edges
   public :: patch_data_range
   public :: reflect_x_y_or_z
   public :: reverse_patch_i
   public :: reverse_patch_j
   public :: scale_and_shift
   public :: surface_connectivity
   public :: surface_vectors
   public :: swap_coordinates
   public :: transpose_patch
   public :: update_patch

!  Private procedures:

   private :: dsq  ! Squared distance between two points

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine centroids_areas (patchin, patchout)

!     For the given surface patch and a companion "patch" with none of its
!     fields allocated yet, assign companion dimensions to represent the number
!     of cells rather than vertices. Allocate new %x, %y, %z fields and a single
!     function field and return them with the cell centroids and areas.
!
!     [Use of the Tecplot_io package with its grid_block_structure module does
!     not permit adding further fields to that derived data type, so the work-
!     around is to declare a second array of this data type at the higher level
!     and construct ancillary fields in the companion array elements.]
!
!     Centroids are simply the average coordinates of the vertices.
!
!     Areas are half the magnitudes of the cross-products of the two diagonals
!     of each quadrilateral cell, ASSUMING THE CELLS ARE (ESSENTIALLY) PLANAR.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent (in)    :: patchin   ! %x,y,z and %ni,nj defined
      type (grid_type), intent (inout) :: patchout  ! Fields are assumed to be
                                                    ! unallocated on input;
                                                    ! returned with cell-related
                                                    ! quantities described above
!     Local constants:

      real, parameter :: fourth = 0.25, half = 0.5

!     Local variables:

      integer :: i, j, ni, nj
      real    :: p(3), q(3), r(3)

!     Execution:

      ni = patchin%ni - 1;  patchout%ni = ni;  patchout%mi = ni ! In case of I/O
      nj = patchin%nj - 1;  patchout%nj = nj;  patchout%mj = nj
                            patchout%nk = 1;   patchout%mk = 1

      allocate (patchout%x(ni,nj,1), patchout%y(ni,nj,1), &
                patchout%z(ni,nj,1), patchout%q(1,ni,nj,1))

      do j = 1, nj
         do i = 1, ni
            patchout%x(i,j,1) = fourth * &
               (patchin%x(i,  j,1) + patchin%x(i+1,  j,1) + &
                patchin%x(i,j+1,1) + patchin%x(i+1,j+1,1))
            patchout%y(i,j,1) = fourth * &
               (patchin%y(i,  j,1) + patchin%y(i+1,  j,1) + &
                patchin%y(i,j+1,1) + patchin%y(i+1,j+1,1))
            patchout%z(i,j,1) = fourth * &
               (patchin%z(i,  j,1) + patchin%z(i+1,  j,1) + &
                patchin%z(i,j+1,1) + patchin%z(i+1,j+1,1))

            p(1) = patchin%x(i+1,j+1,1) - patchin%x(i,j,1)
            p(2) = patchin%y(i+1,j+1,1) - patchin%y(i,j,1)
            p(3) = patchin%z(i+1,j+1,1) - patchin%z(i,j,1)

            q(1) = patchin%x(i,j+1,1) - patchin%x(i+1,j,1)
            q(2) = patchin%y(i,j+1,1) - patchin%y(i+1,j,1)
            q(3) = patchin%z(i,j+1,1) - patchin%z(i+1,j,1)

            r(1) = p(2) * q(3) - p(3) * q(2) ! Vector product
            r(2) = p(3) * q(1) - p(1) * q(3)
            r(3) = p(1) * q(2) - p(2) * q(1)
            
            patchout%q(1,i,j,1) = half * sqrt (r(1)**2 + r(2)**2 + r(3)**2)
         end do
      end do

      end subroutine centroids_areas

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine clone_grid (np, nf, grid1, grid2)

!     Set up a copy of a multiblock grid by allocating the blocks and
!     transcribing the header information (only).

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)          :: np      ! # surface patches
      integer, intent (in)          :: nf      ! nf > 0 means copy %mi, etc. too
      type (grid_type), intent (in) :: grid1(np)
      type (grid_type), pointer     :: grid2(:)

!     Local variables:

      integer :: ib

!     Execution:

      allocate (grid2(np))

      do ib = 1, np
         grid2(ib)%ni = grid1(ib)%ni
         grid2(ib)%nj = grid1(ib)%nj
         grid2(ib)%nk = grid1(ib)%nk
         if (nf > 0) then
            grid2(ib)%mi = grid1(ib)%mi
            grid2(ib)%mj = grid1(ib)%mj
            grid2(ib)%mk = grid1(ib)%mk
         end if
      end do

      end subroutine clone_grid

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real function dsq (x1, y1, z1, x2, y2, z2)  ! Obvious functionality
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, intent (in) :: x1, y1, z1, x2, y2, z2

      dsq = (x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2

      end function dsq

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine match_patch_edges (np, grid, eps)

!     For the given multipatch structured surface grid, construct connectivity
!     information in private internal workspace icon(2,4,np) and apply it to
!     match apparent common edges exactly.  Interior points of interior edges
!     are matched by averaging in pairs.  Corner points are then averaged over
!     all patches sharing each corner.  Finally, points apparently on a y = 0
!     symmetry plane are set to y = 0. exactly.  Other possible symmetry plane
!     points (x= 0., z = 0.) are NOT looked for at present.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer,          intent (in)    :: np        ! # surface patches
      type (grid_type), intent (inout) :: grid(np)  ! Surface grid, updated in-
                                                    ! place with matching edges,
                                                    ! and y = 0 exactly along
                                                    ! apparent symmetry plane
                                                    ! edges
      real,             intent (in)    :: eps       ! Tolerance used to detect
                                                    ! apparent edge matches,
                                                    ! applied to y data range
!     Local constants:

      real, parameter :: half = 0.5, zero = 0.

!     Local variables:

      integer :: i1, i2, iinc, j1, j2, jinc, m, me, mi, mj, n, ne, ni, nj
      integer :: ia, ib, ja, jb, lbuf
      integer, allocatable :: lcon(:,:,:)           ! Needed only locally
      real    :: tol, tolsq
      real,    allocatable :: bufm(:), bufn(:)      ! Shouldn't be needed
      logical :: reverse
      logical, allocatable :: edge_done(:,:)        ! For marking matching edges

!     Execution:

!     We could require determining icon(2,4,np) before calling this routine
!     but we choose to keep the matching self-contained instead, so we use
!     local lcon(2,4,np) to avoid possible clashes with surface_vectors usage.

      call set_up_connectivity ()  ! Internal procedure below

!     Treat interior edge points first (not corner points yet):

      do n = 1, np  ! For each patch n in the connectivity info
         ni = grid(n)%ni
         nj = grid(n)%nj
!!!      write (*, '(a, 3i5)') 'n,ni,nj:', n,ni,nj
         do ne = 1, 4  ! For each edge of patch n
            if (edge_done(ne,n)) cycle
            m = lcon(1,ne,n)   ! m > 0 if a patch m edge matches this edge
!!!         write (*, '(a, 2i5)') '   ne,m:', ne,m
            if (m <= 0) cycle
            select case (ne)   ! Define index range for this edge of patch n
               case (1)
                  i1 = 1;   i2 = 1;   j1 = 2;  j2 = nj - 1;  lbuf = nj - 2
               case (2)
                  i1 = ni;  i2 = ni;  j1 = 2;  j2 = nj - 1;  lbuf = nj - 2
               case (3)
                  i1 = 2;   i2 = ni - 1;  j1 = 1;   j2 = 1;  lbuf = ni - 2
               case (4)
                  i1 = 2;   i2 = ni - 1;  j1 = nj;  j2 = nj; lbuf = ni - 2
            end select
!!!         write (*, '(a, 4i5)') '   i1,i2,j1,j2:', i1,i2,j1,j2
            mi = grid(m)%ni
            mj = grid(m)%nj
            me = lcon(2,ne,n)  ! Edge me of patch m matches edge ne of patch n
!!!         write (*, '(a, 3i5)') '   me,mi,mj:', me,mi,mj
            reverse = me < 0
            me = abs (me)
            iinc = 1;  jinc = 1
            select case (me)   ! Define index range for matching patch m
               case (1)
                  ia = 1;   ib = 1;   ja = 2;  jb = mj - 1
                  if (reverse) then
                     ja = jb;  jb = 2;  jinc = -1
                  end if
               case (2)
                  ia = mi;  ib = mi;  ja = 2;  jb = mj - 1
                  if (reverse) then
                     ja = jb;  jb = 2;  jinc = -1
                  end if
               case (3)
                  ia = 2;   ib = mi - 1;  ja = 1;   jb = 1
                  if (reverse) then
                     ia = ib;  ib = 2;  iinc = -1
                  end if
               case (4)
                  ia = 2;   ib = mi - 1;  ja = mj;  jb = mj
                  if (reverse) then
                     ia = ib;  ib = 2;  iinc = -1
                  end if
            end select

!!!         write (*, '(a, 6i5)') &
!!!            '   ia,ib,iinc,ja,jb,jinc:', ia,ib,iinc,ja,jb,jinc

            allocate (bufn(lbuf), bufm(lbuf))

!           The following tells all, but violates matching shape rules:

!!!         grid(n)%x(i1:i2,j1:j2,1) = (grid(n)%x(i1:i2,j1:j2,1) + &
!!!                                     grid(m)%x(ia:ib:iinc,ja:jb:jinc,1))*half
!!!         grid(m)%x(ia:ib:iinc,ja:jb:jinc,1) = grid(n)%x(i1:i2,j1:j2,1)
!!!         grid(n)%y(i1:i2,j1:j2,1) = (grid(n)%y(i1:i2,j1:j2,1) + &
!!!                                     grid(m)%y(ia:ib:iinc,ja:jb:jinc,1))*half
!!!         grid(m)%y(ia:ib:iinc,ja:jb:jinc,1) = grid(n)%y(i1:i2,j1:j2,1)
!!!         grid(n)%z(i1:i2,j1:j2,1) = (grid(n)%z(i1:i2,j1:j2,1) + &
!!!                                     grid(m)%z(ia:ib:iinc,ja:jb:jinc,1))*half
!!!         grid(m)%z(ia:ib:iinc,ja:jb:jinc,1) = grid(n)%z(i1:i2,j1:j2,1)

            call a2v (ni, nj, grid(n)%x, i1, i2, 1,    j1, j2, 1,    lbuf, bufn)
            call a2v (mi, mj, grid(m)%x, ia, ib, iinc, ja, jb, jinc, lbuf, bufm)
            bufn(:) = half*(bufn(:) + bufm(:))
            call v2a (lbuf, bufn, ni, nj, grid(n)%x, i1, i2, 1,    j1, j2, 1)
            call v2a (lbuf, bufn, mi, mj, grid(m)%x, ia, ib, iinc, ja, jb, jinc)

            call a2v (ni, nj, grid(n)%y, i1, i2, 1,    j1, j2, 1,    lbuf, bufn)
            call a2v (mi, mj, grid(m)%y, ia, ib, iinc, ja, jb, jinc, lbuf, bufm)
            bufn(:) = half*(bufn(:) + bufm(:))
            call v2a (lbuf, bufn, ni, nj, grid(n)%y, i1, i2, 1,    j1, j2, 1)
            call v2a (lbuf, bufn, mi, mj, grid(m)%y, ia, ib, iinc, ja, jb, jinc)

            call a2v (ni, nj, grid(n)%z, i1, i2, 1,    j1, j2, 1,    lbuf, bufn)
            call a2v (mi, mj, grid(m)%z, ia, ib, iinc, ja, jb, jinc, lbuf, bufm)
            bufn(:) = half*(bufn(:) + bufm(:))
            call v2a (lbuf, bufn, ni, nj, grid(n)%z, i1, i2, 1,    j1, j2, 1)
            call v2a (lbuf, bufn, mi, mj, grid(m)%z, ia, ib, iinc, ja, jb, jinc)

            deallocate (bufn, bufm)
            edge_done(ne,n) = .true.
            edge_done(me,m) = .true.
         end do  ! Next edge of patch n
      end do  ! Next patch n

!     Treat corner points:

      if (np > 1) call match_corners ()  ! Internal procedure

!     Treat y = 0 symmetry plane points:

      call y0_edges ()  ! Internal procedure

      deallocate (lcon, edge_done)

!     Internal procedures for subroutine match_patch_edges:

      contains

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine set_up_connectivity ()  ! Determine a sensible tolerance and
                                            ! set up the patch connectivity info

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        Local variables:

         integer :: ip
         real :: ymin, ymax

!        Execution:

         do ip = 1, np
            call patch_data_range (grid(ip))
         end do

         ymin  = minval (grid(:)%ymin)
         ymax  = maxval (grid(:)%ymax)
         tol   = eps * (ymax - ymin)
         tolsq = tol*tol
         write (*, '(a, 3es16.8)') ' ymin, ymax, tol:', ymin, ymax, tol

         allocate (lcon(2,4,np), edge_done(4,np))

         call surface_connectivity (np, grid, tol, lcon)

         edge_done(:,:) = .false.

         end subroutine set_up_connectivity

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine a2v (ni, nj, a, i1, i2, iinc, j1, j2, jinc, lv, v)

!        Copy a(i1:i2:iinc,j1:j2:jinc,1) to v(1:lv).
!        [Array assignments require the shapes to match.]
!        The indicated numbers of elements are assumed to match.
!        For the usages here, one of the subscript ranges of a(:,:,1) is always
!        of size 1.

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        Arguments:

         integer, intent (in)  :: ni, nj       ! Dimensions of surface array a
         real,    intent (in)  :: a(ni,nj,1)   ! May be passed as surf%x, say
         integer, intent (in)  :: i1, i2, iinc ! Ranges of a(:,:,1) to copy
         integer, intent (in)  :: j1, j2, jinc
         integer, intent (in)  :: lv           ! v(1:lv) is the active portion
         real,    intent (out) :: v(lv)        ! Output vector

!        Local variables:

         integer :: i, j, k

!        Execution:

         k = 0
         do j = j1, j2, jinc
            do i = i1, i2, iinc
               k = k + 1
               v(k) = a(i,j,1)
            end do
         end do

         end subroutine a2v

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine v2a (lv, v, ni, nj, a, i1, i2, iinc, j1, j2, jinc)

!        Copy v(1:lv) to a(i1:i2:iinc,j1:j2:jinc,1).
!        [Array assignments require the shapes to match.]
!        The indicated numbers of elements are assumed to match.
!        For the usages here, one of the subscript ranges of a(:,:,1) is always
!        of size 1.

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        Arguments:

         integer, intent (in) :: lv           ! v(1:lv) is the active portion
         real,    intent (in) :: v(lv)        ! Input vector
         integer, intent (in) :: ni, nj       ! Dimensions of 2D array a
         real,    intent (inout) :: a(ni,nj,1)! Part of a(:,:,1) is updated, ...
         integer, intent (in) :: i1, i2, iinc ! ... as indicated by these ranges
         integer, intent (in) :: j1, j2, jinc

!        Local variables:

         integer :: i, j, k

!        Execution:

         k = 0
         do j = j1, j2, jinc
            do i = i1, i2, iinc
               k = k + 1
               a(i,j,1) = v(k)
            end do
         end do

         end subroutine v2a

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine match_corners ()  ! Average corner coords. over apparently
                                      ! common corners

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        Local variables:

         integer :: i, ic, ii, j, jc, jj, m, mi, mj, n, ni, nj, npatch, nsum
         real    :: dsqrd, xc, yc, zc
         real, allocatable :: xyz_average(:,:,:), xyz_corner(:,:,:)

!        Execution:

         allocate (xyz_corner(3,4,np))

!        Store all the unperturbed coordinates of corner points:

         do n = 1, np
            ni = grid(n)%ni;  nj = grid(n)%nj
            xyz_corner(1,1,n) = grid(n)%x( 1, 1,1)
            xyz_corner(2,1,n) = grid(n)%y( 1, 1,1)
            xyz_corner(3,1,n) = grid(n)%z( 1, 1,1)
            xyz_corner(1,2,n) = grid(n)%x(ni, 1,1)
            xyz_corner(2,2,n) = grid(n)%y(ni, 1,1)
            xyz_corner(3,2,n) = grid(n)%z(ni, 1,1)
            xyz_corner(1,3,n) = grid(n)%x( 1,nj,1)
            xyz_corner(2,3,n) = grid(n)%y( 1,nj,1)
            xyz_corner(3,3,n) = grid(n)%z( 1,nj,1)
            xyz_corner(1,4,n) = grid(n)%x(ni,nj,1)
            xyz_corner(2,4,n) = grid(n)%y(ni,nj,1)
            xyz_corner(3,4,n) = grid(n)%z(ni,nj,1)
         end do

         allocate (xyz_average(3,4,np))

!        For each patch, check its four corners against all other corners:

         do npatch = 1, np
            ni = grid(npatch)%ni;  nj = grid(npatch)%nj
            do ic = 1, 4
               select case (ic)
                  case (1)
                     i = 1;  j = 1
                  case (2)
                     i = ni
                  case (3)
                     i = 1;  j = nj
                  case (4)
                     i = ni
               end select
               nsum = 1
               xc = xyz_corner(1,ic,npatch);  xyz_average(1,ic,npatch) = xc
               yc = xyz_corner(2,ic,npatch);  xyz_average(2,ic,npatch) = yc
               zc = xyz_corner(3,ic,npatch);  xyz_average(3,ic,npatch) = zc
               do n = 1, np
                  if (n == npatch) cycle

                  mi = grid(n)%ni;  mj = grid(n)%nj
                  do jc = 1, 4
                     select case (jc)
                        case (1)
                           ii = 1;  jj = 1
                        case (2)
                           ii = mi
                        case (3)
                           ii = 1;  jj = mj
                        case (4)
                           ii = mi
                     end select
                     dsqrd = dsq (grid(n)%x(ii,jj,1), &
                                  grid(n)%y(ii,jj,1), &
                                  grid(n)%z(ii,jj,1), xc, yc, zc)
                     if (dsqrd < tolsq) then
                        nsum = nsum + 1
                        xyz_average(1,ic,npatch) = xyz_average(1,ic,npatch) + &
                                              grid(n)%x(ii,jj,1)
                        xyz_average(2,ic,npatch) = xyz_average(2,ic,npatch) + &
                                              grid(n)%y(ii,jj,1)
                        xyz_average(3,ic,npatch) = xyz_average(3,ic,npatch) + &
                                              grid(n)%z(ii,jj,1)
                     end if
                  end do  ! Next corner
               end do  ! Next patch
               xyz_average(:,ic,npatch) = xyz_average(:,ic,npatch) / real (nsum)
            end do  ! Next corner
         end do  ! Next patch

!        Update corner coordinates with averaged coordinates (or originals):

         do npatch = 1, np
            ni = grid(npatch)%ni;  nj = grid(npatch)%nj
            do ic = 1, 4
               select case (ic)
                  case (1)
                     i = 1;  j = 1
                  case (2)
                     i = ni
                  case (3)
                     i = 1;  j = nj
                  case (4)
                     i = ni
               end select
               grid(npatch)%x(i,j,1) = xyz_average(1,ic,npatch)
               grid(npatch)%y(i,j,1) = xyz_average(2,ic,npatch)
               grid(npatch)%z(i,j,1) = xyz_average(3,ic,npatch)
            end do  ! Next corner
         end do  ! Next patch

         deallocate (xyz_average, xyz_corner)

         end subroutine match_corners

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine y0_edges ()  ! Force y = 0. exactly on apparent symmetry
                                 ! plane edges.

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        Local variables:

         integer :: ie, n, ni, nj
         real    :: ymax, ymin

!        Execution:

         do n = 1, np
            ni = grid(n)%ni;  nj = grid(n)%nj
            do  ie = 1, 4
               select case (ie)
                  case (1)
                     ymin = minval (grid(n)%y(1,:,1))
                     ymax = maxval (grid(n)%y(1,:,1))
                     if (abs (ymin) < tol .and. abs (ymax) < tol) &
                        grid(n)%y(1,:,1) = zero
                  case (2)
                     ymin = minval (grid(n)%y(ni,:,1))
                     ymax = maxval (grid(n)%y(ni,:,1))
                     if (abs (ymin) < tol .and. abs (ymax) < tol) &
                        grid(n)%y(ni,:,1) = zero
                  case (3)
                     ymin = minval (grid(n)%y(:,1,1))
                     ymax = maxval (grid(n)%y(:,1,1))
                     if (abs (ymin) < tol .and. abs (ymax) < tol) &
                        grid(n)%y(:,1,1) = zero
                  case (4)
                     ymin = minval (grid(n)%y(:,nj,1))
                     ymax = maxval (grid(n)%y(:,nj,1))
                     if (abs (ymin) < tol .and. abs (ymax) < tol) &
                        grid(n)%y(:,nj,1) = zero
               end select
            end do
         end do

         end subroutine y0_edges

      end subroutine match_patch_edges

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine patch_data_range (patch)

!     Assign min & max values of x, y, z for the given surface patch (nk = 1).

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent (inout) :: patch

!     Local constants:

      real, parameter :: big = 1.e+32

!     Local variables:

      integer :: i, j

!     Execution:

      patch%xmin = big;  patch%xmax = -big
      patch%ymin = big;  patch%ymax = -big
      patch%zmin = big;  patch%zmax = -big

      do j = 1, patch%nj
         do i = 1, patch%ni
            patch%xmin = min (patch%x(i,j,1), patch%xmin)
            patch%xmax = max (patch%x(i,j,1), patch%xmax)
            patch%ymin = min (patch%y(i,j,1), patch%ymin)
            patch%ymax = max (patch%y(i,j,1), patch%ymax)
            patch%zmin = min (patch%z(i,j,1), patch%zmin)
            patch%zmax = max (patch%z(i,j,1), patch%zmax)
         end do
      end do

      end subroutine patch_data_range

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine reflect_x_y_or_z (patch, n)

!     Reflect the indicated coordinate in the appropriate plane through (0,0,0).

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent (inout) :: patch
      integer, intent (in) :: n           ! n = 1, 2, 3  <->  reflect x, y, z

!     Local variables:

      integer :: i, j

!     Execution:

      select case (n)
         case (1)
            patch%x = -patch%x
         case (2)
            patch%y = -patch%y
         case (3)
            patch%z = -patch%z
      end select

      end subroutine reflect_x_y_or_z

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine reverse_patch_i (patch, nf)

!     Reverse the i indices of a surface patch, in-place.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent (inout) :: patch
      integer,          intent (in)    :: nf   ! # function values; nf >= 0

!     Local variables:

      integer :: i, ii, j, ni, niby2
      real    :: t, tq(nf)

!     Execution:

      ni = patch%ni;  niby2 = (ni + 1) / 2

      do j = 1, patch%nj
         ii = ni
         do i = 1, niby2
            t               = patch%x(i, j,1)
            patch%x(i, j,1) = patch%x(ii,j,1)
            patch%x(ii,j,1) = t
            t               = patch%y(i, j,1)
            patch%y(i, j,1) = patch%y(ii,j,1)
            patch%y(ii,j,1) = t
            t               = patch%z(i, j,1)
            patch%z(i, j,1) = patch%z(ii,j,1)
            patch%z(ii,j,1) = t
            ii = ii - 1 
         end do
      end do

      if (nf > 0) then

         ni = patch%mi;  niby2 = (ni + 1) / 2

         do j = 1, patch%mj
            ii = ni
            do i = 1, niby2
               tq                = patch%q(:,i, j,1)
               patch%q(:,i, j,1) = patch%q(:,ii,j,1)
               patch%q(:,ii,j,1) = tq
               ii = ii - 1  
            end do  
         end do

      end if

      end subroutine reverse_patch_i

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine reverse_patch_j (patch, nf)

!     Reverse the j indices of a surface patch, in-place.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent (inout) :: patch
      integer,          intent (in)    :: nf   ! # function values; nf >= 0

!     Local variables:

      integer :: i, j, jj, nj, njby2
      real    :: t, tq(nf)

!     Execution:

      nj = patch%nj;  njby2 = (nj + 1) / 2

      do i = 1, patch%ni
         jj = nj
         do j = 1, njby2
            t               = patch%x(i,j, 1)
            patch%x(i,j, 1) = patch%x(i,jj,1)
            patch%x(i,jj,1) = t
            t               = patch%y(i,j, 1)
            patch%y(i,j, 1) = patch%y(i,jj,1)
            patch%y(i,jj,1) = t
            t               = patch%z(i,j, 1)
            patch%z(i,j, 1) = patch%z(i,jj,1)
            patch%z(i,jj,1) = t
            jj = jj - 1
         end do
      end do

      if (nf > 0) then

         nj = patch%mj;  njby2 = (nj + 1) / 2

         do i = 1, patch%mi
            jj = nj
            do j = 1, njby2
               tq                = patch%q(:,i,j, 1)
               patch%q(:,i,j, 1) = patch%q(:,i,jj,1)
               patch%q(:,i,jj,1) = tq
               jj = jj - 1
            end do
         end do

      end if

      end subroutine reverse_patch_j

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine scale_and_shift (patch, scale, dx, dy, dz, patch_out) 

!     Scale all of x, y, z and/or shift any of x, y, z, in-place or not:
!
!                  x_out <-- scale * x + dx          etc.
!
!     Unit scaling and zero shifts are not suppressed.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent (inout) :: patch
      real, intent (in)                :: scale, dx, dy, dz
      type (grid_type), intent (inout) :: patch_out  ! May be the input patch

!     Local variables:

      integer :: i, j

!     Execution:

      do j = 1, patch%nj
         do i = 1, patch%ni
            patch_out%x(i,j,1) = scale * patch%x(i,j,1) + dx
            patch_out%y(i,j,1) = scale * patch%y(i,j,1) + dy
            patch_out%z(i,j,1) = scale * patch%z(i,j,1) + dz
         end do
      end do

      end subroutine scale_and_shift

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine surface_connectivity (np, grid, tol, icon)

!     Description:
!
!     For the given multipatch surface grid in 3-space, determine connectivity
!     data that will allow utilities such as surface_vectors to make use of
!     neighboring patch information as opposed to being confined to a single
!     patch the way surface_normal is.
!
!     For now, edges along symmetry planes are not treated because surface_-
!     vectors automatically reflects the grid if a y = 0 edge is detected.
!     If other uses for this utility arise, symmetry planes could be handled by
!     entering say -1/-2/-3 in icon(1,ie,ip) to mean an x/y/z = 0. edge
!     respectively.
!
!     For each patch of the input grid, the following information suffices to
!     determine access to the correct points of neighboring patches:
!
!        o  edges 1, 2, 3, 4 refer to imin, imax, jmin, jmax respectively
!        o  -1, -2, -3, -4 indicate indexing in the reverse direction
!        o  icon(1:2,1:4,1:np) is set up as in these three examples:
!
!        1: For edge point (1,j) of patch ip, icon(1:2,1,ip) = 6, 3 means that
!           patch 6 abuts edge 1 (imin) with edge 3 (jmin) where j is increasing
!           as for patch ip.
!
!        2: For edge pt. (i,jmax) of patch ip, icon(1:2,4,ip) = 9, -2 means that
!           patch 9 abuts edge 4 (jmax) with edge 2 (imax) where i is indexed in
!           the opposite direction from patch ip's jmax edge 4.
!
!        3: For edge point (i,1) of patch ip, icon(1:2,3,ip) = 0, 0 means that
!           its jmin edge 3 does not appear to abut an edge of another patch.
!
!     Common edges are detected by having common data ranges or bounding boxes
!     followed by simple tests to determine indexing direction.
!
!     A certain 10-patch sphere produces connectivity info like this:
!
!                  ip1 ie1 ip2 ie2 ip3 ie3 ip4 ie4
!                    6 -1    2 -3    3 -3    5 -1      (patch 1 info)
!                    5  3    3  1    1 -2    4 -1
!                    2  2    8 -2    1 -3    4  3         "   3   "
!                    2 -4    9 -2    3  4    5  2
!                    1 -4    4  4    2  1   10  3
!                    1 -1    7  4   10  1    8 -4
!                   10 -4    8  1    9  1    6  2
!                    7  2    3 -2    9  4    6 -4         "   8   "
!                    7  3    4 -2   10 -2    8  3
!                    6  3    9 -3    5  4    7 -1         "  10   "
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)          :: np       ! Number of surface patches
      type (grid_type), intent (in) :: grid(np) ! Structured surface patches
      real,    intent (in)          :: tol      ! Distance tolerance for a match
      integer, intent (out)         :: icon(2,4,np)  ! See usage above

!     Local variables:

      integer :: ie, ig, ip, iq, iip, iiq, jip, jiq, ni, nj
      real    :: tolsq
      real    :: drange(2,3,4,np)  ! xmin/max, ymin/max, zmin/max for edges 1:4

!     Execution:

!     Edge data ranges or edge bounding boxes for all patches:

      do ip = 1, np
         ni = grid(ip)%ni
         nj = grid(ip)%nj
         drange(1,1,1,ip) = minval (grid(ip)%x( 1,:,1))
         drange(2,1,1,ip) = maxval (grid(ip)%x( 1,:,1))
         drange(1,1,2,ip) = minval (grid(ip)%x(ni,:,1))
         drange(2,1,2,ip) = maxval (grid(ip)%x(ni,:,1))
         drange(1,1,3,ip) = minval (grid(ip)%x(:, 1,1))
         drange(2,1,3,ip) = maxval (grid(ip)%x(:, 1,1))
         drange(1,1,4,ip) = minval (grid(ip)%x(:,nj,1))
         drange(2,1,4,ip) = maxval (grid(ip)%x(:,nj,1))

         drange(1,2,1,ip) = minval (grid(ip)%y( 1,:,1))
         drange(2,2,1,ip) = maxval (grid(ip)%y( 1,:,1))
         drange(1,2,2,ip) = minval (grid(ip)%y(ni,:,1))
         drange(2,2,2,ip) = maxval (grid(ip)%y(ni,:,1))
         drange(1,2,3,ip) = minval (grid(ip)%y(:, 1,1))
         drange(2,2,3,ip) = maxval (grid(ip)%y(:, 1,1))
         drange(1,2,4,ip) = minval (grid(ip)%y(:,nj,1))
         drange(2,2,4,ip) = maxval (grid(ip)%y(:,nj,1))

         drange(1,3,1,ip) = minval (grid(ip)%z( 1,:,1))
         drange(2,3,1,ip) = maxval (grid(ip)%z( 1,:,1))
         drange(1,3,2,ip) = minval (grid(ip)%z(ni,:,1))
         drange(2,3,2,ip) = maxval (grid(ip)%z(ni,:,1))
         drange(1,3,3,ip) = minval (grid(ip)%z(:, 1,1))
         drange(2,3,3,ip) = maxval (grid(ip)%z(:, 1,1))
         drange(1,3,4,ip) = minval (grid(ip)%z(:,nj,1))
         drange(2,3,4,ip) = maxval (grid(ip)%z(:,nj,1))
      end do

      tolsq = tol*tol  ! No need to do the square roots
      icon  = 0        ! Edges with no neighbors will stay marked with 0s

!     For each edge of each patch, search for a matching edge:

      do ip = 1, np  ! For each patch ip
         do iq = 1, np  ! Check each patch iq
            if (ip == iq) cycle
            do ie = 1, 4  ! For each edge of patch ip
               do ig = 1, 4  ! Compare with each edge of patch iq
!!!               write (40, '(a, 4i3, 2es16.8)') 'ip,iq,ie,ig,dxmin,dxmax:', &
!!!                  ip,iq,ie,ig, abs (drange(1,1,ie,ip) - drange(1,1,ig,iq)),&
!!!                               abs (drange(2,1,ie,ip) - drange(2,1,ig,iq))
!!!               write (40, '(a, 12x, 2es16.8)') '            dymin,dymax:', &
!!!                               abs (drange(1,2,ie,ip) - drange(1,2,ig,iq)),&
!!!                               abs (drange(2,2,ie,ip) - drange(2,2,ig,iq))
!!!               write (40, '(a, 12x, 2es16.8)') '            dzmin,dzmax:', &
!!!                               abs (drange(1,3,ie,ip) - drange(1,3,ig,iq)),&
!!!                               abs (drange(2,3,ie,ip) - drange(2,3,ig,iq))
                  if (abs (drange(1,1,ie,ip) - drange(1,1,ig,iq)) < tol .and. &
                      abs (drange(2,1,ie,ip) - drange(2,1,ig,iq)) < tol .and. &
                      abs (drange(1,2,ie,ip) - drange(1,2,ig,iq)) < tol .and. &
                      abs (drange(2,2,ie,ip) - drange(2,2,ig,iq)) < tol .and. &
                      abs (drange(1,3,ie,ip) - drange(1,3,ig,iq)) < tol .and. &
                      abs (drange(2,3,ie,ip) - drange(2,3,ig,iq)) < tol) then
                         icon(1,ie,ip) = iq
                         icon(2,ie,ip) = ig
                  end if
               end do
            end do
         end do
      end do

!     Check for opposite indexing directions:

      do ip = 1, np
         do ie = 1, 4
            iq = icon(1,ie,ip)
            if (iq /= 0) then
               ig = icon(2,ie,ip)
               select case (ig)  ! Pick left end of relevant edge
                  case (1)
                     iiq = 1;  jiq = 1
                  case (2)
                     iiq = grid(iq)%ni;  jiq = 1
                  case (3)
                     iiq = 1;  jiq = 1
                  case (4)
                     iiq = 1;  jiq = grid(iq)%nj
               end select
               select case (ie)  ! Pick left end of relevant edge
                  case (1)
                     iip = 1;  jip = 1
                  case (2)
                     iip = grid(ip)%ni;  jip = 1
                  case (3)
                     iip = 1;  jip = 1
                  case (4)
                     iip = 1;  jip = grid(ip)%nj
               end select
               if (dsq (grid(ip)%x(iip,jip,1), grid(ip)%y(iip,jip,1), &
                        grid(ip)%z(iip,jip,1), grid(iq)%x(iiq,jiq,1), &
                        grid(iq)%y(iiq,jiq,1), grid(iq)%z(iiq,jiq,1)) > tolsq) &
                  icon(2,ie,ip) = -ig
            end if
         end do
      end do

      write (*, '(a)') '  ip1 ie1 ip2 ie2 ip3 ie3 ip4 ie4    ip'
      write (*, '(i5, i3, i5, i3, i5, i3, i5, i3, i7)') &
         (icon(:,:,ip), ip, ip = 1, np)

      end subroutine surface_connectivity

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine surface_vectors (np, grid, ip, i, j, tani, tanj, unorm)

!     Description:
!
!     This variation of subroutine SURFACE_NORMAL is intended to overcome the
!     weaknesses that arise at surface patch edge points when information from
!     neighboring patches is not made use of.
!
!     For grid point (i,j) of patch ip in a structured surface grid, construct
!     a 3-point row of (x,y,z) points and a 3-point column of (x,y,z) points
!     that are the best possible for determining surface tangent unit vectors in
!     the two index directions, along with a unit normal vector at that (i,j).
!     These should be such that if (i,j) is on an interior edge, the same output
!     vectors are obtained when the common point from a neighboring patch is
!     input. Further, if (i,j) is on an exterior edge and that edge appears to
!     be in a symmetry plane, the outputs correctly take advantage of the
!     symmetry.
!
!     The 3-point line segments suffice to determine partial derivatives at the
!     center point (or at an end-point on an exterior edge that is not also in a
!     symmetry plane) via finite differencing, and these can be used to produce
!     the desired unit tangent vectors and/or unit normal vector, all three of
!     which are returned.
!
!     Note that the specified-patch-only subroutine SURFACE_NORMAL is more
!     general in the sense that it provides a result for any location on the
!     surface, not necessarily a grid point, but it is doubtful that trying to
!     do that here is worth the effort.
!
!     Method:
!
!     This is trivial for an interior point of a patch.  The quality of the
!     outputs may degrade with skewness in the patch cells, but the results
!     from the basic 3-point method in each index direction are unique.  It is
!     the edge points that require careful handling to give identical results
!     for different calls targeting common edge points from other patches.
!
!     Edge Point Issues:
!
!        o  Internal edge/non-corner point results are also unique if we
!           correctly identify the neighboring patch for use by the basic
!           3-point method.
!        o  External edge/non-corner point results depend on whether the edge
!           is along a symmetry plane or not.  Decision: avoid the distinction
!           by looking for symmetry as a one-time first step, and if it is
!           present, reflect the given surface appropriately and make the
!           two-halves reflected grid the working grid, so such external edges
!           become internal edges and the basic internal edge method applies.
!        o  Patch corner points are the most problematic.  Any number of
!           patches may have a common corner: 2, 3, 4, 5, ..., and the only
!           bulletproof way to deal with that appears to be to compute estimates
!           from each patch sharing a corner and return the averaged results.
!           This should benefit corners on a symmetry plane where some patches
!           may be more skewed than others.  It's not clear what averaging means
!           for the tangent vectors, though, so only the last pair is returned.
!           They probably aren't needed anyway.  The unit normal is most likely
!           to be used.  The fissure problem mentioned in the 06/23/15 history
!           item is one example.
!        o  Once any symmetry is found or not, the internal working surface grid
!           can be scanned to establish its connectivity information.  This is
!           much simpler for surfaces than for multiblock volume grids.
!
!     Reasonable Assumptions (Not Checked For):
!
!     o  All surface patches are right-handed: unit normal vectors point out
!        of a convex body.  If this is not true, a RECTIFY_GRID utility is
!        available from the present author.
!     o  All surface patch dimensions are at least 2.
!     o  The patches are point-to-point matched to within a tiny tolerance.
!     o  The y = 0 plane is the only possible symmetry plane handled for now.
!
!     Programmer Notes:
!
!     o  Use of private global variables wgrid and first_unit_normal would
!        preclude application to more than one surface grid in a run.
!        Therefore, first_unit_normal needs to be public and assigned by the
!        application, and wgrid(:) may need deallocating before being allocated.
!        This arrangement still does not allow interleaved accesses to more than
!        one surface grid--only serial application.  If accesses to two surfaces
!        need to be interleaved, a kludge would be to copy the relevant parts of
!        this module and tweak the names of this procedure and of the public
!        logical first_unit_normal.
!
!     History:
!
!     06/20/15  D.A.Saunders  Initial design of nonlinear method for surface
!                             grid points (only), in subroutine form.
!     06/22/15    "     "     Initial implementation of nonlinear method,
!                             still not using neighboring patch data.  Testing
!                             via SURFACE_CURVATURE on a spherical surface
!                             shows better results than from SURFACE_NORMAL
!                             along patch edges, as hoped.
!     06/23/15    "     "     More severe testing of the preliminary version on
!                             an option in SURFACE_PATCHES to impose fissures
!                             along patch edges emphasizes the need to find a
!                             way to use neighboring patch info (next step).
!     09/18/15    "     "     Moved the unfinished subroutine form to this
!                             module.  (The allocatable working grid idea does
!                             not work as an argument.)
!     09/25/15    "     "     The neighboring patch scheme is working, with just
!                             the corner point case to deal with.  First, use
!                             FDCNTR and FD1SID directly instead of using LCSFIT
!                             as originally done.
!     10/06/15    "     "     Averaging of corner point results is implemented.
!
!     Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer,          intent (in)   :: np       ! Number of surface patches
      type (grid_type), intent (inout):: grid(np) ! Surface grid patch array;
                                                  ! %xmax, etc. may be set here
      integer,          intent (in)   :: ip       ! Patch associated with (i,j)
      integer,          intent (in)   :: i, j     ! Indices of target grid point
      real,             intent (out)  :: tani(3)  ! Unit tangent vector in the
                                                  ! i direction, (xu,yu,zu)'
      real,             intent (out)  :: tanj(3)  ! Unit tangent vector in the
                                                  ! j direction, (xv,yv,zv)'
      real,             intent (out)  :: unorm(3) ! Unit normal vector at (i,j),
                                                  ! via the cross product of the
                                                  ! two unit tangent vectors.
                                                  ! This unit normal should be
                                                  ! independent of patch index-
                                                  ! ing at common edge points,
                                                  ! even though the unit tangent
                                                  ! vectors may not be.
!     Local constants:

      real,    parameter :: eps = 5.e-5           ! Data range fraction for an
                                                  ! edge-point match
      logical, parameter :: norm = .false.        ! Don't normalize arc lengths

!     Local variables:

      integer :: i1, i2, ie, iev, inc, j1, j2, jev, ni, nj
      real    :: fpp
      logical :: corner, maxi, mini, maxj, minj

      integer, save :: npw                        ! # patches in working grid
      real,    save :: tol, tolsq                 ! Initialized once per grid

      real :: total, vsize
      real :: rowi(3,3)                           ! rowi(:,1) are x coordinates;
      real :: colj(3,3)                           ! rowi(:,2) are y coordinates;
                                                  ! rowi(:,3) are z coordinates;
                                                  ! colj(:,1) are x coordinates;
                                                  ! colj(:,2) are y coordinates;
                                                  ! colj(:,3) are z coordinates;
      real :: u(3), v(3)                          ! these can be used to compute
                                                  ! 3-point partial derivatives
!     Execution:

!     If indicated, set up a working grid that handles symmetry the easy way
!     using the interior-edge point method appropriate for any interior edge,
!     which requires the patch connectivity information that is also set up on
!     this first call:

      if (first_unit_normal) call set_up_working_grid ()  ! Internal procedure

!     The following, without the calls to seek_neighoring_data, comprises the
!     one-patch/ignore-neighbors method, for non-corner target points.
!     (Corner points need it repeatedly in different form.)
!
!     Indices iev, jev are those of the 3-point arcs at which central or one-
!     sided surface derivatives are evaluated.
!
!     Indexing is in the direction that preserves handedness when we get the
!     unit normal from the cross-product of the tangent vectors.

      ni = wgrid(ip)%ni;  nj = wgrid(ip)%nj
      mini = i == 1;      maxi = i == ni
      minj = j == 1;      maxj = j == nj
      corner = (mini .or. maxi) .and. (minj .or. maxj)

      if (corner) then
         call treat_a_corner ()
         return  ! Avoid indenting
      end if

      if (mini) then
         i1 = 1;    i2 = 3;   iev = 1;   ie = 1  ! For an exterior/non-sym. edge
         inc = 1
      else if (maxi) then
         i1 = ni-2; i2 = ni;  iev = 3;   ie = 2  ! For an exterior/non-sym. edge
         inc = -1
      else
         i1 = i-1;  i2 = i+1; iev = 2            ! The trivial not-an-edge case
      end if
      
      rowi(:,1) = wgrid(ip)%x(i1:i2,j,1)
      rowi(:,2) = wgrid(ip)%y(i1:i2,j,1)
      rowi(:,3) = wgrid(ip)%z(i1:i2,j,1)

!     Instead of this single-patch approach, make use of neighboring patch data
!     if possible (perhaps from reflection):

      if (mini .or. maxi) call seek_neighboring_data ()  ! Internal procedure

!!!   write (60, '(a, 3i4, 3es13.5, 2x, 3es13.5, 2x, 3es13.5)') &
!!!      'ip,i,j,rowx(1:3),rowy(1:3),rowz(1:3): ', &
!!!       ip,i,j,rowi(:,1),rowi(:,2),rowi(:,3)

      call chords3d (3, rowi(:,1), rowi(:,2), rowi(:,3), norm, total, u)

!     3-point finite difference tangent vector in the i direction:

      if (iev == 2) then
         call fdcntr (iev, u, rowi(:,1), tani(1), fpp)  ! See /numodules
         call fdcntr (iev, u, rowi(:,2), tani(2), fpp)
         call fdcntr (iev, u, rowi(:,3), tani(3), fpp)
      else
         call fd1sid (iev, inc, u, rowi(:,1), tani(1), fpp)
         call fd1sid (iev, inc, u, rowi(:,2), tani(2), fpp)
         call fd1sid (iev, inc, u, rowi(:,3), tani(3), fpp)
      end if

      vsize = sqrt (dot_product (tani, tani))
      tani  = tani / vsize

!     Repeat for the second index direction:

      if (minj) then
         j1 = 1;    j2 = 3;   jev = 1;   ie = 3;  inc = 1
      else if (maxj) then
         j1 = nj-2; j2 = nj;  jev = 3;   ie = 4;  inc = -1
      else
         j1 = j-1;  j2 = j+1; jev = 2
      end if

      colj(:,1) = wgrid(ip)%x(i,j1:j2,1)
      colj(:,2) = wgrid(ip)%y(i,j1:j2,1)
      colj(:,3) = wgrid(ip)%z(i,j1:j2,1)

      if (minj .or. maxj) call seek_neighboring_data ()

!!!   write (60, '(a, 3i4, 3es13.5, 2x, 3es13.5, 2x, 3es13.5)') &
!!!      'ip,i,j,colx(1:3),coly(1:3),colz(1:3): ', &
!!!       ip,i,j,colj(:,1),colj(:,2),colj(:,3)

      call chords3d (3, colj(:,1), colj(:,2), colj(:,3), norm, total, v)

      if (jev == 2) then
         call fdcntr (jev, v, colj(:,1), tanj(1), fpp)
         call fdcntr (jev, v, colj(:,2), tanj(2), fpp)
         call fdcntr (jev, v, colj(:,3), tanj(3), fpp)
      else
         call fd1sid (jev, inc, v, colj(:,1), tanj(1), fpp)
         call fd1sid (jev, inc, v, colj(:,2), tanj(2), fpp)
         call fd1sid (jev, inc, v, colj(:,3), tanj(3), fpp)
      end if

      vsize = sqrt (dot_product (tanj, tanj))
      tanj  = tanj / vsize

      call cross (tani, tanj, unorm)

      vsize = sqrt (dot_product (unorm, unorm))  ! Make sure of it
      unorm = unorm / vsize

!!!   write (60, '(a, 15x, 3i4, 3es13.5, 2x, 3es13.5, 2x, 3es13.5)') &
!!!      'ip,i,j,tani,tanj,unorm:', ip,i,j,tani,tanj,unorm

!     Internal procedures for subroutine surface_vectors:

      contains

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine set_up_working_grid ()

!        Subroutine surface_vectors deals with symmetry about the y = 0 plane
!        by reflecting the surface grid to make points on a y = 0 edge appear
!        to be on an interior edge, not an exterior edge, leaving only interior
!        edge points to deal with.  (Corner points till present headaches.)
!
!        This initialization for a given surface grid does any such reflection
!        and also establishes the connectivity information for the working grid
!        whether it was reflected or not.

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        Local variables:

         integer :: ib, ip, ni, nj
         real    :: ymin, ymax
         logical :: reflect

         if (allocated (wgrid)) then
            deallocate (wgrid, icon)
         end if

         do ib = 1, np
            call patch_data_range (grid(ib))
         end do

         ymin  = minval (grid(:)%ymin)
         ymax  = maxval (grid(:)%ymax)
         tol   = eps * (ymax - ymin)
         tolsq = tol*tol  ! Needed for corner points
         write (*, '(a, 3es16.8)') ' ymin, ymax, tol:', ymin, ymax, tol

         reflect = abs (ymin) < tol .or. abs (ymax) < tol

         if (reflect) then
            npw = np + np
         else
            npw = np
         end if

         allocate (wgrid(npw))

         do ib = 1, np
            ni = grid(ib)%ni;  nj = grid(ib)%nj
            wgrid(ib)%ni = ni
            wgrid(ib)%nj = nj
            wgrid(ib)%nk = 1
            allocate (wgrid(ib)%x(ni,nj,1), &
                      wgrid(ib)%y(ni,nj,1), &
                      wgrid(ib)%z(ni,nj,1))
            wgrid(ib)%x  = grid(ib)%x
            wgrid(ib)%y  = grid(ib)%y
            wgrid(ib)%z  = grid(ib)%z
            if (reflect) then
               ip = ib + np
               wgrid(ip)%ni =  ni
               wgrid(ip)%nj =  nj
               wgrid(ip)%nk =  1
               allocate (wgrid(ip)%x(ni,nj,1), &
                         wgrid(ip)%y(ni,nj,1), &
                         wgrid(ip)%z(ni,nj,1))
               wgrid(ip)%x  =  grid(ib)%x
               wgrid(ip)%y  = -grid(ib)%y  ! Reverses handedness
               wgrid(ip)%z  =  grid(ib)%z

               call reverse_patch_i (wgrid(ip), 0)  ! Restore handedness
            end if
         end do

         allocate (icon(2,4,npw))

         call surface_connectivity (npw, wgrid, tol, icon)

         end subroutine set_up_working_grid

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine seek_neighboring_data ()

!        For a target point (i, j) that is at an edge of surface patch ip but
!        NOT AT A CORNER, use the connectivity info established during initial-
!        ization to provide best-possible data from any neighboring patch(es)
!        for the surface tangent vectorss at (i, j) and hence for the unit
!        normal vector.  The tangent vectors come from 3-point differencing.
!        Common edges are assumed to be properly point-to-point matched.
!
!        Any symmetry (about the y = 0 plane only, at present) has been handled
!        by reflecting the working grid to make a y = 0 exterior edge into an
!        interior edge.  Any other exterior edge in the input grid cannot do
!        better than the one-sided 3-point data already set up by the caller.
!
!        Corner points are handled elsewhere via averaging.
!
!        To interpret icon(1:2,1:4,1:np), consider icon(l,m,n):
!
!           n is the target patch number (ip for the current call);
!           m is ie = 1,2,3,4 for edge imin,imax,jmin,jmax of patch n;
!           l = 1 contains either a patch number iq with an edge matching edge
!                 ie of patch ip, or 0 if no such match exists;
!           l = 2 contains (if iq /= 0) the edge number ig of patch iq that
!                 matches edge ie of patch ip; it is negative if the edges are
!                 indexed in opposite directions.
!
!        A certain 10-patch sphere produces connectivity info like this:
!
!                     ip1 ie1 ip2 ie2 ip3 ie3 ip4 ie4
!                       6 -1    2 -3    3 -3    5 -1      (patch 1 info)
!                       5  3    3  1    1 -2    4 -1
!                       2  2    8 -2    1 -3    4  3         "   3   "
!                       2 -4    9 -2    3  4    5  2
!                       1 -4    4  4    2  1   10  3
!                       1 -1    7  4   10  1    8 -4
!                      10 -4    8  1    9  1    6  2
!                       7  2    3 -2    9  4    6 -4         "   8   "
!                       7  3    4 -2   10 -2    8  3
!                       6  3    9 -3    5  4    7 -1         "  10   "

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        Local variables:

         integer :: ig, ii, iq, jj, mi, mj

!        Execution:

!        Using the connectivity data turns out to be surprisingly awkward.
!        There must be a better way to set up the right stencil!

         iq = icon(1,ie,ip)   ! Patch neighboring this edge ie of patch ip
         ig = icon(2,ie,ip)   ! Corresp. edge of patch iq

         if (iq /= 0) then    ! Interior point of an interior edge
            iev = 2;  jev = 2 ! Evaluate at the midpoints of the 3-pt. lines
            mi = wgrid(iq)%ni
            mj = wgrid(iq)%nj
            select case (ig)
               case (1, -1)
                  ii = 2
                  if (ie <= 2) then
                     jj = j
                  else
                     jj = i
                  end if
                  if (ig < 0) jj = mj + 1 - jj
               case (2, -2)
                  ii = mi - 1
                  if (ie <= 2) then
                     jj = j
                  else
                     jj = i
                  end if
                  if (ig < 0) jj = mj + 1 - jj
               case (3, -3)
                  if (ie >= 3) then
                     ii = i
                  else
                     ii = j
                  end if
                  if (ig < 0) ii = mi + 1 - ii
                  jj = 2
               case (4, -4)
                  if (ie >= 3) then
                     ii = i
                  else
                     ii = j
                  end if
                  if (ig < 0) ii = mi + 1 - ii
                  jj = mj - 1
            end select
            select case (ie)  ! Subject edge of patch ip
               case (1, -1)
                  rowi(  1,1) = wgrid(iq)%x(ii,jj,1)
                  rowi(  1,2) = wgrid(iq)%y(ii,jj,1)
                  rowi(  1,3) = wgrid(iq)%z(ii,jj,1)
                  rowi(2:3,1) = wgrid(ip)%x(1:2,j,1)
                  rowi(2:3,2) = wgrid(ip)%y(1:2,j,1)
                  rowi(2:3,3) = wgrid(ip)%z(1:2,j,1)
               case (2, -2)
                  rowi(1:2,1) = wgrid(ip)%x(ni-1:ni,j,1)
                  rowi(1:2,2) = wgrid(ip)%y(ni-1:ni,j,1)
                  rowi(1:2,3) = wgrid(ip)%z(ni-1:ni,j,1)
                  rowi(  3,1) = wgrid(iq)%x(ii,jj,1)
                  rowi(  3,2) = wgrid(iq)%y(ii,jj,1)
                  rowi(  3,3) = wgrid(iq)%z(ii,jj,1)
               case (3, -3)
                  colj(  1,1) = wgrid(iq)%x(ii,jj,1)
                  colj(  1,2) = wgrid(iq)%y(ii,jj,1)
                  colj(  1,3) = wgrid(iq)%z(ii,jj,1)
                  colj(2:3,1) = wgrid(ip)%x(i,1:2,1)
                  colj(2:3,2) = wgrid(ip)%y(i,1:2,1)
                  colj(2:3,3) = wgrid(ip)%z(i,1:2,1)
               case (4, -4)
                  colj(1:2,1) = wgrid(ip)%x(i,nj-1:nj,1)
                  colj(1:2,2) = wgrid(ip)%y(i,nj-1:nj,1)
                  colj(1:2,3) = wgrid(ip)%z(i,nj-1:nj,1)
                  colj(  3,1) = wgrid(iq)%x(ii,jj,1)
                  colj(  3,2) = wgrid(iq)%y(ii,jj,1)
                  colj(  3,3) = wgrid(iq)%z(ii,jj,1)
            end select
         else
!           Interior point on an exterior edge--no neighbor.
!           Use 1-sided differencing as already set up.
         end if

         end subroutine seek_neighboring_data

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine treat_a_corner ()

!        Surface patch corners may be common to any number of patches.  In order
!        to produce the same unit normal vector for all ip/i/j possibilities re-
!        presenting the same corner point, they must all calculate the same av-
!        eraged result.  Tangent vectors at a corner point, on the other hand,
!        cannot very well be averaged.  Therefore, just the last result produced
!        here is used.  This ignores possible neighboring patch data for now, so
!        it may not be as good as it could be, but the tangent vectors probably
!        aren't needed anyway--just the normal.

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        Local variables:

         integer :: ilocal, iplocal, jlocal, nilocal, njlocal, nsum
         real    :: dsqrd, ulocal(3), usum(3), xc, yc, zc

!        Execution:
 
!        Observe that ip, i, j here are those passed to surface_vectors, and
!        they've been found to be at a corner.  We check all corners of all
!        other patches for matching corners.

         call simple_corner_method (ip, i, j, unorm)

!!!      if (ip==1 .or. ip==2 .or. ip==5) &
!!!          write (65, '(3i4, 3es16.8)') ip,i,j, unorm

         if (npw > 1) then
            nsum = 1
            usum(:) = unorm(:)
            xc = wgrid(ip)%x(i,j,1)
            yc = wgrid(ip)%y(i,j,1)
            zc = wgrid(ip)%z(i,j,1)

            do iplocal = 1, npw
               if (iplocal == ip) cycle

               nilocal = wgrid(iplocal)%ni
               njlocal = wgrid(iplocal)%nj

               do jlocal = 1, njlocal, njlocal - 1
                  do ilocal = 1, nilocal, nilocal - 1
                     dsqrd = dsq (wgrid(iplocal)%x(ilocal,jlocal,1), &
                                  wgrid(iplocal)%y(ilocal,jlocal,1), &
                                  wgrid(iplocal)%z(ilocal,jlocal,1), xc, yc, zc)
!!!                  write (65, '(a, 6i5, es16.8)') 'ip,i,j,ip/i/jlocal,dsq:', &
!!!                     ip,i,j,iplocal,ilocal,jlocal,dsqrd
                     if (dsqrd < tolsq) then
                         call simple_corner_method (iplocal, ilocal, jlocal, &
                                                    ulocal)
                         usum(:) = usum(:) + ulocal(:)
                         nsum    = nsum + 1
                     end if
                  end do
               end do
            end do
            unorm = usum / real (nsum)
            vsize = sqrt (dot_product (unorm, unorm))
            unorm = unorm / vsize
         end if

         end subroutine treat_a_corner

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine simple_corner_method (ip, i, j, unorm)

!        This form of the straightforward one-patch method for corner points
!        is arranged to be called for all patches sharing a common corner.
!
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        Arguments:

         integer, intent (in)  :: ip, i, j  ! Patch # and corner indices
         real,    intent (out) :: unorm(3)  ! Unit normal from 1-sided method

!        Local variables:

         integer :: i1, i2, j1, j2, ni, nj
         logical :: maxi, mini, maxj, minj

!        Execution:

         ni = wgrid(ip)%ni;  nj = wgrid(ip)%nj
         mini = i == 1;      maxi = i == ni
         minj = j == 1;      maxj = j == nj

         if (mini) then
            i1 = 1;    i2 = 3;   iev = 1;   ie = 1;   inc = 1
         else if (maxi) then
            i1 = ni-2; i2 = ni;  iev = 3;   ie = 2;   inc = -1
         end if

         rowi(:,1) = wgrid(ip)%x(i1:i2,j,1)
         rowi(:,2) = wgrid(ip)%y(i1:i2,j,1)
         rowi(:,3) = wgrid(ip)%z(i1:i2,j,1)

         call chords3d (3, rowi(:,1), rowi(:,2), rowi(:,3), norm, total, u)

!        1-sided 3-point finite difference tangent vector in the i direction:

         call fd1sid (iev, inc, u, rowi(:,1), tani(1), fpp)
         call fd1sid (iev, inc, u, rowi(:,2), tani(2), fpp)
         call fd1sid (iev, inc, u, rowi(:,3), tani(3), fpp)

         vsize = sqrt (dot_product (tani, tani))
         tani  = tani / vsize

!        Repeat for the second index direction:

         if (minj) then
            j1 = 1;    j2 = 3;   jev = 1;   ie = 3;  inc = 1
         else if (maxj) then
            j1 = nj-2; j2 = nj;  jev = 3;   ie = 4;  inc = -1
         end if

         colj(:,1) = wgrid(ip)%x(i,j1:j2,1)
         colj(:,2) = wgrid(ip)%y(i,j1:j2,1)
         colj(:,3) = wgrid(ip)%z(i,j1:j2,1)

         call chords3d (3, colj(:,1), colj(:,2), colj(:,3), norm, total, v)

         call fd1sid (jev, inc, v, colj(:,1), tanj(1), fpp)
         call fd1sid (jev, inc, v, colj(:,2), tanj(2), fpp)
         call fd1sid (jev, inc, v, colj(:,3), tanj(3), fpp)

         vsize = sqrt (dot_product (tanj, tanj))
         tanj  = tanj / vsize

         call cross (tani, tanj, unorm)

         vsize = sqrt (dot_product (unorm, unorm))  ! Make sure of it
         unorm = unorm / vsize

         end subroutine simple_corner_method

      end subroutine surface_vectors

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine swap_coordinates (patch, m, n)

!     Exchange the indicated pair of coordinates.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent (inout) :: patch
      integer, intent (in) :: m, n        ! 1, 2, 3  <->  x, y, z

!     Local variables:

      integer :: ni, nj

!     Execution:

      ni = patch%ni;  nj = patch%nj

      select case (m)
         case (1)
            select case (n)
               case (2)
                  call swap (ni, nj, patch%x, patch%y)
               case (3)
                  call swap (ni, nj, patch%x, patch%z)
            end select
         case (2)
            select case (n)
               case (1)
                  call swap (ni, nj, patch%y, patch%x)
               case (3)
                  call swap (ni, nj, patch%y, patch%z)
            end select
         case (3)
            select case (n)
               case (1)
                  call swap (ni, nj, patch%z, patch%x)
               case (2)
                  call swap (ni, nj, patch%z, patch%y)
            end select
      end select

      contains

         subroutine swap (ni, nj, cm, cn)

         integer, intent (in)    :: ni, nj
         real,    intent (inout) :: cm(ni,nj,1), cn(ni,nj,1)
         integer :: i, j
         real    :: temp

!        Play safe by not assuming packed 2-D arrays.

         do j = 1, nj
            do i = 1, ni
               temp      = cm(i,j,1)
               cm(i,j,1) = cn(i,j,1)
               cn(i,j,1) = temp
            end do
         end do

         end subroutine swap

      end subroutine swap_coordinates

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine transpose_patch (patch, nf)

!     Reverse the i and j indices of a surface patch, in-place.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent (inout) :: patch
      integer,          intent (in)    :: nf   ! # function values; nf >= 0

!     Local variables:

      integer          :: i, j, mi, mj, ni, nj
      real             :: t, tq(nf)
      type (grid_type) :: tpatch

!     Execution:

      ni = patch%ni;  nj = patch%nj

      allocate (tpatch%x(nj,ni,1), tpatch%y(nj,ni,1), tpatch%z(nj,ni,1))

      do j = 1, nj
         do i = 1, ni
            tpatch%x(j,i,1) = patch%x(i,j,1)
            tpatch%y(j,i,1) = patch%y(i,j,1)
            tpatch%z(j,i,1) = patch%z(i,j,1)
         end do
      end do

      deallocate (patch%x, patch%y, patch%z)
      patch%ni = nj;  patch%nj = ni
      allocate (patch%x(nj,ni,1), patch%y(nj,ni,1), patch%z(nj,ni,1))

      do j = 1, ni
         do i = 1, nj
            patch%x(i,j,1) = tpatch%x(i,j,1)
            patch%y(i,j,1) = tpatch%y(i,j,1)
            patch%z(i,j,1) = tpatch%z(i,j,1)
         end do
      end do

      deallocate (tpatch%x, tpatch%y, tpatch%z)

      if (nf > 0) then

         mi = patch%mi;  mj = patch%mj

         allocate (tpatch%q(nf,mj,mi,1))

         do j = 1, mj
            do i = 1, mi
               tpatch%q(:,j,i,1) = patch%q(:,i,j,1)
            end do
         end do

         deallocate (patch%q)
         patch%mi = mj;  patch%mj = mi
         allocate (patch%q(nf,mj,mi,1))

         do j = 1, mi
            do i = 1, mj
               patch%q(:,i,j,1) = tpatch%q(:,i,j,1)
            end do
         end do

         deallocate (tpatch%q)

      end if

      end subroutine transpose_patch

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine update_patch (patch_in, nf, release_in, release_out, patch_out)

!     Transcribe one surface patch to another.
!     The output patch fields are (re)allocated here before the copy.
!     The input patch fields are deallocated after the copy if indicated.
!     The "release" arguments are adopted because of observed misbehavior with
!     if (allocated(...)) ...
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent (inout) :: patch_in    ! Surface patch to copy
      integer,          intent (in)    :: nf          ! # functions >= 0
      logical,          intent (in)    :: release_in  ! T => deallocate inputs
      logical,          intent (in)    :: release_out ! T => reallocate outputs
      type (grid_type), intent (inout) :: patch_out   ! x

!     Execution:

      if (release_out) then
         deallocate (patch_out%x, patch_out%y, patch_out%z)
         if (nf > 0) deallocate (patch_out%q)
      end if

      patch_out%ni = patch_in%ni
      patch_out%nj = patch_in%nj
      patch_out%nk = 1

      allocate (patch_out%x(patch_in%ni,patch_in%nj,1), &
                patch_out%y(patch_in%ni,patch_in%nj,1), &
                patch_out%z(patch_in%ni,patch_in%nj,1))

      patch_out%x = patch_in%x
      patch_out%y = patch_in%y
      patch_out%z = patch_in%z

      if (nf > 0) then
         patch_out%mi = patch_in%mi
         patch_out%mj = patch_in%mj
         patch_out%mk = 1
         allocate (patch_out%q(nf,patch_in%mi,patch_in%mj,1))
         patch_out%q  = patch_in%q
      end if

      if (release_in) then
         deallocate (patch_in%x, patch_in%y, patch_in%z)
         if (nf > 0) deallocate (patch_in%q)
      end if

      end subroutine update_patch

   end module surface_patch_utilities
