!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program fix_filler
!
!  Description:
!
!     FIX_FILLER is a specialized utility that removes undesirable skewness
!     in the thin grid blocks surrounding a bent protruding gap filler.  Such
!     blocks are generated initially with a Gridgen glf script by Chun Tang.
!
!  Strategy:
!
!     >  Read the input volume grid (4 blocks).
!     >  Construct the wall surface patches and fill in the patch missing
!        where the gap filler is.  We need to determine the gap filler height
!        off this surface at each of its points via ADT searching.
!     >  For each of the 4 blocks:
!        >  Transfer the j = 1 and k = 1 faces to the output block.
!        >  For each point (i,nj,1) on the outer edge:
!           >  Construct a wall normal with length equal to the height of the
!              corresponding point on the gap filler (i,1,nk) from the wall.
!           >  Impose a standard relative k distribution on this line.
!        >  For the i = 1 and i = ni faces:
!           >  Morph the edge (i,1:nj,nk) from the original to the new outer
!              boundary (NULINE3D).
!           >  Morph the interior (WARPQ3D).
!        >  Morph the original k = nk face to the new edges (WARPQ3D).
!        >  Morph the original volume grid (WARP3D).
!     >  Save the adjusted block.
!     >  Save the adjusted offset surface for use by RADIAL_INTERP.  It needs
!        the patch from the top of the gap filler in front of the "collar"
!        patches for PLUG_INTERP reasons.
!
!  Procedures:
!
!     XYZQ_IO package  I/O utilities for PLOT3D grid and function files
!     <Numerous other numerical utilities>
!
!  History:
!
!     08/22/06  D. Saunders  Initial design: intersect j lines with a 2-line
!                            artificial surface off the wall.
!     08/25/06  "     "      Initial implementation.
!     08/26/06  "     "      Realized we need the revised offset surface output.
!     08/28/06  "     "      Unfortunately, some line skewing still appears in
!                            the new outer boundary normal to the wall.  This is
!                            worst opposite the ends of the bent gap filler.
!                            Therefore, give up the line-surface intersection
!                            approach in favor of grid perturbation techniques
!                            by constructing the outer j = nj boundary precisely
!                            with k lines strictly normal to the wall (no more
!                            skewing on the outer boundary, which is important
!                            for the PLUG_INTERP step of the rapid analysis
!                            procedures).
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure
   use xyzq_io_module
   use adt_utilities

   implicit none

!  Constants:

   integer, parameter ::   &
      lunin  =  1,         & ! Skewed volume grid
      lunout =  2,         & ! Unskewed volume grid
      lunsrf =  3,         & ! Corresponding offset surface with middle patch
      lunkbd =  5,         &
      luncrt =  6,         &
      lunq   = -1,         & ! Suppresses function file input
      m      =  4            ! # blocks surrounding the plug

   real, parameter ::      &
      one    = 1.,         &
      zero   = 0.

   logical, parameter ::   &
      false  = .false.,    &
      true   = .true.

!  Variables:

   integer :: &
      i, ib, ic, ios, iquad, j, jc, k, nblocks, ni, nj, nk, npts, nquad, num_q

   real :: &
      dsq, dx, dy, dz, height, p, pint, q, qint, stotal, t, tint, tm1,         &
      un(3), xyzint(3), xyztarget(3)

   logical :: &
      cell_centered, formatted_in, formatted_out

   integer, allocatable :: &
      conn(:,:)                   ! For (patch,i,j) of surface quads. to search

   real, allocatable, dimension (:) :: &
      tline, t1, x1, y1, z1, t2, x2, y2, z2, xline, yline, zline

   real, allocatable, dimension (:,:,:,:) :: &
      s0

   real, allocatable, dimension (:,:,:,:,:) :: &
      dfacei, dfacej, dfacek

   character :: &
      lcs_method * 1

!  Derived data type:

   type (grid_type), pointer, dimension (:) :: &
      xyzin, xyzout  ! Input and output volume grids

   type (grid_type), pointer, dimension (:) :: &
      xyzsurf,    &  ! Adjusted offset surface, with patch 1 at top of gapfiller
      xyzwall        ! Wall patches with missing center patch 1 inserted

!  Execution:
!  ---------

   call file_prompt (lunin, -luncrt, 'input volume grid', 'old', true,         &
                     formatted_in, ios)
   if (ios /= 0) go to 99

   call file_prompt (lunout, -luncrt, 'output grid',  'unknown', true,         &
                     formatted_out, ios)
   if (ios /= 0) go to 99

!  Read all blocks of the input grid.  The "cell_centered" argument isn't
!  relevant: the grid should be vertex-centered. Lunq < 0 suppresses the q file.

   call xyzq_read (lunin, lunq, formatted_in, nblocks, num_q, cell_centered,   &
                   xyzin, ios)
   if (ios /= 0) go to 99

   write (luncrt, '(/, a, /)') ' Block dimensions:'
   write (luncrt, '(4(i6, 2X, 3i5))') &
      (i, xyzin(i)%ni, xyzin(i)%nj, xyzin(i)%nk, i = 1, nblocks)


!                 Assumptions (some easily generalized if necessary):
!                 ---------------------------------------------------

!                 The blocks are right handed.
!                 k = 1 is on the wall.
!                 j is the spoke direction.
!                 j increases towards the outer boundary.
!                 The blocks are in a sensible order around the plug filler, to
!                 simplify constructing the missing inner patches at k = 1 & nk.


!  Set up the output unskewed volume grid file:

   allocate (xyzout(m))   ! m = 4 above

   nj = xyzin(1)%nj;  nk = xyzin(1)%nk  ! Common to all 4 blocks

   do ib = 1, m
      xyzout(ib)%ni = xyzin(ib)%ni
      xyzout(ib)%nj = nj
      xyzout(ib)%nk = nk

      call xyz_allocate (xyzout(ib), ios)
   end do

   call xyz_header_io (2, lunout, formatted_out, nblocks, xyzout, ios)

   if (ios /= 0) go to 99


!  Set up the wall patches as a surface to be searched to determine the
!  plug height at each i.  We need to fill in the missing center patch 1.
!  This was adapted from the code below for the offset modified surface.

   allocate (xyzwall(m+1))

   nquad = 0

   do ib = 2, m + 1
      xyzwall(ib)%ni = xyzin(ib-1)%ni
      xyzwall(ib)%nj = xyzin(ib-1)%nj
      xyzwall(ib)%nk = 1
      nquad = nquad + (xyzwall(ib)%ni - 1) * (xyzwall(ib)%nj - 1)

      call xyz_allocate (xyzwall(ib), ios)

      xyzwall(ib)%x(:,:,1) = xyzin(ib-1)%x(:,:,1)
      xyzwall(ib)%y(:,:,1) = xyzin(ib-1)%y(:,:,1)
      xyzwall(ib)%z(:,:,1) = xyzin(ib-1)%z(:,:,1)
   end do

   xyzwall(1)%ni = xyzwall(2)%ni
   xyzwall(1)%nj = xyzwall(3)%ni
   xyzwall(1)%nk = 1

   call xyz_allocate (xyzwall(1), ios)

   nquad = nquad + (xyzwall(1)%ni - 1) * (xyzwall(1)%nj - 1)

!  Set up the edges for patch 1:

   ni = xyzwall(1)%ni
   nj = xyzwall(1)%nj;  ib = ni

   do i = 1, ni
      xyzwall(1)%x(i,1,1)  = xyzwall(2)%x(ib,1,1)
      xyzwall(1)%y(i,1,1)  = xyzwall(2)%y(ib,1,1)
      xyzwall(1)%z(i,1,1)  = xyzwall(2)%z(ib,1,1)
      ib = ib - 1
      xyzwall(1)%x(i,nj,1) = xyzwall(4)%x(i,1,1)
      xyzwall(1)%y(i,nj,1) = xyzwall(4)%y(i,1,1)
      xyzwall(1)%z(i,nj,1) = xyzwall(4)%z(i,1,1)
   end do

   ib = nj

   do j = 1, nj
      xyzwall(1)%x(ni,j,1) = xyzwall(3)%x(ib,1,1)
      xyzwall(1)%y(ni,j,1) = xyzwall(3)%y(ib,1,1)
      xyzwall(1)%z(ni,j,1) = xyzwall(3)%z(ib,1,1)
      ib = ib - 1
      xyzwall(1)%x(1,j,1)  = xyzwall(5)%x(j,1,1)  ! Confusing!
      xyzwall(1)%y(1,j,1)  = xyzwall(5)%y(j,1,1)
      xyzwall(1)%z(1,j,1)  = xyzwall(5)%z(j,1,1)
   end do

!  Fill the interior of the missing patch.  It may not be precisely on the wall,
!  but that doesn't matter for determining the plug height at each i.

   allocate (s0(2*(ni + nj),1,1,1))

   call tfiq3d (ni, 1, ni, 1, nj, xyzwall(1)%x, xyzwall(1)%y, xyzwall(1)%z, s0)

   deallocate (s0)

!  Build a search tree from the wall patches:

   allocate (conn(3,nquad))  ! For (patch,i,j) of each surface quad.

   call build_adt (nblocks, xyzwall, nquad, conn)


!  The updated offset surface grid needs to be used by RADIAL_INTERP.
!  Set up its dimensions and storage:

   allocate (xyzsurf(m+1))  ! Includes the top of the gap filler

   do ib = 2, m + 1         ! PLUG_INTERP requires the "collar" blocks last
      xyzsurf(ib)%ni = xyzin(ib-1)%ni
      xyzsurf(ib)%nj = xyzin(ib-1)%nj
      xyzsurf(ib)%nk = 1
      call xyz_allocate (xyzsurf(ib), ios)
   end do

   xyzsurf(1)%ni = xyzsurf(2)%ni
   xyzsurf(1)%nj = xyzsurf(3)%ni
   xyzsurf(1)%nk = 1
   call xyz_allocate (xyzsurf(1), ios)


!  Process one volume block at a time:
!  -----------------------------------

   nj = xyzin(1)%nj;  nk = xyzin(1)%nk  ! Common to all blocks
   jc = nj - 1;  q = one                ! Applies to searches at j = jmax

   allocate (xline(nk), yline(nk), zline(nk), tline(nk))
   allocate (x1(nj), y1(nj), z1(nj), t1(nj), x2(nj), y2(nj), z2(nj), t2(nj))

   do ib = 1, nblocks

      ni = xyzin(ib)%ni;  npts = ni * nj * nk

      xyzout(ib)%x(:,1,:) = xyzin(ib)%x(:,1,:)  ! j = 1 face on the gap filler
      xyzout(ib)%y(:,1,:) = xyzin(ib)%y(:,1,:)
      xyzout(ib)%z(:,1,:) = xyzin(ib)%z(:,1,:)

      xyzout(ib)%x(:,:,1) = xyzin(ib)%x(:,:,1)  ! k = 1 face on the wall
      xyzout(ib)%y(:,:,1) = xyzin(ib)%y(:,:,1)
      xyzout(ib)%z(:,:,1) = xyzin(ib)%z(:,:,1)

      p = zero

      do i = 1, ni

!        Establish the normalized k distribution on the gap filler at this i:

         xline(:) = xyzin(1)%x(i,1,:)
         yline(:) = xyzin(1)%y(i,1,:)
         zline(:) = xyzin(1)%z(i,1,:)

         call chords3d (nk, xline, yline, zline, true, stotal, tline)

!        Determine the height of the plug at this i:

         xyztarget(1) = xline(nk)
         xyztarget(2) = yline(nk)
         xyztarget(3) = zline(nk)

         call search_adt (xyztarget, iquad, pint, qint, dsq, true, m,          &
                          xyzwall, nquad, conn, xyzint)

         height = sqrt (dsq)

!        Determine the unit normal to the wall at this (i,nj):

         if (i < ni) then
            ic = i     ! jc  = nj - 1 above
         else    ! Leave ic at ni - 1 (must point to "lower left" cell)
            p = one
         end if

         call surface_normal (ni, nj, xyzwall(ib+1)%x, xyzwall(ib+1)%y,        &
                              xyzwall(ib+1)%z, ic, jc, p, q, un)

         xyzout(ib)%x(i,nj,nk) = xyzout(ib)%x(i,nj,1) + height * un(1)
         xyzout(ib)%y(i,nj,nk) = xyzout(ib)%y(i,nj,1) + height * un(2)
         xyzout(ib)%z(i,nj,nk) = xyzout(ib)%z(i,nj,1) + height * un(3)

!        Straight lines for the new outer boundary in the k direction:

         do k = 2, nk

            t = tline(k);          tm1 = one - t

            xyzout(ib)%x(i,nj,k) = tm1 * xyzout(ib)%x(i,nj,1) +                &
                                     t * xyzout(ib)%x(i,nj,nk)
            xyzout(ib)%y(i,nj,k) = tm1 * xyzout(ib)%y(i,nj,1) +                &
                                     t * xyzout(ib)%y(i,nj,nk)
            xyzout(ib)%z(i,nj,k) = tm1 * xyzout(ib)%z(i,nj,1) +                &
                                     t * xyzout(ib)%z(i,nj,nk)
         end do

      end do  ! Next i

!     Morph the i = 1 and ni faces.  The radial edges at nk still have to be
!     set up.  WARPQ3D and WARP3D need the input block relative arc lengths,
!     so calculate them now.

      allocate      (s0(ni,nj,nk,3))

      call paramxyz (1, ni, 1, nj, 1, nk, 1, ni, 1, nj, 1, nk,                 &
                     xyzin(ib)%x, xyzin(ib)%y, xyzin(ib)%z, s0)

      allocate      (dfacei(3, nj, nk, 2, 4), dfacej(3, ni, nk, 2, 4),         &
                     dfacek(3, ni, nj, 2, 4))

      do i = 1, ni, ni - 1

!        Perturb the radial edge at k = nk:

         x1(:) = xyzin(ib)%x(i,:,nk)
         y1(:) = xyzin(ib)%y(i,:,nk)
         z1(:) = xyzin(ib)%z(i,:,nk)

         x2(1) = x1(1);  x2(nj) = xyzout(ib)%x(i,nj,nk)
         y2(1) = y1(1);  y2(nj) = xyzout(ib)%y(i,nj,nk)
         z2(1) = z1(1);  z2(nj) = xyzout(ib)%z(i,nj,nk)

         call nuline3d (1, nj, x1, y1, z1, x2, y2, z2)

         xyzout(ib)%x(i,:,nk) = x2(:)
         xyzout(ib)%y(i,:,nk) = y2(:)
         xyzout(ib)%z(i,:,nk) = z2(:)

!        Morph the interior points of this i face.  WARPQ3D is TFI-like, but
!        it uses the original interior relative spacings.

         call warpq3d (1, ni, 1, nj, 1, nk, i, i, 1, nj, 1, nk,                &
                       xyzin(ib)%x, xyzin(ib)%y, xyzin(ib)%z,                  &
                       s0, dfacei, dfacej, dfacek,                             &
                       xyzout(ib)%x, xyzout(ib)%y, xyzout(ib)%z)
      end do

!     Morph the interior points of the k = nk face.

      call warpq3d (1, ni, 1, nj, 1, nk, 1, ni, 1, nj, nk, nk,                 &
                    xyzin(ib)%x, xyzin(ib)%y, xyzin(ib)%z,                     &
                    s0, dfacei, dfacej, dfacek,                                &
                    xyzout(ib)%x, xyzout(ib)%y, xyzout(ib)%z)

!     The perturbed block boundaries are now in place.
!     Morph the volume interior points from the originals:

      call warp3d   (1, ni, 1, nj, 1, nk, 1, ni, 1, nj, 1, nk,                 &
                     xyzin(ib)%x, xyzin(ib)%y, xyzin(ib)%z, s0,                &
                     dfacei, dfacej, dfacek,                                   &
                     xyzout(ib)%x, xyzout(ib)%y, xyzout(ib)%z)

      deallocate (dfacei, dfacej, dfacek, s0)

!     Save the perturbed block:

      call xyz_block_io (2, lunout, formatted_out, npts, xyzout(ib)%x,         &
                            xyzout(ib)%y, xyzout(ib)%z, ios)
      if (ios /= 0) go to 99

!     Transfer its k = nk face to the offset surface:

      xyzsurf(ib+1)%x(:,:,1) = xyzout(ib)%x(:,:,nk)
      xyzsurf(ib+1)%y(:,:,1) = xyzout(ib)%y(:,:,nk)
      xyzsurf(ib+1)%z(:,:,1) = xyzout(ib)%z(:,:,nk)

      deallocate (xyzout(ib)%x, xyzout(ib)%y, xyzout(ib)%z)

   end do ! Next block

   close (lunout)

   deallocate (conn, xline, yline, zline, tline, x1, y1, z1, t1, x2, y2, z2, t2)

   do ib = 1, nblocks
      deallocate (xyzin  (ib)%x, xyzin  (ib)%y, xyzin  (ib)%z)
      deallocate (xyzwall(ib)%x, xyzwall(ib)%y, xyzwall(ib)%z)
   end do

   deallocate (xyzin, xyzout, xyzwall)

!  Finally, construct the surface patch 1 at the top of the gap filler
!  for use by the RADIAL_INTERP procedure.

   ni = xyzsurf(1)%ni
   nj = xyzsurf(1)%nj;  ib = ni

   do i = 1, ni
      xyzsurf(1)%x(i,1,1)  = xyzsurf(2)%x(ib,1,1)
      xyzsurf(1)%y(i,1,1)  = xyzsurf(2)%y(ib,1,1)
      xyzsurf(1)%z(i,1,1)  = xyzsurf(2)%z(ib,1,1)
      ib = ib - 1
      xyzsurf(1)%x(i,nj,1) = xyzsurf(4)%x(i,1,1)
      xyzsurf(1)%y(i,nj,1) = xyzsurf(4)%y(i,1,1)
      xyzsurf(1)%z(i,nj,1) = xyzsurf(4)%z(i,1,1)
   end do

   ib = nj

   do j = 1, nj
      xyzsurf(1)%x(ni,j,1) = xyzsurf(3)%x(ib,1,1)
      xyzsurf(1)%y(ni,j,1) = xyzsurf(3)%y(ib,1,1)
      xyzsurf(1)%z(ni,j,1) = xyzsurf(3)%z(ib,1,1)
      ib = ib - 1
      xyzsurf(1)%x(1,j,1)  = xyzsurf(5)%x(j,1,1)  ! Note the j,1,1
      xyzsurf(1)%y(1,j,1)  = xyzsurf(5)%y(j,1,1)
      xyzsurf(1)%z(1,j,1)  = xyzsurf(5)%z(j,1,1)
   end do

   allocate (s0(2*(ni + nj),1,1,1))

   call tfiq3d (ni, 1, ni, 1, nj, xyzsurf(1)%x, xyzsurf(1)%y, xyzsurf(1)%z, s0)

   deallocate (s0)

   if (formatted_out) then
      open (lunsrf, file='adjusted-surf.grd', status='unknown')
   else
      open (lunsrf, file='adjusted-surf.gu', form='unformatted',               &
            status='unknown')
   end if

   nblocks = m + 1

   call xyz_write (lunsrf, formatted_out, nblocks, xyzsurf, ios)

   do ib = 1, nblocks
      deallocate (xyzsurf(ib)%x, xyzsurf(ib)%y, xyzsurf(ib)%z)
   end do

   deallocate (xyzsurf)

99 continue

! *** stop ! Avoid system dependencies.

   end program fix_filler
