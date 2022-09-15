!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program gsmooth
!
!  Description:
!
!     GSMOOTH applies the elliptic grid smoothing utility, ELLIP3D, to the
!  specified block of a multiblock grid file (or to a single-block grid).
!  The smoothing is CPU-intensive enough that treating more than one block
!  can reasonably be expected to be done with separate runs.  The standard
!  PLOT3D formats are supported (formatted or unformatted, with or without
!  IBLANK).  (Actually, translation from the CFD_IO_PACKAGE to the XYZQ_IO
!  package has lost the IBLANK option, though it could easily be restored
!  if the need arises since XYZQ_IO has iblank analogues.)
!
!     This version can be used to check cell volumes of the target block
!  (and nothing else except perhaps extract the block) by using ITMAX = 0.
!  ITMAX = -1 checks the cell volumes for ALL blocks (and nothing else).
!
!     Only Euler-type grid blocks are likely to be smoothed successfully--
!  Navier-Stokes-type blocks should be derived from smoothed Euler blocks.
!     
!  Sample Control File (gsmooth.inp):
!
!     ITMAX OMEGA CONV   EXPI1/I2/J1/J2/K1/K2    BGPHI  BGPSI  BGOMG  FGMODE
!      100   1.3  4.0  .45 .45 .45 .45 .45 .45  'YYYY' 'YYYY' 'YYYY' 'YYYYYY'
!     PRINT FGLIM URFG
!        F   1.0  0.1
!
!  Control Descriptions:
!
!     (Here, "foreground" refers to control of orthogonality at a boundary;
!      "background" refers to control of interior spacing from a boundary.)
!
!     ITMAX     Number of smoothing iterations:
!               0 means just check cell volumes for the target block (and
!                 maybe extract the block);
!              -1 means check the cell volumes of ALL blocks (& nothing else)
!     OMEGA     Over-relaxation factor for SLOR smoothing iterations
!     CONV      Number of orders of magnitude reduction in max. dX, dY, or dZ
!     EXPI1     Exponential decay factors for the foreground terms at the
!     EXPI2     6 faces
!     EXPJ1
!     EXPK1
!     EXPK2
!     BGPHI     'YYYY' turns on background control of interior spacing from
!               the I faces (Phi effects) at edges J1, J2, K1, K2 resp.
!     BGPSI     'YYYY' turns on background control of interior spacing from
!               the J faces (Psi effects) at edges I1, I2, K1, K2 resp.
!     BGOMG     'YYYY' turns on background control of interior spacing from
!               the K faces (Omega effects) at edges I1, I2, J1, J2, resp.
!     FGMODE    'NNNNNN' turns off foreground effects at all 6 faces;
!               'YYYYYY' turns on index-based orthogonality control;
!               any mixture of controls is permitted (e.g., 'NYYYNN')
!
!     PRINT     T displays 3-D smoothing iterations in gsmooth.out
!     FGLIM     Growth limiter for foreground terms
!     URFG      Under-relaxation factors for changes in foreground terms
!
!  History:
!
!     06/30/99  D.A.Saunders  Initial implementation
!     07/03/99    "      "    ELLIP3D's writes to unit 6 must come to the
!                             screen, so forget about a printable output file.
!     02/22/00    "      "    Suppress the full-size output file if the number
!                             of smoothing iterations is 0.  The optional copy
!                             of the target block remains an option.  Also,
!                             check ALL block cell volumes if ITMAX = -1.
!     02/10/14    "      "    Replaced CFD_IO_PACKAGE with xyzq_io package.
!                             Omitted the iblank option for now: it's easily
!                             incorporated if the need arises.
!     03/12/22    "      "    Handy test program for verifying pyramid_volume
!                             and hex_volume (alternatives to Antony Jameson's
!                             approach in subroutine cellvol, intended for
!                             unstructured grids).
!     03/15/22    "      "    Test alternative tet_vol and hex_vol similarly.
!     03/16/22    "      "    The 03/15/22 approach is preferable.  Retain it
!                             as corroboration of the Jameson scheme.  Monitor
!                             the volume summed over all blocks as a further
!                             comparison.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, CA.
!           Now with:       AMA, Inc. at ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure   ! Derived data type for one grid block
   use xyzq_io_module         ! PLOT3D file I/O package

   implicit none

!  Constants:

   integer, parameter :: &
      luninp = 1,        &    ! Control file
      lunin  = 2,        &    ! Input grid
      lunout = 3,        &    ! Output grid
      lunblk = 4,        &    ! Smoothed block if requested
      lunkbd = 5,        &
      luncrt = 6

   real, parameter :: &
      one = 1., zero = 0.

!  Variables:

   integer :: &
      iblock, ios, itmax, ivmin, jvmin, kvmin, l, n, nblock, ni, nj, nk, npts, &
      numneg

   real :: &
      conv,  dmax,  fglimit, omega, urfg, volmin, voltot1, voltot2, &
      voltotal1, voltotal2, expi1, expi2, expj1, expj2, expk1, expk2

   logical :: &
      all_block_volumes, extracopy, print, smoothing, &
      cell_centered, formatted_in, formatted_out

   character :: &
      answer*1, bgphi*4, bgpsi*4, bgomg*4, fgmode*6, text*132

   real, allocatable :: &
      s(:,:,:,:), xyz(:,:,:,:)

!  Composite data types:

   type (grid_type), pointer :: &
      blockgrid(:), ingrid(:)

!  Execution:

   open (luninp, file='gsmooth.inp', status='old', iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Unable to open file gsmooth.inp.'
      go to 99
   end if

   read (luninp, *)
   read (luninp, *) &
      itmax, omega, conv, expi1, expi2, expj1, expj2, expk1, expk2, &
      bgphi, bgpsi, bgomg, fgmode
   read (luninp, *)
   read (luninp, *) print, fglimit, urfg

   close (luninp)

   smoothing         = itmax > 0
   all_block_volumes = itmax < 0

!  Get the input grid name and formatting:

   call file_prompt (lunin, -luncrt, 'input grid', 'old', .true., &
                     formatted_in, ios)

!  Read its header:

   call xyz_header_io (1, lunin, formatted_in, nblock, ingrid, ios)
   if (ios /= 0) go to 99

!  Get the output grid name and formatting?

   if (smoothing) then
      call file_prompt (lunout, -luncrt, 'output grid', 'unknown', .true., &
                        formatted_out, ios)
      if (ios /= 0) go to 99

!     The output grid header info is the same as the input header:

      call xyz_header_io (2, lunout, formatted_out, nblock, ingrid, ios)
      if (ios /= 0) go to 99
   end if

!  Display the block dimensions to help selection of one:

   write (luncrt, '(/, (i6, 2x, 3i5))') &
      (n, ingrid(n)%ni, ingrid(n)%nj, ingrid(n)%nk, n = 1, nblock)

   if (smoothing) then
      write (luncrt, '(/, a)', advance='no') &
         ' Number of the block to be smoothed: '
      read  (lunkbd, *) n
   else if (.not. all_block_volumes) then
      write (luncrt, '(/, a)', advance='no') &
         ' Number of the block to be processed: '
      read  (lunkbd, *) n
   end if

   if (.not. all_block_volumes) then
      iblock = n
      ni = ingrid(n)%ni
      nj = ingrid(n)%nj
      nk = ingrid(n)%nk

      extracopy = nblock > 1

      if (extracopy) then
         write (luncrt, '(/, a)', advance='no') &
            ' Write the processed block to an extra file for plotting? (y/n): '
         read  (lunkbd, '(a)') answer
         extracopy = answer == 'y' .or. answer == 'Y'
      end if
   else
      extracopy = .false.
   end if

   if (extracopy) then ! Write the header before smoothing in case of trouble

      allocate (blockgrid(1))

      blockgrid(1)%ni = ingrid(iblock)%ni
      blockgrid(1)%nj = ingrid(iblock)%nj
      blockgrid(1)%nk = ingrid(iblock)%nk

      open (lunblk, file='modified_block.g', status='unknown')

      n = 1  ! # blocks has to be a variable

      call xyz_header_io (2, lunblk, .true., n, blockgrid, ios)

   end if

!  Transfer grid blocks from input to output, smoothing the specified one:

   voltotal1 = zero;  voltotal2 = zero

   do n = 1, nblock

      ni = ingrid(n)%ni
      nj = ingrid(n)%nj
      nk = ingrid(n)%nk;   npts = ni*nj*nk

      call xyz_allocate (ingrid(n), ios)

      call xyz_block_io (1, lunin, formatted_in, npts, &
                         ingrid(n)%x, ingrid(n)%y, ingrid(n)%z, ios)

      if (all_block_volumes) iblock = n

      if (n == iblock) then

         if (smoothing) then

            allocate (xyz(ni,nj,nk,3))  ! Too much trouble to change ELLIP3D
            allocate (s(ni,nj,nk,3))    ! For arc lengths in each index dir.

            xyz(:,:,:,1) = ingrid(n)%x
            xyz(:,:,:,2) = ingrid(n)%y
            xyz(:,:,:,3) = ingrid(n)%z

            call ellip3d (ni, nj, nk, 1, ni, 1, nj, 1, nk, xyz, s,    &
                          print, itmax, conv,  dmax,   omega,         &
                          bgphi, bgpsi, bgomg, fgmode, fglimit, urfg, &
                          expi1, expi2, expj1, expj2, expk1, expk2,   &
                          0, 999, 999)

            ingrid(n)%x = xyz(:,:,:,1)
            ingrid(n)%y = xyz(:,:,:,2)
            ingrid(n)%z = xyz(:,:,:,3)

            deallocate (xyz, s)

         end if

!        Check the cell volumes:

         allocate (s(ni-1,nj-1,nk-1,2))  ! For two sets of cell volume calcs.

         call check_vols (n, ni, nj, nk, one, &
                          ingrid(n)%x, ingrid(n)%y, ingrid(n)%z, &
                          s(:,:,:,1), luncrt, voltot1, volmin,    &
                          ivmin, jvmin, kvmin, numneg)

         voltotal1 = voltotal1 + voltot1

         call check_vol2 (n, ni, nj, nk, one, &
                          ingrid(n)%x, ingrid(n)%y, ingrid(n)%z, &
                          s(:,:,:,2), luncrt, voltot2, volmin,    &
                          ivmin, jvmin, kvmin, numneg)

         voltotal2 = voltotal2 + voltot2

         deallocate (s)

         if (extracopy) then ! Save (x,y,z)s of extra copy for plotting

            call xyz_block_io (2, lunblk, .true., npts, &
                               ingrid(n)%x, ingrid(n)%y, ingrid(n)%z, ios)
            close (lunblk)

            deallocate (blockgrid)

         end if

      end if

      if (smoothing) then  ! Save the smoothed block

         call xyz_block_io (2, lunout, formatted_out, npts, &
                            ingrid(n)%x, ingrid(n)%y, ingrid(n)%z, ios)
      end if

      deallocate (ingrid(n)%x, ingrid(n)%y, ingrid(n)%z)

   end do ! Next block

   write (luncrt, '(10x, a, es16.8)') &
      'Total volume of all blocks, Jameson formulation:   ', voltotal1, &
      'Total volume of all blocks, hex_vol/tet_vol method:', voltotal2

   deallocate (ingrid)
   close (lunout)

99 continue

   end program gsmooth

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine check_vols (n, il, jl, kl, hand, x, y, z, vol, lunwrt, &
                          voltot, volmin, ivmin, jvmin, kvmin, numneg)
!
!  Compute cell volumes; count volumes < 0. & print the least volume.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: n                 ! Block number
   integer, intent (in)  :: il, jl, kl, lunwrt
   real,    intent (in)  :: hand
   real,    intent (in)  :: x(il,jl,kl), y(il,jl,kl), z(il,jl,kl)
   integer, intent (out) :: ivmin, jvmin, kvmin, numneg
   real,    intent (out) :: vol(2:il,2:jl,2:kl), voltot, volmin

!  Local constants:

   real, parameter :: zero = 0.

!  Local variables:

   integer :: i, j, k, ncells
   logical, save :: first = .true.

!  Execution:
      
   call cellvol (il, jl, kl, x, y, z, vol, hand)

   volmin = vol(2,2,2)
   voltot = zero
   ivmin  = 2
   jvmin  = 2
   kvmin  = 2
   numneg = 0

   do k = 2, kl
      do j = 2, jl
         do i = 2, il
            voltot = voltot + vol(i,j,k)
            if (vol(i,j,k) < zero) numneg = numneg + 1
            if (vol(i,j,k) < volmin) then
               volmin = vol(i,j,k)
               ivmin = i - 1     ! "Lower left" indices
               jvmin = j - 1
               kvmin = k - 1
            end if
         end do
      end do
   end do

   ncells = (il - 1) * (jl - 1) * (kl - 1)

   if (first) then
      first = .false.
      write (lunwrt, '(a)')
   end if

   write (lunwrt, '(a, i4, a, es16.8, a, 3i4, a, es16.8, a, i9, a, i9, i5)') &
      ' Block:', n, '   Total vol.:', voltot, '   Cell', &
      ivmin, jvmin, kvmin, '  has least volume:', volmin, &
      '   # cells with negative volume:', numneg, '  out of', ncells, n

   end subroutine check_vols

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine check_vol2 (n, il, jl, kl, hand, x, y, z, vol, lunwrt, &
                          voltot, volmin, ivmin, jvmin, kvmin, numneg)
!
!  Compute cell volumes by the second method;
!  count volumes < 0. & print the least volume.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: n                 ! Block number
   integer, intent (in)  :: il, jl, kl, lunwrt
   real,    intent (in)  :: hand
   real,    intent (in)  :: x(il,jl,kl), y(il,jl,kl), z(il,jl,kl)
   integer, intent (out) :: ivmin, jvmin, kvmin, numneg
   real,    intent (out) :: vol(2:il,2:jl,2:kl), voltot, volmin

!  Local constants:

   real, parameter :: zero = 0.

!  Local variables:

   integer :: i, j, k, ncells

!  Execution:

   call cellvol2 (il, jl, kl, x, y, z, vol, hand)

   volmin = vol(2,2,2)
   voltot = zero
   ivmin  = 2
   jvmin  = 2
   kvmin  = 2
   numneg = 0

   do k = 2, kl
      do j = 2, jl
         do i = 2, il
            voltot = voltot + vol(i,j,k)
            if (vol(i,j,k) < zero) numneg = numneg + 1
            if (vol(i,j,k) < volmin) then
               volmin = vol(i,j,k)
               ivmin = i - 1     ! "Lower left" indices
               jvmin = j - 1
               kvmin = k - 1
            end if
         end do
      end do
   end do

   ncells = (il - 1) * (jl - 1) * (kl - 1)

   write (lunwrt, '(a, i4, a, es16.8, a, 3i4, a, es16.8, a, i9, a, i9, i5)') &
      ' Block:', n, '   Total vol.:', voltot, '   Cell', &
      ivmin, jvmin, kvmin, '  has least volume:', volmin, &
      '   # cells with negative volume:', numneg, '  out of', ncells, n

   end subroutine check_vol2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine cellvol (il, jl, kl, x, y, z, vol, hand)
!
!  Argument-driven portion of the Antony Jameson FLO87 METRIC routine for
!  calculating cell volumes, as needed for a grid quality check. HAND = +/-1
!  for RH/LH xyz resp. VP1, VP3, & VP5 negative signs added empirically.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in) :: &
      il, jl, kl

   real, intent (in) :: &
      x(il,jl,kl), y(il,jl,kl), z(il,jl,kl)

   real, intent (in) :: &
      hand

   real, intent (out) :: &
      vol(2:il,2:jl,2:kl)

!  Local constants:

   real, parameter :: &
      fourth = 1./ 4., sixth = 1./ 6., r8th = 1./ 8.

!  Local variables:

   integer :: &
      i, j, k, l, m, n
   real :: &
      factor, xp, yp, zp, vp1, vp2, vp3, vp4, vp5, vp6, &
      volpym, xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd

!  Statement function (more efficient than internal procedure):

   volpym (xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd) = &
         (xp - fourth * (xa + xb + xc + xd))                 &
      * ((ya - yc) * (zb - zd) - (za - zc) * (yb - yd))      &
      +  (yp - fourth * (ya + yb + yc + yd))                 &
      * ((za - zc) * (xb - xd) - (xa - xc) * (zb - zd))      &
      +  (zp - fourth * (za + zb + zc + zd))                 &
      * ((xa - xc) * (yb - yd) - (ya - yc) * (xb - xd))

!  VOLPYM = Volume of a pyramid with four-sided base, times 6
!  (adapted from Jameson code in FLO87).

!  Execution:

   factor = hand*sixth

   do k = 2, kl
      n = k - 1
      do j = 2, jl
         m = j - 1
         do i = 2, il
            l = i - 1
            xp = (x(i,j,k) +x(i,m,k) +x(i,m,n) +x(i,j,n)     &
                + x(l,j,k) +x(l,m,k) +x(l,m,n) +x(l,j,n))*r8th
            yp = (y(i,j,k) +y(i,m,k) +y(i,m,n) +y(i,j,n)     &
                + y(l,j,k) +y(l,m,k) +y(l,m,n) +y(l,j,n))*r8th
            zp = (z(i,j,k) +z(i,m,k) +z(i,m,n) +z(i,j,n)     &
                + z(l,j,k) +z(l,m,k) +z(l,m,n) +z(l,j,n))*r8th
            vp1 = -volpym (x(i,j,k), y(i,j,k), z(i,j,k),     &
                           x(i,m,k), y(i,m,k), z(i,m,k),     &
                           x(i,m,n), y(i,m,n), z(i,m,n),     &
                           x(i,j,n), y(i,j,n), z(i,j,n))
            vp2 =  volpym (x(l,j,k), y(l,j,k), z(l,j,k),     &
                           x(l,m,k), y(l,m,k), z(l,m,k),     &
                           x(l,m,n), y(l,m,n), z(l,m,n),     &
                           x(l,j,n), y(l,j,n), z(l,j,n))
            vp3 = -volpym (x(i,j,k), y(i,j,k), z(i,j,k),     &
                           x(i,j,n), y(i,j,n), z(i,j,n),     &
                           x(l,j,n), y(l,j,n), z(l,j,n),     &
                           x(l,j,k), y(l,j,k), z(l,j,k))
            vp4 =  volpym (x(i,m,k), y(i,m,k), z(i,m,k),     &
                           x(i,m,n), y(i,m,n), z(i,m,n),     &
                           x(l,m,n), y(l,m,n), z(l,m,n),     &
                           x(l,m,k), y(l,m,k), z(l,m,k))
            vp5 = -volpym (x(i,j,k), y(i,j,k), z(i,j,k),     &
                           x(l,j,k), y(l,j,k), z(l,j,k),     &
                           x(l,m,k), y(l,m,k), z(l,m,k),     &
                           x(i,m,k), y(i,m,k), z(i,m,k))
            vp6 =  volpym (x(i,j,n), y(i,j,n), z(i,j,n),     &
                           x(l,j,n), y(l,j,n), z(l,j,n),     &
                           x(l,m,n), y(l,m,n), z(l,m,n),     &
                           x(i,m,n), y(i,m,n), z(i,m,n))
            vol(i,j,k) = (vp1 + vp2 + vp3 + vp4 + vp5 + vp6)*factor
         end do
      end do
   end do

   end subroutine cellvol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine cellvol2 (il, jl, kl, x, y, z, vol, hand)
!
!  Adaptation of cellvol routine to test hex_vol (five tetrahedra) against
!  Antony Jameson's six-pyramid method.  This one can't detect a negative vol.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in) :: &
      il, jl, kl

   real, intent (in) :: &
      x(il,jl,kl), y(il,jl,kl), z(il,jl,kl)

   real, intent (in) :: &
      hand

   real, intent (out) :: &
      vol(2:il,2:jl,2:kl)

!  Local variables:

   integer :: &
      i, j, k, l, m, n
   real, dimension (3) :: &
      p1, p2, p3, p4, p5, p6, p7, p8

!  Execution:

   do k = 2, kl
      n = k - 1
      do j = 2, jl
         m = j - 1
         do i = 2, il
            l = i - 1
            p1(1) = x(l,m,n);  p1(2) = y(l,m,n);  p1(3) = z(l,m,n)
            p2(1) = x(i,m,n);  p2(2) = y(i,m,n);  p2(3) = z(i,m,n)
            p3(1) = x(i,j,n);  p3(2) = y(i,j,n);  p3(3) = z(i,j,n)
            p4(1) = x(l,j,n);  p4(2) = y(l,j,n);  p4(3) = z(l,j,n)
            p5(1) = x(l,m,k);  p5(2) = y(l,m,k);  p5(3) = z(l,m,k)
            p6(1) = x(i,m,k);  p6(2) = y(i,m,k);  p6(3) = z(i,m,k)
            p7(1) = x(i,j,k);  p7(2) = y(i,j,k);  p7(3) = z(i,j,k)
            p8(1) = x(l,j,k);  p8(2) = y(l,j,k);  p8(3) = z(l,j,k)

            call hex_vol (p1, p2, p3, p4, p5, p6, p7, p8, vol(i,j,k))
         end do
      end do
   end do

   end subroutine cellvol2
