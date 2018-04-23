!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      program plug_interp
!
!     Description:
!
!        PLUG_INTERP performs a specialized interpolation of the flow field on
!     the outer boundary of the blocks around a Shuttle wing leading edge plug
!     repair.  The same program proves handy for a protruding gap filler also.
!
!        This version applies what was originally done on the outer boundary
!     only to the interior of the "collar" blocks as well, to overcome what
!     seems to be the cause of slow convergence (especially for the gap filler
!     case).  The surface of the plug/gap filler is NOT touched.
!
!        The topology (via Gridgen scripting followed by the RADIAL_INTERP and
!     COMBINE_BLOCKS procedures) is assumed to be as follows, with the thin
!     "collar" blocks always appearing last as listed in template.inp.2 used
!     by the PREPARE_RAPID_ANALYSIS scheme.
!
!     Square-cornered plug:
!
!        Blocks 1:5  are on the plug (circular before projecting to the surface)
!           "   6:9  surround 1:5 but are offset from the OML by the plug height
!           "  10:13 surround the plug on the OML with about 1 mm height
!
!        Originally, blocks 6:9 corresponded to 10:13 with consistent (i,k)
!        indexing in the outer boundary faces.  However, Chun Tang changed his
!        Gridgen script in a way that makes corresponding pairs uncertain.
!        Therefore: determine the connectivity from the template.con file that
!        program TEMPLATE generates for the combined blocks.
!
!     Rounded-corner plug:
!
!        This has an extra set of 4 blocks off half of the rounded corner, so
!        the collar blocks are 14:17.
!
!     Protruding gap filler:
!
!        The collar blocks are 6:9.  Use 'plug' in template.inp.2.
!        The surface j direction is assumed to be 1 at the gap filler.
!
!     The following explanation applies to the original square corner plug case.
!     No changes other than the control inputs are needed for the other cases.
!
!        RADIAL_INTERP is given the new surface that contains the plug patches
!     and the "tops" of the thin surrounding patches.  This means the baseline
!     flow solution is interpolated into blocks 6:9.  COMBINE_BLOCKS uses its
!     specialized option to initialize the flow for the thin blocks 10:13 (as
!     implemented originally for cavity blocks).  The only places that aren't
!     good enough are on the outer faces of the thin blocks.  This program
!     carefully redistributes the slightly shifted flow on the outer faces of
!     blocks 6:9 so that it correctly covers the outer faces of blocks 10:13 as
!     well.
!
!        (Later:) It also reverses the scaling of species densities, but cannot
!                 recover the zeroed velocity components.
!        (Later still:)  The same redistribution is now applied for every j
!                        except j = 1 on the protuberance.)
!        (Even later:)   Include j = 1 so that for plugs the i and j directions
!                        are not significant.
!
!        Cell-centered grid and flow files should be specified at the prompts.
!
!     Procedures:
!
!        XYZQ_IO package  I/O utilities for PLOT3D grid and function files
!        <various numerics utilities also>
!
!     History:
!
!        02/24/06  D. Saunders  Initial adaptation of EXTRACT_BLOCKS.
!        03/10/06   "      "    Chun reordered the blocks, so make use of
!                               template.con to find corresponding pairs.
!        03/15/06   "      "    Undo COMBINE_BLOCKS' reduction of pressure for
!                               the last set of blocks (suited to cavities, not
!                               plug repairs).  We cannot recover the zeroed
!                               velocity components, though.
!        08/28/06   "      "    Apply the redistribution at every (i,j) except
!                               j = 1 (vs. just (i,nj) originally).
!        09/12/06   "      "    High protruding gap fillers tend to diverge.
!                               Try varying the velocity components with j.
!        03/13/07   "      "    Todd White's plug case with j not necessarily in
!                               the spoke direction (out from the plug) failed
!                               because of omitting j = 1, assumed to be at the
!                               protuberance wall.  Including j = 1 should not
!                               hurt, and allows arbitrary i/j convention.
!        02/05/10   "      "    The decay of velocity components towards the
!                               gap filler was also decaying temperature!
!
!     Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Modules:

      use grid_block_structure
      use xyzq_io_module

      implicit none

!     Constants:

      integer, parameter :: &
         luncon   = 1,      & ! For template.con
         lunkbd   = 5,      &
         luncrt   = 6,      &
         lunxyz   = 7,      &
         lunq_in  = 8,      &
         lunq_out = 9,      &
         m        = 4       ! # blocks surrounding the plug

      real, parameter ::    &
         eps    = 1.e-4,    & ! To check for indexing mismatch
         one    = 1.0,      &
         pscale = 10.         ! Inverse of what COMBINE_BLOCKS uses

      logical, parameter :: &
         false = .false.,   &
         true  = .true.

      character, parameter :: &
         mono * 1 = 'M'     ! Monotonic spline interpolation

!     Variables:

      integer :: &
         i, ibp, ibt, ibv, ios, j, k, n, nblocks, &
         nip, niv, nj, njp, njv, nkp, nkv, nspecies, num_q, &
         ib_thin(m)

      real :: &
         dx, dy, dz, r, stotal

      logical :: &
         allocate_radial, cell_centered, formatted_in, formatted_out,          &
         gap_filler, new

      real, allocatable, dimension (:) :: &
         fd, sd, xd, yd, zd, &
         fp, sp, xp, yp, zp, &
         fv, sv, xv, yv, zv, unused

!     Derived data type:

      type (grid_type), pointer, dimension (:) :: &
         xyzq

!     Execution:

!     Determine the block pairs that are to have their outer faces matched up.

      open (luncon, file='template.con', status='old', iostat=ios)

      if (ios /= 0) then
         write (luncrt, '(/, 2a)') ' Cannot open template.con file',           &
                                   ' to determine corresponding block pairs.'
         go to 999
      end if

      read (luncon, *) nblocks
      do ibv = 1, nblocks - 2*m
         read (luncon, '(a)') ! Skip a line
      end do

      k = 1
      do ibv = nblocks - 2*m + 1, nblocks - m ! Initially ibv = 6, 7, 8, 9
         read (luncon, '(i5,i3,3i2,3(1x,i3,3i2),i4)') &
            j,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,ib_thin(k) 
         if (j /= ibv) then
            write (luncrt, '(a, 2i4)') &
               ' Apparent connectivity error.  ibv, j:', ibv, j
            go to 999
         end if
         k = k + 1
      end do

      close (luncon)

      call file_prompt (lunxyz,  -luncrt, 'input cell-centered grid', 'old',   &
                        true, formatted_in, ios)
      if (ios /= 0) go to 999

      call file_prompt (lunq_in, -luncrt, 'input cell-centered flow', 'old',   &
                        false, formatted_in, ios)
      if (ios /= 0) go to 999

      call file_prompt (lunq_out, -luncrt, 'output flow', 'unknown', true,     &
                        formatted_out, ios)
      if (ios /= 0) go to 999

!     Read all blocks of the input grid and flow field.  The "cell_centered"
!     argument isn't relevant: the grid and flow should both be cell-centered.

      call xyzq_read (lunxyz, lunq_in, formatted_in, nblocks, num_q,           &
                      cell_centered, xyzq, ios)
      if (ios /= 0) go to 999

!     Display the block dimensions to reassure the user:

      write (luncrt, '(/, a, i3, a, /)') &
         ' # flow variables found:', num_q, '    Block dimensions:'
      write (luncrt, '(4(i6, 2X, 3i5))') &
         (i, xyzq(i)%ni, xyzq(i)%nj, xyzq(i)%nk, i = 1, nblocks)

!     Process the last m block pairs, where m is 4 initially and nblocks is 13.

!     First, set up for undoing some of COMBINE_BLOCKS' cavity-related changes:

      if (num_q <= 9) then  ! Assume only 1 temperature
         nspecies = num_q - 4
      else                  ! Assume 2 temperatures
         nspecies = num_q - 5
      end if

      gap_filler = nblocks == 9 .and. num_q == 9

      allocate_radial = true
      ibt = 1  ! Pointer into ib_thin(*)

      do ibv = nblocks - 2*m + 1, nblocks - m ! Initially ibv = 6, 7, 8, 9 for
                                              ! plugs; 2, 3, 4, 5 (gap fillers)
         niv = xyzq(ibv)%ni
         njv = xyzq(ibv)%nj
         nkv = xyzq(ibv)%nk

         ibp = ib_thin(ibt)  ! ibp = around-plug block <-> main volume block ibv
         ibt = ibt + 1

         nip = xyzq(ibp)%ni
         njp = xyzq(ibp)%nj
         nkp = xyzq(ibp)%nk

         write (luncrt, '(a, i2, a, i3)') ' Processing blocks', ibv, ' &', ibp

         if (nip /= niv .or. njp /= njv) then
            write (luncrt, '(/, a, 2i5, /, a, /, (2i5))') &
               ' Mismatched (i,j) dimensions for blocks', ibv, ibp, &
               '   ni   nj', niv, njv, nip, njp
            go to 999
         end if

         nj = njv  ! = njp
         j  = nj / 2

         if (abs (xyzq(ibp)%x(2,j,nkp) - xyzq(ibv)%x(2,j,1)) > eps .or.        &
             abs (xyzq(ibp)%y(2,j,nkp) - xyzq(ibv)%y(2,j,1)) > eps .or.        &
             abs (xyzq(ibp)%z(2,j,nkp) - xyzq(ibv)%z(2,j,1)) > eps) then
            write (luncrt, '(/, a, 2i3, /, (1p, 3e16.8))') &
               ' Likely index direction mismatch for blocks ', ibp, ibv,       &
               xyzq(ibp)%x(2,j,nkp),xyzq(ibp)%y(2,j,nkp),xyzq(ibp)%z(2,j,nkp), &
               xyzq(ibv)%x(2,j,1),  xyzq(ibv)%y(2,j,1),  xyzq(ibv)%z(2,j,1)   
            write (luncrt, '(a)') ' template.con should tell more.'
            go to 999
         end if

!        Undo COMBINE_BLOCKS' scaling of species densities in the thin block:

         xyzq(ibp)%q(1:nspecies,:,:,:) = pscale * xyzq(ibp)%q(1:nspecies,:,:,:)
  
!        Work-space for one radial line:

         if (allocate_radial) then

            allocate_radial = false

            allocate (xv(nkv), yv(nkv), zv(nkv), sv(nkv), &  ! Main vol. coords.
                      xp(nkp), yp(nkp), zp(nkp), sp(nkp), &  ! Plug blk. coords.
                      xd(nkv), yd(nkv), zd(nkv), sd(nkv))    ! "Data" (shifted)

            allocate (fd(nkv), fp(nkp), fv(nkv))  ! "Data" for one flow variable

            k = max (nkp, nkv)
            allocate (unused(k))                  ! Unused derivatives
                                                  ! and for interpolated values
         end if

!        For each (i,j) of the (i,nj,k) face pair.  Include j = 1 in case it
!        is not at the plug wall (but it has to be below for a gap filler).

         do j = 1, nj

            do i = 1, nip

               do k = 1, nkp
                  xp(k) = xyzq(ibp)%x(i,j,k)
                  yp(k) = xyzq(ibp)%y(i,j,k)
                  zp(k) = xyzq(ibp)%z(i,j,k)
               end do

               call chords3d (nkp, xp, yp, zp, false, stotal, sp)

               dx = xp(1) - xp(nkp)  ! Undo RADIAL_INTERP's shift
               dy = yp(1) - yp(nkp)
               dz = zp(1) - zp(nkp)

               do k = 1, nkv
                  xv(k) = xyzq(ibv)%x(i,j,k);  xd(k) = xv(k) + dx
                  yv(k) = xyzq(ibv)%y(i,j,k);  yd(k) = yv(k) + dy
                  zv(k) = xyzq(ibv)%z(i,j,k);  zd(k) = zv(k) + dz
               end do

               call chords3d (nkv, xd, yd, zd, false, stotal, sd)

               sv(:) = sd(:) + sp(nkp)  ! Offset target arcs in the main block

!!!            write (6, '(2a)') '           xp           yp           zp',    &
!!!                              '           sp'
!!!            write (6, '(1p, 4e13.5)') &
!!!               (xp(k), yp(k), zp(k), sp(k), k = 1, nkp)

!!!            write (6, '(2a)') '           xv           yv           zv',    &
!!!                              '           sv           sd'
!!!            write (6, '(1p, 5e13.5)') &
!!!               (xv(k), yv(k), zv(k), sv(k), sd(k), k = 1, nkv)

!              Interpolate the flow for each k of the plug block:

               do n = 1, num_q

                  fd(:) = xyzq(ibv)%q(n,i,j,:)

                  call lcsfit (nkv, sd, fd, true, mono, nkp, sp, fp, unused)

                  xyzq(ibp)%q(n,i,j,:) = fp(:)

               end do

!              Interpolate the flow for each k of the main block:

               k = nkv - 1

               do n = 1, num_q

                  fd(:) = xyzq(ibv)%q(n,i,j,:)

                  call lcsfit (nkv, sd, fd, true, mono, k, sv, fv, unused)

                  xyzq(ibv)%q(n,i,j,1:k) = fv(1:k)

               end do

            end do ! Next i

         end do ! Next j

         if (gap_filler) then

!           Decay velocities towards a gap filler:

            do j = 1, nj - 1
               r = real (j - 1) / real (nj - 1)
               do k = 1, nkp
                  do i = 1, nip
                     xyzq(ibp)%q(6:8,i,j,k) = r *xyzq(ibp)%q(6:8,i,nj,k)
                  end do
               end do
            end do

         end if

      end do ! Next pair of blocks

      deallocate (xv, yv, zv, sv, fv, xp, yp, zp, sp, fp, xd, yd, zd, sd, fd,  &
                  unused)

      if (gap_filler) then

!        Special treatment of the narrow block above the gap filler, to remove
!        the starting guess discontinuity at the shock caused by its height:

         do k = 1, xyzq(1)%nk
            do j = 1, xyzq(1)%nj
               do i = 1, xyzq(1)%ni
                  xyzq(1)%q(:,i,j,k) = xyzq(2)%q(:,2,2,k)
               end do
            end do
         end do

         xyzq(1)%q(6:8,:,:,1) = 0.  ! Top-of-gapfiller velocities [redundant?]

      end if

!     Write the updated flow field:

      call q_write (lunq_out, formatted_out, nblocks, num_q, xyzq, ios)

      if (ios /= 0) go to 999

      do i = 1, nblocks
         deallocate (xyzq(i)%x, xyzq(i)%y, xyzq(i)%z, xyzq(i)%q)
      end do

      deallocate (xyzq)

  999 continue

! *** stop ! Avoid system dependencies.

      end program plug_interp
