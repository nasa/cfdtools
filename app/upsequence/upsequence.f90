!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      program upsequence
!
!     Description:
!
!        UPSEQUENCE is a specialized utility for treating coarse-grid forms of
!     the multiblock flow solutions for cavity or plug configurations produced
!     by NASA Ames procedures for rapid local analysis of damage or repair with
!     known grid topology.  It simply copies the coarse grid solution to the
!     finer grid cells indicated by iseq, jseq, kseq (same for all blocks is
!     assumed), and then transcribes the outer boundary flow from the original
!     fine-grid solution starting guess as part of speeding the calculation on
!     the finer grid.
!
!        All flow solution files should contain cell-centered data with halo
!     cells included.
!
!     Assumptions for cavity cases:
!
!        o  6 blocks, with block 6 in the cavity (no outer boundary)
!
!     Assumptions for plug cases:
!
!        o  13 blocks, with blocks 1:5 above the plug (kmax boundary to treat)
!
!        o  Blocks 6:13 have their jmax faces at outer local boundaries
!
!     Assumptions for gap filler cases:
!
!        o  9 blocks, with blocks 1:5 above the filler (kmax boundary to treat)
!
!        o  Blocks 2:9 have their jmax faces at outer local boundaries
!
!     Procedures:
!
!        XYZQ_IO package  I/O utilities for PLOT3D grid and function files
!
!     History:
!
!        06/06/2006  D. Saunders  Adaptation of THIN_FLOW.
!        06/07/2006   "       "   Testing via THIN_GRID was misleading: more
!                                 care is needed, treating interior cells the
!                                 expected way and halo cells separately.
!        09/13/2006   "       "   Added the gap-filler case.
!
!     Author:  David Saunders, ELORET/NASA Ames Research Ctr., Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Modules:

      use grid_block_structure
      use xyzq_io_module

      implicit none

!     Constants:

      integer, parameter :: &
         luninc = 1,        &  ! For the coarse flow input
         luninf = 2,        &  ! For the fine original flow input
         lunout = 3,        &  ! For the upsequenced flow
         lunkbd = 5,        &
         luncrt = 6,        &
         m      = 4            ! # blocks surrounding a plug or gap filler

      logical, parameter :: &
         true   = .true.

!     Variables:

      integer :: &
         ib, ic, if, im1, inc, ios, iq, jc, jf, jm1, jnc, kc, kf, km1, knc,    &
         nblocks, nf, mi, mj, mk, ni, nj, nk

      logical :: &
         cavity, formatted_inc, formatted_inf, formatted_out, gapfiller, plug

!     Composite data types:

      type (grid_type), pointer, dimension (:) :: &
         gridinc, gridinf, gridout

!     Execution:
!     ----------

      call file_prompt (luninc, -luncrt, 'input coarse flow solution', 'old',  &
                        true, formatted_inc, ios)
      if (ios /= 0) go to 999

!     Read the number of blocks, allocate them, and read dimensions:

      if (formatted_inc) then
         read (luninc, *, iostat=ios) nblocks
      else
         read (luninc,    iostat=ios) nblocks
      end if

      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble reading the number of blocks.'
         go to 999
      end if

      cavity    = nblocks == 6
      plug      = nblocks == 13
      gapfiller = nblocks == 9

      rewind (luninc)

      allocate (gridinc(nblocks))

      call q_header_io (1, luninc, formatted_inc, nblocks, nf, gridinc, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble reading coarse flow header.'
         go to 999
      end if

      call file_prompt (luninf, -luncrt, 'input original fine flow', 'old',  &
                        true, formatted_inf, ios)
      if (ios /= 0) go to 999

      allocate (gridinf(nblocks))

      call q_header_io (1, luninf, formatted_inf, nblocks, nf, gridinf, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble reading fine flow header.'
         go to 999
      end if

      call file_prompt (lunout, -luncrt, 'output upsequenced flow',         &
                        'unknown', true, formatted_out, ios)
      if (ios /= 0) go to 999

      allocate (gridout(nblocks))

      do ib = 1, nblocks
         gridout(ib)%mi = gridinf(ib)%mi
         gridout(ib)%mj = gridinf(ib)%mj
         gridout(ib)%mk = gridinf(ib)%mk
      end do

      call q_header_io (2, lunout, formatted_out, nblocks, nf, gridout, ios)

      if (ios /= 0) go to 999

!     Upsequence the blocks one at a time:

      do ib = 1, nblocks

!        Allocate and read a coarse block:

         call q_allocate (gridinc(ib), nf, ios)
         if (ios /= 0) go to 999

         mi = gridinc(ib)%mi;  mj = gridinc(ib)%mj;  mk = gridinc(ib)%mk

         call q_block_io (1, luninc, formatted_inc, nf, mi, mj, mk,            &
                          gridinc(ib)%q, ios)
         if (ios /= 0) go to 999

!        Allocate and read a fine block:

         call q_allocate (gridinf(ib), nf, ios)
         if (ios /= 0) go to 999

         ni = gridinf(ib)%mi;  nj = gridinf(ib)%mj;  nk = gridinf(ib)%mk

         call q_block_io (1, luninf, formatted_inf, nf, ni, nj, nk,            &
                          gridinf(ib)%q, ios)
         if (ios /= 0) go to 999

!        The increments are probably 2, but maybe not in the k direction?
!        The halo cells always make 2 extra cells treated differently.

         inc = (ni - 2) / (mi - 2);  im1 = inc - 1
         jnc = (nj - 2) / (mj - 2);  jm1 = jnc - 1
         knc = (nk - 2) / (mk - 2);  km1 = knc - 1

         write (6, '(a, i3, a, 3i4)') &
            '   Block', ib, ':  apparent sequencing:', inc, jnc, knc

!        Allocate the fine output block and copy coarse cells the way FCONVERT
!        (presumably) does:

         call q_allocate (gridout(ib), nf, ios)
         if (ios /= 0) go to 999

!        Block corner halo cells:

         kf = 1
         do kc = 1, mk, mk - 1
            jf = 1
            do jc = 1, mj, mj - 1
               if = 1
               do ic = 1, mi, mi - 1
                  gridout(ib)%q(:,if,jf,kf) = gridinc(ib)%q(:,ic,jc,kc)
                  if = ni
               end do
               jf = nj
            end do
            kf = nk
         end do

!        Interior edge cells in the i direction:

         kf = 1
         do kc = 1, mk, mk - 1
            jf = 1
            do jc = 1, mj, mj - 1
               do ic = 2, mi - 1
                  if = ic * inc - 2*im1
                  do iq = 1, nf
                     gridout(ib)%q(iq,if:if+im1,jf,kf) = &
                        gridinc(ib)%q(iq,ic,jc,kc)
                  end do
               end do
               jf = nj
            end do
            kf = nk
         end do

!        Interior edge cells in the j direction:

         kf = 1
         do kc = 1, mk, mk - 1
            if = 1
            do ic = 1, mi, mi - 1
               do jc = 2, mj - 1
                  jf = jc * jnc - 2*jm1
                  do iq = 1, nf
                     gridout(ib)%q(iq,if,jf:jf+jm1,kf) = &
                        gridinc(ib)%q(iq,ic,jc,kc)
                  end do
               end do
               if = ni
            end do
            kf = nk
         end do

!        Interior edge cells in the k direction:

         jf = 1
         do jc = 1, mj, mj - 1
            if = 1
            do ic = 1, mi, mi - 1
               do kc = 2, mk - 1
                  kf = kc * knc - 2*km1
                  do iq = 1, nf
                     gridout(ib)%q(iq,if,jf,kf:kf+km1) = &
                        gridinc(ib)%q(iq,ic,jc,kc)
                  end do
               end do
               if = ni
            end do
            jf = nj
         end do

!        Interior face cells at k = 1 and kmax:

         kf = 1
         do kc = 1, mk, mk - 1
            do jc = 2, mj - 1
               jf = jc * jnc - 2*jm1
               do ic = 2, mi - 1
                  if = ic * inc - 2*im1
                  do iq = 1, nf
                     gridout(ib)%q(iq,if:if+im1,jf:jf+jm1,kf) = &
                        gridinc(ib)%q(iq,ic,jc,kc)
                  end do
               end do
            end do
            kf = nk
         end do

!        Interior face cells at j = 1 and jmax:

         jf = 1
         do jc = 1, mj, mj - 1
            do kc = 2, mk - 1
               kf = kc * knc - 2*km1
               do ic = 2, mi - 1
                  if = ic * inc - 2*im1
                  do iq = 1, nf
                     gridout(ib)%q(iq,if:if+im1,jf,kf:kf+km1) = &
                        gridinc(ib)%q(iq,ic,jc,kc)
                  end do
               end do
            end do
            jf = nj
         end do

!        Interior face cells at i = 1 and imax:

         if = 1
         do ic = 1, mi, mi - 1
            do kc = 2, mk - 1
               kf = kc * knc - 2*km1
               do jc = 2, mj - 1
                  jf = jc * jnc - 2*jm1
                  do iq = 1, nf
                     gridout(ib)%q(iq,if,jf:jf+jm1,kf:kf+km1) = &
                        gridinc(ib)%q(iq,ic,jc,kc)
                  end do
               end do
            end do
            if = ni
         end do

!        Interior volume cells:

         do kc = 2, mk - 1
            kf = kc * knc - 2*km1
            do jc = 2, mj - 1
               jf = jc * jnc - 2*jm1
               do ic = 2, mi - 1
                  if = ic * inc - 2*im1
                  do iq = 1, nf
                     gridout(ib)%q(iq,if:if+im1,jf:jf+jm1,kf:kf+km1) = &
                        gridinc(ib)%q(iq,ic,jc,kc)
                  end do
               end do
            end do
         end do

!        Recover the outer boundary flow from the input original fine grid:

         if (cavity) then

            if (ib < 5) then ! jmin boundary

               do kf = 1, nk
                  do if = 1, ni
                     gridout(ib)%q(:,if,1,kf) = gridinf(ib)%q(:,if,1,kf)
                  end do
               end do

            end if

            if (ib < 6) then ! kmax boundary

               do jf = 1, nj
                  do if = 1, ni
                     gridout(ib)%q(:,if,jf,nk) = gridinf(ib)%q(:,if,jf,nk)
                  end do
               end do

            end if

         else if (plug) then

            if (ib > 5)  then ! jmax boundary

               do kf = 1, nk
                  do if = 1, ni
                     gridout(ib)%q(:,if,nj,kf) = gridinf(ib)%q(:,if,nj,kf)
                  end do
               end do

            end if

            if (ib < 10) then ! kmax boundary

               do jf = 1, nj
                  do if = 1, ni
                     gridout(ib)%q(:,if,jf,nk) = gridinf(ib)%q(:,if,jf,nk)
                  end do
               end do

            end if

         else ! Gap filler case

            if (ib > 1)  then ! jmax boundary

               do kf = 1, nk
                  do if = 1, ni
                     gridout(ib)%q(:,if,nj,kf) = gridinf(ib)%q(:,if,nj,kf)
                  end do
               end do

            end if

            if (ib < 6) then ! kmax boundary

               do jf = 1, nj
                  do if = 1, ni
                     gridout(ib)%q(:,if,jf,nk) = gridinf(ib)%q(:,if,jf,nk)
                  end do
               end do

            end if

         end if

!        Save this output block:

         call q_block_io (2, lunout, formatted_out, nf, ni, nj, nk,           &
                          gridout(ib)%q, ios)
         if (ios /= 0) go to 999
         
         deallocate (gridinc(ib)%q, gridinf(ib)%q, gridout(ib)%q)

      end do

  999 continue

! *** stop ! Avoid system dependencies.

      end program upsequence
