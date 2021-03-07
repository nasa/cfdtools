!*******************************************************************************
!
      program compare_flows
!
!     COMPARE_FLOWS reads two PLOT3D function files and locates the
!     largest point difference in the specified subscript range(s).
!     The grids may be 2- or 3D and are assumed to be multiblock, formatted
!     or not.
!
!     07/21/03  D.A.Saunders   Original COMPARE_BLOCKS (grid coordinates).
!     09/21/09    "     "      COMPARE_FLOWS derived from COMPARE_BLOCKS.
!     05/20/11    "     "      The reading has been wrong all along!
!                              Storing all functions of all blocks contiguously
!                              cannot also put all functions of a given block
!                              together, so subscripting during I/O is required.
!
!     David Saunders, ELORET Corp./NASA Ames Research Center, Moffett Field, CA
!
!*******************************************************************************

      implicit none

!     Constants:

      integer, parameter :: &
         lun1   = 1,        &
         lun2   = 2,        &
         luncrt = 6,        &
         lunkbd = 5,        &
         lunout = 3

      real, parameter ::    &
         zero = 0.

      character, parameter :: &
         form * 11 = 'unformatted'

!     Allocatable variables:

      integer, allocatable, dimension (:) ::   &
         ip, nf, ni, nj, nk, nf2, ni2, nj2, nk2

      real,    allocatable, dimension (:,:) :: &
         f1, f2

!     Local variables:

      integer :: &
         i, j, k, l, m, i1, i2, j1, j2, k1, k2, l1, l2, &
         iostat, imax, jmax, kmax, lmax, &
         n, nblocks, nblocks2, ncompare, np, npoints, numf

      real :: &
         d, dmax, drms

      logical :: &
         formatted, print, twod

      character :: &
         answer * 1, filename * 96

!     Execution:
!     ----------

      write (luncrt, '(a)', advance='no') &
         ' 2D or 3D files? [2|3]: '
      read  (lunkbd, *) n
      twod = n == 2

      write (luncrt, '(a)', advance='no') &
         ' Formatted or unformatted? [f|u]: '
      read  (lunkbd, '(a)') answer
      formatted = answer == 'f' .or. answer == 'F'

      write (luncrt, '(a)', advance='no') ' 1st function file: '
      read  (lunkbd, '(a)') filename

      k1 = 1
      if (formatted) k1 = 3

      open (unit=lun1, file=filename, status='old', form=form(k1:11))

      write (luncrt, '(a)', advance='no') ' 2nd function file: '
      read  (lunkbd, '(a)') filename

      if (formatted) then
         read (lun1, *, iostat=iostat) nblocks
      else
         read (lun1,    iostat=iostat) nblocks
      end if

      if (iostat /= 0) then
         write (luncrt, '(/, a, i6)') &
            ' Error reading # blocks for file 1.  iostat: ', iostat
         go to 999
      end if

      allocate (nf(nblocks),  ni(nblocks),  nj(nblocks),  nk(nblocks), &
               nf2(nblocks), ni2(nblocks), nj2(nblocks), nk2(nblocks))

      if (twod) then
         if (formatted) then
            read (lun1, *, iostat=iostat) &
               (ni(l), nj(l), nf(l), l = 1, nblocks)
         else
            read (lun1,    iostat=iostat) &
               (ni(l), nj(l), nf(l), l = 1, nblocks)
         end if
         nk(:) = 1
      else
         if (formatted) then
            read (lun1, *, iostat=iostat) &
               (ni(l), nj(l), nk(l), nf(l), l = 1, nblocks)
         else
            read (lun1,    iostat=iostat) &
               (ni(l), nj(l), nk(l), nf(l), l = 1, nblocks)
         end if
      end if

      if (iostat /= 0) then
         write (luncrt, '(/, a, i6)') &
            ' Error reading block dimensions, file 1.  iostat: ', iostat
         go to 999
      end if

      write (luncrt, '(/, a)') ' Block #, dimensions & nf for file 1: '
      if (twod) then
         write (luncrt, '(1x, i4, 2i7, i4)') &
            (l, ni(l), nj(l), nf(l), l = 1, nblocks)
      else
         write (luncrt, '(1x, i4, 3i7, i4)') &
            (l, ni(l), nj(l), nk(l), nf(l), l = 1, nblocks)
      end if

!     Pack the blocks in memory.  Determine the start of each block:

      allocate (ip(nblocks + 1))

      ip(1) = 1
      do l = 1, nblocks
         np = ni(l) * nj(l) * nk(l)
         ip(l+1) = ip(l) + np
      end do

      numf    = nf(1)           ! Assumed constant
      npoints = ip(nblocks+1) - 1

      allocate (f1(npoints,numf), f2(npoints,numf), stat=iostat)

      if (iostat /= 0) then
         write (luncrt, '(/, a, i4, i12)') &
            'Trouble allocating function storage. # functions & # points:', &
            numf, npoints
         go to 999
      end if

      do l = 1, nblocks

         i1 = ip(l)
         np = ip(l+1) - i1

         call rd_blk (lun1, formatted, npoints, numf, i1, np, f1, iostat)

         if (iostat /= 0) then
            write (luncrt, '(/, a, i5, a, i5)') &
               ' Trouble reading file 1.  Block #:', l, '.  iostat:', iostat
            go to 999
         end if

      end do

      close (lun1)

      open (unit=lun2, file=filename, status='old', form=form(k1:11))

      if (formatted) then
         read (lun2, *, iostat=iostat) nblocks2
      else
         read (lun2,    iostat=iostat) nblocks2
      end if

      if (iostat /= 0) then
         write (luncrt, '(/, a, i6)') &
            ' Error reading nblocks for file 2.  iostat: ', iostat
         go to 999
      end if

      if (nblocks2 /= nblocks) then
         write (luncrt, '(a, /, i8)') ' # blocks in file 2:', nblocks2
         go to 999
      end if

      if (twod) then
         if (formatted) then
            read (lun2, *, iostat=iostat) &
               (ni2(l), nj2(l), nf2(l), l = 1, nblocks)
         else
            read (lun2,    iostat=iostat) &
               (ni2(l), nj2(l), nf2(l), l = 1, nblocks)
         end if
         nk2(:) = 1
      else
         if (formatted) then
            read (lun2, *, iostat=iostat) &
               (ni2(l), nj2(l), nk2(l), nf2(l), l = 1, nblocks)
         else
            read (lun2,    iostat=iostat) &
               (ni2(l), nj2(l), nk2(l), nf2(l), l = 1, nblocks)
         end if
      end if

      if (iostat /= 0) then
         write (luncrt, '(/, a, i6)') &
            ' Error reading block sizes for file 2.  iostat: ', iostat
         go to 999
      end if

      if (nf2(1) /= numf) then
         write (luncrt, '(/, a, 2i5)') &
            ' Mismatched function counts:', numf, nf2(1)
         go to 999
      end if

      do l = 1, nblocks

         if (ni2(l) /= ni(l) .or. nj2(l) /= nj(l) .or. nk2(l) /= nk(l)) then
            write (luncrt, '(a, i3, a, (/, 3i8))') &
               ' Mismatched dimensions for block # ', l, ':', &
               ni(l), nj(l), nk(l), ni2(l), nj2(l), nk2(l) 
            go to 999
         end if

         i1 = ip(l)
         np = ip(l+1) - i1

         call rd_blk (lun2, formatted, npoints, numf, i1, np, f2, iostat)

         if (iostat /= 0) theN
            write (luncrt, '(/, a, i5, a, i5)') &
               ' Trouble reading file 2.  Block #:', l, '.  iostat:', iostat
            go to 999
         end if

      end do

      close (lun2)

      open (unit=lunout, file='compare.out', status='unknown', form='formatted')

!     Loop over blocks and functions to compare:

  100 continue

         write (luncrt, '(a)', advance='no') &
            ' Function number to be compared: '

         read  (lunkbd, *, end=999) m
         if (m < 1 .or. m > numf) go to 100

  101    l1 = 1;  l2 = 1
         if (nblocks > 1) then
            write (luncrt, '(a)', advance='no') &
               ' First and last blocks to be compared: '

            read  (lunkbd, *, end=999) l1, l2
            if (l1 < 1 .or. l1 > nblocks) go to 101
            if (l2 < 1 .or. l2 > nblocks) go to 101
         end if

         if (l1 == l2) then

            l = l1
            if (twod) then
               write (luncrt, '(a, 2i5)') ' Dimensions: ', ni(l), nj(l)
            else
               write (luncrt, '(a, 3i5)') ' Dimensions: ', ni(l), nj(l), nk(l)
            end if

  102       continue
            if (twod) then
               write (luncrt, '(a)', advance='no') &
                  ' (i,j) range to compare (2 pairs): '
               read  (lunkbd, *, end=999) i1, i2, j1, j2
               k1 = 1;  k2 = 1
            else
               write (luncrt, '(a)', advance='no') &
                  ' (i,j,k) range to compare (3 pairs): '
               read  (lunkbd, *, end=999) i1, i2, j1, j2, k1, k2
            end if

            if (i1 < 1 .or. i1 > ni(l)) go to 102
            if (i2 < 1 .or. i2 > ni(l)) go to 102
            if (j1 < 1 .or. j1 > nj(l)) go to 102
            if (j2 < 1 .or. j2 > nj(l)) go to 102
            if (k1 < 1 .or. k1 > nk(l)) go to 102
            if (k2 < 1 .or. k2 > nk(l)) go to 102

            ncompare = (i2 - i1 + 1) * (j2 - j1 + 1) * (k2 - k1 + 1)
            print = ncompare <= 1000

            if (ncompare == 1) then
               n = ip(l) + i1 - 1 + (j1 - 1) * ni(l) + (k1 - 1) * ni(l) * nj(l)
               write (luncrt, '(a, 1p, e25.15)') &
                  ' f1: ', f1(n,m), ' f2: ', f2(n,m)
            end if

         else
            print = .false.
         end if

         if (print) write (lunout, '(//, a, /)') &
            '  blk    i    j    k    f         df'

         dmax = -1.
         drms = zero
         ncompare = 0

         do l = l1, l2

            if (l1 /= l2) then
               i1 = 1
               i2 = ni(l)
               j1 = 1
               j2 = nj(l)
               k1 = 1
               k2 = nk(l)
            end if

            do k = k1, k2
               do j = j1, j2
                  n = ip(l) + i1 - 1 + (j - 1) * ni(l) + (k - 1) * ni(l) * nj(l)
                  do i = i1, i2
                     d = f2(n,m) - f1(n,m)
                     n = n + 1

                     if (print) write (lunout, '(5i5, 1p, e25.15)') &
                        l, i, j, k, m, d

                     d = abs (d)
                     drms = drms + d*d
                     ncompare = ncompare + 1

                     if (d > dmax) then
                        dmax = d
                        imax = i
                        jmax = j
                        kmax = k
                        lmax = l
                     end if
                  end do
               end do
            end do
         end do

         drms = sqrt (drms / real (ncompare))

         if (twod) then
            write (luncrt, '(a,1p,e25.15,a,2i5,a,i5,a,i4)') &
               ' Max. difference:', dmax, '  (i,j): ', imax, jmax, &
               '  block:', lmax, '  fn.:', m, ' RMS  difference:', drms
         else
            write (luncrt, '(a,1p,e25.15,a,3i5,a,i5,a,i4)') &
               ' Max. difference:', dmax, '  (i,j,k): ', imax, jmax, kmax, &
               '  block:', lmax, '  fn.:', m, ' RMS  difference:', drms
         end if
         n = ip(lmax) + imax - 1 + (jmax - 1) * ni(lmax) + &
                                   (kmax - 1) * ni(lmax) * nj(lmax)
         write (luncrt, '(a, 1p, e25.15)') &
            ' Corresp. f1:    ', f1(n,m),  &
            '          f2:    ', f2(n,m)
      go to 100

  999 continue

      contains ! Internal procedure for COMPARE_FLOWS

!       ------------------------------------------------------------------
        subroutine rd_blk (lun, formatted, npoints, nf, ip, np, f, iostat)
!       ------------------------------------------------------------------

!       RD_BLK reads all functions for one block of a PLOT3D function file,
!       formatted or not.  Trying to avoid subscripting while packing all
!       functions of all blocks proves impossible, because at the higher
!       level, the only way to get a rectangular 2-D array is to dimension
!       it f(npoints,nf) where npoints is the total number of points in all
!       blocks.  This does not have the functions of a given block contiguous
!       the way they are in the input file.
!       The logical unit is assumed to be positioned at the start of the block.

!       Arguments:

        integer, intent (in)  :: lun            ! Logical unit to read from
        logical, intent (in)  :: formatted      ! T or F
        integer, intent (in)  :: npoints        ! # points in this block
        integer, intent (in)  :: nf             ! # functions in the file
        integer, intent (in)  :: ip             ! f(ip:ip+np-1,1:nf) are ..
        integer, intent (in)  :: np             ! .. filled (np = # in the blk.)
        real,    intent (out) :: f(npoints,nf)  ! Packed block function values
        integer, intent (out) :: iostat         ! 0 means no read error

!       Execution:

        if (formatted) then
           read (lun, *, iostat=iostat) f(ip:ip+np-1,:)
        else
           read (lun,    iostat=iostat) f(ip:ip+np-1,:)
        end if

        end subroutine rd_blk

      end program compare_flows
