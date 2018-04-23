c-----------------------------------------------------------------------
c
      program JOIN_GRID
c
c  Description:
c
c     Utility for joining or rearranging multiblock grid blocks.
c
c     This version assumes block faces match neighboring faces exactly.
c     (Use the SPLIT_GRID program first if this is not true.)
c
c     Blocks are joined in one, two, or three directions at once,
c     first in the i direction, then j, then k (although the ijk
c     convention within each input block can be permuted, and each
c     direction can be reversed, via +/-1, 2, 3 input controls).
c
c     The first block entered for a given joined block corresponds
c     to the (1,1,1) element of a 3-D array of blocks being joined.
c     The relative positions of further blocks follow the Fortran
c     convention of varying the left-most index first, etc.
c
c     Multiple joined blocks may be produced in one run.
c
c     Reordering input blocks or simply transmitting them is achieved
c     by entering 1 1 1 for the i/j/kjoin counts.
c
c  Sample control file (join.inp):
c
c     Block joining control inputs:
c     Number of blocks after joining
c        2
c     -----------------------------------------------
c     Joined block ID
c        1
c     ijoin jjoin kjoin (# joins)
c        3     2     2
c     nblock(1:ijoin,1:jjoin,1:kjoin) & permuting indices
c        37          1  2  3  ! Base block
c        38          3 -2 -1  ! 1st block joined in i direction
c        30          1  2  3  ! Block (3,1,1)
c        40          1  2  3  ! (1,2,1)
c        41          1  2  3  ! (2,2,1)
c        50          1  2  3  ! (3,2,1)
c        51          1  2  3  ! (1,1,2)
c        60          1  2  3  ! (2,1,2)
c        :           :  :  :
c        99          1  2  3  ! (3,2,2)
c     -----------------------------------------------
c     Joined block ID
c        2
c     ijoin jjoin kjoin (# joins)
c        1     1     1
c     Input block
c       101          1  2  3  ! Transcribe this block
c     -----------------------------------------------
c
c  History:
c
c     ??/??/??  James Reuther   Original split_xyz.
c     11/12/99  David Saunders  Overhauled for use by mere mortals;
c                               no attempt to handle other than PLOT3D
c                               unformatted files; x(1:3,i,j,k) ordering
c                               in memory is inefficient.
c     11/16/99    "      "      Abandoned obscure original control scheme.
c                               Functionality (probably different from the
c                               original, however that worked) is spelled
c                               out above.  Avoid unnecessary rewinds.
c
c  Origin:  NASA Ames Research Center, Moffett Field, CA.
c
c-----------------------------------------------------------------------

      implicit none

c     Constants:

      integer, parameter ::
     .  lunctl = 1, lunin = 2, lunout = 3, lunkbd = 5, luncrt = 6

c     Variables:

      integer ::
     .  i, j, k, is, js, ks, isum, jsum, ksum, iperm, jperm, kperm,
     .  ioffset, joffset, koffset, ios, m0, ngin, ngout,
     .  l, m, n, lll, mmm, nnn, mxl, mxm, mxn, njoini, njoinj, njoink

      integer, target  ::
     .  ii, jj, kk

      integer, pointer ::
     .  ll, mm, nn

      integer, allocatable, dimension (:) ::
     .  ib, jb, kb, idimout, jdimout, kdimout, idm, jdm, kdm,
     .  ijoin, jjoin, kjoin, nbase, njoin

      integer, allocatable, dimension (:,:) ::
     .  imax, jmax, kmax, in, jn, kn

      integer, allocatable, dimension (:,:,:) ::
     .  mblockin, ip, jp, kp

      real, allocatable, dimension (:,:,:,:) ::
     .  xi, xo

      logical
     .  first

      character ::
     .  filename * 48

c     Execution:
c     ----------

      open (lunctl, file='join.inp', status='OLD', iostat=ios)

      if (ios /= 0) then
        write (luncrt, '(/, a)')
     .    ' Unable to open control file join.inp.'
        go to 999
      end if

      write (luncrt, '(/, a)', advance='NO')
     .  ' Input grid  (PLOT3D /mg/unf): '
      read  (lunkbd, '(a)') filename
      open (lunin, file=filename, form='unformatted', status='OLD')

      write (luncrt, '(/, a)', advance='NO')
     .  ' Output grid (PLOT3D /mg/unf): '
      read  (lunkbd, '(a)') filename
      open (lunout, file=filename, form='unformatted', status='NEW')

      read (lunin, iostat=ios) ngin

      if (ios /= 0) then
        write (luncrt, '(/, a)') ' Error reading # input blocks.'
        go to 999
      end if

      allocate (idm(ngin), jdm(ngin), kdm(ngin))

      read (lunin, iostat=ios) (idm(m), jdm(m), kdm(m), m = 1, ngin)

      if (ios /= 0) then
        write (luncrt, '(//, a)') ' Error reading input block sizes.'
        go to 999
      else
        write (luncrt, '(//, a)') ' Block    ni    nj    nk   (input)'
      end if

      write (luncrt, '(i4, 2x, 3i6)')
     .  (m, idm(m), jdm(m), kdm(m), m = 1, ngin)


c     Scan the control file to determine all output block dimensions:

      call SKIP (2) ! Internal procedure skips lines

      read (lunctl, *, iostat=ios) ngout

      if (ios /= 0) then
        write (luncrt, '(/, a)') ' Error reading # output blocks.'
        go to 999
      end if

      allocate (idimout(ngout), jdimout(ngout), kdimout(ngout))

      do n = 1, ngout ! For each group of blocks to be joined

        call SKIP (4)

        read (lunctl, *, iostat=ios) njoini, njoinj, njoink

        if (ios /= 0) then
          write (luncrt, '(/, a, i4)')
     .      ' Error reading # joins.  Output block ID:', n
          go to 999
        end if

        call SKIP (1)

c       Store the block #s & their permuting indices for this output block:

        allocate
     .   (mblockin(njoini,njoinj,njoink), ip(njoini,njoinj,njoink),
     .          jp(njoini,njoinj,njoink), kp(njoini,njoinj,njoink))

        do k = 1, njoink
          do j = 1, njoinj
            do i = 1, njoini

              read (lunctl, *, iostat=ios) m, iperm, jperm, kperm

              if (ios /= 0) then
                write (luncrt, '(/, a, i4, /, a, 3i4)')
     .            ' Error reading blocks to join.  Output block #:', n,
     .            ' Current (i,j,k) counts: ', i, j, k
                go to 999
              end if

              mblockin(i,j,k) = m
              ip(i,j,k) = iperm
              jp(i,j,k) = jperm
              kp(i,j,k) = kperm

            end do
          end do
        end do

c       Determine output block dimensions and check that other sums match:

        idimout(n) = 1
        first = .true.

        do k = 1, njoink
          do j = 1, njoinj
            isum = 1
            do i = 1, njoini
              m = mblockin(i,j,k)

              call PERMUTE (ip(i,j,k), jp(i,j,k), kp(i,j,k),
     .                      idm(m),    jdm(m),    kdm(m))

              if (first) then
                idimout(n) = idimout(n) + mxl - 1 ! Assume an overlap
              else
                isum = isum + mxl - 1
              end if

            end do

            if (first) then
              first = .false.
            else
              if (isum /= idimout(n)) then
                 write (luncrt, '(/, (a, i6))')
     .             ' Mismatched input i dimensions.  Output block:', n,
     .             ' Expected i dimension:', idimout(n),
     .             ' Differing dimension: ', isum,
     .             ' Input block j:       ', j,
     .             ' Input block k:       ', k
                 go to 999
              end if
            end if

          end do ! Next j
        end do ! Next k

        jdimout(n) = 1
        first = .true.

        do i = 1, njoini
          do k = 1, njoink
            jsum = 1
            do j = 1, njoinj
              m = mblockin(i,j,k)

              call PERMUTE (ip(i,j,k), jp(i,j,k), kp(i,j,k),
     .                      idm(m),    jdm(m),    kdm(m))

              if (first) then
                jdimout(n) = jdimout(n) + mxm - 1
              else
                jsum = jsum + mxm - 1
              end if

            end do

            if (first) then
              first = .false.
            else
              if (jsum /= jdimout(n)) then
                 write (luncrt, '(/, (a, i6))')
     .             ' Mismatched input j dimensions.  Output block:', n,
     .             ' Expected j dimension:', jdimout(n),
     .             ' Differing dimension: ', jsum,
     .             ' Input block i:       ', i,
     .             ' Input block k:       ', k
                 go to 999
              end if
            end if

          end do ! Next k
        end do ! Next i

        kdimout(n) = 1
        first = .true.

        do j = 1, njoinj
          do i = 1, njoini
            ksum = 1
            do k = 1, njoink
              m = mblockin(i,j,k)

              call PERMUTE (ip(i,j,k), jp(i,j,k), kp(i,j,k),
     .                      idm(m),    jdm(m),    kdm(m))

              if (first) then
                kdimout(n) = kdimout(n) + mxn - 1
              else
                ksum = ksum + mxn - 1
              end if

            end do

            if (first) then
              first = .false.
            else
              if (ksum /= kdimout(n)) then
                 write (luncrt, '(/, (a, i6))')
     .             ' Mismatched input k dimensions.  Output block:', n,
     .             ' Expected k dimension:', kdimout(n),
     .             ' Differing dimension: ', ksum,
     .             ' Input block i:       ', i,
     .             ' Input block j:       ', j
                 go to 999
              end if
            end if

          end do ! Next i
        end do ! Next j

        deallocate (mblockin, ip, jp, kp)

      end do ! Next output block controls

      write (luncrt, '(/, a)') ' Block    ni    nj    nk   (output)'

      write (luncrt, '(i4, 2x, 3i6)')
     .  (n, idimout(n), jdimout(n), kdimout(n), n = 1, ngout)

      write (luncrt, '(/, a)') ' Begin construction of joined blocks.'

      write (lunout) ngout
      write (lunout) (idimout(n), jdimout(n), kdimout(n), n = 1, ngout)

      rewind (lunctl)
      call SKIP (3)


c     Main loop over output blocks:
c     -----------------------------

      do n = 1, ngout

        write (luncrt, '(/, a, i4, a)', advance='NO') '   Block', n, ':'

        allocate (xo(idimout(n),jdimout(n),kdimout(n),3))

        call SKIP (4)

        read (lunctl, *) njoini, njoinj, njoink

        call SKIP (1)

        m0 = 999999 ! Previous input block read
        koffset = 0

        do k = 1, njoink

          joffset = 0

          do j = 1, njoinj

            ioffset = 0

            do i = 1, njoini

              read (lunctl, *) m, iperm, jperm, kperm

              write (luncrt, '(i5)', advance='NO') m

              if (m > m0) then ! No need to rewind

                do is = m0 + 1, m - 1
                  read (lunin)
                end do

              else

                rewind (lunin)

                do is = 1, m + 1
                  read (lunin) ! Skip nblocks, dimens, & blocks 1 : m - 1
                end do

              end if

              m0 = m
 
              allocate (xi(idm(m),jdm(m),kdm(m),3))

              read (lunin) xi

              call PERMUTE (iperm,  jperm,  kperm,
     .                      idm(m), jdm(m), kdm(m))

c             Swap and flip:            ! Tricky formula from FLO107-MB

              do kk = 1, kdm(m)
                do jj = 1, jdm(m)
                  do ii = 1, idm(m)
                    lll = ll * (is + 1)/2 - (mxl - (ll-1)) * (is - 1)/2
     .                    + ioffset
                    mmm = mm * (js + 1)/2 - (mxm - (mm-1)) * (js - 1)/2
     .                    + joffset
                    nnn = nn * (ks + 1)/2 - (mxn - (nn-1)) * (ks - 1)/2
     .                    + koffset
                    xo(lll,mmm,nnn,1) = xi(ii,jj,kk,1)
                    xo(lll,mmm,nnn,2) = xi(ii,jj,kk,2)
                    xo(lll,mmm,nnn,3) = xi(ii,jj,kk,3)
                  end do
                end do
              end do

              deallocate (xi)

              ioffset = ioffset + mxl - 1

            end do ! Next join in the i direction

            joffset = joffset + mxm - 1

          end do ! Next join in the j direction

          koffset = koffset + mxn - 1

        end do ! Next join in the k direction

        write (lunout) xo

        deallocate (xo)

      end do ! Next output block

      close (lunout)
      close (lunctl)


  999 continue ! Avoid system-dependent STOP behavior


!     Internal procedures for JOIN_GRID:

      contains

!       ----------------------------------------------------------------------
!
        subroutine PERMUTE (ip, jp, kp, idm, jdm, kdm)
!
!       Handle permuting/reversing i,j,k for an input block.
!       ----------------------------------------------------------------------

!       Arguments:

        integer, intent (in)  :: ip, jp, kp ! Permutation indices (+/- 1, 2, 3)
        integer, intent (in)  :: idm, jdm, kdm ! Current input block dimensions

!       Host variables used:

!       integer, intent (in)  :: ii, jj, kk ! Loop indices over input block
!       integer, intent (out) :: ll, mm, nn ! Pointers assigned to ii, jj, or kk
!       integer, intent (out) :: mxl, mxm, mxn ! Permutations of idm, jdm, kdm
!       integer, intent (out) :: is, js, ks ! Directn. increments (i/j/kp signs)

!       Local variables:

        integer ijk

c       Permute block i,j,k indices?

        ijk = abs (ip)
        if      (ijk == 1) then
          ll  => ii
          mxl =  idm
        else if (ijk == 2) then
          ll  => jj
          mxl =  jdm
        else if (ijk == 3) then
          ll  => kk
          mxl =  kdm
        end if

        ijk = abs (jp)
        if      (ijk == 1) then
          mm  => ii
          mxm =  idm
        else if (ijk == 2) then
          mm  => jj
          mxm =  jdm
        else if (ijk == 3) then
          mm  => kk
          mxm =  kdm
        end if

        ijk = abs (kp)
        if      (ijk == 1) then
          nn  => ii
          mxn =  idm
        else if (ijk == 2) then
          nn  => jj
          mxn =  jdm
        else if (ijk == 3) then
          nn  => kk
          mxn =  kdm
        end if

c       Reverse index directions?

        is = sign (1, ip)
        js = sign (1, jp)
        ks = sign (1, kp)

        end subroutine PERMUTE

!       ----------------------------------------------------------
        subroutine SKIP (nlines) ! Skip lines of the control file.
!       ----------------------------------------------------------

        integer l, nlines

        do l = 1, nlines
          read (lunctl, *, iostat=ios)
        end do

        end subroutine SKIP

      end program JOIN_GRID
