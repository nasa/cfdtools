!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      program MERGE_BLOCKS

!  Description:
!
!     Utility for joining or rearranging multiblock grid blocks.
!     This adaptation of program join_grid allows for an optional function file.
!
!     The input block faces are assumed to match neighboring faces exactly.
!     (Use the SPLIT_GRID program first if this is not true.)
!
!     Blocks are joined in one, two, or three directions at once, first in the
!     i direction, then j, then k (although the ijk convention within each input
!     block can be permuted, and each direction can be reversed, via +/-1, 2, 3
!     input controls).
!
!     The first block entered for a given joined block corresponds to the
!     (1,1,1) element of a 3-D array of blocks being joined.  The relative
!     positions of further blocks follow the Fortran convention of varying the
!     left-most index first, etc.
!
!     Multiple joined blocks may be produced in one run.
!
!     Reordering input blocks or simply transmitting them is achieved by
!     entering 1 1 1 for the i/j/kjoin counts.
!
!  Sample control file (merge.inp):
!
!     Block merging control inputs:
!     Number of blocks after merging
!        2
!     -----------------------------------------------
!     Joined block ID
!        1
!     ijoin jjoin kjoin (# joins)
!        3     2     2
!     nblock(1:ijoin,1:jjoin,1:kjoin) & permuting indices
!        37          1  2  3  ! Base block
!        38          3 -2 -1  ! 1st block joined in i direction
!        30          1  2  3  ! Block (3,1,1)
!        40          1  2  3  ! (1,2,1)
!        41          1  2  3  ! (2,2,1)
!        50          1  2  3  ! (3,2,1)
!        51          1  2  3  ! (1,1,2)
!        60          1  2  3  ! (2,1,2)
!        :           :  :  :
!        99          1  2  3  ! (3,2,2)
!     -----------------------------------------------
!     Joined block ID
!        2
!     ijoin jjoin kjoin (# joins)
!        1     1     1
!     Input block
!       101          1  2  3  ! Transcribe this block
!     -----------------------------------------------
!
!  History:
!
!     ??/??/??  James Reuther   Original SPLIT_XYZ, renamed JOIN_GRID here.
!     11/12/99  David Saunders  Overhauled for use by mere mortals;
!                               no attempt to handle other than PLOT3D
!                               unformatted files; x(1:3,i,j,k) ordering
!                               in memory is inefficient.
!     11/16/99    "      "      Abandoned obscure original control scheme.
!                               Functionality (probably different from the
!                               original, however that worked) is spelled
!                               out above.  Avoid unnecessary rewinds.
!     10/24/06    "      "      MERGE_BLOCKS adapted from JOIN_GRID in order
!                               to deal with optional function files too.
!
!  Origin:  NASA Ames Research Center, Moffett Field, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Constants:

      integer, parameter :: &
        lunctl = 1, lunfi = 2, lunfo = 3, lungi = 4, lungo = 7, lunkbd = 5,    &
        luncrt = 6

!     Variables:

      integer :: &
        i, j, k, is, js, ks, isum, jsum, ksum, iperm, jperm, kperm,            &
        ioffset, joffset, koffset, ios, m0, nbin, nbout, nf,                   &
        l, m, n, lll, mmm, nnn, mxl, mxm, mxn, njoini, njoinj, njoink

      integer, target  :: &
        ii, jj, kk

      integer, pointer :: &
        ll, mm, nn

      integer, allocatable, dimension (:) :: &
        ib, jb, kb, idimout, jdimout, kdimout, idm, jdm, kdm,                  &
        ijoin, jjoin, kjoin, nbase, njoin

      integer, allocatable, dimension (:,:) :: &
        imax, jmax, kmax, in, jn, kn

      integer, allocatable, dimension (:,:,:) :: &
        mblockin, ip, jp, kp

      real, allocatable, dimension (:,:,:,:) :: &
        fi, fo, xi, xo

      logical :: &
        first, function

      character :: &
        filename * 64

!     Execution:
!     ----------

      open (lunctl, file='merge.inp', status='OLD', iostat=ios)
      if (ios /= 0) then
        write (luncrt, '(/, a)') ' Unable to open control file merge.inp.'
        go to 999
      end if

      write (luncrt, '(a)', advance='NO') ' Input grid  (PLOT3D /mg/unf): '
      read  (lunkbd, '(a)') filename
      open  (lungi, file=filename, form='unformatted', status='OLD')
      if (ios /= 0) then
        write (luncrt, '(/, a)') ' Unable to open input grid file.'
        go to 999
      end if

      write (luncrt, '(a)', advance='NO') ' Associated function file | none): '
      read  (lunkbd, '(a)') filename
      function = filename(1:4) /= 'none'

      if (function) then
        open  (lunfi, file=filename, form='unformatted', status='OLD')
        if (ios /= 0) then
          write (luncrt, '(/, a)') ' Unable to open input function file.'
          go to 999
        end if
      end if

      write (luncrt, '(a)', advance='NO') ' Output grid (PLOT3D /mg/unf): '
      read  (lunkbd, '(a)') filename
      open  (lungo, file=filename, form='unformatted', status='UNKNOWN')

      if (function) then
        write (luncrt, '(a)', advance='NO') ' Output function file: '
        read  (lunkbd, '(a)') filename
        open  (lunfo, file=filename, form='unformatted', status='UNKNOWN')
      end if

!     Read the number of input blocks:

      read (lungi, iostat=ios) nbin

      if (ios /= 0) then
        write (luncrt, '(/, a)') ' Error reading # input grid blocks.'
        go to 999
      end if

      if (function) then
         read (lunfi, iostat=ios) m
         if (ios /= 0) then
           write (luncrt, '(/, a)') ' Error reading # input function blocks.'
           go to 999
        end if
        if (m /= nbin) then
          write (luncrt, '(/, a, 2i5)') ' Mismatched # blocks: ', nbin, m
          go to 999
        end if
      end if

      allocate (idm(nbin), jdm(nbin), kdm(nbin))

      read (lungi, iostat=ios) (idm(m), jdm(m), kdm(m), m = 1, nbin)

      if (ios /= 0) then
        write (luncrt, '(//, a)') ' Error reading input grid block sizes.'
        go to 999
      end if

      write (luncrt, '(/, a)') ' Block    ni    nj    nk   (input)'
      write (luncrt, '(i4, 2x, 3i6)') (m, idm(m), jdm(m), kdm(m), m = 1, nbin)

      if (function) then
         read (lunfi, iostat=ios) (idm(m), jdm(m), kdm(m), nf, m = 1, nbin)
         if (ios /= 0) then
           write (luncrt, '(/, a)') ' Error reading input function block sizes.'
           go to 999
        end if
      end if

!     Scan the control file to determine all output block dimensions:

      call SKIP (2) ! Internal procedure skips lines

      read (lunctl, *, iostat=ios) nbout

      if (ios /= 0) then
        write (luncrt, '(/, a)') ' Error reading # output blocks.'
        go to 999
      end if

      allocate (idimout(nbout), jdimout(nbout), kdimout(nbout))

!     For each group of blocks to be joined:

      do n = 1, nbout

        call SKIP (4)

        read (lunctl, *, iostat=ios) njoini, njoinj, njoink

        if (ios /= 0) then
          write (luncrt, '(/, a, i4)') &
            ' Error reading # joins.  Output block ID:', n
          go to 999
        end if

        call SKIP (1)

!       Store the block #s & their permuting indices for this output block:

        allocate (mblockin(njoini,njoinj,njoink), ip(njoini,njoinj,njoink),    &
                  jp(njoini,njoinj,njoink), kp(njoini,njoinj,njoink))

        do k = 1, njoink
          do j = 1, njoinj
            do i = 1, njoini

              read (lunctl, *, iostat=ios) m, iperm, jperm, kperm

              if (ios /= 0) then
                write (luncrt, '(/, a, i4, /, a, 3i4)') &
                  ' Error reading blocks to join.  Output block #:', n, &
                  ' Current (i,j,k) counts: ', i, j, k
                go to 999
              end if

              mblockin(i,j,k) = m
              ip(i,j,k) = iperm
              jp(i,j,k) = jperm
              kp(i,j,k) = kperm

            end do
          end do
        end do

!       Determine output block dimensions and check that other sums match:

        idimout(n) = 1
        first = .true.

        do k = 1, njoink
          do j = 1, njoinj
            isum = 1
            do i = 1, njoini
              m = mblockin(i,j,k)

              call PERMUTE (ip(i,j,k), jp(i,j,k), kp(i,j,k), &
                            idm(m),    jdm(m),    kdm(m))

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
                 write (luncrt, '(/, (a, i6))') &
                   ' Mismatched input i dimensions.  Output block:', n, &
                   ' Expected i dimension:', idimout(n),                &
                   ' Differing dimension: ', isum,                      &
                   ' Input block j:       ', j,                         &
                   ' Input block k:       ', k
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

              call PERMUTE (ip(i,j,k), jp(i,j,k), kp(i,j,k), &
                            idm(m),    jdm(m),    kdm(m))

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
                 write (luncrt, '(/, (a, i6))') &
                   ' Mismatched input j dimensions.  Output block:', n, &
                   ' Expected j dimension:', jdimout(n),                &
                   ' Differing dimension: ', jsum,                      &
                   ' Input block i:       ', i,                         &
                   ' Input block k:       ', k
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

              call PERMUTE (ip(i,j,k), jp(i,j,k), kp(i,j,k), &
                            idm(m),    jdm(m),    kdm(m))

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
                 write (luncrt, '(/, (a, i6))') &
                   ' Mismatched input k dimensions.  Output block:', n, &
                   ' Expected k dimension:', kdimout(n),                &
                   ' Differing dimension: ', ksum,                      &
                   ' Input block i:       ', i,                         &
                   ' Input block j:       ', j
                 go to 999
              end if
            end if

          end do ! Next i
        end do ! Next j

        deallocate (mblockin, ip, jp, kp)

      end do ! Next output block controls

      write (luncrt, '(/, a)') ' Block    ni    nj    nk   (output)'

      write (luncrt, '(i4, 2x, 3i6)') &
        (n, idimout(n), jdimout(n), kdimout(n), n = 1, nbout)

      write (luncrt, '(/, a)') ' Begin construction of joined blocks.'

      write (lungo) nbout
      write (lungo) (idimout(n), jdimout(n), kdimout(n), n = 1, nbout)

      if (function) then
        write (lunfo) nbout
        write (lunfo) (idimout(n), jdimout(n), kdimout(n), nf, n = 1, nbout)
      end if

      rewind (lunctl)
      call SKIP (3)


!     Main loop over output blocks:
!     -----------------------------

      do n = 1, nbout

        write (luncrt, '(/, a, i4, a)', advance='NO') '   Block', n, ':'

        allocate (xo(idimout(n),jdimout(n),kdimout(n),3))

        if (function) allocate (fo(idimout(n),jdimout(n),kdimout(n),nf))

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
                  read (lungi)
                end do

                if (function) then
                  do is = m0 + 1, m - 1
                    read (lunfi)
                  end do
                end if

              else

                rewind (lungi)

                do is = 1, m + 1
                  read (lungi) ! Skip nblocks, dimens, & blocks 1 : m - 1
                end do

                if (function) then
                  rewind (lunfi)
                  do is = 1, m + 1
                    read (lunfi) ! Skip nblocks, dimens, & blocks 1 : m - 1
                  end do
                end if

              end if

              m0 = m
 
              allocate (xi(idm(m),jdm(m),kdm(m),3))

              read (lungi) xi

              call PERMUTE (iperm,  jperm,  kperm, idm(m), jdm(m), kdm(m))

!             Swap and flip:            ! Tricky formula from FLO107-MB

              do kk = 1, kdm(m)
                do jj = 1, jdm(m)
                  do ii = 1, idm(m)
                    lll = ll * (is + 1)/2 - (mxl - (ll-1)) * (is - 1)/2 +      &
                          ioffset
                    mmm = mm * (js + 1)/2 - (mxm - (mm-1)) * (js - 1)/2 +      &
                          joffset
                    nnn = nn * (ks + 1)/2 - (mxn - (nn-1)) * (ks - 1)/2 +      &
                          koffset
                    xo(lll,mmm,nnn,1) = xi(ii,jj,kk,1)
                    xo(lll,mmm,nnn,2) = xi(ii,jj,kk,2)
                    xo(lll,mmm,nnn,3) = xi(ii,jj,kk,3)
                  end do
                end do
              end do

              deallocate (xi)

              if (function) then

                allocate (fi(idm(m),jdm(m),kdm(m),nf))

                read (lunfi) fi

                do kk = 1, kdm(m)
                  do jj = 1, jdm(m)
                    do ii = 1, idm(m)
                      lll = ll * (is + 1)/2 - (mxl - (ll-1)) * (is - 1)/2 +    &
                            ioffset
                      mmm = mm * (js + 1)/2 - (mxm - (mm-1)) * (js - 1)/2 +    &
                            joffset
                      nnn = nn * (ks + 1)/2 - (mxn - (nn-1)) * (ks - 1)/2 +    &
                            koffset
                      fo(lll,mmm,nnn,:) = fi(ii,jj,kk,:)
                    end do
                  end do
                end do

                deallocate (fi)

              end if

              ioffset = ioffset + mxl - 1

            end do ! Next join in the i direction

            joffset = joffset + mxm - 1

          end do ! Next join in the j direction

          koffset = koffset + mxn - 1

        end do ! Next join in the k direction

        write (lungo) xo
        deallocate   (xo)

        if (function) then
          write (lunfo) fo
          deallocate   (fo)
        end if

      end do ! Next output block

      write (luncrt, '(/, a)')
      close (lunctl)
      close (lungo)
      if (function) close (lunfo)


  999 continue ! Avoid system-dependent STOP behavior


!     Internal procedures for MERGE_BLOCKS:

      contains

!       ----------------------------------------------------------------------
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

!       Permute block i,j,k indices?

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

!       Reverse index directions?

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

      end program MERGE_BLOCKS
