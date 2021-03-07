!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      program check_boundary
!
!     Description:
!
!        CHECK_BOUNDARY is intended to help catch hypersonic grids that do not
!     completely contain a bow shock - probably through overzealous tailoring
!     of the outer boundary, itself likely to be the result of excessive use of
!     smoothing.
!
!        PLOT3D grid and function files are expected.
!
!        Grids are assumed to be of the type handled by OUTBOUND: one layer of
!     grid blocks, with k = kmax at the outer boundary.  For each block, the
!     the number of radial grid lines with the fewest number of points in the
!     free stream is reported.  Any value less than 2 means trouble, and 3 is
!     recommended.
!
!        The function file should contain Mach number as function # 1.  (Any
!     further functions are ignored.  Obviously, this could be generalized.)
!     The free stream value is taken from block 1, point (1,1,nk).
!
!        Use standard extensions (.g or .gu, .f or .fu) to save some prompting.
!
!     History:
!
!        12/20/06  D. Saunders  Initial implementation.
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
         luncrt = 6,        &
         lunkbd = 5,        &
         lungrd = 1,        &
         lunfun = 2

      character, parameter :: &
         format * 11 = 'unformatted'

!     Variables:

      integer :: &
         i, i1, ib, ios, iworst, j, jworst, k, kbelow, kworst, kworstp1, l,    &
         nblocks, ni, nj, nk, nfree, nfun, nkfree, npts, nworst

      real :: &
         Minf, Mtest

      logical :: &
         formatted

      character :: &
         filename * 64, flag * 5

!     Derived data types:

      type (grid_type), pointer, dimension (:) :: &
         grid

!     Execution:

      write (luncrt, '(/, a)', advance='no') ' Input grid [*.g | *.gu]:    '
      read  (lunkbd, *) filename
      l = len_trim (filename)
      formatted = filename(l:l) /= 'u'
      i1 = 1;  if (formatted) i1 = 3

      open (lungrd, file=filename(1:l), form=format(i1:11), status='old',      &
            iostat=ios)

      if (ios /= 0) then
         write (luncrt, '(/, 2a)') &
            ' Unable to open input grid: ', filename(1:l)
         go to 999
      end if

      call xyz_header_io (1, lungrd, formatted, nblocks, grid, ios)
      if (ios /= 0) go to 999

      write (luncrt, '(a)', advance='no') ' Function file [*.f | *.fu]: '
      read  (lunkbd, *) filename
      l = len_trim (filename)
      formatted = filename(l:l) /= 'u'
      i1 = 1;  if (formatted) i1 = 3

      open (lunfun, file=filename(1:l), form=format(i1:11), status='old',      &
            iostat=ios)
      if (ios /= 0) then
         write (luncrt, '(/, 2a)') &
            ' Unable to open input function file: ', filename(1:l)
         go to 999
      end if

      call q_header_io (1, lunfun, formatted, nblocks, nfun, grid, ios)
      if (ios /= 0) go to 999

!     Process one block at a time:

      write (luncrt, '(/, 2a)') &
         '  ib nMinf  i    j              x              y              z',    &
         '         M(kbelow+1)           M(kbelow)         M(kbelow-1)'

      do ib = 1, nblocks

         call xyz_allocate (grid(ib), ios)

         ni = grid(ib)%ni
         nj = grid(ib)%nj
         nk = grid(ib)%nk

         npts = ni * nj * nk

         call xyz_block_io (1, lungrd, formatted, npts,                        &
                            grid(ib)%x, grid(ib)%y, grid(ib)%z, ios)
         if (ios /= 0) go to 999

         call q_allocate (grid(ib), nfun, ios)

         call q_block_io (1, lunfun, formatted, nfun, ni, nj, nk, grid(ib)%q,  &
                          ios)
         if (ios /= 0) go to 999

         if (ib == 1) then
            Minf  = grid(ib)%q(1,1,1,nk)
            Mtest = Minf - 0.000001
         end if

!        Find the radial line in this block with fewest points in the free strm:

         kworst = -1

         do j = 1, nj

            do i = 1, ni

               do k = nk, 2, -1
                  if (grid(ib)%q(1,i,j,k) < Mtest) then
                     kbelow = k
                     exit
                  end if
               end do

               if (kbelow > kworst) then
                  kworst = kbelow
                  iworst = i
                  jworst = j
               end if

            end do

         end do

         nkfree   = nk - kworst
         kworstp1 = min (kworst+1, nk)
         flag = '     ';  if (nkfree < 3) flag = '  ***'

         write (luncrt, '(2i4, 2i5, 3f15.6, 3f20.12, a5)') &
            ib, nkfree, iworst, jworst,                    &
            grid(ib)%x(iworst,jworst,kworst), grid(ib)%y(iworst,jworst,kworst),&
            grid(ib)%z(iworst,jworst,kworst),                                  &
            grid(ib)%q(1,iworst,jworst,kworstp1),                              &
            grid(ib)%q(1,iworst,jworst,kworst),                                &
            grid(ib)%q(1,iworst,jworst,kworst-1), flag

         deallocate (grid(ib)%x, grid(ib)%y, grid(ib)%z, grid(ib)%q)

      end do ! Next block

      deallocate (grid)

      close (lungrd)
      close (lunfun)

  999 continue

! *** stop ! Avoid system dependencies.

      end program check_boundary
