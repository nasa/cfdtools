!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program iblank_move
!
!  Purpose:
!
!     The overset option of flow solver DPLR allows extraction of blanking
!     information as one function in a PLOT3D-type function file.  Utilities
!     associated with the OVERFLOW solver expect the blanking information to be
!     contained in the grid file.  This utility simply reads output from DPLR's
!     POSTFLOW utility that is in a PLOT3D function file with the blanking
!     variable as one of the functions (the number is prompted for), and
!     transfers it as the integer equivalent to an output iblanked form of the
!     input grid file.  The function data proper are output in PLOT3D's Q file
!     form, not the function file form originally implemented:
!
!        ngrid
!        (ni(n),nj(n),nk(n),n=1,ngrid),nq,nspecies
!        for n=1:ngrid
!           refmach,alpha,rey,time,gaminf,beta,tinf,igam,htinf,ht1,ht2,
!           rgas(1:nspecies),fsmach,tvref,dtvref
!           ((((q(i,j,k,l),i=1,ni(n)),j=1:nj(n)),k=1,nk(n)),l=1,nq)
!
!        where nq = 6 + nspecies + 0|1|2 (# turbulent terms)
!
!     See online documentation for more on the variable descriptions. Here,
!     the extras are written as place holders, because they are not required
!     for the expected usage by the PLOT3D-related overset utilities to provide
!     input to the BLAYER utility.
!
!     The formatting of the input grid is determined automatically, and the
!     input function file is assumed to be formatted the same ASCII|unformatted
!     way.  Output files preserve the formatting type.
!
!     A few prompts suffice for file names, the iblank function number, the
!     number of species, and the number of turbulence field quantities.  The
!     grid blocks are processed one at a time for efficiency reasons.
!
!  History:
!     07/19/2019  DAS  Initial implementation at the request of Chun Tang.
!     07/22/2019   "   Perplexing failure on unformatted input files traced to
!                      forgetting about the form= keyword.
!     07/23/2019   "   Suppress embedded .g[u]-type extensions in the output
!                      file names.
!     08/14/2019   "   The output was supposed to include a Q file, not a
!                      function file. Make the file name *.q[u].
!  Author:
!     David Saunders, AMA, Inc. at NASA Ames Research Center.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure  ! Module for a data type representing one block
   use xyzq_io_module        ! PLOT3D-type I/O pacakge, including iblank options

   implicit none

!  Constants:

   integer, parameter :: lungin  = 1   ! Input grid file
   integer, parameter :: lunfin  = 2   ! Input function file with iblank info
   integer, parameter :: lungout = 3   ! Output iblanked grid file
   integer, parameter :: lunqout = 4   ! Output Q file (no iblanks)
   integer, parameter :: lunkbd  = 5   ! Keyboard inputs
   integer, parameter :: luncrt  = 6   ! Screen prompts & diagnostics

   integer, parameter :: iz      = 0   ! Place holder only
   real,    parameter :: rz      = 0.  !   "     "     "

   character (11), parameter :: form = 'unformatted'

!  Variables:

   integer :: ib, ibnum, i, ios, j, k, lf, lg, mi, mj, mk, n, nblocks, &
              nfin, nfout, npts, nq, nspecies, nturb
   logical :: formatted
   character (132) :: filenamef, filenameg
   type (grid_type), pointer, dimension (:) :: xyzf_in, xyzf_out

!  Execution:

!  Open the grid file and determine its formatting:

   write (luncrt, '(/, a)', advance='no') 'Input grid file name:     '
   read  (lunkbd, '(a)') filenameg
   lg = len_trim (filenameg)
   call determine_grid_form (filenameg(1:lg), lungin, formatted, ios)
   if (ios /= 0) then
      write (luncrt, '(2a)') '*** Trouble opening file ', filenameg(1:lg), &
         '    Bad file name?'
      go to 99
   end if
   i = 1;  if (formatted) i = 3  ! index into the form string

   open (lungin, file=filenameg(1:lg), status='old', form=form(i:), iostat=ios)
   if (ios /= 0) go to 99  ! Shouldn't be possible

!  Repeat for the input function file; assume the formatting is the same:

   write (luncrt, '(a)', advance='no') 'Input function file name: '
   read  (lunkbd, '(a)') filenamef
   lf = len_trim (filenamef)
   open (lunfin, file=filenamef(1:lf), status='old', form=form(i:), iostat=ios)
   if (ios /= 0) go to 99

   write (luncrt, '(a)', advance='no') 'iblank function number: '
   read  (lunkbd, *) ibnum

!! write (luncrt, '(a)', advance='no') '# species: '
!! read  (lunkbd, *) nspecies
   nspecies = 0  ! We just need to output nfin - 1 functions

!! write (luncrt, '(a)', advance='no') '# turbulent terms [0|1|2]: '
!! read  (lunkbd, *) nturb
   nturb = 0  ! See nspecies

!! nq = 6 + nspecies + nturb  ! Use nfout = nfin - 1

!  Read the input grid and function file headers.  The grid header read
!  allocates the array of blocks, xyzf_in(:):

   call xyz_header_io (1, lungin, formatted, nblocks, xyzf_in, ios)
   write (luncrt, '(a, i5)') 'Number of blocks found:   ', nblocks
   if (ios /= 0) go to 99

   call q_header_io (1, lunfin, formatted, nblocks, nfin, xyzf_in, ios)
   if (ios /= 0) go to 99

   write (luncrt, '(a, i5)') 'Number of functions found:', nfin
   if (ibnum < 1 .or.ibnum > nfin) then
      write (luncrt, '(a, i5)') '*** Bad iblank function number:', ibnum
      go to 99
   end if

!  Display the block dimensions to reassure the user:

   write (luncrt, '(a)')
   write (luncrt, '(4(i6, 2x, 3i5))') &
      (ib, xyzf_in(ib)%ni, xyzf_in(ib)%nj, xyzf_in(ib)%nk, ib = 1, nblocks)

!  Set up the output files:

   lg = index (filenameg(1:lg), '.', .true.)  ! Backward search for .g[u]
   filenameg(lg:) = '.iblank.gu'
   lf = index (filenamef(1:lf), '.', .true.)
   filenamef(lf:) = '.iblank.qu'
   lg = lg + 9;  lf = lf + 9
   if (formatted) then
      lg = lg - 1
      lf = lf - 1
   end if

   open (lungout, file=filenameg(1:lg), status='unknown', form=form(i:), &
         iostat=ios)
   if (ios /= 0) then
      write (luncrt, '(2a)') '*** Unable to open output grid file: ', &
         filenameg(1:lg)
      go to 99
   end if

   open (lunqout, file=filenamef(1:lf), status='unknown', form=form(i:), &
         iostat=ios)
   if (ios /= 0) then
      write (luncrt, '(2a)') '*** Unable to open output function file: ', &
         filenamef(1:lf)
      go to 99
   end if

   allocate (xyzf_out(nblocks))
   do ib = 1, nblocks
      xyzf_out(ib)%ni = xyzf_in(ib)%ni
      xyzf_out(ib)%nj = xyzf_in(ib)%nj
      xyzf_out(ib)%nk = xyzf_in(ib)%nk
      xyzf_out(ib)%mi = xyzf_in(ib)%mi
      xyzf_out(ib)%mj = xyzf_in(ib)%mj
      xyzf_out(ib)%mk = xyzf_in(ib)%mk
   end do

   call xyz_header_io (2, lungout, formatted, nblocks, xyzf_out, ios)
   if (ios /= 0) go to 99

   nfout = nfin - 1  ! = nq if nspecies and nturb are both set to 0

!!!call q_header_io (2, lunqout, formatted, nblocks, nfout, xyzf_out, ios)
!!!if (ios /= 0) go to 99

   if (formatted) then
      write (lunqout, '(i6)') nblocks
      write (lunqout, '(3i6)') &
         (xyzf_out(ib)%ni,xyzf_out(ib)%nj,xyzf_out(ib)%nk, ib = 1, nblocks), &
         nfout, nspecies
   else
      write (lunqout) nblocks
      write (lunqout) &
         (xyzf_out(ib)%ni,xyzf_out(ib)%nj,xyzf_out(ib)%nk, ib = 1, nblocks), &
         nfout, nspecies
   end if

!  Process one block at a time:

   do ib = 1, nblocks

      call xyz_allocate (xyzf_in(ib), ios)
      if (ios /= 0) then
         write (luncrt, '(a, i5)') '*** Trouble allocating input grid block:', &
            ib
         go to 99
      end if

      npts = xyzf_in(ib)%ni*xyzf_in(ib)%nj*xyzf_in(ib)%nk

      call xyz_block_io (1, lungin, formatted, npts, xyzf_in(ib)%x, &
                         xyzf_in(ib)%y, xyzf_in(ib)%z, ios)
      if (ios /= 0) then
         write (luncrt, '(a, i5)') '*** Trouble reading grid block:', ib
         go to 99
      end if

      call q_allocate (xyzf_in(ib), nfin, ios)
      if (ios /= 0) then
         write (luncrt, '(a, i5)') '*** Trouble allocating input fn. block:', &
            ib
         go to 99
      end if

      mi = xyzf_in(ib)%mi;  mj = xyzf_in(ib)%mj;  mk = xyzf_in(ib)%mk
      call q_block_io (1, lunfin, formatted, nfin, mi, mj, mk, xyzf_in(ib)%q, &
                       ios)
      if (ios /= 0) then
         write (luncrt, '(a, i5)') '*** Trouble reading function block:', ib
         go to 99
      end if

      call xyzi_allocate (xyzf_out(ib), ios)

      xyzf_out(ib)%x = xyzf_in(ib)%x
      xyzf_out(ib)%y = xyzf_in(ib)%y
      xyzf_out(ib)%z = xyzf_in(ib)%z
      xyzf_out(ib)%iblank(:,:,:) = xyzf_in(ib)%q(ibnum,:,:,:)

      deallocate (xyzf_in(ib)%x, xyzf_in(ib)%y, xyzf_in(ib)%z)

      call xyzi_block_io (2, lungout, formatted, npts, xyzf_out(ib)%x, &
                       xyzf_out(ib)%y, xyzf_out(ib)%z, xyzf_out(ib)%iblank, ios)
      if (ios /= 0) then
         write (luncrt, '(a, i5)') '*** Trouble writing grid block:', ib
         go to 99
      end if

      deallocate (xyzf_out(ib)%x, xyzf_out(ib)%y, xyzf_out(ib)%z, &
                  xyzf_out(ib)%iblank)

      call q_allocate (xyzf_out(ib), nfout, ios)

      do n = 1, ibnum - 1
         xyzf_out(ib)%q(n,:,:,:) = xyzf_in(ib)%q(n,:,:,:)
      end do

      do n = ibnum + 1, nfin
         xyzf_out(ib)%q(n-1,:,:,:) = xyzf_in(ib)%q(n,:,:,:)
      end do

!!!   call q_block_io (2, lunqout, formatted, nfout, mi, mj, mk, &
!!!                    xyzf_out(ib)%q, ios)

      if (formatted) then
         write (lunqout, '(7f3.0, i2, 50f3.0)', iostat=ios) &
            rz,rz,rz,rz,rz,rz,rz,iz,rz,rz,rz,(rz,n=1,max(2,nspecies)),rz,rz,rz
         write (lunqout, '(6es18.10)', iostat=ios) &
            ((((xyzf_out(ib)%q(n,i,j,k),i=1,mi),j=1,mj),k=1,mk),n=1,nfout)
      else
         write (lunqout, iostat=ios) &
            rz,rz,rz,rz,rz,rz,rz,iz,rz,rz,rz,(rz,n=1,max(2,nspecies)),rz,rz,rz
         write (lunqout, iostat=ios) &
            ((((xyzf_out(ib)%q(n,i,j,k),i=1,mi),j=1,mj),k=1,mk),n=1,nfout)
      end if

      if (ios /= 0) then
         write (luncrt, '(a, i5)') '*** Trouble writing function block:', ib
         go to 99
      end if

      deallocate (xyzf_in(ib)%q, xyzf_out(ib)%q)

      write (luncrt, '(a, i5)') '   Block completed:', ib

   end do  ! Next block

   deallocate (xyzf_in, xyzf_out)

99 continue

   end program iblank_move
