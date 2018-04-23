!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program surface_peaks
!
!  Description:
!
!     For a structured multiblock surface dataset, calculate two sets of peak
!  values and their locations: one set along grid lines in the i direction, and
!  one set in the j direction.  Initially, this is intended to help display peak
!  wing leading edge heating, but it may have further uses.
!
!     The input file should be in Tecplot format (ASCII); output is a Tecplot
!  ASCII or binary file according to the indicating file name extension.
!
!     In order to group the two sets of results properly, the entire surface is
!  read.  All zones are processed for peaks in the i direction then in the j
!  direction before results are written as a full dataset.  (Reading and writing
!  one block at a time would be inappropriate here.)
!
!  History:
!
!     12/10/06  D. A. Saunders  Initial adaptation of SURFACE_CURVATURE.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Ctr, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_header_structure  ! See Tecplot_io.f90 for this module
   use grid_block_structure   !  "     "     "
   use tecplot_io_module      !  "     "     "

   implicit none

!  Constants:

   integer, parameter :: &
      lunkbd = 5,        &
      luncrt = 6,        &
      lunin  = 1,        &
      lunout = 2,        &
      ndim   = 3

   real, parameter :: &
      zero = 0.

!  Variables:

   integer :: &
      i, ib, ifun, ios, ipk, j, jb, jpk, len, nblocks, nbout, nf, ni, nj

   real :: &
      peak

   logical :: &
      formatted

   character :: &
      filename * 80

!  Composite data types:

   type (grid_header) :: &
      header_in, header_out

   type (grid_type), pointer, dimension (:) :: &
      surf_in, surf_out

!  Tecplot function:

   integer :: TecEnd110

!  Execution:

   write (luncrt, '(/, a)', advance='no') ' Input surface dataset:  '
   read  (lunkbd, *) filename;  len = len_trim (filename)
   header_in%filename     = filename;  formatted = filename(len-2:len) == 'dat'
   header_in%formatted    = formatted
   header_in%ndim         = ndim

   write (luncrt, '(a)', advance='no') ' Function number to locate peaks of: '
   read  (lunkbd, *) ifun

   write (luncrt, '(a)', advance='no') ' Output surface dataset: '
   read  (lunkbd, *) filename;  len = len_trim (filename)
   header_out%filename    = filename;  formatted = filename(len-2:len) /= 'plt'
   header_out%formatted   = formatted
   header_out%ndim        = ndim
   header_out%numq        = 1
   header_out%datapacking = 0  ! Point order

!  Read the input dataset:

   call Tecplot_read (lunin, header_in, surf_in, ios)

   if (ios /= 0) then
      len = len_trim (header_in%filename)
      write (luncrt, '(/, 2a)') &
         ' Unable to read input dataset: ', header_in%filename(1:len)
      go to 99
   end if

   nblocks = header_in%nblocks
   nbout   = 2 * nblocks
   nf      = header_in%numq
   if (nf < ifun) then
      write (luncrt, '(/, a, 2i4)') &
         ' Infeasible function choice.  ifun, nf: ', ifun, nf
      go to 99
   end if

   header_out%nblocks     = nbout
   header_out%datapacking = 0  ! POINT order
   header_out%title       = header_in%title
   header_out%ndatasetaux = 0

   allocate (header_out%varname(4))

   header_out%varname(1:ndim) = header_in%varname(1:ndim)
   header_out%varname(ndim+1) = header_in%varname(ndim+ifun)

   allocate (surf_out(nbout))

!  Process one block at a time, but the desired output order prevents writing:

   do ib = 1, nblocks

      jb = nblocks + ib
      surf_out(ib)%zone_title   = surf_in(ib)%zone_title
      surf_out(jb)%zone_title   = surf_in(ib)%zone_title
      surf_out(ib)%nzoneaux     = 0
      surf_out(jb)%nzoneaux     = 0
      surf_out(ib)%solutiontime = -999.  ! Undefined
      surf_out(jb)%solutiontime = -999.

      ni =  surf_in(ib)%ni;  nj =  surf_in(ib)%nj
      surf_out(ib)%ni = nj;  surf_out(ib)%mi = nj
      surf_out(ib)%nj = 1;   surf_out(ib)%mj = 1
      surf_out(ib)%nk = 1;   surf_out(ib)%mk = 1
      surf_out(jb)%ni = ni;  surf_out(jb)%mi = ni
      surf_out(jb)%nj = 1;   surf_out(jb)%mj = 1
      surf_out(jb)%nk = 1;   surf_out(jb)%mk = 1

!     Peaks in the i direction:

      call Tec_block_allocate (surf_out(ib), ndim, 1, ios)
      if (ios /= 0) go to 99

      do j = 1, nj
         peak = -1.e+30
         do i = 1, ni
            if (peak < surf_in(ib)%q(ifun,i,j,1)) then
                peak = surf_in(ib)%q(ifun,i,j,1)
                ipk  = i
            end if
         end do
         surf_out(ib)%x(  j,1,1) = surf_in(ib)%x(ipk,j,1)
         surf_out(ib)%y(  j,1,1) = surf_in(ib)%y(ipk,j,1)
         surf_out(ib)%z(  j,1,1) = surf_in(ib)%z(ipk,j,1)
         surf_out(ib)%q(1,j,1,1) = peak
      end do

!     Peaks in the j direction:

      call Tec_block_allocate (surf_out(jb), ndim, 1, ios)
      if (ios /= 0) go to 99

      do i = 1, ni
         peak = -1.e+30
         do j = 1, nj
            if (peak < surf_in(ib)%q(ifun,i,j,1)) then
                peak = surf_in(ib)%q(ifun,i,j,1)
                jpk  = j
            end if
         end do
         surf_out(jb)%x(  i,1,1) = surf_in(ib)%x(i,jpk,1)
         surf_out(jb)%y(  i,1,1) = surf_in(ib)%y(i,jpk,1)
         surf_out(jb)%z(  i,1,1) = surf_in(ib)%z(i,jpk,1)
         surf_out(jb)%q(1,i,1,1) = peak
      end do

      call deallocate_blocks (ib, ib, ndim, nf, surf_in, ios)

   end do ! Next block

   call Tecplot_write (lunout, header_out, surf_out, ios)
   if (ios /= 0) go to 99

   if (header_out%formatted) then
      close (lunout)
   else
      ios = TecEnd110 ()
   end if

   call deallocate_blocks (1, nbout, ndim, 1, surf_out, ios)

   deallocate (header_in%varname, header_out%varname)

99 continue

!! stop ! Avoid system dependencies

   end program surface_peaks
