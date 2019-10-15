!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program surface_peaks
!
!  Description:
!
!     For a structured multiblock surface dataset, calculate two sets of peak
!  values and their locations: one set along grid lines in the i direction, and
!  one set in the j direction. Two more sets of mininum values are also computed
!  now, for tabulation only.  Initially, this was intended to help display peak
!  Shuttle wing leading edge heating, but it may have further uses.
!
!     The input file should be in Tecplot format (ASCII); output is a Tecplot
!  ASCII or binary file according to the indicated file name extension, in the
!  form of 2xnzones 1-dimensional zones, the first nzones for the i-direction
!  and the second for the j direction.  The locations of the maxima (only) are
!  stored as the output zone coordinates, because the minima locations are not
!  likely to be needed and would overwrite the peak locations in this scheme.
!  Global maximum and minimum function values are also tabulated with their
!  locations on standard output.
!
!     In order to group the four sets of results properly, the entire surface is
!  read.  All zones are processed for peaks in the i direction then in the j
!  direction before results are written as a full dataset.  (Reading and writing
!  one block at a time would be inappropriate here.)  [?? I'm not sure why,
!  twelve years later.  Probably to do with doubling the number of zones and
!  duplicating x,y,z as a way of keeping the index directions separate.]
!
!     Note that the maxima at each index line of each surface patch could be
!  hard to interpret if the patches were not indexed consistently.
!
!  History:
!
!     12/10/06  D.A.Saunders  Initial adaptation of SURFACE_CURVATURE.
!                             The number of output zones was doubled in order to
!                             distinguish the i and j directions, and only the
!                             locations of the maxima in each zone were written
!                             as the zone x/y/z coordinates.
!     05/06/19    "     "     Added tabulation of zone peaks on standard output,
!                             along with global maximum and minimum.
!                             Avoiding doubling of the zone count in the output
!                             dataset proved misguided, but it may be hard to
!                             make usefule plots from this dataset.
!     06/13/19    "     "     Belated completion of prints to standard output
!                             and i/jdirection.dat after a hiatus.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Ctr, Moffett Field, CA
!           Now with AMA, Inc.  at NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_header_structure  ! See Tecplot_io.f90 for this module
   use grid_block_structure   !  "     "     "
   use tecplot_io_module      !  "     "     "

   implicit none

!  Constants:

   integer, parameter :: &
      lunkbd  = 5,       &    ! Keyboard inputs
      luncrt  = 6,       &    ! Screen  outputs
      lunin   = 1,       &    ! Input surface dataset (Tecplot ASCII)
      lunout  = 2,       &    ! Output dataset (potentially plottable, 2x zones)
      lunidir = 3,       &    ! Tabulated mins/maxes in the i direction
      lunjdir = 4,       &    ! Tabulated mins/maxes in the j direction
      lunglbl = 7,       &    ! Tabulated global mins/maxes
      ndim    = 3

   real, parameter :: &
      zero = 0.

!  Variables:

   integer :: &
      i, ib, ibmax, ibmin, ifmax, ifmin, igmax, igmin, ifun, ios, &
      j, jb, jbmax, jbmin, jfmax, jfmin, jgmax, jgmin, l, nblocks, nbout, &
      nf, ni, nj, numf

   real :: &
      fmax, fmin, gmax, gmin

   logical :: &
      formatted

   character (132) :: &
      filename

!  Composite data types:

   type (grid_header) :: &
      header_in, header_out

   type (grid_type), pointer, dimension (:) :: &
      surf_in, surf_out

!  Tecplot function:

   integer :: TecEnd110

!  Execution:

   write (luncrt, '(/, a)', advance='no') ' Input surface dataset:  '
   read  (lunkbd, *, iostat=ios) filename;  l = len_trim (filename)
   if (ios /= 0) go to 99
   header_in%filename     = filename;  formatted = filename(l-2:l) == 'dat'
   header_in%formatted    = formatted
   header_in%ndim         = ndim

   write (luncrt, '(a)', advance='no') &
      ' Function number to locate extremes of: '
   read  (lunkbd, *, iostat=ios) ifun
   if (ios /= 0) go to 99

   write (luncrt, '(a)', advance='no') &
      ' Output surface dataset (*.dat|*.plt): '
   read  (lunkbd, *, iostat=ios) filename;  l = len_trim (filename)
   if (ios /= 0) go to 99
   header_out%filename    = filename;  formatted = filename(l-2:l) /= 'plt'
   header_out%formatted   = formatted
   header_out%ndim        = ndim
   header_out%numq        = 1
   header_out%datapacking = 0  ! Point order

!  Read the input dataset:

   call Tecplot_read (lunin, header_in, surf_in, ios)
   if (ios /= 0) then
      l = len_trim (header_in%filename)
      write (luncrt, '(/, 2a)') &
         ' Unable to read input dataset: ', header_in%filename(1:l)
      go to 99
   end if

   nf = header_in%numq
   if (nf < ifun) then
      write (luncrt, '(/, a, 2i4)') &
         ' Infeasible function choice.  ifun, nf: ', ifun, nf
      go to 99
   end if
   write (luncrt, '(2a)') ' Variable name: ', header_in%varname(ndim+ifun)

   nblocks = header_in%nblocks
   nbout   = nblocks + nblocks
   header_out%nblocks     = nbout
   header_out%datapacking = 0  ! POINT order
   header_out%title       = header_in%title
   header_out%ndatasetaux = 0

   allocate (header_out%varname(ndim+1))

   header_out%varname(1:ndim) = header_in%varname(1:ndim)
   header_out%varname(ndim+1) = header_in%varname(ndim+ifun)

   allocate (surf_out(nbout))

!  Read and scan one block at a time:

   open (lunidir, file='idirection.dat', status='unknown')
   open (lunjdir, file='jdirection.dat', status='unknown')

   write (lunidir, '(a)') 'Zone  j          fmin imin          fmax imax'
   write (lunjdir, '(a)') 'Zone  i          fmin jmin          fmax jmax'

   gmin = huge (gmin);  gmax = -gmin  ! Global extremes

   do ib = 1, nblocks
      jb = ib + nblocks
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
         fmin = huge (fmin);  fmax = -fmin
         do i = 1, ni
            if (surf_in(ib)%q(ifun,i,j,1) < fmin) then
                fmin  = surf_in(ib)%q(ifun,i,j,1)
                ifmin = i;  jfmin = j
            end if
            if (surf_in(ib)%q(ifun,i,j,1) > fmax) then
                fmax  = surf_in(ib)%q(ifun,i,j,1)
                ifmax = i;  jfmax = j
            end if
         end do
         write (lunidir, '(i3, i4, 2(es14.6, i5))') &
            ib, j, fmin, ifmin, fmax, ifmax

!        Save x/y/z for max only:

         surf_out(ib)%x(  j,1,1) = surf_in(ib)%x(ifmax,j,1)
         surf_out(ib)%y(  j,1,1) = surf_in(ib)%y(ifmax,j,1)
         surf_out(ib)%z(  j,1,1) = surf_in(ib)%z(ifmax,j,1)
         surf_out(ib)%q(1,j,1,1) = fmax

!        Global max and min over all blocks:

         if (fmax  > gmax) then
             gmax  = fmax
             igmax = ifmax;  jgmax = jfmax;  ibmax = ib
         end if
         if (fmin  < gmin) then
             gmin  = fmin
             igmin = ifmin;  jgmin = jfmin;  ibmin = ib
         end if
      end do

!     Peaks in the j direction:

      call Tec_block_allocate (surf_out(jb), ndim, 1, ios)
      if (ios /= 0) go to 99

      do i = 1, ni
         fmin = huge (fmin);  fmax = -fmin
         do j = 1, nj
            if (surf_in(ib)%q(ifun,i,j,1) < fmin) then
                fmin  = surf_in(ib)%q(ifun,i,j,1)
                ifmin = i;  jfmin = j
            end if
            if (surf_in(ib)%q(ifun,i,j,1) > fmax) then
                fmax  = surf_in(ib)%q(ifun,i,j,1)
                ifmax = i;  jfmax = j
            end if
         end do
         write (lunjdir, '(i3, i4, 2(es14.6, i5))') &
            ib, i, fmin, ifmin, fmax, ifmax

!        Save x/y/z for j-direction max (only):

         surf_out(jb)%x(  i,1,1) = surf_in(ib)%x(i,jfmax,1)
         surf_out(jb)%y(  i,1,1) = surf_in(ib)%y(i,jfmax,1)
         surf_out(jb)%z(  i,1,1) = surf_in(ib)%z(i,jfmax,1)
         surf_out(jb)%q(1,i,1,1) = fmax
      end do

   end do ! Next block

   open (lunglbl, file='extrema.dat', status='unknown')
   do l = luncrt, lunglbl, lunglbl - luncrt
      write (l, '(es14.6, 2x, 3es14.6, 3i5, a)') gmin, &
         surf_in(ibmin)%x(igmin,jgmin,1), surf_in(ibmin)%y(igmin,jgmin,1), &
         surf_in(ibmin)%z(igmin,jgmin,1), ibmin, igmin, jgmin, '   # minimum', &
         gmax, &
         surf_in(ibmax)%x(igmax,jgmax,1), surf_in(ibmax)%y(igmax,jgmax,1), &
         surf_in(ibmax)%z(igmax,jgmax,1), ibmax, igmax, jgmax, '   # maximum'
   end do

   call Tecplot_write (lunout, header_out, surf_out, ios)
   if (ios /= 0) go to 99

   if (header_out%formatted) then
      close (lunout)
   else
      ios = TecEnd110 ()
   end if

   deallocate (header_in%varname, header_out%varname)
   call deallocate_blocks (1, nblocks, ndim, nf, surf_in,  ios)
   call deallocate_blocks (1, nbout,   ndim, 1,  surf_out, ios)

99 continue

!! stop ! Avoid system dependencies

   end program surface_peaks
