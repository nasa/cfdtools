!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program surface_peaks
!
!  Original Description:
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
!  This Version:
!
!     An option has been retrofitted that could be used to locate the stagnation
!  point in a single dataset or (as for what prompted the extension) to show how
!  the stag. point varies with angle of attack.
!
!     This option is invoked by entering 0 for the function number.  It assumes
!  that the shear stress components are in the input dataset, but the magnitude
!  of shear stress is not (as is the case for output from BLAYER).  If this mag-
!  nitude is actually present, simply enter its function number.  Then the min.
!  line in extrema.dat is relevant. If |tauw| is not present, and the requested
!  function number is 0, a function name starting with tau is located and taken
!  to be the x component, followed by the y and z components. The result is
!  written to tauw.dat rather than extrema.dat to avoid overwriting original
!  usages.
!
!     Since the minimum |tauw| grid point may actually be on the aft body rather
!  than very close to the true forebody stag. point, this option looks for an
!  x (presumably near that of the max. diameter) beyond which surface points are
!  ignored in the search.  Enter this x in place of the output dataset line in
!  the control file, line 3.  [No: enter it on line 4 now; none is an option for
!  line 3 as well. This x applies to the extrema.dat file only.]]
!
!     Ancillary files are also written:  idirection.dat and jdirection.dat
!  These indicate peaks along each index line for each surface zone.  E.g.,:
!
!        Zone  j          fmin imin          fmax imax
!          1   1  2.584626E+05   19  2.826368E+05   13
!          1   2  2.581957E+05   19  2.819879E+05   13
!          1   3  2.574232E+05   19  2.805266E+05   13
!          :   :   :       :      :   :       :      :
!
!  Control:
!
!     Either answer a few prompts, or enter something like this on standard in:
!
!        blayer.dat       ! Input surface dataset (Tecplot ASCII)
!        17               ! Function number to treat; 0 => |tauw|
!        peaks.dat        ! Output Tecplotable file name, or none
!        0.9999           ! x beyond which to ignore surface points
!
!     If the function number is not 0, the output file contains (x,y,z,f) for
!  each surface grid line in the i direction followed by likewise for the j
!  direction.  No such file is written if the input file name is 'none'.
!
!     File extrema.dat is written to show the global max. and min., except this
!  is written as tauw.dat instead if ifun = 0 to avoid clobbering existing files
!  named extrema.dat.
!
!     [Later:]  Now, the global min/max file is written as extrema_nnn.dat for
!  function number nnn < 1000 as a better way to avoid overwriting.
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
!                             make useful plots from this dataset.
!     06/13/19    "     "     Belated completion of prints to standard output
!                             and i/jdirection.dat after a hiatus.
!     01/27/2023  "     "     Retrofitted treatment of the stag. point location
!                             as an option (ifun = 0 on line 2 and xignore on
!                             line 3 of the control file).  tauw.dat is written.
!     01/31/2023  "     "     File extrema.dat is now named extrema_n.dat to
!                             avoid overwriting an output for a different fn. n.
!                             The xignore input is now entered on line 4 so it
!                             can apply to all cases.  The output file can be
!                             suppressed by entering the name as none.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
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
      j, jb, jbmax, jbmin, jfmax, jfmin, jgmax, jgmin, l, length, &
      nblocks, nbout, nf, ni, nj, numf

   real :: &
      fmax, fmin, gmax, gmin, xignore

   logical :: &
      formatted, ij_line_output

   character (132) :: &
      extrema_name, filename

!  Derived data types:

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
   header_in%filename  = filename;  formatted = filename(l-2:l) == 'dat'
   header_in%formatted = formatted
   header_in%ndim      = ndim

   write (luncrt, '(a)', advance='no') &
      ' Function number to locate extremes of [0 = stag. pt. location]: '
   read  (lunkbd, *, iostat=ios) ifun
   if (ios /= 0) go to 99

   if (ifun > 0) write (luncrt, '(a)', advance='no') &
      ' Output surface dataset (*.dat|*.plt): '
   read  (lunkbd, *, iostat=ios) filename;  l = len_trim (filename)
   ij_line_output = filename(1:l) /= 'none'

   write (luncrt, '(a)', advance='no') &
      ' x beyond which to suppress surface point search: '
   read  (lunkbd, *, iostat=ios) xignore
   if (ios /= 0) go to 99

!  It is easier to leave the output file handling intact; just don't write it
!  if the file name is 'none.

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

   nblocks = header_in%nblocks
   nf      = header_in%numq

   if (ifun == 0) then  ! Treat the retrofitted stag. point case separately
      call stag_pt_case ()
      go to 99
   end if

   if (nf < ifun) then
      write (luncrt, '(/, a, 2i4)') &
         ' Infeasible function choice.  ifun, nf: ', ifun, nf
      go to 99
   end if
   write (luncrt, '(2a)') ' Variable name: ', header_in%varname(ndim+ifun)

   nbout                  = nblocks + nblocks
   header_out%nblocks     = nbout
   header_out%datapacking = 0  ! POINT order
   header_out%title       = header_in%title
   header_out%ndatasetaux = 0

   allocate (header_out%varname(ndim+1))

   header_out%varname(1:ndim) = header_in%varname(1:ndim)
   header_out%varname(ndim+1) = header_in%varname(ndim+ifun)

   allocate (surf_out(nbout))

!  Scan one block at a time:

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
            if (surf_in(ib)%x(i,j,1) > xignore) cycle
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
            if (surf_in(ib)%x(i,j,1) > xignore) cycle
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

   call numbered_name ('extrema_', ifun, extrema_name, length)
   extrema_name(length+1:length+4) = '.dat';  length = length + 4

   open (lunglbl, file=extrema_name(1:length), status='unknown')

   do l = luncrt, lunglbl, lunglbl - luncrt
      write (l, '(es14.6, 2x, 3es14.6, 3i5, a)') &
         gmin, &
         surf_in(ibmin)%x(igmin,jgmin,1), surf_in(ibmin)%y(igmin,jgmin,1),    &
         surf_in(ibmin)%z(igmin,jgmin,1), ibmin, igmin, jgmin, '  # minimum', &
         gmax, &
         surf_in(ibmax)%x(igmax,jgmax,1), surf_in(ibmax)%y(igmax,jgmax,1),    &
         surf_in(ibmax)%z(igmax,jgmax,1), ibmax, igmax, jgmax, '  # maximum'
   end do
   close (lunglbl)

   if (ij_line_output) then
      call Tecplot_write (lunout, header_out, surf_out, ios)
      if (ios /= 0) go to 99

      if (header_out%formatted) then
         close (lunout)
      else
         ios = TecEnd110 ()
      end if
   end if

   deallocate (header_in%varname, header_out%varname)
   call deallocate_blocks (1, nblocks, ndim, nf, surf_in,  ios)
   call deallocate_blocks (1, nbout,   ndim, 1,  surf_out, ios)

99 continue

!! stop ! Avoid system dependencies

!  Internal procedure for program surface_peaks:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine stag_pt_case ()  ! ifun = 0 has been read; look for tauw com-
                                  ! ponents and compute |tauw| min/max.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i, i1, j
      real    :: tauw
      
!     Look for three function names with tau in them:

      do i = 1, nf
         i1 = index (header_in%varname(ndim+i), 'tau')
         if (i1 > 0 .and. i < nf) then
             i1 = index (header_in%varname(ndim+i+1), 'tau')
             if (i1 > 0 .and. i < nf) then
                 i1 = index (header_in%varname(ndim+i+2), 'tau')
                 if (i1 > 0) ifun = i
                 exit
             end if
         end if
      end do

      if (i1 <= 0) then
         write (luncrt, '(a)') &
            '*** Can''t find contiguous tauwx, tauwy, tauwz. Aborting.'
         go to 99
      end if

!     Find the max and min |tauw|:

      gmin = huge (gmin);  gmax = -gmin  ! Global extremes

      do ib = 1, nblocks
         do j = 1, surf_in(ib)%nj
            do i = 1, surf_in(ib)%ni
               if (surf_in(ib)%x(i,j,1) > xignore) cycle

               tauw = sqrt (surf_in(ib)%q(ifun,  i,j,1)**2 + &
                            surf_in(ib)%q(ifun+1,i,j,1)**2 + &
                            surf_in(ib)%q(ifun+2,i,j,1)**2)
               if (tauw < gmin) then
                  ibmin = ib;  igmin = i;  jgmin = j; gmin = tauw
               end if
               if (tauw > gmax) then
                  ibmax = ib;  igmax = i;  jgmax = j; gmax = tauw
               end if
            end do
         end do
      end do

      open (lunglbl, file='tauw.dat', status='unknown')
      do l = luncrt, lunglbl, lunglbl - luncrt
         write (l, '(es14.6, 2x, 3es14.6, 3i5, a)') &
            gmin, &
            surf_in(ibmin)%x(igmin,jgmin,1), surf_in(ibmin)%y(igmin,jgmin,1),  &
            surf_in(ibmin)%z(igmin,jgmin,1), ibmin, igmin, jgmin, ' # minimum',&
            gmax, &
            surf_in(ibmax)%x(igmax,jgmax,1), surf_in(ibmax)%y(igmax,jgmax,1),  &
            surf_in(ibmax)%z(igmax,jgmax,1), ibmax, igmax, jgmax, ' # maximum'
      end do
      close (lunglbl)

 99   return

      end subroutine stag_pt_case

   end program surface_peaks
