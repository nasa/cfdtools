!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program p3d2tec
!
!  Description:
!
!     Convert between PLOT3D and Tecplot file formats, in either direction.
!     Structured multiblock grids in 3-space with optional function files are
!     assumed.  File names are prompted for, and standard file extensions are
!     used to distinguish between formatted and unformatted files.  Note that
!     reading unformatted Tecplot files is not an option, though:
!
!        *.g[u] [and *.f[u]]  <-->  *.dat | *.plt   in BLOCK or POINT order
!
!     Blocks are processed one at a time (no need to store whole files).
!
!     This version looks for an optional 'p3d2tec.inp.2' file if Tecplot output
!     is implied.  This should contain the desired variable names "" delimited:
!
!        "x, m"            [Just one per line]
!        "y, m"
!        "z, m"
!        "Cp"
!        "qw, W/cm^2"
!        "T, K"
!        "..."
!
!     If this file is not found, or the wrong number of names is found, the
!     function names default to "function 1", "function 2", ...
!
!  Option to Insert Arc Lengths:
!
!     For the case of a flow solution interpolated along one or more lines of
!     sight normal to the geometry (see LINES_OF_SIGHT and FLOW_INTERP by the
!     present author, which work with PLOT3D files), it can be convenient to
!     plot such profiles against distance from the wall.  Therefore, if all
!     input blocks are found to be lines (only one dimension > 1), the user is
!     prompted about an option to insert cumulative arc lengths after the "z"
!     coordinate and before the first flow function variable.  The option is
!     suppressed by entering an empty name for arc length.  Otherwise, enter
!     something like "s, m" or "s, in".
!
!  History:
!
!     11/22/06  D. A. Saunders  Initial implementation.
!     11/27/06  "  "  "  "  "   Added the option to supply variable names.
!                               The *.inp.2 file follows other utility usage.
!     08/22/08  "  "  "  "  "   Added the option to insert arc lengths if the
!                               input grid blocks all have only one dimension
!                               greater than 1.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Ctr, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_header_structure  ! See Tecplot_io.f90 for this module
   use grid_block_structure   !  "     "     "
   use tecplot_io_module      !  "     "     "
   use xyzq_io_module

   implicit none

!  Constants:

   integer, parameter :: &
      lun1a  = 1,        &    ! Input grid or Tecplot file
      lun1b  = 2,        &    ! Input function file, if any
      lun2a  = 3,        &    ! Output grid or Tecplot file
      lun2b  = 4,        &    ! Output function file, if any
      lunkbd = 5,        &    ! Keyboard inputs
      luncrt = 6,        &    ! Screen
      lunvar = 7,        &    ! Optional p3d2tec.inp.2 for variable names
      ndim   = 3

   character, parameter :: &
      blank * 1   = ' ',   &
      format * 11 = 'unformatted'

!  Variables:

   integer :: &
      i, i1, ib, iext, ios, j, k, last, lunf, n, nblocks, nf, nfin,            &
      ni, nj, nk, npts, nvar

   real :: &
      total_arc

   real, allocatable :: &
      arcs(:)

   logical :: &
      formatted_in, formatted_out, insert_arcs, not_lines, plot3d_in

   character :: &
      answer * 1, arc_length_name * 32, filename * 64

!  Composite data types:

   type (grid_header) :: &
      header

   type (grid_type) :: one_line

   type (grid_type), pointer, dimension (:) :: &
      grid

!  Tecplot function:

   integer :: TecEnd110

!  Execution:

   write (luncrt, '(/, a)') &
      ' Conversion between PLOT3D *.g[u] [and *.f[u]] and Tecplot *.dat | *.plt'

!  Prompt for the [first?] input file:

   write (luncrt, '(/, a)', advance='no') ' Input file name: '
   read  (lunkbd, *) filename
   last =  len_trim (filename)
   iext = 1 + index (filename(1:last), '.', BACK=.true.)
   plot3d_in    = filename(iext:iext) == 'g'
   formatted_in = filename(last:last) /= 'u' .or. filename(iext:last) == 'dat'
   i1 = 1;  if (formatted_in) i1 = 3

!  Open it, read the header info., and look for a second file if any.
!  Also, set up the output file(s) and write the header(s):

   header%ndim = ndim  ! All cases

   if (plot3d_in) then

      open (lun1a, file=filename, form=format(i1:11), status='old', iostat=ios)

      if (ios /= 0) then
         write (luncrt, '(/, 2a)') ' Unable to open file ', filename(1:last)
         go to 99
      end if

      call xyz_header_io (1, lun1a, formatted_in, nblocks, grid, ios)
      if (ios /= 0) go to 99

!     Are all blocks just lines?

      not_lines = .false.
      do ib = 1, nblocks
         ni = grid(ib)%ni;  nj = grid(ib)%nj;  nk = grid(ib)%nk
         if (ni * nj * nk > max (ni, nj, nk)) not_lines = .true.
      end do
      insert_arcs = .not. not_lines

      write (luncrt, '(a)', advance='no') ' Input function file name | none: '
      read  (lunkbd, *) filename
      last =  len_trim (filename)

      if (filename(1:last) == 'none') then
         nfin = 0
         lunf = -lun1b
      else
         lunf =  lun1b
         open (lun1b, file=filename, form=format(i1:11), status='old',         &
               iostat=ios)
         if (ios /= 0) then
            write (luncrt, '(/, 2a)') ' Unable to open file ', filename(1:last)
            go to 99
         end if
         call q_header_io (1, lun1b, formatted_in, nblocks, nfin, grid, ios)
         if (ios /= 0) go to 99
      end if

      nf = nfin  ! For output, unless ...

!     Option to insert arc lengths for flow profile plots:

      if (insert_arcs) then

         write (luncrt, '(2a)') &
            ' All target input blocks are lines, so you have the option of ',  &
            'inserting arc lengths.'
         write (luncrt, '(2a)', advance='no') &
            ' Arc length variable name [e.g., "s, m" including quotes; ',      &
            '"" => do not insert]: '
         read  (lunkbd, *) arc_length_name
         if (len_trim (arc_length_name) > 2) then
            nf = nf + 1
         else
            insert_arcs = .false.
         end if

      end if

      header%numq        = nf
      header%ndim        = ndim
      header%nblocks     = nblocks
      header%ndatasetaux = 0

!     Prompt for the output Tecplot file details:

      write (luncrt, '(a)', advance='no') ' Output Tecplot file name: '
      read  (lunkbd, *)  header%filename
      last  =  len_trim (header%filename)
      formatted_out    = header%filename(last-2:last) == 'dat'
      header%formatted = formatted_out
      header%title     = header%filename(1:last)

      write (luncrt, '(a)', advance='no') ' POINT or BLOCK output? [p | b]: '
      read  (lunkbd, *) answer
      header%datapacking = 0
      if (answer == 'b') header%datapacking = 1

      nvar = ndim + nf

      allocate (header%varname(nvar))

!     Look for variable names in the optional file p3d2tec.inp.2, or set up
!     generic names:

      call set_var_names (lunvar, nvar, header%varname) ! Internal procedure

      call Tec_header_write (lun2a, header, ios)

      if (ios /= 0) go to 99

   else ! Tecplot file input

      insert_arcs = .false.  ! For now

      header%filename  = filename(1:last)
      header%formatted = formatted_in

      call Tec_header_read (lun1a, header, grid, ios)

      if (ios /= 0) go to 99

      nf      = header%numq
      nblocks = header%nblocks

!     Prompt for the output PLOT3D file details:

      write (luncrt, '(a)', advance='no') ' Output PLOT3D grid name: '
      read  (lunkbd, *) filename
      last  = len_trim (filename)
      formatted_out = filename(last:last) /= 'u'
      i1 = 1;  if (formatted_out) i1 = 3

      open (lun2a, file=filename, form=format(i1:11), status='unknown')

      call xyz_header_io (2, lun2a, formatted_out, nblocks, grid, ios)
      if (ios /= 0) go to 99

      if (nf > 0) then
         i = last - 1
         if (formatted_out) i = last 
         filename(i:i) = 'f'

         open (lun2b, file=filename, form=format(i1:11), status='unknown')

         call q_header_io (2, lun2b, formatted_out, nblocks, nf, grid, ios)
         if (ios /= 0) go to 99
      end if

   end if

!  Process one block at a time:

   do ib = 1, nblocks 

      call Tec_block_allocate (grid(ib), ndim, nf, ios)
      if (ios /= 0) go to 99

      ni = grid(ib)%ni;  nj = grid(ib)%nj;  nk = grid(ib)%nk
      npts = ni * nj *nk

      if (plot3d_in) then

         call xyz_block_io (1, lun1a, formatted_in, npts, grid(ib)%x,          &
                            grid(ib)%y, grid(ib)%z, ios)
         if (ios /= 0) then
            write (luncrt, '(/, a, i4)') ' Trouble reading block #', ib
            go to 99
         end if

         grid(ib)%zone_title   = blank
         write (grid(ib)%zone_title, '(a, i4)') 'Zone:', ib
         grid(ib)%nzoneaux     = 0
         grid(ib)%solutiontime = -999.

         if (insert_arcs) then

            if (lunf > 0) then  ! We have to buffer the input flow data

               allocate (one_line%q(nfin,ni,nj,nk))

               call q_block_io (1, lunf, formatted_in, nfin, grid(ib)%mi,      &
                                grid(ib)%mj, grid(ib)%mk, one_line%q, ios)
               if (ios /= 0) go to 99 

               grid(ib)%q(2:nf,:,:,:) = one_line%q(1:nfin,:,:,:)

               deallocate (one_line%q)

            end if

         else  ! No need to buffer the flow data

            if (lunf > 0) then

               call q_block_io (1, lunf, formatted_in, nfin, grid(ib)%mi,      &
                                grid(ib)%mj, grid(ib)%mk, grid(ib)%q, ios)
               if (ios /= 0) go to 99 

            end if

         end if

      else  ! Tecplot input

         call Tec_block_read (lun1a, header, grid(ib), ios)

         if (ios /= 0) then
            write (luncrt, '(/, a, i4)') ' Trouble reading block #', ib
            go to 99
         end if

      end if

!     Insert arc lengths for blocks that are really lines?

      if (insert_arcs) then

         allocate (arcs(npts))

         call chords3d (npts, grid(ib)%x, grid(ib)%y, grid(ib)%z, .false.,     &
                        total_arc, arcs)
         n = 0
         do k = 1, nk
            do j = 1, nj
               do i = 1, ni
                  n = n + 1
                  grid(ib)%q(1,i,j,k) = arcs(n)
               end do
            end do
         end do

         deallocate (arcs)

      end if

      if (plot3d_in) then

         header%numq = nf  ! In case arc lengths were inserted

         call Tec_block_write (lun2a, header, grid(ib), ios)

      else

         call xyz_block_io (2, lun2a, formatted_out, npts, grid(ib)%x,         &
                            grid(ib)%y, grid(ib)%z, ios)
         if (nf > 0) then
            if (ios /= 0) go to 99
            call q_block_io (2, lun2b, formatted_out, nf, grid(ib)%mi,         &
                             grid(ib)%mj, grid(ib)%mk, grid(ib)%q, ios)
         end if
                                      
      end if

      if (ios /= 0) go to 99

      call deallocate_blocks (ib, ib, ndim, nf, grid, ios)
      if (ios /= 0) go to 99

   end do ! Next block

!  Clean up:

   if (plot3d_in) then
      if (formatted_out) then
         close (lun2a)
      else
         ios = TecEnd110 ()
      end if
   else
      close (lun2a)
      if (nf > 0) close (lun2b)
   end if

99 continue

!! stop ! Avoid system dependencies

!  Internal procedure for program p3d2tec:

   contains

!     --------------------------------------------------------------------------

      subroutine set_var_names (lunvar, nvar, varname)

!     Assign generic variable names but also look for an optional file named
!     'p3d2tec.inp.2' containing one variable name per line, including x, y, and
!     z, all bracketed by double quotes.  Also, handle the option to insert arc
!     lengths after z, even if no PLOT3D function file was input.

!     --------------------------------------------------------------------------

      implicit none

!     Arguments:

      integer,   intent (in)  :: lunvar  ! Logical unit, opened & closed here
      integer,   intent (in)  :: nvar    ! # names, including x, y, z names
      character, intent (out) :: varname (nvar) * (*)  ! Variable names to use

!     Local constants:

      integer,   parameter :: name_limit = 32  ! Longest within Tecplot_io.f90
      character, parameter :: quotes * 1 = '"'

!     Local variables:

      integer   :: first, last, mark, n
      character :: buffer * (name_limit + 2) ! Allow for ""

!     Execution:

      open (lunvar, file='p3d2tec.inp.2', status='old', iostat=ios)

      if (ios == 0) then

         n = 0
         if (insert_arcs) n = 1  ! Put the function names in the right place;
                                 ! shift the x, y, z, names backwards later
         do ! Until EOF

            read (lunvar, '(a)', iostat=ios) buffer
            if (ios /= 0) exit

            last = len_trim (buffer);  first = 1

            call scan4 (buffer(1:last), quotes, first, last, mark)

            if (mark < 0) exit

            n = n + 1
            if (n <= nvar) then
               varname(n) = buffer(first:mark)
            else
               write (luncrt, '(2a)') &
                  ' *** Excess variable name: ', buffer(first:mark)
            end if

         end do

         close (lunvar)

         if (n == nvar) ios = 0

         if (insert_arcs .and. ios == 0) then
            varname(1) = varname(2)
            varname(2) = varname(3)
            varname(3) = varname(4)
            varname(4) = arc_length_name
         end if

      end if

      if (ios /= 0) then  ! Use generic names
         varname(1) = 'x'
         varname(2) = 'y'
         varname(3) = 'z'
         n = ndim
         if (insert_arcs) then
            varname(4) = arc_length_name
            n = n + 1
         end if

         do i = 1, nfin
            write (varname(n + i), '(a, i3)') 'function', i
         end do
      end if

      end subroutine set_var_names

   end program p3d2tec
