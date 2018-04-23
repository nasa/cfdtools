!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program v2c
!
!  Description:
!
!     V2C converts a grid and optional flow variables from cell vertices to cell
!     centers.  It treats 2-space and 3-space multiblock grids, automatically
!     determining whether the input grid is 2D or 3D and formatted or not.  The
!     output format must still be prompted for.
!
!     DPLR-type halo cells may be included or suppressed in the output file(s).
!     Halo cells are on the grid block boundaries.
!
!     The grid may represent a volume, a surface, or just a line.  There can be
!     no halo cells in degenerate directions.  Thus an input dimension n becomes
!     (in the output) either n + 1 (halos included), n - 1 (halos suppressed) or
!     1 if n = 1.
!
!     Original intentions of providing cell-centers-to-vertices options have
!     been abandoned as redundant.
!
!  Procedures:
!
!     DETERMINE_GRID_DIM    Determines whether a grid is 2D or 3D
!     DETERMINE_GRID_FORM   Distinguishes between formatted & unformatted grids
!     READER                Prompting utility
!     XYZQ_IO               I/O utilities for PLOT3D-type multiblock grid files
!     VERTICES_TO_CENTERS   Subroutine that does most of the work
!
!  History:
!
!     02/05/08  D.A.Saunders  Initial implementation (3-space only).
!     01/06/10    "     "     Todd White noticed that outputs with no halos
!                             weren't right at the high end in each index.
!                             Only subroutine vertices_to_centers needed a fix.
!     07/30/15    "     "     Added handling of 2-space files.
!     11/06/15    "     "     A 2D output function file was using the wrong
!                             logical formatting variable.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA.
!           Now with AMA, Inc. at NASA/ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure   ! Derived data type for one grid block
   use xyq_io_module          ! 2-space multiblock file I/O package
   use xyzq_io_module         ! 3-space multiblock file I/O package

!  Constants:

   implicit none

   integer, parameter :: &
      lungi  = 1,        &
      lunfi  = 2,        &
      lungo  = 3,        &
      lunfo  = 4,        &
      lunkbd = 5,        &
      luncrt = 6

   character, parameter :: &
      format * 11 = 'unformatted'

!  Variables:

   integer :: &
      i, i1, ib, ios, mode, nblock, ndim, nf, npts

   logical :: &
      cell_centered, cr, eof, formattedgi, formattedfi, &
      formattedgo, formattedfo, halos, initialize, threeD, twoD

   character (1) :: &
      answer

   character (64) :: &
      filegi, filefi, filego, filefo

!  Composite data types:

   type (grid_type), pointer, dimension (:) :: &
      blocki, blocko                           ! Grid block arrays

!  Execution:

!  Prompt for the input grid and check for existence and apparent form:

   ios = 1
   do while (ios /= 0)
      call reads (luncrt, 'Input grid file name (multiblock, 2D|3D):   ', &
                  lunkbd, filegi, cr, eof)
      if (eof) go to 999  ! Quit
      if (cr) cycle

      call determine_grid_form (filegi, lungi, formattedgi, ios)
   end do

   call determine_grid_dim (filegi, lungi, formattedgi, ndim, ios)

   twoD   = ndim == 2
   threeD = ndim == 3

   i1 = 1;  if (formattedgi) i1 = 3

   open (lungi, file=filegi, form=format(i1:11), status='OLD')

!  Accompanying function file?

   ios = 1
   do while (ios /= 0)
      call reads (luncrt, 'Accompanying function file? [<cr> = none]:  ', &
                  lunkbd, filefi, cr, eof)
      if (eof) go to 999  ! Quit
      if (cr) exit

      call determine_grid_form (filefi, lunfi, formattedfi, ios)
   end do

   if (cr) then
      nf = 0
   else
      nf = 1;  i1 = 1;  if (formattedfi) i1 = 3
      open (lunfi, file=filefi, form=format(i1:11), status='OLD')
   end if

!  Output grid name and formatting:

   ios = 1
   do while (ios /= 0)
      call reads (luncrt, 'Output cell-centered grid file name:        ', &
                  lunkbd, filego, cr, eof)
      if (eof) go to 999
      if (cr) cycle

      formattedgo = .true.  ! Bias in favor of DPLR pbca files
      call ready (luncrt, 'Formatted (y) or unformatted (n); <cr>=yes: ', &
                  lunkbd, formattedgo, cr, eof)
      if (eof) go to 999

      ios = 0
   end do

   i1 = 1;  if (formattedgo) i1 = 3

   open (lungo, file=filego, form=format(i1:11), status='UNKNOWN')

!  Output function file?  (Same format as the grid.)

   if (nf /= 0) then
      cr = .true.
      do while (cr)
         call reads (luncrt, 'Output function file name:                  ', &
                     lunkbd, filefo, cr, eof)
         if (eof) go to 999
      end do

      open (lunfo, file=filefo, form=format(i1:11), status='UNKNOWN')
      formattedfo = formattedgo
   end if

   halos = .true.
   call ready (luncrt, 'Include halo cells? (y|n; <cr> = yes): ', &
               lunkbd, halos, cr, eof)
   if (eof) go to 999

   mode = 2;  if (halos) mode = 1

!  Read the input grid header.
!  This finds nblock, allocates blocki(1:nblock), and sets block dimensions.

   if (twoD)   call  xy_header_io (1, lungi, formattedgi, nblock, blocki, ios)
   if (threeD) call xyz_header_io (1, lungi, formattedgi, nblock, blocki, ios)
   if (ios /= 0) go to 999

!  Read the input function file header if indicated:

   if (nf /= 0) then
      if (twoD) then
         call q_header_io_2d (1, lunfi, formattedfi, nblock, nf, blocki, ios)
      else
         call q_header_io    (1, lunfi, formattedfi, nblock, nf, blocki, ios)
      end if
      if (ios /= 0) go to 999
   end if

   allocate (blocko(nblock))

!  Determine the output block dimensions, which depend on halos = T or F:

   initialize = .true.

   do ib = 1, nblock
      if (twoD) then
         call vertices_to_centers_2d (initialize,mode,nf,blocki(ib), blocko(ib))
      else
         call vertices_to_centers (initialize, mode, nf, blocki(ib), blocko(ib))
      end if
   end do

   initialize = .false.

!  Write the output file header(s):

   if (twoD) then
      call xy_header_io (2, lungo, formattedgo, nblock, blocko, ios)
      if (ios /= 0) go to 999

      if (nf /= 0) then
         call q_header_io_2d (2, lunfo, formattedfo, nblock, nf, blocko, ios)
         if (ios /= 0) go to 999
      end if
   else
      call xyz_header_io (2, lungo, formattedgo, nblock, blocko, ios)
      if (ios /= 0) go to 999

      if (nf /= 0) then
         call q_header_io (2, lunfo, formattedfo, nblock, nf, blocko, ios)
         if (ios /= 0) go to 999
      end if
   end if

!  Process one block at a time:
!  ----------------------------

   if (twoD) then

      do ib = 1, nblock
         call xy_allocate (blocki(ib), ios)
         call xy_allocate (blocko(ib), ios)

         npts = blocki(ib)%ni * blocki(ib)%nj

         call xy_block_io (1, lungi, formattedgi, npts, &
                           blocki(ib)%x, blocki(ib)%y, ios)
         if (ios /= 0) go to 999

         if (nf /= 0) then
            call q_allocate_2d (blocki(ib), nf, ios)
            call q_allocate_2d (blocko(ib), nf, ios)
            call q_block_io_2d (1, lunfi, formattedfi, nf, &
                                blocki(ib)%mi, blocki(ib)%mj, blocki(ib)%q, ios)
            if (ios /= 0) go to 999
         end if

         call vertices_to_centers_2d (initialize,mode,nf,blocki(ib), blocko(ib))

         deallocate (blocki(ib)%x, blocki(ib)%y)

         npts = blocko(ib)%ni * blocko(ib)%nj

         call xy_block_io (2, lungo, formattedgo, npts, &
                           blocko(ib)%x, blocko(ib)%y, ios)

         deallocate (blocko(ib)%x, blocko(ib)%y)

         if (nf /= 0) then
            deallocate (blocki(ib)%q)
            call q_block_io_2d (2, lunfo, formattedfo, nf, &
                                blocko(ib)%mi, blocko(ib)%mj, blocko(ib)%q, ios)
            if (ios /= 0) go to 999
            deallocate (blocko(ib)%q)
         end if

      end do ! Next block

   else  ! 3D

      do ib = 1, nblock
         call xyz_allocate (blocki(ib), ios)
         call xyz_allocate (blocko(ib), ios)

         npts = blocki(ib)%ni * blocki(ib)%nj * blocki(ib)%nk

         call xyz_block_io (1, lungi, formattedgi, npts, &
                            blocki(ib)%x, blocki(ib)%y, blocki(ib)%z, ios)
         if (ios /= 0) go to 999

         if (nf /= 0) then
            call q_allocate (blocki(ib), nf, ios)
            call q_allocate (blocko(ib), nf, ios)
            call q_block_io (1, lunfi, formattedfi, nf, &
                             blocki(ib)%mi, blocki(ib)%mj, blocki(ib)%mk, &
                             blocki(ib)%q, ios)
            if (ios /= 0) go to 999
         end if

         call vertices_to_centers (initialize, mode, nf, blocki(ib), blocko(ib))

         deallocate (blocki(ib)%x, blocki(ib)%y, blocki(ib)%z)

         npts = blocko(ib)%ni * blocko(ib)%nj * blocko(ib)%nk

         call xyz_block_io (2, lungo, formattedgo, npts, &
                            blocko(ib)%x, blocko(ib)%y, blocko(ib)%z, ios)

         deallocate (blocko(ib)%x, blocko(ib)%y, blocko(ib)%z)

         if (nf /= 0) then
            deallocate (blocki(ib)%q)
            call q_block_io (2, lunfo, formattedfo, nf, &
                             blocko(ib)%mi, blocko(ib)%mj, blocko(ib)%mk, &
                             blocko(ib)%q, ios)
            if (ios /= 0) go to 999
            deallocate (blocko(ib)%q)
         end if

      end do ! Next block

   end if

   deallocate (blocki, blocko)
 
999 continue

! *** STOP ! Avoid system dependencies.

   end program v2c
