!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module xyq_io_module

!  This is the 2D analogue of xyzq_io_module, from which the following description is retained with slight adjustments.

!  This module packages the I/O for a 2-space PLOT2D-type multiblock grid and/or associated flow field "q" file, with optional
!  blanking control.
!
!  The grid file is separate from the PLOT2D-type function file, which contains a variable number of flow field quantities.
!
!  However, the "grid_type" data structure adopted for LaRC compatibility obliges us to assume reading of a "q" file will follow
!  reading of a related grid file, because reading the grid requires allocation of the grid blocks (here), and these blocks have
!  the requisite "q" fields (which are allocated only if a "q" file is requested).
!
!  On the other hand, writing of an xy file need not precede writing of a q file.  These functions are separated.
!
!  Three pairs of lower level utilities used by the whole-grid I/O utilities are made public to improve flexibility, as in cases
!  where only one block of a grid needs to be processed at a time.  If blanking is desired, use the appropriate variant of the
!  grid routines, since the PLOT2D convention includes "iblank" in the same record as x, y.
!
!  Note that, in order to pass unallocated arrays as arguments, the calling program must declare the x/y/iblank/q data type with
!  the "pointer" attribute as opposed to "allocatable".
!
!  Prompting for input and output file names, and opening of the files, are provided here for completeness, but some applications
!  may choose to read file names from a control file and open the files explicitly.
!
!  History (3-space):
!
!  04/02/2004   D. A. Saunders   Initial implementation, for grid-morphing use at NASA LaRC (and other CFD utilities).
!  01/15/2010    "    "    "     After years of regretting the decision to avoid possible machine dependencies with advance='no',
!                                the file prompting now leaves the cursor at the end of the prompt.
!  03/03/2010    "    "    "     Added "iblank" variants of the grid utilities.
!  10/16/2014    "    "    "     Added xyq_2d_to_3d to ease treatment of a 2-D plane as plane j = 1 of a 3-D grid.
!  07/30/2015    "    "    "     Header output for function files needed 3i5 format, not 4i5.
!
!  History (2-space):
!
!  08/29/2012   D. A. Saunders   2D analogue, prompted by THIN_GRID_2D.
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, Moffett Field, CA
! 
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure ! Or equivalent to define "grid_type"

   implicit none

   private

   public :: file_prompt_2d ! Prompts for a file name and open it

   public :: xyq_read       ! Reads a multiblock grid and (optionally) an associated flow solution (PLOT2D-style)
   public :: xyq_2d_to_3d   ! Converts xy(ni,nj) data to xyz(ni,1,nk) data with y -> 0. and nk = nj
   public :: xy_write       ! Writes a multiblock grid (PLOT2D-style)
   public :: q_write_2d     ! Writes a flow solution   (  "     "   )

   public :: xy_header_io   ! Reads or writes grid header records
   public :: q_header_io_2d ! Reads or writes flow header records
   public :: xy_allocate    ! Allocates one grid block (deallocates can be done explicitly)
   public :: q_allocate_2d  ! Allocates one flow block (  "     "     "     "     "     " )
   public :: xy_block_io    ! Reads or writes one grid block
   public :: q_block_io_2d  ! Reads or writes one flow block

   public :: xyiq_read      ! "iblank" variants of the xy utilities above
   public :: xyi_write
   public :: xyi_allocate
   public :: xyi_block_io

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine file_prompt_2d (lun_file, file_type, status, format_prompt, formatted, iostatus)

!     Prompt for a file of the indicated type and status.  Standard input and output are used (logical unit *) for the prompting.
!     The journaling option from the 3D package has been dispensed with.  The indicated file is opened here.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer,   intent (in)    :: lun_file      ! Logical unit for the file; opened here
      character, intent (in)    :: file_type*(*) ! Description of the file type, lower case; e.g., 'input grid' or 'output q'
      character, intent (in)    :: status*(*)    ! File status: 'old' for input; 'new' or 'unknown' for output
      logical,   intent (in)    :: format_prompt ! True means prompt for and output the "formatted" argument, else use it as input
      logical,   intent (inout) :: formatted     ! True if an ASCII file is indicated; false means unformatted
      integer,   intent (out)   :: iostatus      ! 0 means no error detected;
                                                 ! diagnostics are written to standard output
!     Local constants:

      character, parameter :: format * 11 = 'unformatted'

!     Local variables:

      integer   :: i1, ntries
      character :: answer * 1, filename * 80

!     Execution:
!     ----------

      write (*, '(/, 3a)', advance='no') ' Enter the ', file_type, ' file name: '
      read  (*, '(a)') filename

      if (format_prompt) then
         write (*, '(/, a)', advance='no') ' Is the file formatted (y) or unformatted (n)? '
         read  (*, *) answer
         formatted = answer == 'y' .or. answer == 'Y'
      end if

      i1 = 1;  if (formatted) i1 = 3

      do ntries = 1, 3
         open (lun_file, file=filename, form=format(i1:11), status=status, iostat=iostatus)
         if (iostatus == 0) exit
         if (ntries < 3) then
            write (*, '(a)') ' Bad file name? Try again:'
            read  (*, '(a)') filename
         else
            write (*, '(a)') ' Unable to open the file.'
         end if
      end do

      end subroutine file_prompt_2d

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine xyq_read (lun_xy, lun_q, formatted, nblocks, num_q, cell_centered, xyq, iostatus)

!     Read a PLOT2D-type multiblock 2-space grid and (optionally) a flow-field in PLOT2D-type function file format.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)      :: lun_xy        ! Logical unit for the xy file;
                                                 ! assumed open upon entry; closed here before the return
      integer, intent (in)      :: lun_q         ! Logical unit for the optional "q" file;
                                                 ! assumed open upon entry; closed here before the return;
                                                 ! lun_q < 0 means suppress the reading of a q file;
                                                 ! lun_q = lun_xy is permissible
      logical, intent (in)      :: formatted     ! True for ASCII file(s); false means unformatted reads
      integer, intent (out)     :: nblocks       ! Number of grid blocks found in the grid file
      integer, intent (out)     :: num_q         ! Number of "q" file variables found (if lun_q > 0)
      logical, intent (out)     :: cell_centered ! True if the "q" file is found to be cell-centered
      type (grid_type), pointer :: xyq(:)        ! Array of grid blocks containing x,y and optional q data;
                                                 ! allocated here prior to the reading
      integer, intent (out)     :: iostatus      ! 0 means no error detected, else diagnostics are written to
                                                 ! standard output; early return follows
!     Local variables:

      integer :: ib, ios, mi, mj, mk, mq, ni, nj, nk, npts
      logical :: centers, vertices

!     Execution:
!     ----------

!     The input grid file is not optional:
!     ------------------------------------

      call xy_header_io (1, lun_xy, formatted, nblocks, xyq, ios)

      if (ios /= 0) then
         write (*, '(a)') ' Error allocating the input grid or reading its dimensions.'
         go to 999
      end if

      do ib = 1, nblocks

         call xy_allocate (xyq(ib), ios)

         if (ios /= 0) then
            write (*, '(a, i4)') ' Trouble allocating input grid block #', ib
            go to 999
         end if

         npts = xyq(ib)%ni * xyq(ib)%nj

         call xy_block_io (1, lun_xy, formatted, npts, xyq(ib)%x, xyq(ib)%y, ios)

         if (ios /= 0) then
            write (*, '(a, i4)')  ' Trouble reading input grid.  Block #:', ib
            write (*, '(a, 2i5)') ' Dimensions: ', xyq(ib)%ni, xyq(ib)%nj
            go to 999
         end if

      end do

      close (lun_xy)

      if (lun_q < 0) go to 999 ! Avoid the indenting


!     Optional input "q" file:
!     ------------------------

      call q_header_io_2d (1, lun_q, formatted, nblocks, mq, xyq, ios)

      if (ios /= 0) then
         write (*, '(a)') ' Error reading the input flow solution dimensions.'
         go to 999
      end if

      vertices = .true.;  centers = .true.

      do ib = 1, nblocks

         mi = xyq(ib)%mi;  mj = xyq(ib)%mj
         ni = xyq(ib)%ni;  nj = xyq(ib)%nj

         vertices = mi == ni   .and. mj == nj
         centers  = mi == ni-1 .and. mj == nj-1

         if (.not. vertices .and. .not. centers) then
            write (*, '(/, a, i4, /, a, 2i5)') ' Grid/flow field dimension mismatch.  Block:', ib, &
               ' Grid: ', ni, nj, ' Flow: ', mi, mj
            ios = 1
            go to 999
         end if

         call q_allocate_2d (xyq(ib), mq, ios)  ! Put "q" variables in 1st index (poor for I/O but better for interpolations)

         if (ios /= 0) then
            write (*, '(a, i4)') ' Trouble allocating input flow solution block #', ib
            go to 999
         end if

         call q_block_io_2d (1, lun_q, formatted, mq, mi, mj, xyq(ib)%q, ios)

         if (ios /= 0) then
            write (*, '(a, i4, a, i4)') ' Trouble reading input flow solution.  Block #:', ib, '  I/O status:', ios
            write (*, '(a, 3i5)') ' Dimensions and # functions: ', mi, mj, mq
            go to 999
         end if

      end do

      close (lun_q)

      num_q = mq
      cell_centered = centers

  999 iostatus = ios

      end subroutine xyq_read

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine xyq_2d_to_3d (nblocks, nf, xyzq, iostatus)

!     Convert an xy(ni,nj) 2-D plane to plane j = 1 of a 3-D grid.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)      :: nblocks       ! Number of input blocks; nblocks >= 1
      integer, intent (in)      :: nf            ! Number of input functions;   nf >= 0
      type (grid_type), pointer :: xyzq(:)       ! Array of grid blocks containing x,y[,q] data (input); x,z[,q] and y = 0. (output)
      integer, intent (out)     :: iostatus      ! 0 means no error detected, else diagnostics are written to
                                                 ! standard output; early return follows
!     Local variables:

      integer :: ib, ni, nj, mi, mj
      real, allocatable :: xt(:,:,:), ft(:,:,:)

!     Execution:

      do ib = 1, nblocks

         ni = xyzq(ib)%ni
         nj = xyzq(ib)%nj
         xyzq(ib)%nj = 1
         xyzq(ib)%nk = nj
         allocate (xt(ni,nj,2))
         xt(:,:,1) = xyzq(ib)%x(:,:,1)
         xt(:,:,2) = xyzq(ib)%y(:,:,1)
         deallocate (xyzq(ib)%x, xyzq(ib)%y)
         allocate   (xyzq(ib)%x(ni,1,nj), xyzq(ib)%y(ni,1,nj), xyzq(ib)%z(ni,1,nj), stat=iostatus)
         if (iostatus /= 0) then
            write (*, '(a, i5)') ' xyq_2d_to_3d trouble reallocating %x/y/z; block #:', ib
            go to 99
         end if
         xyzq(ib)%x(:,1,:) = xt(:,:,1)
         xyzq(ib)%y(:,1,:) = 0.
         xyzq(ib)%z(:,1,:) = xt(:,:,2)
         deallocate (xt)

         if (nf > 0) then
            mi = xyzq(ib)%mi
            mj = xyzq(ib)%mj
            xyzq(ib)%mj = 1
            xyzq(ib)%mk = mj
            allocate (ft(nf,mi,mj))
            ft(:,:,:) = xyzq(ib)%q(:,:,:,1)
            deallocate (xyzq(ib)%q)
            allocate   (xyzq(ib)%q(nf,mi,1,mj), stat=iostatus)
            if (iostatus /= 0) then
               write (*, '(a, i5)') ' xyq_2d_to_3d trouble reallocating %q; block #:', ib
               go to 99
            end if
            xyzq(ib)%q(:,:,1,:) = ft(:,:,:)
            deallocate (ft)
         end if

      end do

 99   return

      end subroutine xyq_2d_to_3d

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine xy_write (lun_xy, formatted, nblocks, xy, iostatus)

!     Write a PLOT2D-type multiblock 2-space grid to the indicated unit.
!     The file is assumed to be open upon entry and is closed before the return here.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)      :: lun_xy        ! Logical unit for the xy file;
                                                 ! assumed open upon entry; closed here before the return
      logical, intent (in)      :: formatted     ! True for ASCII file(s); false means unformatted writes
      integer, intent (inout)   :: nblocks       ! Number of grid blocks (an input, but inout because of the xy_header_io call)
      type (grid_type), pointer :: xy(:)         ! Array of grid blocks containing x,y data
      integer, intent (out)     :: iostatus      ! 0 means no error detected, else diagnostics are written to
                                                 ! standard output; early return follows
!     Local variables:

      integer   :: ib, ios, npts

!     Execution:
!     ----------

      call xy_header_io (2, lun_xy, formatted, nblocks, xy, ios)

      if (ios /= 0) then
         write (*, '(a)') ' Error writing the output grid dimensions.'
         go to 999
      end if

      do ib = 1, nblocks

         npts = xy(ib)%ni * xy(ib)%nj

         call xy_block_io (2, lun_xy, formatted, npts, xy(ib)%x, xy(ib)%y, ios)

         if (ios /= 0) then
            write (*, '(a, i4)')  ' Trouble writing output grid.  Block #:', ib
            write (*, '(a, 2i5)') ' Dimensions: ', xy(ib)%ni, xy(ib)%nj
            go to 999
         end if

      end do

      close (lun_xy)

  999 iostatus = ios

      end subroutine xy_write

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine q_write_2d (lun_q, formatted, nblocks, num_q, q, iostatus)

!     Write a multiblock flow-field in PLOT2D-type function file format.
!     The file is assumed to be open upon entry and is closed before the return here.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)      :: lun_q         ! Logical unit for the output q file;
                                                 ! assumed open upon entry; closed here before the return
      logical, intent (in)      :: formatted     ! True for ASCII file(s); false means unformatted writes
      integer, intent (inout)   :: nblocks       ! Number of grid blocks (an input, but inout because of the q_header_io call)
      integer, intent (inout)   :: num_q         ! Number of q variables ( "   "  ... )
      type (grid_type), pointer :: q(:)          ! Array of blocks containing q data
      integer, intent (out)     :: iostatus      ! 0 means no error detected, else diagnostics are written to
                                                 ! standard output; early return follows
!     Local variables:

      integer :: ib, ios

!     Execution:
!     ----------

      call q_header_io_2d (2, lun_q, formatted, nblocks, num_q, q, ios)

      if (ios /= 0) then
         write (*, '(a)') ' Error writing the output flow solution dimensions.'
         go to 999
      end if

      do ib = 1, nblocks

         call q_block_io_2d (2, lun_q, formatted, num_q, q(ib)%mi, q(ib)%mj, q(ib)%q, ios)

         if (ios /= 0) then
            write (*, '(a, i4, a, i4)') ' Trouble writing output flow solution.  Block #:', ib, '  I/O status:', ios
            write (*, '(a, 4i5)') ' Dimensions and # functions: ', q(ib)%mi, q(ib)%mj, num_q
            go to 999
         end if

      end do

      close (lun_q)

  999 iostatus = ios

      end subroutine q_write_2d

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine xy_header_io (mode, lun, formatted, nblocks, xy, ios)
!
!     Read or write header records for a 2-space PLOT2D-type grid.  Blocks are allocated here before reading their dimension fields.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)       :: mode       ! 1 = read, 2 = write
      integer, intent (in)       :: lun        ! Logical unit number
      logical, intent (in)       :: formatted  ! T|F
      integer, intent (inout)    :: nblocks    ! Number of blocks: output for a read; input for a write
      type (grid_type), pointer  :: xy(:)      ! Array of grid blocks: allocated here for a read
      integer, intent (out)      :: ios        ! 0 = no error

!     Local variables:

      integer :: ib

!     Execution:

      select case (mode)

         case (1) ! Allocate blocks and read header records

            if (formatted) then
               read (lun, *, iostat=ios) nblocks
            else
               read (lun,    iostat=ios) nblocks
            end if

            if (ios /= 0) then
               write (*, '(a)') ' Error reading the number of grid blocks.'
               go to 999
            end if

            allocate (xy(nblocks), stat=ios)

            if (ios /= 0) then
               write (*, '(a)') ' Error allocating the grid blocks.'
               go to 999
            end if

            if (formatted) then
               do ib = 1, nblocks  ! Skip any 1s for nk
                  read (lun, *, iostat=ios) xy(ib)%ni, xy(ib)%nj
               end do
            else
                  read (lun,    iostat=ios) (xy(ib)%ni, xy(ib)%nj, ib = 1, nblocks)
            end if

            if (ios /= 0) then
               write (*, '(a)') ' Error reading the input grid dimensions.'
               go to 999
            end if

         case (2) ! Write header records

            if (formatted) then
               write (lun, '(i5)', iostat=ios) nblocks
            else
               write (lun,         iostat=ios) nblocks
            end if

            if (ios /= 0) then
               write (*, '(a)') ' Error writing the number of grid blocks.'
               go to 999
            end if

            if (formatted) then
               write (lun, '(2i5)', iostat=ios) (xy(ib)%ni, xy(ib)%nj, ib = 1, nblocks)
            else
               write (lun,          iostat=ios) (xy(ib)%ni, xy(ib)%nj, ib = 1, nblocks)
            end if

            if (ios /= 0) then
               write (*, '(a)') ' Error writing the output grid dimensions.'
               go to 999
            end if

         case default

            ios = 999
            write (*, '(a, i6)') ' Bad mode argument passed to xy_header_io:', mode

      end select

  999 return

      end subroutine xy_header_io

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine q_header_io_2d (mode, lun, formatted, nblocks, mq, q, ios)
!
!     Read or write header records for a 2-space PLOT2D-type function file.  Blocks are assumed to have been allocated already
!     for reading of the associated grid.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)       :: mode        ! 1 = read, 2 = write
      integer, intent (in)       :: lun         ! Logical unit number
      logical, intent (in)       :: formatted   ! T|F
      integer, intent (in)       :: nblocks     ! Number of grid blocks: if a different number is found, this is an error
      type (grid_type), pointer  :: q(:)        ! Array of grid blocks with fields for the dimensions of the function ("q") data
      integer, intent (inout)    :: mq          ! Number of flow variables at each cell or vertex: output from a read, else input
      integer, intent (out)      :: ios         ! 0 = no error

!     Local variables:

      integer :: ib, nb

!     Execution:

      select case (mode)

         case (1) ! Read q file header records

            if (formatted) then
               read (lun, *, iostat=ios) nb
            else
               read (lun,    iostat=ios) nb
            end if

            if (ios /= 0) then
               write (*, '(a)') ' Error reading the number of flow solution blocks.'
               go to 999
            end if

            if (nb /= nblocks) then
               write (*, '(a, i4, a, i4)') ' Block count mismatch.  Grid:', nblocks, '   Flow solution:', nb
               ios = 1
               go to 999
            end if

            if (formatted) then
               do ib = 1, nb  ! Skip and 1s for nk
                  read (lun, *, iostat=ios) q(ib)%mi, q(ib)%mj, mq
               end do
            else
                  read (lun,    iostat=ios) (q(ib)%mi, q(ib)%mj, mq, ib = 1, nb)
            end if

            if (ios /= 0) then
               write (*, '(a)') ' Error reading the input flow solution dimensions.'
               go to 999
            end if

         case (2) ! Write q file header records

            if (formatted) then
               write (lun, '(i5)', iostat=ios) nblocks
            else
               write (lun,         iostat=ios) nblocks
            end if

            if (ios /= 0) then
               write (*, '(a)') ' Error writing the number of flow solution blocks.'
               go to 999
            end if

            if (formatted) then
               write (lun, '(3i5)', iostat=ios) (q(ib)%mi, q(ib)%mj, mq, ib = 1, nblocks)
            else
               write (lun,          iostat=ios) (q(ib)%mi, q(ib)%mj, mq, ib = 1, nblocks)
            end if

            if (ios /= 0) then
               write (*, '(a)') ' Error writing the output flow solution dimensions.'
               go to 999
            end if

         case default

            ios = 999
            write (*, '(a, i6)') ' Bad mode argument passed to q_header_io_2d:', mode

      end select

  999 return

      end subroutine q_header_io_2d

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine xy_allocate (block, ios)
!
!     Allocate a single grid block.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent(inout) :: block
      integer, intent (out)           :: ios

!     Local variables:

      integer :: ni, nj

!     Execution:

      ni = block%ni;  nj = block%nj

      allocate (block%x(ni,nj,1), block%y(ni,nj,1), stat=ios)

      end subroutine xy_allocate

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine q_allocate_2d (block, mq, ios)
!
!     Allocate a single flow solution block.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent(inout) :: block
      integer, intent (in)            :: mq
      integer, intent (out)           :: ios

!     Execution:

      allocate (block%q(mq, block%mi, block%mj, 1), stat=ios)

      end subroutine q_allocate_2d

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine xy_block_io (mode, lun, formatted, npts, x, y, ios)
!
!     Read or write one 2-space PLOT2D-type grid block efficiently (all coordinates packed as one record of length 2 * npts).
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)                   :: mode       ! 1 = read, 2 = write
      integer, intent (in)                   :: lun        ! Logical unit number
      logical, intent (in)                   :: formatted  ! T|F
      integer, intent (in)                   :: npts       ! Number of (x,y) points in the block
      real, dimension (npts), intent (inout) :: x, y       ! Grid block coordinates
      integer, intent (out)                  :: ios        ! 0 = no error

!     Execution:

      select case (mode)
         case (1)
            if (formatted) then
               read (lun, *, iostat=ios) x, y
            else
               read (lun,    iostat=ios) x, y
            end if
         case (2)
           if (formatted) then
               write (lun, '(6es19.11)', iostat=ios) x, y
            else
               write (lun,   iostat=ios) x, y
            end if
         case default
            ios = 999
      end select

      end subroutine xy_block_io

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine q_block_io_2d (mode, lun, formatted, mq, mi, mj, q, ios)
!
!     Read or write one 2-space PLOT2D-type function file.  The different functions correspond to the first index in memory.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)                          :: mode        ! 1 = read, 2 = write
      integer, intent (in)                          :: lun         ! Logical unit number
      logical, intent (in)                          :: formatted   ! T|F
      integer, intent (in)                          :: mq          ! Number of flow variables at each cell or vertex
      integer, intent (in)                          :: mi, mj      ! Number of cells or vertices in each direction
      real, dimension (mq,mi,mj,1), intent (inout)  :: q           ! Flow function values
      integer, intent (out)                         :: ios         ! 0 = no error

!     Local variables:

      integer :: i, j, m

!     Execution:

      select case (mode)
         case (1)
            if (formatted) then
               read (lun, *,  iostat=ios) (((q(m,i,j,1), i = 1, mi), j = 1, mj), m = 1, mq)
            else
               read (lun,     iostat=ios) (((q(m,i,j,1), i = 1, mi), j = 1, mj), m = 1, mq)
            end if
         case (2)
           if (formatted) then
               write (lun, '(6es18.10)', iostat=ios) (((q(m,i,j,1), i = 1, mi), j = 1, mj), m = 1, mq)
            else
               write (lun,               iostat=ios) (((q(m,i,j,1), i = 1, mi), j = 1, mj), m = 1, mq)
            end if
         case default
            ios = 999
      end select

      end subroutine q_block_io_2d
!                                             +--------------------------------------+
!                                             | "iblank" variants are grouped below. |
!                                             +--------------------------------------+

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine xyiq_read (lun_xy, lun_q, formatted, nblocks, num_q, cell_centered, xyiq, iostatus)

!     Read an iblanked PLOT2D-type multiblock 2-space grid and (optionally) a flow-field in PLOT2D-type function file format.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)      :: lun_xy        ! Logical unit for the xy file;
                                                 ! assumed open upon entry; closed here before the return
      integer, intent (in)      :: lun_q         ! Logical unit for the optional "q" file;
                                                 ! assumed open upon entry; closed here before the return;
                                                 ! lun_q < 0 means suppress the reading of a q file;
                                                 ! lun_q = lun_xy is permissible
      logical, intent (in)      :: formatted     ! True for ASCII file(s); false means unformatted reads
      integer, intent (out)     :: nblocks       ! Number of grid blocks found in the grid file
      integer, intent (out)     :: num_q         ! Number of "q" file variables found (if lun_q > 0)
      logical, intent (out)     :: cell_centered ! True if the "q" file is found to be cell-centered
      type (grid_type), pointer :: xyiq(:)       ! Array of grid blocks containing x,y,iblank and optional q data;
                                                 ! allocated here prior to the reading
      integer, intent (out)     :: iostatus      ! 0 means no error detected, else diagnostics are written to
                                                 ! standard output; early return follows
!     Local variables:

      integer :: ib, ios, mi, mj, mq, ni, nj, npts
      logical :: centers, vertices

!     Execution:
!     ----------

!     The input grid file is not optional:
!     ------------------------------------

      call xy_header_io (1, lun_xy, formatted, nblocks, xyiq, ios)

      if (ios /= 0) then
         write (*, '(a)') ' Error allocating the input grid or reading its dimensions.'
         go to 999
      end if

      do ib = 1, nblocks

         call xyi_allocate (xyiq(ib), ios)

         if (ios /= 0) then
            write (*, '(a, i4)') ' Trouble allocating input grid block #', ib
            go to 999
         end if

         npts = xyiq(ib)%ni * xyiq(ib)%nj

         call xyi_block_io (1, lun_xy, formatted, npts, xyiq(ib)%x, xyiq(ib)%y, xyiq(ib)%iblank, ios)

         if (ios /= 0) then
            write (*, '(a, i4)')  ' Trouble reading input grid.  Block #:', ib
            write (*, '(a, 2i5)') ' Dimensions: ', xyiq(ib)%ni, xyiq(ib)%nj
            go to 999
         end if

      end do

      close (lun_xy)

      if (lun_q < 0) go to 999 ! Avoid the indenting


!     Optional input "q" file:
!     ------------------------

      call q_header_io_2d (1, lun_q, formatted, nblocks, mq, xyiq, ios)

      if (ios /= 0) then
         write (*, '(a)') ' Error reading the input flow solution dimensions.'
         go to 999
      end if

      vertices = .true.;  centers = .true.

      do ib = 1, nblocks

         mi = xyiq(ib)%mi;  mj = xyiq(ib)%mj
         ni = xyiq(ib)%ni;  nj = xyiq(ib)%nj

         vertices = mi == ni   .and. mj == nj   .and. vertices
         centers  = mi == ni-1 .and. mj == nj-1 .and. centers

         if (.not. vertices .and. .not. centers) then
            write (*, '(/, a, i4, /, a, 2i5)') ' Grid/flow field dimension mismatch.  Block:', ib, &
               ' Grid: ', ni, nj, ' Flow: ', mi, mj
            ios = 1
            go to 999
         end if

         call q_allocate_2d (xyiq(ib), mq, ios)  ! Put "q" variables in 1st index (poor for I/O but better for interpolations)

         if (ios /= 0) then
            write (*, '(a, i4)') ' Trouble allocating input flow solution block #', ib
            go to 999
         end if

         call q_block_io_2d (1, lun_q, formatted, mq, mi, mj, xyiq(ib)%q, ios)

         if (ios /= 0) then
            write (*, '(a, i4, a, i4)') ' Trouble reading input flow solution.  Block #:', ib, '  I/O status:', ios
            write (*, '(a, 3i5)') ' Dimensions and # functions: ', mi, mj, mq
            go to 999
         end if

      end do

      close (lun_q)

      num_q = mq
      cell_centered = centers

  999 iostatus = ios

      end subroutine xyiq_read

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine xyi_write (lun_xy, formatted, nblocks, xyi, iostatus)

!     Write an iblanked PLOT2D-type multiblock 2-space grid to the indicated unit.
!     The file is assumed to be open upon entry and is closed before the return here.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)      :: lun_xy        ! Logical unit for the xyi file;
                                                 ! assumed open upon entry; closed here before the return
      logical, intent (in)      :: formatted     ! True for ASCII file(s); false means unformatted writes
      integer, intent (inout)   :: nblocks       ! Number of grid blocks (an input, but inout because of the xy_header_io call)
      type (grid_type), pointer :: xyi(:)        ! Array of grid blocks containing x,y,iblank data
      integer, intent (out)     :: iostatus      ! 0 means no error detected, else diagnostics are written to
                                                 ! standard output; early return follows
!     Local variables:

      integer :: ib, ios, npts

!     Execution:
!     ----------

      call xy_header_io (2, lun_xy, formatted, nblocks, xyi, ios)

      if (ios /= 0) then
         write (*, '(a)') ' Error writing the output grid dimensions.'
         go to 999
      end if

      do ib = 1, nblocks

         npts = xyi(ib)%ni * xyi(ib)%nj

         call xyi_block_io (2, lun_xy, formatted, npts, xyi(ib)%x, xyi(ib)%y, xyi(ib)%iblank, ios)

         if (ios /= 0) then
            write (*, '(a, i4)')  ' Trouble writing output grid.  Block #:', ib
            write (*, '(a, 2i5)') ' Dimensions: ', xyi(ib)%ni, xyi(ib)%nj
            go to 999
         end if

      end do

      close (lun_xy)

  999 iostatus = ios

      end subroutine xyi_write

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine xyi_allocate (block, ios)
!
!     Allocate a single grid block.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent(inout) :: block
      integer, intent (out)           :: ios

!     Local variables:

      integer :: ni, nj, nk

!     Execution:

      ni = block%ni;  nj = block%nj

      allocate (block%x(ni,nj,1), block%y(ni,nj,1), block%iblank(ni,nj,1), stat=ios)

      end subroutine xyi_allocate

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine xyi_block_io (mode, lun, formatted, npts, x, y, iblank, ios)
!
!     Read or write one 2-space/iblank PLOT2D-type grid block efficiently (all coordinates and blanking array packed as one record).
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)                      :: mode       ! 1 = read, 2 = write
      integer, intent (in)                      :: lun        ! Logical unit number
      logical, intent (in)                      :: formatted  ! T|F
      integer, intent (in)                      :: npts       ! Number of (x,y) points in the block
      real,    dimension (npts), intent (inout) :: x, y       ! Grid block coordinates
      integer, dimension (npts), intent (inout) :: iblank     ! 0, 1, or -block number;
      integer, intent (out)                     :: ios        ! 0 = no error

!     Execution:

      select case (mode)
         case (1)
            if (formatted) then
               read (lun, *, iostat=ios) x, y, iblank
            else
               read (lun,    iostat=ios) x, y, iblank
            end if
         case (2)
           if (formatted) then
               write (lun, '(6es19.11)', iostat=ios) x, y
               if (ios == 0) write (lun, '(10i6)', iostat=ios) iblank  ! Format allows for up to (-)9999 blocks
            else
               write (lun, iostat=ios) x, y, iblank
            end if
         case default
            ios = 999
      end select

      end subroutine xyi_block_io

   end module xyq_io_module
