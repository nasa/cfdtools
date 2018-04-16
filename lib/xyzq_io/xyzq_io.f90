!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module xyzq_io_module

!  This module packages the I/O for a 3-space PLOT3D-type multiblock grid and/or associated flow field "q" file, with optional
!  blanking control.
!
!  The grid file is separate from the PLOT3D-type function file, which contains a variable number of flow field quantities.
!
!  However, the "grid_type" data structure adopted for LaRC compatibility obliges us to assume reading of a "q" file will follow
!  reading of a related grid file, because reading the grid requires allocation of the grid blocks (here), and these blocks have
!  the requisite "q" fields (which are allocated only if a "q" file is requested).
!
!  On the other hand, writing of an xyz file need not precede writing of a q file.  These functions are separated.
!
!  Three pairs of lower level utilities used by the whole-grid I/O utilities are made public to improve flexibility, as in cases
!  where only one block of a grid needs to be processed at a time.  If blanking is desired, use the appropriate variant of the
!  grid routines, since the PLOT3D convention includes "iblank" in the same record as x, y, z.
!
!  Note that, in order to pass unallocated arrays as arguments, the calling program must declare the x/y/z/iblank/q data type with
!  the "pointer" attribute as opposed to "allocatable".
!
!  Prompting for input and output file names, and opening of the files, are provided here for completeness, but some applications
!  may choose to read file names from a control file and open the files explicitly.
!
!  04/02/2004   D. A. Saunders   Initial implementation, for grid-morphing use at NASA LaRC (and other CFD utilities).
!  01/15/2010    "    "    "     After years of regretting the decision to avoid possible machine dependencies with advance='no',
!                                the file prompting now leaves the cursor at the end of the prompt.
!  03/03/2010    "    "    "     Added "iblank" variants of the grid utilities.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, Moffett Field, CA
! 
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure ! Or equivalent to define "grid_type"

   implicit none

   private

   public :: file_prompt    ! Prompts for a file name and open it

   public :: xyzq_read      ! Reads a multiblock grid and (optionally) an associated flow solution (PLOT3D-style)
   public :: xyz_write      ! Writes a multiblock grid (PLOT3D-style)
   public :: q_write        ! Writes a flow solution   (  "     "   )

   public :: xyz_header_io  ! Reads or writes grid header records
   public :: q_header_io    ! Reads or writes flow header records
   public :: xyz_allocate   ! Allocates one grid block (deallocates can be done explicitly)
   public :: q_allocate     ! Allocates one flow block (  "     "     "     "     "     " )
   public :: xyz_block_io   ! Reads or writes one grid block
   public :: q_block_io     ! Reads or writes one flow block

   public :: xyziq_read     ! "iblank" variants of the xyz utilities above
   public :: xyzi_write
   public :: xyzi_allocate
   public :: xyzi_block_io

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine file_prompt (lun_file, lun_journal, file_type, status, format_prompt, formatted, iostatus)

!     Prompt for a file of the indicated type and status.  Standard input and output are used (logical unit *) for the prompting.
!     Responses are echoed to a journal file unless lun_journal <= 0.  The indicated file is opened here.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)    :: lun_file      ! Logical unit for the file; opened here

      integer, intent (in)    :: lun_journal   ! Logical unit for journal file, if > 0

      character, intent (in)  :: file_type*(*) ! Description of the file type, lower case; e.g., 'input grid' or 'output q'

      character, intent (in)  :: status*(*)    ! File status: 'old' for input; 'new' or 'unknown' for output

      logical, intent (in)    :: format_prompt ! True means prompt for and output the "formatted" argument, else use it as input

      logical, intent (inout) :: formatted     ! True if an ASCII file is indicated; false means unformatted

      integer, intent (out)   :: iostatus      ! 0 means no error detected;
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
      if (lun_journal > 0) write (lun_journal, '(a)') filename(1:len_trim(filename))

      if (format_prompt) then
         write (*, '(/, a)', advance='no') ' Is the file formatted (y) or unformatted (n)? '
         read  (*, *) answer
         formatted = answer == 'y' .or. answer == 'Y'
         if (lun_journal > 0) write (lun_journal, '(a)') answer
      end if

      i1 = 1;  if (formatted) i1 = 3

      do ntries = 1, 3
         open (lun_file, file=filename, form=format(i1:11), status=status, iostat=iostatus)
         if (iostatus == 0) exit
         if (ntries < 3) then
            write (*, '(a)') ' Bad file name? Try again:'
            read (*, '(a)') filename
            if (lun_journal > 0) write (lun_journal, '(a)') filename(1:len_trim(filename))
         else
            write (*, '(a)') ' Unable to open the file.'
         end if
      end do

      end subroutine file_prompt

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine xyzq_read (lun_xyz, lun_q, formatted, nblocks, num_q, cell_centered, xyzq, iostatus)

!     Read a PLOT3D-type multiblock 3-space grid and (optionally) a flow-field in PLOT3D-type function file format.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)      :: lun_xyz       ! Logical unit for the xyz file;
                                                 ! assumed open upon entry; closed here before the return

      integer, intent (in)      :: lun_q         ! Logical unit for the optional "q" file;
                                                 ! assumed open upon entry; closed here before the return;
                                                 ! lun_q < 0 means suppress the reading of a q file;
                                                 ! lun_q = lun_xyz is permissible

      logical, intent (in)      :: formatted     ! True for ASCII file(s); false means unformatted reads

      integer, intent (out)     :: nblocks       ! Number of grid blocks found in the grid file

      integer, intent (out)     :: num_q         ! Number of "q" file variables found (if lun_q > 0)

      logical, intent (out)     :: cell_centered ! True if the "q" file is found to be cell-centered

      type (grid_type), pointer :: xyzq(:)       ! Array of grid blocks containing x,y,z and optional q data;
                                                 ! allocated here prior to the reading

      integer, intent (out)     :: iostatus      ! 0 means no error detected, else diagnostics are written to
                                                 ! standard output; early return follows
!     Local variables:

      integer   :: ib, ios, mi, mj, mk, mq, ni, nj, nk, npts

      logical   :: centers, vertices

!     Execution:
!     ----------

!     The input grid file is not optional:
!     ------------------------------------

      call xyz_header_io (1, lun_xyz, formatted, nblocks, xyzq, ios)

      if (ios /= 0) then
         write (*, '(a)') ' Error allocating the input grid or reading its dimensions.'
         go to 999
      end if

      do ib = 1, nblocks

         call xyz_allocate (xyzq(ib), ios)

         if (ios /= 0) then
            write (*, '(a, i4)') ' Trouble allocating input grid block #', ib
            go to 999
         end if

         npts = xyzq(ib)%ni * xyzq(ib)%nj * xyzq(ib)%nk

         call xyz_block_io (1, lun_xyz, formatted, npts, xyzq(ib)%x, xyzq(ib)%y, xyzq(ib)%z, ios)

         if (ios /= 0) then
            write (*, '(a, i4)')  ' Trouble reading input grid.  Block #:', ib
            write (*, '(a, 4i5)') ' Dimensions: ', xyzq(ib)%ni, xyzq(ib)%nj, xyzq(ib)%nk
            go to 999
         end if

      end do

      close (lun_xyz)

      if (lun_q < 0) go to 999 ! Avoid the indenting


!     Optional input "q" file:
!     ------------------------

      call q_header_io (1, lun_q, formatted, nblocks, mq, xyzq, ios)

      if (ios /= 0) then
         write (*, '(a)') ' Error reading the input flow solution dimensions.'
         go to 999
      end if

      vertices = .true.;  centers = .true.

      do ib = 1, nblocks

         mi = xyzq(ib)%mi;  mj = xyzq(ib)%mj;  mk = xyzq(ib)%mk
         ni = xyzq(ib)%ni;  nj = xyzq(ib)%nj;  nk = xyzq(ib)%nk

         vertices = mi == ni   .and. mj == nj   .and. mk == nk   .and. vertices
         centers  = mi == ni-1 .and. mj == nj-1 .and. mk == nk-1 .and. centers

         if (.not. vertices .and. .not. centers) then
            write (*, '(/, a, i4, /, a, 3i5)') ' Grid/flow field dimension mismatch.  Block:', ib, &
               ' Grid: ', ni, nj, nk, ' Flow: ', mi, mj, mk
            ios = 1
            go to 999
         end if

         call q_allocate (xyzq(ib), mq, ios)  ! Put "q" variables in 1st index (poor for I/O but better for interpolations)

         if (ios /= 0) then
            write (*, '(a, i4)') ' Trouble allocating input flow solution block #', ib
            go to 999
         end if

         call q_block_io (1, lun_q, formatted, mq, mi, mj, mk, xyzq(ib)%q, ios)

         if (ios /= 0) then
            write (*, '(a, i4, a, i4)') ' Trouble reading input flow solution.  Block #:', ib, '  I/O status:', ios
            write (*, '(a, 4i5)') ' Dimensions and # functions: ', mi, mj, mk, mq
            go to 999
         end if

      end do

      close (lun_q)

      num_q = mq
      cell_centered = centers

  999 iostatus = ios

      end subroutine xyzq_read

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine xyz_write (lun_xyz, formatted, nblocks, xyz, iostatus)

!     Write a PLOT3D-type multiblock 3-space grid to the indicated unit.
!     The file is assumed to be open upon entry and is closed before the return here.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)      :: lun_xyz       ! Logical unit for the xyz file;
                                                 ! assumed open upon entry; closed here before the return

      logical, intent (in)      :: formatted     ! True for ASCII file(s); false means unformatted writes

      integer, intent (inout)   :: nblocks       ! Number of grid blocks (an input, but inout because of the xyz_header_io call)

      type (grid_type), pointer :: xyz(:)        ! Array of grid blocks containing x,y,z data

      integer, intent (out)     :: iostatus      ! 0 means no error detected, else diagnostics are written to
                                                 ! standard output; early return follows
!     Local variables:

      integer   :: ib, ios, npts

!     Execution:
!     ----------

      call xyz_header_io (2, lun_xyz, formatted, nblocks, xyz, ios)

      if (ios /= 0) then
         write (*, '(a)') ' Error writing the output grid dimensions.'
         go to 999
      end if

      do ib = 1, nblocks

         npts = xyz(ib)%ni * xyz(ib)%nj * xyz(ib)%nk

         call xyz_block_io (2, lun_xyz, formatted, npts, xyz(ib)%x, xyz(ib)%y, xyz(ib)%z, ios)

         if (ios /= 0) then
            write (*, '(a, i4)')  ' Trouble writing output grid.  Block #:', ib
            write (*, '(a, 4i5)') ' Dimensions: ', xyz(ib)%ni, xyz(ib)%nj, xyz(ib)%nk
            go to 999
         end if

      end do

      close (lun_xyz)

  999 iostatus = ios

      end subroutine xyz_write

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine q_write (lun_q, formatted, nblocks, num_q, q, iostatus)

!     Write a multiblock flow-field in PLOT3D-type function file format.
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

      integer   :: ib, ios

!     Execution:
!     ----------

      call q_header_io (2, lun_q, formatted, nblocks, num_q, q, ios)

      if (ios /= 0) then
         write (*, '(a)') ' Error writing the output flow solution dimensions.'
         go to 999
      end if

      do ib = 1, nblocks

         call q_block_io (2, lun_q, formatted, num_q, q(ib)%mi, q(ib)%mj, q(ib)%mk, q(ib)%q, ios)

         if (ios /= 0) then
            write (*, '(a, i4, a, i4)') ' Trouble writing output flow solution.  Block #:', ib, '  I/O status:', ios
            write (*, '(a, 4i5)') ' Dimensions and # functions: ', q(ib)%mi, q(ib)%mj, q(ib)%mk, num_q
            go to 999
         end if

      end do

      close (lun_q)

  999 iostatus = ios

      end subroutine q_write

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine xyz_header_io (mode, lun, formatted, nblocks, xyz, ios)
!
!     Read or write header records for a 3-space PLOT3D-type grid.  Blocks are allocated here before reading their dimension fields.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)       :: mode       ! 1 = read, 2 = write
      integer, intent (in)       :: lun        ! Logical unit number
      logical, intent (in)       :: formatted  ! T|F
      integer, intent (inout)    :: nblocks    ! Number of blocks: output for a read; input for a write
      type (grid_type), pointer  :: xyz(:)     ! Array of grid blocks: allocated here for a read
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

            allocate (xyz(nblocks), stat=ios)

            if (ios /= 0) then
               write (*, '(a)') ' Error allocating the grid blocks.'
               go to 999
            end if

            if (formatted) then
               read (lun, *, iostat=ios) (xyz(ib)%ni, xyz(ib)%nj, xyz(ib)%nk, ib = 1, nblocks)
            else
               read (lun,    iostat=ios) (xyz(ib)%ni, xyz(ib)%nj, xyz(ib)%nk, ib = 1, nblocks)
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
               write (lun, '(3i5)', iostat=ios) (xyz(ib)%ni, xyz(ib)%nj, xyz(ib)%nk, ib = 1, nblocks)
            else
               write (lun,          iostat=ios) (xyz(ib)%ni, xyz(ib)%nj, xyz(ib)%nk, ib = 1, nblocks)
            end if

            if (ios /= 0) then
               write (*, '(a)') ' Error writing the output grid dimensions.'
               go to 999
            end if

         case default

            ios = 999
            write (*, '(a, i6)') ' Bad mode argument passed to xyz_header_io:', mode

      end select

  999 return

      end subroutine xyz_header_io

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine q_header_io (mode, lun, formatted, nblocks, mq, q, ios)
!
!     Read or write header records for a 3-space PLOT3D-type function file.  Blocks are assumed to have been allocated already
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
               read (lun, *, iostat=ios) (q(ib)%mi, q(ib)%mj, q(ib)%mk, mq, ib = 1, nb)
            else
               read (lun,    iostat=ios) (q(ib)%mi, q(ib)%mj, q(ib)%mk, mq, ib = 1, nb)
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
               write (lun, '(4i5)', iostat=ios) (q(ib)%mi, q(ib)%mj, q(ib)%mk, mq, ib = 1, nblocks)
            else
               write (lun,          iostat=ios) (q(ib)%mi, q(ib)%mj, q(ib)%mk, mq, ib = 1, nblocks)
            end if

            if (ios /= 0) then
               write (*, '(a)') ' Error writing the output flow solution dimensions.'
               go to 999
            end if

         case default

            ios = 999
            write (*, '(a, i6)') ' Bad mode argument passed to q_header_io:', mode

      end select

  999 return

      end subroutine q_header_io

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine xyz_allocate (block, ios)
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

      ni = block%ni;  nj = block%nj;  nk = block%nk

      allocate (block%x(ni,nj,nk), block%y(ni,nj,nk), block%z(ni,nj,nk), stat=ios)

      end subroutine xyz_allocate

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine q_allocate (block, mq, ios)
!
!     Allocate a single flow solution block.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent(inout) :: block
      integer, intent (in)            :: mq
      integer, intent (out)           :: ios

!     Execution:

      allocate (block%q(mq, block%mi, block%mj, block%mk), stat=ios)

      end subroutine q_allocate

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine xyz_block_io (mode, lun, formatted, npts, x, y, z, ios)
!
!     Read or write one 3-space PLOT3D-type grid block efficiently (all coordinates packed as one record of length 3 * npts).
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)                   :: mode       ! 1 = read, 2 = write
      integer, intent (in)                   :: lun        ! Logical unit number
      logical, intent (in)                   :: formatted  ! T|F
      integer, intent (in)                   :: npts       ! Number of (x,y,z) points in the block
      real, dimension (npts), intent (inout) :: x, y, z    ! Grid block coordinates
      integer, intent (out)                  :: ios        ! 0 = no error

!     Execution:

      select case (mode)
         case (1)
            if (formatted) then
               read (lun, *, iostat=ios) x, y, z
            else
               read (lun, iostat=ios) x, y, z
            end if
         case (2)
           if (formatted) then
               write (lun, '(1p, 6e19.11)', iostat=ios) x, y, z
            else
               write (lun, iostat=ios) x, y, z
            end if
         case default
            ios = 999
      end select

      end subroutine xyz_block_io

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine q_block_io (mode, lun, formatted, mq, mi, mj, mk, q, ios)
!
!     Read or write one 3-space PLOT3D-type function file.  The different functions correspond to the first index in memory.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)                          :: mode        ! 1 = read, 2 = write
      integer, intent (in)                          :: lun         ! Logical unit number
      logical, intent (in)                          :: formatted   ! T|F
      integer, intent (in)                          :: mq          ! Number of flow variables at each cell or vertex
      integer, intent (in)                          :: mi, mj, mk  ! Number of cells or vertices in each direction
      real, dimension (mq,mi,mj,mk), intent (inout) :: q           ! Flow function values
      integer, intent (out)                         :: ios         ! 0 = no error

!     Local variables:

      integer :: i, j, k, m

!     Execution:

      select case (mode)
         case (1)
            if (formatted) then
               read (lun, *,  iostat=ios) ((((q(m,i,j,k), i = 1, mi), j = 1, mj), k = 1, mk), m = 1, mq)
            else
               read (lun,     iostat=ios) ((((q(m,i,j,k), i = 1, mi), j = 1, mj), k = 1, mk), m = 1, mq)
            end if
         case (2)
           if (formatted) then
               write (lun, '(1p, 6e18.10)', iostat=ios) ((((q(m,i,j,k), i = 1, mi), j = 1, mj), k = 1, mk), m = 1, mq)
            else
               write (lun,                  iostat=ios) ((((q(m,i,j,k), i = 1, mi), j = 1, mj), k = 1, mk), m = 1, mq)
            end if
         case default
            ios = 999
      end select

      end subroutine q_block_io
!                                             +--------------------------------------+
!                                             | "iblank" variants are grouped below. |
!                                             +--------------------------------------+

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine xyziq_read (lun_xyz, lun_q, formatted, nblocks, num_q, cell_centered, xyziq, iostatus)

!     Read an iblanked PLOT3D-type multiblock 3-space grid and (optionally) a flow-field in PLOT3D-type function file format.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)      :: lun_xyz       ! Logical unit for the xyz file;
                                                 ! assumed open upon entry; closed here before the return

      integer, intent (in)      :: lun_q         ! Logical unit for the optional "q" file;
                                                 ! assumed open upon entry; closed here before the return;
                                                 ! lun_q < 0 means suppress the reading of a q file;
                                                 ! lun_q = lun_xyz is permissible

      logical, intent (in)      :: formatted     ! True for ASCII file(s); false means unformatted reads

      integer, intent (out)     :: nblocks       ! Number of grid blocks found in the grid file

      integer, intent (out)     :: num_q         ! Number of "q" file variables found (if lun_q > 0)

      logical, intent (out)     :: cell_centered ! True if the "q" file is found to be cell-centered

      type (grid_type), pointer :: xyziq(:)      ! Array of grid blocks containing x,y,z,iblank and optional q data;
                                                 ! allocated here prior to the reading

      integer, intent (out)     :: iostatus      ! 0 means no error detected, else diagnostics are written to
                                                 ! standard output; early return follows
!     Local variables:

      integer   :: ib, ios, mi, mj, mk, mq, ni, nj, nk, npts

      logical   :: centers, vertices

!     Execution:
!     ----------

!     The input grid file is not optional:
!     ------------------------------------

      call xyz_header_io (1, lun_xyz, formatted, nblocks, xyziq, ios)

      if (ios /= 0) then
         write (*, '(a)') ' Error allocating the input grid or reading its dimensions.'
         go to 999
      end if

      do ib = 1, nblocks

         call xyzi_allocate (xyziq(ib), ios)

         if (ios /= 0) then
            write (*, '(a, i4)') ' Trouble allocating input grid block #', ib
            go to 999
         end if

         npts = xyziq(ib)%ni * xyziq(ib)%nj * xyziq(ib)%nk

         call xyzi_block_io (1, lun_xyz, formatted, npts, xyziq(ib)%x, xyziq(ib)%y, xyziq(ib)%z, xyziq(ib)%iblank, ios)

         if (ios /= 0) then
            write (*, '(a, i4)')  ' Trouble reading input grid.  Block #:', ib
            write (*, '(a, 4i5)') ' Dimensions: ', xyziq(ib)%ni, xyziq(ib)%nj, xyziq(ib)%nk
            go to 999
         end if

      end do

      close (lun_xyz)

      if (lun_q < 0) go to 999 ! Avoid the indenting


!     Optional input "q" file:
!     ------------------------

      call q_header_io (1, lun_q, formatted, nblocks, mq, xyziq, ios)

      if (ios /= 0) then
         write (*, '(a)') ' Error reading the input flow solution dimensions.'
         go to 999
      end if

      vertices = .true.;  centers = .true.

      do ib = 1, nblocks

         mi = xyziq(ib)%mi;  mj = xyziq(ib)%mj;  mk = xyziq(ib)%mk
         ni = xyziq(ib)%ni;  nj = xyziq(ib)%nj;  nk = xyziq(ib)%nk

         vertices = mi == ni   .and. mj == nj   .and. mk == nk   .and. vertices
         centers  = mi == ni-1 .and. mj == nj-1 .and. mk == nk-1 .and. centers

         if (.not. vertices .and. .not. centers) then
            write (*, '(/, a, i4, /, a, 3i5)') ' Grid/flow field dimension mismatch.  Block:', ib, &
               ' Grid: ', ni, nj, nk, ' Flow: ', mi, mj, mk
            ios = 1
            go to 999
         end if

         call q_allocate (xyziq(ib), mq, ios)  ! Put "q" variables in 1st index (poor for I/O but better for interpolations)

         if (ios /= 0) then
            write (*, '(a, i4)') ' Trouble allocating input flow solution block #', ib
            go to 999
         end if

         call q_block_io (1, lun_q, formatted, mq, mi, mj, mk, xyziq(ib)%q, ios)

         if (ios /= 0) then
            write (*, '(a, i4, a, i4)') ' Trouble reading input flow solution.  Block #:', ib, '  I/O status:', ios
            write (*, '(a, 4i5)') ' Dimensions and # functions: ', mi, mj, mk, mq
            go to 999
         end if

      end do

      close (lun_q)

      num_q = mq
      cell_centered = centers

  999 iostatus = ios

      end subroutine xyziq_read

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine xyzi_write (lun_xyz, formatted, nblocks, xyzi, iostatus)

!     Write an iblanked PLOT3D-type multiblock 3-space grid to the indicated unit.
!     The file is assumed to be open upon entry and is closed before the return here.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)      :: lun_xyz       ! Logical unit for the xyzi file;
                                                 ! assumed open upon entry; closed here before the return

      logical, intent (in)      :: formatted     ! True for ASCII file(s); false means unformatted writes

      integer, intent (inout)   :: nblocks       ! Number of grid blocks (an input, but inout because of the xyz_header_io call)

      type (grid_type), pointer :: xyzi(:)       ! Array of grid blocks containing x,y,z,iblank data

      integer, intent (out)     :: iostatus      ! 0 means no error detected, else diagnostics are written to
                                                 ! standard output; early return follows
!     Local variables:

      integer   :: ib, ios, npts

!     Execution:
!     ----------

      call xyz_header_io (2, lun_xyz, formatted, nblocks, xyzi, ios)

      if (ios /= 0) then
         write (*, '(a)') ' Error writing the output grid dimensions.'
         go to 999
      end if

      do ib = 1, nblocks

         npts = xyzi(ib)%ni * xyzi(ib)%nj * xyzi(ib)%nk

         call xyzi_block_io (2, lun_xyz, formatted, npts, xyzi(ib)%x, xyzi(ib)%y, xyzi(ib)%z, xyzi(ib)%iblank, ios)

         if (ios /= 0) then
            write (*, '(a, i4)')  ' Trouble writing output grid.  Block #:', ib
            write (*, '(a, 4i5)') ' Dimensions: ', xyzi(ib)%ni, xyzi(ib)%nj, xyzi(ib)%nk
            go to 999
         end if

      end do

      close (lun_xyz)

  999 iostatus = ios

      end subroutine xyzi_write

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine xyzi_allocate (block, ios)
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

      ni = block%ni;  nj = block%nj;  nk = block%nk

      allocate (block%x(ni,nj,nk), block%y(ni,nj,nk), block%z(ni,nj,nk), block%iblank(ni,nj,nk), stat=ios)

      end subroutine xyzi_allocate

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine xyzi_block_io (mode, lun, formatted, npts, x, y, z, iblank, ios)
!
!     Read or write one 3-space/iblank PLOT3D-type grid block efficiently (all coordinates and blanking array packed as one record).
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)                      :: mode       ! 1 = read, 2 = write
      integer, intent (in)                      :: lun        ! Logical unit number
      logical, intent (in)                      :: formatted  ! T|F
      integer, intent (in)                      :: npts       ! Number of (x,y,z) points in the block
      real,    dimension (npts), intent (inout) :: x, y, z    ! Grid block coordinates
      integer, dimension (npts), intent (inout) :: iblank     ! 0, 1, or -block number;
      integer, intent (out)                     :: ios        ! 0 = no error

!     Execution:

      select case (mode)
         case (1)
            if (formatted) then
               read (lun, *, iostat=ios) x, y, z, iblank
            else
               read (lun, iostat=ios) x, y, z, iblank
            end if
         case (2)
           if (formatted) then
               write (lun, '(1p, 6e19.11)', iostat=ios) x, y, z
               if (ios == 0) write (lun, '(10i6)', iostat=ios) iblank  ! Format allows for up to (-)9999 blocks
            else
               write (lun, iostat=ios) x, y, z, iblank
            end if
         case default
            ios = 999
      end select

      end subroutine xyzi_block_io

   end module xyzq_io_module
