!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module grid_block_structure

!  Module for a grid block or zone with associated flow field variables and a zone name (Tecplot_io_module version).
!  This version is a superset of that originally used with the xyzq_io package, and is intended to be compatible with that package.
!  If both packages are used by an application, this version of the module is the one that should be compiled.

   type grid_type
      real,    dimension(:,:,:),    pointer :: x,y,z  ! grid coordinates
      real,    dimension(:,:,:,:),  pointer :: q      ! dependent variables
      integer, dimension(:,:,:),    pointer :: iblank ! 0 => suppress the point
      real    :: xmin, xmax         ! Data range for the grid block in the x direction
      real    :: ymin, ymax         ! ... and the y direction ...
      real    :: zmin, zmax         ! ... and the z direction
      integer :: ni                 ! number of nodes in i direction
      integer :: nj                 ! number of nodes in j direction
      integer :: nk                 ! number of nodes in k direction
      integer :: mi                 ! number of dependent variables in i direction
      integer :: mj                 ! number of dependent variables in j direction
      integer :: mk                 ! number of dependent variables in k direction
      character :: zone_title * 32  ! Zone or block title; may contain embedded blanks, but its length is limited
   end type grid_type

   end module grid_block_structure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module Tecplot_io_module

!  This module packages the I/O for Tecplot binary (and ASCII) files representing 3-space structured multiblock grids and associated
!  function data.  Both BLOCK and POINT data packing formats are supported.  The function values are assumed to be vertex-centered,
!  with (x,y,z) grid coordinates as the first 3 variables. 
!
!  Several lower-level utilities used by the whole-grid I/O utilities are made public to permit processing of one block at a time
!  (as opposed to reading or writing the whole file in one call).
!
!  Note that, in order to pass unallocated arrays as arguments, the calling program must declare the x/y/z/q data type with the
!  "pointer" attribute as opposed to "allocatable".
!
!  04/02/04  David Saunders  Initial implementation of XYZQ_IO package for PLOT3D files.
!  07/29/04    "       "     Adaptation for Tecplot files, but reading of Tecplot binaries is unsolved.
!  09/24/04    "       "     Automated handling of -r4 and -r8 compilation via new isdouble internal variable.
!  05/04/05    "       "     Variable names are now returned from a formatted read (originally overlooked).
!  05/06/05    "       "     Verbose mode can be activated by entering ios = 1 to Tec_header_read.
!  05/13/05    "       "     The datapacking (1|2 = POINT|BLOCK) should be output for a header read, not input.
!  05/14/05    "       "     Binary outputs cannot handle embedded commas in variable names.  Substitute colons.
!  06/16/05    "       "     Don't assume 'ZONE T' starts in column 1 in count_blocks; handle trailing control Ms.
!  11/30/05    "       "     Added deallocate_blocks.
!  12/04/05    "       "     Fixed Tec_block_read where num_q == 3 should have been num_q == 0.
!  12/05/05    "       "     If num_q = 0 in Tec_block_allocate, allocate %q(1,1,1,1) to avoid undefined arguments elsewhere.
!  12/06/05    "       "     Extended to handle simpler input of the type written by flow solver DPLR's post-processor, as for 2-D.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
! 
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure ! Or equivalent to define "grid_type"

   implicit none

   private

!  Constants used by the package:

   integer,   parameter :: len_string = 64        ! Room and then some for a ' I=iii, J=jjj, K=kkk' string
   integer,   parameter :: max_blocks = 1000      ! Limit on number of blocks in a file being read
   integer,   parameter :: max_length = 1000      ! Limit on length of string into which variable names are packed
   integer,   parameter :: name_limit = 32        ! Limit on the length of variable names; increase it if necessary

   character, parameter :: blank * 1  = ' ',      &
                           null * 1   = char (0), &
                           quotes * 1 = '"'

!  Awkward variables used by the package:

   integer   :: isdouble                               ! 0 = single precision data (-r4); 1 = double precision data (-r8)
   real      :: eps                                    ! epsilon (eps) allows isdouble to be assigned for binary reads & writes
   logical   :: verbose                                ! Pass ios = 1 to Tec_header_read to activate printing of header/zone info.
   character :: dimensions (max_blocks) * (len_string) ! For temporary storage of block dimensions
   character :: names_buffer * (max_length)            ! For reading an unknown number of variables names from one or more lines
   character :: packed_names * (max_length)            ! For binary output; reused for some of the reading
                                                       ! Fortran 90 doesn't support dynamically variable string lengths
!  Parsing utility needed from another package:

   external     scan2, scan4

!  Utilities provided for applications:

   public :: Tecplot_read       ! Reads  an entire structured Tecplot file
   public :: Tecplot_write      ! Writes an entire structured Tecplot file

   public :: Tec_header_read    ! Reads  Tecplot file header records
   public :: Tec_header_write   ! Writes Tecplot file header records
   public :: Tec_block_allocate ! Allocate the array fields for one block or zone of a Tecplot file
   public :: Tec_block_read     ! Reads  one block of a structured Tecplot file
   public :: Tec_block_write    ! Writes one block of a structured Tecplot file

   public :: deallocate_blocks  ! Deallocates indicated blocks of a grid

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Tecplot_read (lun, filename, formatted, datapacking, title, nblocks, num_q, names, xyzq, ios)

!     Read a Tecplot file of 3-space structured data, binary or ASCII, with BLOCK or POINT data packing.
!     The first 3 variables are returned as the x, y, z fields of each block in xyzq(:).
!     Remaining variables become the "q" field packed in "point" order (n-tuples, where n = num_q = # variables found minus 3).
!     The input file title is also returned, along with the variable names (up to name_limit characters each).
!     Titles for the blocks or zones are returned in the appropriate field of each element of the xyzq array.
!     The file is opened and closed here.
!
!     05/13/05   DAS   The "datapacking" argument has been changed from intent (in) to intent (out).
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer,   intent (in)    :: lun           ! Logical unit for the Tecplot file;
                                                 ! opened here upon entry; closed here before the return
      character, intent (in)    :: filename*(*)  ! Name of the Tecplot file; a terminating null character is appended here for the
                                                 ! Tecplot binary case
      logical,   intent (in)    :: formatted     ! True for reading ASCII data; false for reading a Tecplot binary file

      integer,   intent (out)   :: datapacking   ! 1 means DATAPACKING=POINT; 2 means DATAPACKING=BLOCK

      character, intent (out)   :: title * (*)   ! Title found as the first line of the Tecplot file

      integer,   intent (out)   :: nblocks       ! Number of blocks or zones found

      integer,   intent (out)   :: num_q         ! Number of variables found minus 3; num_q = 0 for a pure 3-space grid file

      character,        pointer :: names (:) * (name_limit) ! Variable names, left-justified; array size will be 3 + num_q

      type (grid_type), pointer :: xyzq(:)       ! Array of structured blocks containing x,y,z and optional q data;
                                                 ! allocated here prior to the reading

      integer,   intent (inout) :: ios           ! 1 on input activates writing of header/zone info. during the read;
                                                 ! 0 on output means no error detected, else diagnostics are written to
                                                 ! standard output; early return follows
!     Local variables:

      integer    :: ib, mi, mj, mk, mq, ni, nj, nk

!     Execution:
!     ----------

!     Open the file, read the header records, determine the number of blocks or zones, allocate an array of grid structures,
!     and set up the block dimensions:

      call Tec_header_read (lun, filename, formatted, datapacking, title, nblocks, num_q, names, xyzq, ios)

      if (ios /= 0) go to 999

      do ib = 1, nblocks

         call Tec_block_allocate (xyzq(ib), num_q, ios)

         if (ios /= 0) then
            write (*, '(a, i4)') ' Tecplot_read: Trouble allocating x,y,z,q arrays for block #', ib
            write (*, '(2a)')    ' File name: ', trim (filename)
            go to 999
         end if

         call Tec_block_read (lun, formatted, datapacking, xyzq(ib)%zone_title, xyzq(ib)%ni, xyzq(ib)%nj, xyzq(ib)%nk, num_q, &
                              xyzq(ib)%x, xyzq(ib)%y, xyzq(ib)%z, xyzq(ib)%q, ios)

         if (ios /= 0) then
            write (*, '(a, i4)')  ' Tecplot_read:  Trouble reading block #', ib
            write (*, '(2a)')     ' File name: ', trim (filename)
            write (*, '(a, 4i5)') ' Block dimensions & num_q: ', xyzq(ib)%ni, xyzq(ib)%nj, xyzq(ib)%nk, num_q
            go to 999
         end if

      end do

  999 close (lun)

      end subroutine Tecplot_read

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Tecplot_write (lun, filename, formatted, datapacking, title, nblocks, num_q, names, xyzq, ios)

!     Write a Tecplot file of 3-space structured data, binary or ASCII, with BLOCK or POINT data packing.
!     The first 3 variables written for each block are the x, y, z fields of xyzq.
!     Any "q" variables are appended, giving a total of num_q + 3 variables.
!     Titles for the blocks or zones are written from the appropriate field of each element of the xyzq array.
!     The file is opened and closed here.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer,          intent (in)    :: lun          ! Logical unit for the Tecplot file;
                                                       ! assumed open upon entry; closed here before the return
      character,        intent (in)    :: filename*(*) ! Name of the Tecplot file; a terminating null character is appended here
                                                       ! for the Tecplot binary case
      logical,          intent (in)    :: formatted    ! True for writing ASCII data; false for writing a Tecplot binary file

      integer,          intent (in)    :: datapacking  ! 1 means DATAPACKING=POINT; 2 means DATAPACKING=BLOCK (for all zones)

      character,        intent (inout) :: title * (*)  ! Title written as the first line of the Tecplot file; really intent (in)

      integer,          intent (inout) :: nblocks      ! Number of blocks or zones; really intent (in)

      integer,          intent (inout) :: num_q        ! Number of variables other than x, y, z; num_q = 0 for a pure 3-space grid;
                                                       ! really intent (in) [but the compiler demands inout]
      character, pointer :: names (:) * (name_limit)   ! Variable names, left-justified, dimension (3 + num_q)

      type (grid_type), pointer        :: xyzq(:)      ! Array of structured blocks containing x,y,z and optional q data

      integer,          intent (out)   :: ios          ! 0 means no error detected, else diagnostics are written to standard output;
                                                       ! early return follows
!     Local variables:

      integer :: ib

!     Execution:
!     ----------

!     Open the file and write the header records:

      call Tec_header_write (lun, filename, formatted, datapacking, title, num_q, names, ios)

      if (ios /= 0) go to 999

      do ib = 1, nblocks

         call Tec_block_write (lun, formatted, datapacking, xyzq(ib)%zone_title, xyzq(ib)%ni, xyzq(ib)%nj, xyzq(ib)%nk, num_q, &
                               xyzq(ib)%x, xyzq(ib)%y, xyzq(ib)%z, xyzq(ib)%q, ios)

         if (ios /= 0) then
            write (*, '(a, i4)')  ' Tecplot_write:  Trouble writing block #', ib
            write (*, '(2a)')     ' File name: ', trim (filename)
            write (*, '(a, 4i5)') ' Block dimensions & num_q: ', xyzq(ib)%ni, xyzq(ib)%nj, xyzq(ib)%nk, num_q
            go to 999
         end if

      end do

  999 close (lun)

      end subroutine Tecplot_write

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine Tec_header_read (lun, filename, formatted, datapacking, title, nblocks, num_q, names, xyzq, ios)
!
!     Open a Tecplot data file and read the header records for 3-space structured Tecplot data.
!     Before returning, the file is rewound and advanced ready for reading the first block.
!
!     The simplest type of header has no title and a single list of variable names like this from DPLR's Postflow:
!
!      variables=x,y,z,rho,p,T,Tv,n2,o2,no,no+,n2+,o2+,n,o,n+,o+,e,u,v,w,h,M,mu,kap,Vel,ent
!      zone t="flow3d" F=point, i=  97 j=  33 k= 129   <start of first block>
!
!     The format originally handled looks like this:
!
!     TITLE     = "CFD11L_DPLR_TCG_PhyMod1"
!     VARIABLES = "x, m"
!     "y, m"
!     "z, m"
!     "rho, kg/m^3"
!     "p, Pa"
!     "T, K"
!     "c_N_2"
!     "c_O_2"
!     "c_N_O"
!     :::::::
!     ZONE T="G1"  <start of first block>
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer,   intent (in)    :: lun          ! Logical unit number
      character, intent (in)    :: filename*(*) ! Name of the Tecplot file; a terminating null character is appended here for the
                                                ! Tecplot binary case
      logical,   intent (in)    :: formatted    ! T|F
      integer,   intent (out)   :: datapacking  ! 1 means DATAPACKING=POINT; 2 means DATAPACKING=BLOCK (for all zones);
                                                ! output here for Tec_block_read purposes (all zones assumed the same order)
      character, intent (out)   :: title*(*)    ! Title found on the first line of the Tecplot file (blank if absent)
      integer,   intent (out)   :: nblocks      ! Number of blocks found
      integer,   intent (out)   :: num_q        ! Number of variables found, not counting x,y,z
      character,        pointer :: names (:) * (name_limit) ! Variable names, left-justified; allocated here (3 + num_q)
      type (grid_type), pointer :: xyzq(:)      ! Array of derived data types for structured grid blocks and optional function data,
                                                ! allocated here
      integer,   intent (inout) :: ios          ! 1 on input means write header/zone info. during reading (verbose mode);
                                                ! 0 on output means no error

!     Local constants:

      character, parameter :: routine * 18 = ' Tec_header_read: '

!     Local variables:

      integer :: first, ib, iend, istart, last, mark, nvar ! Shared by count_blocks and pack_names
      logical :: double_quotes                             !    "    "    "    "    "    "    "

!     Execution:

      verbose = ios == 1

      if (formatted) then

         open (lun, file=filename, status='old', iostat=ios)

      else ! Binary inputs are not implemented yet

         eps = epsilon (eps)
         if (eps < 1.e-10) then
            isdouble = 1
         else
            isdouble = 0
         end if

!        ios = Tec_xxx (...) ! Some utility not available yet in the Tec_io library
         write (*, '(/, a)') ' Reading of Tecplot binaries is not an option yet.'
         ios = 999
         go to 999
      end if

      if (ios /= 0) then
         write (*, '(3a)') routine, 'Trouble opening file ', trim (filename)
         go to 999
      end if

      if (formatted) then

!        Count the blocks, saving dimension records as strings, then allocate work space and decode the block sizes.
!        Variable names are also counted and returned.

         call count_blocks () ! Local procedure below

!        The file is left ready for reading the first block.

      else ! Not implemented yet; what's here applies to PLOT3D files
         read (lun, iostat=ios) nblocks
      end if

      if (ios /= 0) then
         write (*, '(a)') ' Error determining the number of grid blocks.'
         go to 999
      end if

      allocate (xyzq(nblocks), stat=ios)

      if (ios /= 0) then
         write (*, '(a)') ' Error allocating the array of blocks.'
         go to 999
      end if

      if (formatted) then

         call decode_dimensions () ! Local procedure below

      else ! Not implemented
         read (lun, iostat=ios) (xyzq(ib)%ni, xyzq(ib)%nj, xyzq(ib)%nk, ib = 1, nblocks)
      end if

      if (ios /= 0) then
         write (*, '(2a)') routine, ' Error determining the input block dimensions.'
      end if

  999 return

!     Local procedures for subroutine Tec_header_read:

      contains

!        -----------------------
         subroutine count_blocks ()
!        -----------------------

!        Count the blocks of a Tecplot ASCII file, and save block dimensions in a string array.
!        Also, return the title, the number of functions other than (x,y,z), and all variable names.
!        Assume the file starts with TITLE, followed by the variable names on one or more lines.
!
!        [TITLE = "..."]
!        VARIABLES = <any reasonable list, delimited either by commas or by ".." possibly with embedded commas, on 1 or more lines.
!
!        Assume the blocks (zones) start with strings like this:
!
!        ZONE T="G1"
!         I=iii, J=jjj, K=kkk, ...
!         DATAPACKING=.....
!         DT=(...... ...... )
!
!        The file is left ready for reading the first block or zone.

!        Local constants:

         character, parameter :: zone_t_l * 6 = 'zone t', zone_t_u * 6 = 'ZONE T'

!        Local variables:

         integer   :: i, i1, l, line, nb, nheader
         logical   :: looking_for_datapacking, looking_for_variables, title_present, variables_present
         character :: buffer * 80 ! Enough for more than just ZONE T;
                                  ! DPLR's Postflow can put datapacking & block dimensions on the same line as the zone title  :(

!        Execution:

!        Count the header lines before the first occurrence of ZONE:

         title_present = .false.;  variables_present = .false.
         nheader = 0
         buffer  = blank

         do while (index (buffer, 'ZONE') == 0) ! Up-cased below
            read (lun, '(a)') buffer
            l = len_trim (buffer)
            call upcase (buffer(1:l))
            nheader = nheader + 1
            if (index (buffer(1:l), 'TITLE') > 0) title_present = .true.
            if (index (buffer(1:l), 'VAR'  ) > 0) variables_present = .true.
         end do

         nheader = nheader - 1
         rewind (lun)

         if (verbose) then
            write (*, '(a, i4)') ' # header lines found:', nheader
            write (*, '(a, l3)') ' Title present?         ', title_present, &
                                 ' Variable names present?', variables_present
         end if

         title = blank

         if (title_present) then

            do line = 1, nheader
               read (lun, '(a)') packed_names ! Handy buffer used for the binary write option
               i = index (packed_names, 'TITLE')
               if (i == 0) i = index (packed_names, 'title')  ! Avoid upcasing variable names

               if (i > 0) then
                  last = len_trim (packed_names)
                  i = index (packed_names(1:last), quotes)
                  if (i == 0) then ! The title string must be on a separate line
                     read (lun, '(a)') packed_names
                     last = len_trim (packed_names)
                     i = index (packed_names(1:last), quotes)
                  end if
                  if (i < last - 1) title = packed_names(i+1:last-1)
                  if (verbose) write (*, '(1x, a)') packed_names(1:last)
                  exit
               end if
            end do

            rewind (lun)

         end if

!        Determine the number of variables from their names if possible.
!        If no variable names are found, they will be counted the original way from the DT=(SINGLE ... or (DOUBLE ...) line,
!        and defaulted to "V1", "V2", ... (and assumed the same for all zones).

         nvar = 0;  istart = 1 ! See pack_names

         if (variables_present) then

            looking_for_variables = .true.
            line = 0

            do while (line < nheader)

               line = line + 1
               read (lun, '(a)') names_buffer
               last =  len_trim (names_buffer)

               if (looking_for_variables) then

                  i = index (names_buffer(1:last), 'VARIABLES')
                  if (i == 0) i = index (names_buffer(1:last), 'variables')  ! Avoid up-casing variable names

                  if (i > 0) then
                     looking_for_variables = .false.
                     first = index (names_buffer(1:last), quotes)
                     double_quotes = first > 0

                     if (double_quotes) then

                        call pack_names ()  ! Variable name(s) on same line as VARIABLES; see local procedure below

                     else ! Handle names delimited by commas or spaces instead of double quote pairs

                        first = i + 10  ! Character after 'VARIABLES'

                        if (first <= last) call pack_names ()
                     end if

                  end if

                  cycle

               end if

               first = 1;  call pack_names () ! Any further lines in the header are assumed to be variable names

            end do

            if (verbose) write (*, '(a, i4)') ' Number of variables found:', nvar

            rewind (lun)

         end if

!        Initialize scanning for the number of zones by setting up to read a ZONE line:

         do line = 1, nheader
            read (lun, '(a)')
         end do

         line = nheader
         nb   = 0
         looking_for_datapacking = .true.

         do ! Until EOF

            line = line + 1
            read (lun, '(a)', iostat=ios) buffer

            if (ios < 0) then ! EOF
               ios = 0
               exit
            end if

            i = index (buffer(1:8), zone_t_l)           ! Assume ZONE T or zone t is at the beginning, to save some searching
            if (i == 0) i = index (buffer(1:8), zone_t_u)

            if (i > 0) then ! Start of a new zone

               nb = nb + 1
               if (verbose) write (*, '(a, i4)') ' Block:', nb

               if (nb > max_blocks) then
                  write (*, '(3a, i6)') routine, 'Too many blocks in file ', trim (filename), '.  Limit:', max_blocks
                  ios = 999
                  go to 999
               end if

               l = len_trim (buffer)
               if (buffer(l:l) /= quotes) then    ! There must be more info. on the zone title line
                  dimensions(nb) = buffer(i+7:l)  ! As though it were read
               else
                  line = line + 1
                  read (lun, '(a)', iostat=ios) dimensions(nb) ! Leading part of an I=iii, J=jjj, K=kkk, ... string

                  if (ios /= 0) then
                     write (*, '(3a, i6)') routine, 'Trouble reading dimension string in file ', trim (filename), '.  Line #:', line
                     ios = 999
                     go to 999
                  end if
               end if

               l = len_trim (dimensions(nb))
               if (verbose) write (*, '(1x, a)') dimensions(nb)(1:l)

               call upcase (dimensions(nb)(1:l)) ! In case datapacking is also on the dimensions line

               i = index (dimensions(nb)(1:l), 'POINT')
               if (i > 0) then
                  if (looking_for_datapacking) then
                     datapacking = 1
                     looking_for_datapacking = .false.
                  end if
               else
                  i = index (dimensions(nb)(1:l), 'BLOCK')
                  if (i > 0) then
                     datapacking = 2
                     looking_for_datapacking = .false.
                  end if
               end if

               if (i > 0) then  ! For decode_dimensions, we need to avoid the datapacking info.in dimensions(nb)
                  i = index (dimensions(nb)(1:l), 'I=')
                  if (i == 0) i = index (dimensions(nb)(1:l), 'I ')
                  dimensions(nb) = dimensions(nb)(i:)
               end if

               if (nb == 1) then ! Determine the datapacking type, and count the number of variables if no names were present

                  if (looking_for_datapacking) then

                     line = line + 1
                     read (lun, '(a)') names_buffer ! Handy buffer
                     last = len_trim  (names_buffer)

                     if (verbose) write (*, '(1x, a)') names_buffer(1:last)

                     if (index (names_buffer(1:last), 'POINT') > 0) then
                        datapacking = 1
                     else if (index (names_buffer(1:last), 'BLOCK') > 0) then
                        datapacking = 2
                     else
                        write (*, '(3a)') routine, 'Bad data packing found in zone 1: ', names_buffer(1:last), &
                                          ' File name: ', trim (filename)
                        ios = 999
                        go to 999
                     end if

                     looking_for_datapacking = .false.

                  end if

                  if (nvar == 0) then ! Names were evidently missing; count them from a DT=(....LE ....LE ) line

                     line = line + 1
                     read (lun, '(a)') names_buffer
                     last =  len_trim (names_buffer)

                     if (verbose) write (*, '(a)')  names_buffer(1:last)
                     i1 = index (names_buffer(1:last), 'DT')

                     if (i1 == 0) then
                        write (*, '(2a)') routine, 'DT record not found where expected in zone 1:', &
                                          ' Trying to count the variables in: ', names_buffer(1:last), &
                                          ' File name: ', trim (filename)
                        ios = 999
                        go to 999
                     end if

                     do i = i1, last
                        if (names_buffer(i:i) == 'E') nvar = nvar + 1 ! Works for SINGLE or DOUBLE
                     end do

                  end if

                  allocate (names(nvar))

                  num_q = nvar - 3
                  if (verbose) write (*, '(a, i4)') ' # variables found apart from (x,y,z):', num_q

                  if (variables_present) then ! Unpack them into the output array

                     first = 1;  last = iend - 1 ! See pack_names

                     do i = 1, nvar

                        call scan4 (packed_names, quotes, first, last, mark)

                        names(i) = packed_names(first:mark)
                        first = mark + 2
                     end do

                  else ! Default the variable names to "V1", "V2", ...

                     do i = 1, nvar
                        names(i) = blank
                        if (i < 10) then
                           write (names(i)(1:2), '(a1, i1)') 'V', i
                        else ! Assume no more than 99 variables
                           write (names(i)(1:3), '(a1, i2)') 'V', i
                        end if
                     end do

                  end if

                  if (verbose) write (*, '(a, (10(3x, a)))') ' Variable names:', (trim (names(i)), i = 1, nvar)

               end if ! End check of zone 1 (only)

            end if

         end do

!        Leave the file ready for reading the first block:

         rewind (lun)

         do line = 1, nheader
            read (lun, '(a)')
         end do

  999    nblocks = nb

         end subroutine count_blocks

!        ---------------------
         subroutine pack_names ()  ! Transfer names from names_buffer to next portion (istart) of packed_names
!        ---------------------

!        Pack_names ought to be internal to count_blocks, but internal subroutines cannot have internal subroutines.
!        Allowing for embedded blanks requires (new) scan4; the standard scan2 won't work.
!        The packed names are returned like this:  "X, m" "Y, m" "Z, m", "T, K" ...

         mark = 0

         do ! Until mark = 0

            if (double_quotes) then
               call scan4 (names_buffer, quotes, first, last, mark)
            else
               call scan2 (names_buffer,  '=, ', first, last, mark)
            end if

            if (mark == 0) exit

            nvar = nvar + 1

            packed_names(istart:istart) = quotes
            istart = istart + 1
            iend = istart + mark - first
            packed_names(istart:iend) = names_buffer(first:mark)
            iend = iend + 1
            packed_names(iend:iend) = quotes
            iend = iend + 1
            packed_names(iend:iend) = blank
            istart = iend + 1
            first  = mark + 2

         end do

         end subroutine pack_names

!        ----------------------------
         subroutine decode_dimensions ()
!        ----------------------------

!        Convert stored block dimension strings to integers.
!        The strings look like this:
!
!         I=iii, J=jjj, K=kkk, ...

!        Local constants:

         character, parameter :: delimiters * 6 = ' ,=IJK'

!        Local variables:

         integer   :: first, ib, last, mark, n, ndim(3)
         character :: format_string * 4

!        Execution:

         format_string = '(i?)'

         do ib = 1, nblocks

            first = 1; last = len_trim (dimensions(ib))

            do n = 1, 3
               call scan2 (dimensions(ib), delimiters, first, last, mark)

               write (format_string(3:3), '(i1)') mark - first + 1
               read (dimensions(ib)(first:mark), format_string, IOSTAT=ios) ndim(n)

               if (ios /= 0) then
                  write (*, '(3a, i5)') routine, 'Trouble decoding block size from file ', trim (filename), '.  Block #:', ib
                  go to 999
               end if

               first = mark + 2
            end do ! Next i,j,k dimension

            xyzq(ib)%ni = ndim(1);  xyzq(ib)%mi = ndim(1) ! These may be redundant
            xyzq(ib)%nj = ndim(2);  xyzq(ib)%mj = ndim(2)
            xyzq(ib)%nk = ndim(3);  xyzq(ib)%mk = ndim(3)

         end do ! Next block

   999   return

         end subroutine decode_dimensions

      end subroutine Tec_header_read

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine Tec_header_write (lun, filename, formatted, datapacking, title, num_q, names, ios)
!
!     Open a Tecplot data file and write the header records for 3-space structured Tecplot data.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer,   intent (in)  :: lun          ! Logical unit number
      character, intent (in)  :: filename*(*) ! Name of the Tecplot file
      logical,   intent (in)  :: formatted    ! T|F
      integer,   intent (in)  :: datapacking  ! 1 means DATAPACKING=POINT; 2 means DATAPACKING=BLOCK (for all zones)
      character, intent (in)  :: title*(*)    ! First line of the Tecplot file
      integer,   intent (in)  :: num_q        ! Number of variables, not counting x,y,z
      character, intent (in)  :: names (num_q+3) * (name_limit) ! Variable names, left-justified
      integer,   intent (out) :: ios          ! 0 = no error

!     Local constants:

      character, parameter :: colon * 1 = ':', comma * 1 = ',', routine * 19 = ' Tec_header_write: '

!     Local variables:

      integer    :: k, l, m, n, nvars

!     Tecplot I/O library function:

      integer    :: TecIni

!     Execution:

      nvars = 3 + num_q

      if (formatted) then

         open (lun, file=filename, status='unknown', iostat=ios)

         if (ios /= 0) then
            write (*, '(3a)') routine, 'Trouble opening file ', trim (filename)
            go to 999
         end if

         write (lun, '(4a)') 'TITLE = ', quotes, trim (title), quotes, 'VARIABLES ='
         write (lun, '(3a)') (quotes, trim (names(n)), quotes, n = 1, nvars)

      else ! Use the Tecplot I/O function for initializing writing to a binary file

         l = 0
         do n = 1, nvars
            l = len_trim (names(n)) + l
         end do
         l = l + nvars + 1 ! For the commas and the final null

         if (l > max_length) then
            write (*, '(2a, i5)') routine, 'Packed variable names are too long.  Length required:', l
            ios = 999
            go to 999
         end if

         m = 0
         do n = 1, nvars
            l = len_trim (names(n))
            packed_names(m+1:m+l) = names(n)(1:l);  l = m + l + 1
            do k = m+1, m+l
               if (packed_names(k:k) == comma) packed_names(k:k) = colon  ! Double quotes still don't allow embedded commas
            end do
            packed_names(l:l)     = comma;          m = l
         end do
         l = m + 1
         packed_names(l:l) = null

!!!      write (*, '(a, i6)') ' Packed names(1:l) passed to TecIni, with l =', l
!!!      write (*, '(a)') packed_names(1:l)

         eps = epsilon (eps)
         if (eps < 1.e-10) then
            isdouble = 1
         else
            isdouble = 0
         end if

         ios = TecIni (trim (title) // null, packed_names(1:l),       &
                       trim (filename) // null, '.' // null, 0, isdouble)    ! 0 = No debug

         if (ios /= 0) then
            write (*, '(3a)') routine, 'Trouble initializing Tecplot file ', trim (filename)
         end if

      end if

  999 return

      end subroutine Tec_header_write

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine Tec_block_allocate (block, num_q, ios)
!
!     Allocate the x, y, z arrays and optional q array for one block or zone.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent(inout) :: block
      integer,          intent (in)   :: num_q
      integer,          intent (out)  :: ios

!     Local constants:

      character, parameter :: routine_etc * 33 = ' Tec_block_allocate trouble with '

!     Local variables:

      integer :: ni, nj, nk

!     Execution:

      ni = block%ni;  nj = block%nj;  nk = block%nk

      allocate (block%x(ni,nj,nk), block%y(ni,nj,nk), block%z(ni,nj,nk), stat=ios)

      if (ios /= 0) then
         write (*, '(2a)') routine_etc, 'x, y, z.'
         go to 999
      end if

      if (num_q > 0) then

         allocate (block%q(num_q,ni,nj,nk), stat=ios)

      else ! Avoid undefined arguments elsewhere

         allocate (block%q(1,1,1,1), stat=ios)

      end if

      if (ios /= 0) write (*, '(2a)') routine_etc, 'q array.'

  999 if (ios /= 0) write (*, '(a, 4i5)') '  ni, nj, nk, num_q:', ni, nj, nk, num_q

      end subroutine Tec_block_allocate

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine Tec_block_read (lun, formatted, datapacking, zone_title, ni, nj, nk, num_q, x, y, z, q, ios)
!
!     Read one 3-space block of structured data from a Tecplot file (binary or ASCII, POINT or BLOCK packing order).
!     The block dimensions have already been determined as part of reading the file header.
!     The file is presumed to be positioned ready to read a zone title as the first record.
!
!     For now, assume the contents of a Tecplot ASCII file look like this (with a leading blank or two allowed in header records):
!
!      ZONE T="xxx"
!      I=iii, J=jjj, K=kkk, ZONETYPE=Ordered
!      DATAPACKING=xxxxx
!      DT=(.....E .....E .... )
!      6.647333145E+00 ....
!
!     Alternatively, output from DPLR's Postflow has zone headers that look like this and these are now handled:
!
!      zone t="flow3d" F=point, i=  97 j=  33 k= 129
!
!     Reading of a Tecplot binary file awaits possible addtions to Amtec's tecio.a library.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer,   intent (in)                              :: lun              ! Logical unit number
      logical,   intent (in)                              :: formatted        ! T|F
      integer,   intent (in)                              :: datapacking      ! 1 means DATAPACKING=POINT; 2 means DATAPACKING=BLOCK
      character, intent (out)                             :: zone_title * (*) ! Zone title
      integer,   intent (in)                              :: ni, nj, nk       ! Grid block dimensions
      integer,   intent (in)                              :: num_q            ! Number of additional functions; num_q >= 0
      real,      intent (out), dimension (ni,nj,nk)       :: x, y, z          ! Grid block coordinates
      real,      intent (out), dimension (num_q,ni,nj,nk) :: q                ! Further functions, if num_q > 0
      integer,   intent (out)                             :: ios              ! 0 = no error

!     Local variables:

      integer   :: i, j, k, l, n

      character :: buffer * 80

!     Execution:

      read (lun, '(a)') buffer ! ZONE T="..."
 
      l = len_trim (buffer)
      if (index (buffer(1:l), 'ZONE T') <= 0) then
         if (index (buffer(1:l), 'zone t') <= 0) then
            write (*, '(/, (a))') ' Tec_block_read: Zone title not found in first line:', buffer(1:l)
            ios = 1
            go to 99
         end if
      end if

      i = index (buffer(1:l), quotes) + 1
      j = index (buffer(i:l), quotes) + i - 1
      zone_title = buffer(i:j-1)

      if (j == l) then     ! Assume the following 3 lines follow a zone title-only line
!!!      write (*,  '(a)') buffer(1:j+1)
         read (lun, '(a)') buffer ! I=iii, ...
!!!      write (*,  '(a)') trim (buffer)
         read (lun, '(a)') buffer ! DATAPACKING=...
!!!      write (*,  '(a)') trim (buffer)
         read (lun, '(a)') buffer ! DT=(...
!!!      write (*,  '(a)') trim (buffer)
      else                 ! Assume only one zone header line
      end if

      select case (datapacking)

         case (1) ! POINT

            if (formatted) then

               if (num_q == 0) then
                  read (lun, *, iostat=ios) (((x(i,j,k), y(i,j,k), z(i,j,k), i = 1, ni), j = 1, nj), k = 1, nk)
               else
                  read (lun, *, iostat=ios) (((x(i,j,k), y(i,j,k), z(i,j,k), q(:,i,j,k), i = 1, ni), j = 1, nj), k = 1, nk)
               end if

               if (ios /= 0) then
                  write (*, '(/, a, i4)') ' Trouble reading data in POINT order.  Logical unit:', lun
                  go to 99
               end if

            else ! Not implemented yet
               if (num_q == 0) then
                  read (lun, iostat=ios) (((x(i,j,k), y(i,j,k), z(i,j,k), i = 1, ni), j = 1, nj), k = 1, nk)
               else
                  read (lun, iostat=ios) (((x(i,j,k), y(i,j,k), z(i,j,k), q(:,i,j,k), i = 1, ni), j = 1, nj), k = 1, nk)
               end if
            end if

         case (2) ! BLOCK

            if (formatted) then

               read (lun, *, iostat=ios) x, y, z
               if (ios /= 0) then
                  write (*, '(/, a, i4)') ' Trouble reading (x,y,z) data in BLOCK order.  Logical unit:', lun
                  go to 99
               end if

               if (num_q > 0) then
                  read (lun, *, iostat=ios) (q(n,1:ni,1:nj,1:nk), n = 1, num_q)
                  if (ios /= 0) then
                     write (*, '(/, a, 2i4)') ' Trouble reading q data in BLOCK order.  n and lun:', n, lun
                     go to 99
                  end if
               end if

            else ! Not implemented yet
               read (lun, iostat=ios) x, y, z
               if (num_q > 0) read (lun, iostat=ios) q
            end if

         case default

            ios = 1

      end select

   99 return

      end subroutine Tec_block_read

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine Tec_block_write (lun, formatted, datapacking, zone_title, ni, nj, nk, num_q, x, y, z, q, ios)
!
!     Write one 3-space block of structured data from a Tecplot file (binary or ASCII, POINT or BLOCK packing order).
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer,   intent (in)                             :: lun              ! Logical unit number
      logical,   intent (in)                             :: formatted        ! T|F
      integer,   intent (in)                             :: datapacking      ! 1 means DATAPACKING=POINT; 2 means DATAPACKING=BLOCK
      character, intent (in)                             :: zone_title * (*) ! Zone title
      integer,   intent (in)                             :: ni, nj, nk       ! Grid block dimensions
      integer,   intent (in)                             :: num_q            ! Number of additional functions; num_q >= 0
      real,      intent (in), dimension (ni,nj,nk)       :: x, y, z          ! Grid block coordinates
      real,      intent (in), dimension (num_q,ni,nj,nk) :: q                ! Further functions, if num_q > 0
      integer,   intent (out)                            :: ios              ! 0 = no error

!     Local constants:

      character, parameter :: routine * 18 = ' Tec_block_write: '

!     Local variables:

      integer    :: i, j, k, n, npts, nvars
      character  :: aformat*6, eformat*14

      real,      allocatable, dimension (:)     :: one_point
      real,      allocatable, dimension (:,:,:) :: one_var
      character, allocatable, dimension (:)     :: chars * 7

!     Amtec's Tecplot I/O library functions:

      integer :: TecDat, TecZne

!     Execution:

      npts  = ni * nj * nk
      nvars = 3 + num_q

      if (formatted) then

         write (lun, '(4a)') 'ZONE T=', quotes, trim (zone_title), quotes
         write (lun, '(a, i3, a, i3, a, i3, a)') ' I=', ni, ', J=', nj, ', K=', nk, ', ZONETYPE=Ordered'

         aformat = '(nnna)'
         write (aformat(2:4), '(i3)') nvars + 2

         allocate (chars(nvars))

         if (isdouble == 1) then
            chars(:) = 'DOUBLE '
            eformat = '(1p, nnne17.9)'
         else ! 0
            chars(:) = 'SINGLE '
            eformat = '(1p, nnne15.7)'
         end if

         write (eformat(6:8), '(i3)') nvars

         select case (datapacking)

         case (1) ! POINT

            write (lun, '(a)') ' DATAPACKING=POINT'
            write (lun, aformat) ' DT=(', chars, ')'

            allocate (one_point(nvars))

            do k = 1, nk
               do j = 1, nj
                  do i = 1, ni
                     one_point(1) = x(i,j,k)
                     one_point(2) = y(i,j,k)
                     one_point(3) = z(i,j,k)
                     if (nvars > 3) one_point(4:nvars) = q(:,i,j,k)
                     write (lun, eformat) one_point
                  end do
               end do
            end do

            deallocate (one_point)

         case (2) ! BLOCK

            write (lun, '(a)') ' DATAPACKING=BLOCK'
            write (lun, aformat) ' DT=(', chars, ')'
            write (lun, '(1p, 6e17.9)') x
            write (lun, '(1p, 6e17.9)') y
            write (lun, '(1p, 6e17.9)') z

            if (nvars > 3) then
               do n = 1, num_q
                  write (lun, '(1p, 6e17.9)') q(n,1:ni,1:nj,1:nk)
               end do
            end if

         case default

            ios = 999

         end select

         deallocate (chars)

      else ! Binary output

         select case (datapacking)

         case (1) ! POINT

            ios = TecZne (trim (zone_title) // null, ni, nj, nk, 'POINT' // null, null)

            if (ios /= 0) then
               write (*, '(/, 2a, i4)') routine, 'Trouble writing zone info.  Logical unit:', lun
               go to 999
            end if

            allocate (one_point(nvars))

            do k = 1, nk
               do j = 1, nj
                  do i = 1, ni
                     one_point(1) = x(i,j,k)
                     one_point(2) = y(i,j,k)
                     one_point(3) = z(i,j,k)
                     if (nvars > 3) one_point(4:nvars) = q(:,i,j,k)
                     ios = TecDat (nvars, one_point, isdouble)
                  end do
               end do
            end do

            deallocate (one_point)

         case (2) ! BLOCK

            ios = TecZne (trim (zone_title) // null, ni, nj, nk, 'BLOCK' // null, null)

            if (ios /= 0) then
               write (*, '(/, 2a, i4)') routine, 'Trouble writing zone info.  Logical unit:', lun
               go to 999
            end if

            ios = TecDat (npts, x, isdouble)  ! Assume x(:,:,:) has been dynamically allocated to fit the data
            ios = TecDat (npts, y, isdouble)
            ios = TecDat (npts, z, isdouble)

            if (nvars > 3) then
               allocate (one_var(ni,nj,nk))
               do n = 1, num_q
                  one_var = q(n,1:ni,1:nj,1:nk)
                  ios = TecDat (npts, one_var, isdouble)
               end do
               deallocate (one_var)
            end if

         case default

            ios = 999

         end select

      end if

  999 continue

      end subroutine Tec_block_write

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine deallocate_blocks (ib1, ib2, nf, xyzq, ios)

!     Deallocate the indicated blocks.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Arguments:

      integer, intent (in)  :: ib1, ib2              ! Block range to deallocate
      integer, intent (in)  :: nf                    ! > 0 means deallocate %q too
      type (grid_type), intent (inout) :: xyzq(ib2)  ! At least this many
      integer, intent (out) :: ios                   ! 0 means no problem

!     Local variables:

      integer :: ib

!     Execution:

      do ib = ib1, ib2
         deallocate (xyzq(ib)%x, xyzq(ib)%y, xyzq(ib)%z, stat=ios)
         if (ios /= 0) go to 90
         if (nf > 0) then
            deallocate (xyzq(ib)%q, stat=ios)
            if (ios /= 0) go to 90
         end if
      end do

      go to 99

   90 write (*,*) &
        ' Trouble deallocating grid_type variable. ib1, ib2, nf: ', ib1, ib2, nf
   99 return

      end subroutine deallocate_blocks

   end module Tecplot_io_module
