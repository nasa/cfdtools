!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module grid_header_structure

!  Module prompted by Tecplot 360 dataset auxiliaries.  Storing them with block 1 (only) was attempted but is not a good idea.
!  The dataset file name, title, and variable names, etc., belong here too.

   type grid_header
      character          :: filename * 80           ! Name of file if the dataset is to be read or written
      logical            :: formatted               ! T | F for ASCII | binary
      integer            :: ndim                    ! 2 | 3 for 2D | 3D data (x and y only, or x, y, and z)
      integer            :: numq                    ! # additional variables beyond ndim; numq >= 0
      integer            :: nblocks                 ! # blocks in dataset (not always known?)
      integer            :: datapacking             ! Same for all zones;  0 => POINT order; 1 => BLOCK order (more efficient)
      character          :: title * 80              ! Dataset title
      character, pointer :: varname (:) * 32        ! Variable names; embedded blanks are permitted; same length as in Tecplot_io

      integer            :: ndatasetaux             ! # auxiliary data strings associated with dataset; ndatasetaux >= 0
      character, pointer :: datasetauxname (:) * 32 ! Names of dataset auxiliary items, if any (no embedded blanks)
      character, pointer :: datasetauxvalue(:) * 32 ! Corresponding values (excluding leading and trailing double quotes)
   end type grid_header

   end module grid_header_structure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module grid_block_structure

!  Module for a structured grid block or zone with associated flow field variables and a zone name (Tecplot_io_module version).
!  This version is a superset of that originally used with the xyzq_io package, and is intended to be compatible with that package.
!  If both packages are used by an application, this version of the module is the one that should be compiled.

   type grid_type
      integer            :: nzoneheaderlines        ! # lines in zone header, to facilitate reading zones after counting them
      character          :: zone_title * 32         ! Zone/block title; embedded blanks are permitted; same length as in Tecplot_io
      integer            :: nzoneaux                ! # auxiliary data strings associated with zone; nzoneaux >= 0
      character, pointer :: zoneauxname (:) * 32    ! Names of zone auxiliary items, if any (no embedded blanks)
      character, pointer :: zoneauxvalue(:) * 32    ! Corresponding values (excluding leading and trailing double quotes)
      real               :: solutiontime            ! Time (or some other useful real number) associated with zone
      integer            :: ni                      ! # points in i direction
      integer            :: nj                      ! # points in j direction
      integer            :: nk                      ! # points in k direction
      integer            :: mi                      ! # dependent variables in i direction
      integer            :: mj                      ! # dependent variables in j direction
      integer            :: mk                      ! # dependent variables in k direction
      real               :: xmin, xmax              ! Data range for the grid block in the x direction
      real               :: ymin, ymax              ! ... and the y direction ...
      real               :: zmin, zmax              ! ... and the z direction
      real, dimension(:,:,:),    pointer :: x,y,z   ! Grid coordinates
      real, dimension(:,:,:,:),  pointer :: q       ! Dependent variables
      integer, dimension(:,:,:), pointer :: iblank  ! 0 => suppress the point
   end type grid_type

   end module grid_block_structure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module Tecplot_io_module

!  This module packages I/O for Tecplot ASCII and binary files representing 2|3-space structured multiblock grids and associated
!  function data.  Both BLOCK and POINT data packing formats are supported.  The function values are assumed to be vertex-centered,
!  with (x,y) or (x,y,z) grid coordinates as the first 2 or 3 variables.
!
!  Several lower-level utilities used by the whole-grid I/O utilities are made public to permit processing of one block at a time
!  (as opposed to reading or writing the whole file in one call).
!
!  Note that, in order to pass unallocated arrays as arguments, the calling program must declare the x/y/z/q data type with the
!  "pointer" attribute as opposed to "allocatable".
!
!  04/02/04  David Saunders  Initial implementation of XYZQ_IO package for PLOT3D files.
!  07/29/04    "       "     Adaptation for Tecplot files, but reading of Tecplot binaries is unsolved.
!  09/24/04    "       "     Automated handling of -r4 and -r8 compilation via new IsDouble internal variable.
!  05/04/05    "       "     Variable names are now returned from a formatted read (originally overlooked).
!  05/06/05    "       "     Verbose mode can be activated by entering ios = 1 to Tec_header_read.
!  05/13/05    "       "     The datapacking (1|2 = POINT|BLOCK) should be output for a header read, not input.
!  05/14/05    "       "     Binary outputs cannot handle embedded commas in variable names.  Substitute colons.
!  06/16/05    "       "     Don't assume 'ZONE T' starts in column 1 in count_blocks; handle trailing control-Ms.
!  11/30/05    "       "     Added deallocate_blocks.
!  12/04/05    "       "     Fixed Tec_block_read where numq == 3 should have been numq == 0.
!  12/05/05    "       "     If numq = 0 in Tec_block_allocate, allocate %q(1,1,1,1) to avoid undefined arguments elsewhere.
!  12/06/05    "       "     Extended to handle simpler input of the type written by flow solver DPLR's post-processor.
!  07/27/06 -  "       "     Tecplot 360 extensions via new fields in grid_type data structure for auxiliary data and solution time.
!  08/10/06                  Handling of keywords has been made more systematic.  The calling sequences for Tec_header_write and
!                            Tec_block_read/_write have changed.  Tec_block_allocate knows the block number now.
!  08/10/06    "       "     Merged the 2D case into the 3D package via public variable io_dim = 2 or 3 (now moved to the header).
!  08/11/06    "       "     Introduced a separate derived data type for dataset header information to deal with auxiliaries,
!                            after storing them with zone 1 (only) was found too awkward.  File name, title, and variable names,
!                            etc., naturally belong in the same data structure.
!  09/07/06    "       "     Enter header%ndim = -1 to Tecplot_read or Tec_header_read to automate distinguishing 2D from 3D.
!  10/30/06    "       "     A 58-variable DT=(SINGLE SINGLE ... ) line required a bigger buffer and a test for exceeding it.
!  11/03/06    "       "     Blank zone titles need to be handled by checking for mark < 0 returned by scan4.
!  11/27/06    "       "     Realized that one global implicit none does NOT cover the internal procedures (where all variables were
!                            declared anyway).   Now, implicit none is present in each routine as it should be.
!  03/02/07    "       "     Reading BLOCK data x, y[, z] first then reading q if numq > 0 can lose the first part of q if the
!                            coordinates don't fill the line.  Read coordinates and functions in the same read.
!  09/30/12    "       "     Limiting zone dimensions to 3 digits was short-sighted.  5 is now the limit.
!  10/01/12    "       "     Dinesh needed 70 variables; raised string lengths to be more than enough.
!  12/21/15    "       "     Raised the zone limit from 1000 to 3000.
!
!  Author:  David Saunders, ELORET Corporation (originally; now with ERC, Inc.), NASA Ames Research Center, Moffett Field, CA
!
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_header_structure ! Defines "grid_header" derived data type for one grid
   use grid_block_structure  ! Tecplot-specific version to define "grid_type" derived data type for one block

   implicit none

   private

!  Constants used by the package:

   integer,   parameter :: len_string = 64         ! Room and then some for a ' I=iii, J=jjj, K=kkk' string
   integer,   parameter :: len_buffer = 1500       ! Limit on length of an input line, most likely the DT=(SINGLE SINGLE ... ) line
   integer,   parameter :: max_blocks = 3000       ! Limit on number of blocks in a file being read
   integer,   parameter :: max_length = 3000       ! Limit on length of string into which variable names are packed
   integer,   parameter :: name_limit = 32         ! Limit on the length of variable names; must match value in modules above

   real,      parameter :: undefined  = -999.

   logical,   parameter :: false = .false.,        &
                           true  = .true.
   character, parameter :: blank  * 1  = ' ',      &
                           null   * 1  = char (0), &
                           quotes * 1  = '"'

!  Internal variables used by the package:

   integer   :: IsDouble                           ! 0 = single precision data (-r4); 1 = double precision data (-r8)
   integer   :: nzoneheaderlines (max_blocks)      ! Temporary storage for # zone header lines while # zones is being counted
   integer   :: nzoneauxlines (max_blocks)         !     "         "             auxiliary  "
   integer   :: isize (max_blocks)                 ! Temporary storage for block dimensions
   integer   :: jsize (max_blocks)
   integer   :: ksize (max_blocks)
   integer   :: ndim                               ! Internal equivalent of header variable ndim
   integer   :: numq                               ! Internal equivalent of header variable numq
   integer   :: nvar                               ! Internal variable = numq + ndim
   integer   :: ndatasetaux                        ! Internal equivalent of header variable ndatasetaux
   integer   :: nblocks                            ! Internal equivalent of header variable nblocks
   real      :: eps                                ! epsilon (eps) allows IsDouble to be assigned for binary reads & writes
   logical   :: formatted                          ! Internal equivalent of header variable formatted
   logical   :: two_d                              ! Internal equivalent of header variable ndim
   logical   :: valid, verbose                     ! Pass ios = 1 to Tec_header_read to activate printing of header/zone info.
   character :: buffer * (len_buffer)              ! For a line of input; its length is dominated by DT=(SINGLE SINGLE ... ) lines
   character :: packed_names * (max_length)        ! For binary output; reused for some of the reading
                                                   ! Fortran 90 doesn't support dynamically variable string lengths
!  Parsing utilities used by the package:

   logical   :: number
   external  :: number, scan2, scan4, upcase

!  Tecplot procedures:

   integer   :: TecIni110, TecAuxStr110, TecZne110, TecZAuxStr110, TecDat110, TecEnd110
   external  :: TecIni110, TecAuxStr110, TecZne110, TecZAuxStr110, TecDat110, TecEnd110

!  Utilities provided for applications:

   public :: Tecplot_read       ! Reads  an entire structured Tecplot file
   public :: Tecplot_write      ! Writes an entire structured Tecplot file

   public :: Tec_header_read    ! Reads  Tecplot file header records
   public :: Tec_header_write   ! Writes Tecplot file header records
   public :: Tec_block_allocate ! Allocate the array fields for one block or zone of a Tecplot file
   public :: Tec_block_read     ! Reads  one block of a structured Tecplot file
   public :: Tec_block_write    ! Writes one block of a structured Tecplot file

   public :: clone_header       ! Derives one dataset header from another
   public :: clone_zone         ! Derives one zone (block) from another and allocates array fields
   public :: deallocate_header  ! Deallocates any dataset auxiliaries in the given grid header
   public :: deallocate_blocks  ! Deallocates all array fields of the indicated blocks of a grid

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Tecplot_read (lun, header, grid, ios)

!     Read a Tecplot file of 3-space structured data, binary or ASCII, with BLOCK or POINT data packing.
!     The first ndim variables are returned as the x, y[, z] fields of each block in grid(:).
!     Remaining variables become the "q" field packed in "point" order (n-tuples, where n = numq = # variables found minus ndim).
!     The input file title is also returned, along with the variable names (up to name_limit characters each).
!     Titles for the blocks or zones are returned in the appropriate field of each element of the grid array.
!     Any auxiliary dataset and auxiliary zone records are returned in the header structure and zone array structures respectively.
!     The file is opened and closed here.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Arguments:

      integer, intent (in)               :: lun      ! Logical unit for the Tecplot file;
                                                     ! opened here upon entry; closed here before the return
      type (grid_header), intent (inout) :: header   ! Data structure containing dataset header information; input with file name

      type (grid_type), pointer          :: grid(:)  ! Array of structured blocks containing x,y[,z] and optional q data;
                                                     ! allocated here prior to the reading

      integer, intent (inout)            :: ios      ! 1 on input activates writing of header/zone info. during the read;
                                                     ! 0 on output means no error detected, else diagnostics are written to
                                                     ! standard output; early return follows
!     Local variables:

      integer :: ib

!     Execution:
!     ----------

!     Open the file, read the header records, determine the number of blocks or zones, allocate an array of grid structures,
!     and set up the block dimensions:

      call Tec_header_read (lun, header, grid, ios)

      if (ios /= 0) go to 999

      ndim = header%ndim
      numq = header%numq

      do ib = 1, header%nblocks

         call Tec_block_allocate (grid(ib), ndim, numq, ios)

         if (ios /= 0) then
            write (*, '(a, i5)') ' Tecplot_read: Trouble allocating block #', ib
            write (*, '(2a)')    ' File name: ', trim (header%filename)
            go to 999
         end if

         call Tec_block_read (lun, header, grid(ib), ios)

         if (ios /= 0) then
            write (*, '(a, i4)') ' Tecplot_read:  Trouble reading block #', ib
            write (*, '(2a)')    ' File name: ', trim (header%filename)
            if (header%ndim == 2) then
               write (*, '(a, 3i5)') ' Block dimensions & numq: ', grid(ib)%ni, grid(ib)%nj, numq
            else
               write (*, '(a, 4i5)') ' Block dimensions & numq: ', grid(ib)%ni, grid(ib)%nj, grid(ib)%nk, numq
            end if
            go to 999
         end if

      end do

  999 if (header%formatted) then
         close (lun)
      else
         ios = TecEnd110 ()
         if (ios /= 0) write (*, '(2a)') ' Trouble closing binary file ', trim (header%filename)
      end if

      end subroutine Tecplot_read

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Tecplot_write (lun, header, grid, ios)

!     Write a Tecplot file of 2- or 3-space structured data, binary or ASCII, with BLOCK or POINT data packing.
!     The first ndim variables written for each block are the x, y[, z] fields of the grid blocks.
!     Any further variables are appended, giving a total of ndim + numq variables.
!     Titles for the blocks or zones are written from the appropriate field of each element of the grid array.
!     Likewise for any auxiliary data records.
!     The file is opened and closed here.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Arguments:

      integer,            intent (in)  :: lun      ! Logical unit for the Tecplot file if it's ASCII;
                                                   ! opened here upon entry; closed here before the return
      type (grid_header), intent (in)  :: header   ! Data structure containing dataset header information; input with file name
      type (grid_type),   pointer      :: grid(:)  ! Array of structured blocks containing x, y[, z] and optional q data
      integer,            intent (out) :: ios      ! 0 means no error detected, else diagnostics are written to standard output;
                                                   ! early return follows
!     Local variables:

      integer :: ib

!     Execution:
!     ----------

!     Open the file and write the header records:

      call Tec_header_write (lun, header, ios)

      if (ios /= 0) go to 999

      numq = header%numq

      do ib = 1, header%nblocks

         call Tec_block_write (lun, header, grid(ib), ios)

         if (ios /= 0) then
            write (*, '(a, i4)')  ' Tecplot_write:  Trouble writing block #', ib
            write (*, '(2a)')     ' File name: ', trim (header%filename)
            if (header%ndim == 2) then
               write (*, '(a, 3i5)') ' Block dimensions & numq: ', grid(ib)%ni, grid(ib)%nj, numq
            else
               write (*, '(a, 4i5)') ' Block dimensions & numq: ', grid(ib)%ni, grid(ib)%nj, grid(ib)%nk, numq
            end if
            go to 999
         end if

      end do

  999 if (header%formatted) then
         close (lun)
      else
         ios = TecEnd110 ()
         if (ios /= 0) write (*, '(2a)') ' Trouble closing binary file ', trim (header%filename)
      end if

      end subroutine Tecplot_write

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine Tec_header_read (lun, header, grid, ios)
!
!     Open a Tecplot data file and read the file and zone header records for 2- or 3-space structured Tecplot data.
!     The grid(*) array is allocated and zone header information is entered.
!     Before returning, the file is rewound and advanced ready for reading the first block.
!
!     The simplest type of header has no title and a single list of variable names like this from DPLR's Postflow:
!
!      variables=x,y,z,rho,p,T,Tv,n2,o2,no,no+,n2+,o2+,n,o,n+,o+,e,u,v,w,h,M,mu,kap,Vel,ent
!      zone t="flow3d" F=point, i=  97 j=  33 k= 129   <start of first block>
!
!     If no variable names are present, they are defaulted to "X", "Y", ["Z", ]"V1", "V2", ..., but a DT=(... ) record must then
!     accompany zone 1 so that the number of variables can be determined.
!
!     If header%ndim is input as 2 or 3, that is assumed to be correct, but entering -1 means this routine will look for z or Z
!     as the name of the third variable and set header%ndim appropriately.  The latter option requires names in the data file.
!
!     The format originally handled looked like this:          Tecplot 360 extensions allow for the format below:
!
!     TITLE     = "CFD11L_DPLR_TCG_PhyMod1"                    TITLE = "Dataset title"                 [optional]
!     VARIABLES = "x, m"                                       VARIABLES = "varname1", "varname2", ... [0 or more lines]
!     "y, m"                                                   DATASETAUXDATA name="value"             [0 or more, 1 per line]
!     "z, m"                                                   ZONE T="Zone 1 title"                   [optional zone 1 title]
!     "rho, kg/m^3"                                            ZONETYPE=Ordered                        [understood]
!     "p, Pa"                                                  I=idim, J=jdim, K=kdim                  [commas optional, all lines]
!     "T, K"                                                   DATAPACKING=BLOCK                       [or POINT; also: F=POINT]
!     "c_N_2"                                                  DT=(SINGLE SINGLE .... )                [or (DOUBLE DOUBLE ... )]
!     "c_O_2"                                                  SOLUTIONTIME=time                       [optional real number]
!     "c_N_O"                                                  AUXDATA name="value"                    [0+ per zone, 1 per line]
!     :::::::                                                  xxx.xx  xxx.xx                          [Data proper]
!     ZONE T="G1"  [start of first zone]                       ::::::  ::::::
!                                                              ZONE T="Zone 2 title"                   [Repeat for further zones]
!
!     Blank header lines and comment lines starting with # are also handled now.
!     Reading of binary files is incomplete.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Arguments:

      integer, intent (in)               :: lun     ! Logical unit number
      type (grid_header), intent (inout) :: header  ! Data structure containing dataset header information; input with file name
      type (grid_type), pointer          :: grid(:) ! Array of derived data types for structured grid blocks, allocated here
      integer, intent (inout)            :: ios     ! 1 on input means write header/zone info. during reading (verbose mode);
                                                    ! 0 on output means no error detected
!     Local constants:

      character, parameter :: routine * 18 = ' Tec_header_read: '

!     Local variables:

      integer   :: ib, first, iend, istart, last, lentrim, line, mark, nb, nfileheaderlines
      logical   :: double_quotes, names_in_header

!     Execution:

      verbose = ios == 1;  ndim = header%ndim;  two_d = ndim == 2  ! This covers the absent-names case, but see if ndim = -1 below.

      formatted = header%formatted

      if (formatted) then

         open (lun, file=trim (header%filename), status='old', iostat=ios)

      else ! Binary inputs are not implemented yet

         eps = epsilon (eps)
         if (eps < 1.e-10) then
            IsDouble = 1
         else
            IsDouble = 0
         end if

!        ios = Tec_xxx (...) ! Some utility not available yet in the Tec_io library
         write (*, '(/, a)') ' Reading of Tecplot binaries is not an option yet.'
         ios = 999
         go to 999
      end if

      if (ios /= 0) then
         write (*, '(3a)') routine, 'Trouble opening file ', trim (header%filename)
         go to 999
      end if

      if (formatted) then

!        Read the file header lines up to and including the first 'ZONE' line.
!        Any auxiliary zone data cannot be stored yet, though, till the grid blocks have been allocated.
!        Variable names are not assigned here either, because they may not be in the file header; count_blocks assigns them.

         call read_file_header () ! Local procedure below

!        Count the blocks, saving dimension records as strings for decoding after the right number of blocks has been allocated.
!        If variable names weren't in the file header, they'll be counted from zone 1 and defaulted as X, Y, [Z, ]V1, V2, ...

         call count_blocks ()     ! Local procedure below

      else ! Not implemented yet; what's here applies to PLOT3D files

         read (lun, iostat=ios) nblocks

      end if

      if (ios /= 0) then
         write (*, '(2a)') ' Error determining the number of grid zones in ', trim (header%filename)
         go to 999
      end if

!     Option to detect the number of dimensions:

      if (ndim == -1) then ! The variable names were present (else count_blocks sets it)

         if (header%varname(3)(1:1) == 'z' .or. header%varname(3)(1:1) == 'Z') then
            ndim = 3
         else
            ndim = 2
         end if

         header%ndim = ndim;  two_d = ndim == 2

      end if

      numq        = nvar - ndim
      header%numq = numq

      if (verbose) write (*, '(a, i2)') ' # dimensions:', ndim, ' # functions: ', numq

!     Allocate an array of grid block derived data types:

      allocate (grid(nblocks), stat=ios)

      if (ios /= 0) then
         write (*, '(a)') ' Error allocating array of grid blocks.'
         go to 999
      end if

!     Now we can assign block dimensions and the numbers of header and auxiliary data records to each zone.

      do ib = 1, nblocks
         grid(ib)%ni = isize(ib);  grid(ib)%mi = isize(ib)
         grid(ib)%nj = jsize(ib);  grid(ib)%mj = jsize(ib)
         if (two_d) then
            grid(ib)%nk = 1;          grid(ib)%mk = 1       ! These may be used, but z coordinates are not read
         else
            grid(ib)%nk = ksize(ib);  grid(ib)%mk = ksize(ib)
         end if
         grid(ib)%nzoneheaderlines = nzoneheaderlines(ib)
         grid(ib)%nzoneaux         = nzoneauxlines(ib)
         grid(ib)%solutiontime     = undefined
      end do

!     Also, any dataset auxiliary info. can be reread and stored:

      if (formatted) then

         call reread_header ()  ! Local procedure below

      else ! Not implemented

!!!      read (lun, iostat=ios) (grid(ib)%ni, grid(ib)%nj, grid(ib)%nk, ib = 1, nblocks)

!!!      if (ios /= 0) then
!!!         write (*, '(2a, i4)') routine, ' Error rereading the file header.  Unit number:', lun
!!!      end if

      end if

  999 return

!     Local procedures for subroutine Tec_header_read, in the order they're used:

      contains

!        ---------------------------
         subroutine read_file_header ()
!        ---------------------------

!        Read lines up to the first 'ZONE' string, looking for optional title, variable names, and auxiliary dataset records.
!        Return with the file ready to process the first zone header line (already read).

         implicit none

!        Local variables:

         integer :: i, first_local, last_local, mark_local
         logical :: read_a_line

!        Execution:

         header%title      = blank
         names_in_header   = false
         read_a_line       = true
         nfileheaderlines  = 0
         nvar              = 0  ! Recall that the number of variables may have to be found from the first zone DT=(xxxxLE ... )
         istart            = 1  ! See pack_names
         ndatasetaux       = 0

         do ! While the first token on a line is not 'ZONE'

!!!         write (6, '(a, l2, i4)') ' read_file_header: read_a_line, nfileheaderlines =', read_a_line, nfileheaderlines

            if (read_a_line) then
               read (lun, '(a)') buffer
               lentrim = len_trim (buffer);  last = lentrim  ! SCAN2 expects last as input, and can update it
               if (lentrim == 0) then
                  nfileheaderlines = nfileheaderlines + 1
                  cycle  ! Skip a blank line
               end if
               first = 1
            else
               read_a_line = true  ! Because no other keyword is expected on the same file header line
            end if

            call scan2  (buffer(1:lentrim), ' =', first, last, mark)
            call upcase (buffer(first:mark))

!!!         write (6, '(a, 3i5, 2a)') ' read_file_header: first, mark, last =', first, mark, last, '  b(f:m): ', buffer(first:mark)

            select case (buffer(first:first))

            case ('T') ! TITLE, which could follow the keyword on the next line

               nfileheaderlines = nfileheaderlines + 1
               first = mark + 2

               if (first > lentrim) then
                  read (lun, '(a)') buffer
                  nfileheaderlines = nfileheaderlines + 1
                  first = 1;  lentrim = len_trim (buffer);  last = lentrim
               end if

               call scan4 (buffer(1:lentrim), quotes, first, last, mark)

               if (mark > 0) header%title = buffer(first:mark)  ! Omit the double quotes

               if (verbose) then
                  if (mark > 0) then
                     write (*, '(1x, 3a)') 'TITLE = "', buffer(first:mark), quotes
                  else
                     write (*, '(1x,  a)') 'TITLE = ""'
                  end if
               end if

            case ('V') ! VARIABLES

               nfileheaderlines = nfileheaderlines + 1
               names_in_header  = true
               double_quotes    = true  ! Maybe
               first = mark + 2

               if (first < lentrim) then   ! More text on this line; could be "X" or just x,y,...

                  i = index (buffer(first:lentrim), quotes)

                  if (i == 0) then ! Simple list of names on same line (as for DPLR's Postflow)
                     double_quotes = false

                     call pack_names ()   ! Ready for transfer to header%varname(*) after that array has been allocated

                     cycle  ! Next input line; no more variables
                  else
                     first = i
                  end if

                  call pack_names ()  ! For 1 or more names on current line

               end if

!              Now we need to read any remaining lines containing double-quoted variables.

               do ! Until another header keyword or ZONE is encountered

                  read (lun, '(a)') buffer
                  first = 1;  lentrim = len_trim (buffer);  last = lentrim

                  call scan2 (buffer(1:lentrim), ' =', first, last, mark)

                  if (buffer(first:first) == quotes) then
                     nfileheaderlines = nfileheaderlines + 1

                     call pack_names ()

                  else ! No more variables

                     read_a_line = false
                     exit

                  end if

               end do

            case ('D') ! DATASETAUXDATA, to be reparsed later once the total number is known and storage is allocated

               ndatasetaux      = ndatasetaux      + 1
               nfileheaderlines = nfileheaderlines + 1

            case ('#') ! Comment - skip it

               nfileheaderlines = nfileheaderlines + 1

            case ('Z') ! ZONE means end of file header lines, with first zone line already read

               exit

            case default ! Unknown keyword

               nfileheaderlines = nfileheaderlines + 1
               write (*, '(3a)') ' *** Unknown file header keyword: ', buffer(first:mark), '.  Proceeding.'

            end select

         end do ! Next file header line


!        Unpack the variable names if they've been found in the header:

         if (nvar > 0) then  ! Else they're counted and set up from the first DT = (... ) record in count_blocks

            allocate (header%varname(nvar))

            first_local = 1;  last_local = iend - 1 ! See pack_names

            do i = 1, nvar

               call scan4 (packed_names, quotes, first_local, last_local, mark_local)

               header%varname(i) = packed_names(first_local:mark_local)
               first_local = mark_local + 2
            end do

            if (verbose) then
               write (*, '(a, i4)') ' # variables found, including (x,y[,z]):', nvar
               write (*, '(a, (10(3x, a)))') ' Variable names:', (trim (header%varname(i)), i = 1, nvar)
            end if

         else ! Protect against failure by the application to indicate # dimensions:

            ndim = header%ndim
            if (ndim /= 2 .and. ndim /= 3) then
                ndim  = 3;  two_d = false;  header%ndim = ndim
            end if

         end if

         end subroutine read_file_header

!        ---------------------
         subroutine pack_names ()  ! Transfer names from buffer to next portion (istart) of packed_names
!        ---------------------

!        Allowing for embedded blanks requires (new) scan4; the standard scan2 won't work.
!        The packed names are returned like this:  "X, m" "Y, m" ["Z, m", ]"T, K" ...

         implicit none

         mark = 0

         do ! Until mark = 0

            if (double_quotes) then
               call scan4 (buffer, quotes, first, last, mark)
            else
               call scan2 (buffer,  '=, ', first, last, mark)
            end if

            if (mark == 0) exit

            nvar = nvar + 1

            packed_names(istart:istart) = quotes
            istart = istart + 1
            iend = istart + mark - first
            packed_names(istart:iend) = buffer(first:mark)
            iend = iend + 1
            packed_names(iend:iend) = quotes
            iend = iend + 1
            packed_names(iend:iend) = blank
            istart = iend + 1
            first  = mark + 2

         end do

         end subroutine pack_names

!        -----------------------
         subroutine count_blocks ()  ! Count the blocks or zones of a Tecplot ASCII file; save block dimensions in a string array.
!        -----------------------

!        Any file header records have been read.  The first line of the first zone header is in the buffer.
!        DPLR's Postflow can put datapacking & block dimensions on the same line as the zone title.  :(
!        If the number of variables is not known from names in the file header, count them from the first zone.
!
!        Some programming considerations:
!
!        The case statement can't easily work with more than the first character of a keyword.
!        Using IF ... THEN ... ELSE ... would allow strings with different lengths to be compared, but is eschewed nonetheless.
!        Trapping of unknown keywords requires precise matches (after upcasing).
!
         implicit none

!        Local variables:

         integer :: i

!        Execution:

         line = nfileheaderlines + 1;  nb = 1;  nzoneheaderlines(1) = 1

!        Always (re)start the loop over possible keywords with a keyword in hand, as when read_file_header is done.
!        The default case (first numerical value of zonal data proper) necessarily reads on until a non-numeric token (or EOF)
!        is found, so all other cases must do likewise.

         do ! Until the numeric case encounters EOF and exits; allow for more than one keyword and value on a line

!!!         write (6, '(a, i5, 2a)') ' Zone:', nb, '  keyword: ', buffer(first:mark)
!!!         write (6, '(a, 6i9)') 'count_blocks: nb, line, first, mark, last, lentrim =', nb, line, first, mark, last, lentrim

            select case (buffer(first:first))

            case ('Z')  ! ZONE or ZONETYPE

               if (mark - first == 3) then ! ZONE (start of a new zone);  zone # is incremented after reading data proper

                  if (buffer(first:mark) == 'ZONE') then
                     nzoneauxlines(nb) = 0
                  else
                     call unknown_keyword ()
                  end if

               else if (mark - first == 7) then ! ZONETYPE is assumed to be ordered

                  if (buffer(first:mark) == 'ZONETYPE') then
                     call next_token ()  ! Skip the "ordered" value
                  else
                     call unknown_keyword ()
                  end if

               else
                  call unknown_keyword ()
               end if

               call next_token ()

            case ('T')  ! Title keyword for this zone; its value may have embedded blanks; assume it's on the same line as T

               first = mark + 2

               call scan4 (buffer(1:lentrim), quotes, first, last, mark)

               if (mark < 0) mark = -mark  ! Scan4 signals a null token this way

               call next_token ()

            case ('D', 'F')  ! Data types or DATAPACKING (or F, used by Postflow to mean DATAPACKING)

               if (mark - first == 1) then  ! DT = (....LE ....LE ...LE )

                  if (buffer(first:mark) == 'DT') then

                     call scan4 (buffer(1:lentrim), '(', first, last, mark)

                     if (mark <= first) then
                        write (*, '(/, a, 2i5)') &
                           ' Missing '')''. Tecplot_io.f90 limit exceeded? lentrim, len_buffer: ', lentrim, len_buffer
                        ios = -len_buffer
                        go to 999
                     end if

                     if (nb == 1) then

                        if (.not. names_in_header) then  ! Default the variable names to "X", "Y", ["Z", ]"V1", "V2", ...

                           do i = first, mark
                              if (buffer(i:i) == 'E' .or. buffer(i:i) == 'e') nvar = nvar + 1 ! Works for SINGLE or DOUBLE
                           end do

                           allocate (header%varname(nvar))

                           header%varname(1) = 'X';  header%varname(2) = 'Y';  if (.not. two_d) header%varname(3) = 'Z'

                           do i = ndim + 1, nvar
                              header%varname(i) = blank
                              numq = i - ndim
                              if (numq < 10) then
                                 write (header%varname(i)(1:2), '(a1, i1)') 'V', numq
                              else ! Assume no more than 99 variables
                                 write (header%varname(i)(1:3), '(a1, i2)') 'V', numq
                              end if
                           end do

                           if (verbose) write (*, '(a, i4)') ' # variables found apart from (x,y,z):', numq

                        end if

                     end if

                  else
                     call unknown_keyword ()
                  end if

               else ! DATAPACKING or F[ORM?]

                  valid = false
                  if (mark - first == 10) valid = buffer(first:mark) == 'DATAPACKING'
                  if (mark - first ==  0) valid = buffer(first:mark) == 'F'

                  if (valid) then

                     call next_token ()

                     if (nb == 1) then
                        header%datapacking = 0   ! POINT
                        if (buffer(first:mark) == 'BLOCK') header%datapacking = 1
                     end if

                  else
                     call unknown_keyword ()
                  end if

               end if

               call next_token ()

            case ('I')  ! I = idim

               if (mark == first) then

                  call next_token ()
                  call char_to_integer (buffer, first, mark, isize(nb))

               else
                  call unknown_keyword ()
               end if

               call next_token ()

            case ('J')  ! J = jdim

               if (mark == first) then

                  call next_token ()
                  call char_to_integer (buffer, first, mark, jsize(nb))

               else
                  call unknown_keyword ()
               end if

               call next_token ()

            case ('K')  ! K = kdim

               if (mark == first) then

                  call next_token ()
                  call char_to_integer (buffer, first, mark, ksize(nb))

               else
                  call unknown_keyword ()
               end if

               call next_token ()

            case ('S')  ! SOLUTIONTIME = xxx

               if (mark - first == 11) then

                  if (buffer(first:mark) == 'SOLUTIONTIME') then
                     call next_token ()
!!!                  call char_to_real (buffer, first, mark), grid(nb)%solutiontime)  ! Can't do it till all blocks are allocated
                  else
                     call unknown_keyword ()
                  end if

               else
                  call unknown_keyword ()
               end if

               call next_token ()

            case ('A')  ! AUXDATA

               if (mark - first == 6) then

                  if (buffer(first:mark) == 'AUXDATA') then
                     nzoneauxlines(nb) = nzoneauxlines(nb) + 1
                     mark = last  ! Force a new line to be read - we can't store the auxiliary data yet
                  else
                     call unknown_keyword ()
                  end if

               else
                  call unknown_keyword ()
               end if

               call next_token ()

            case ('#') ! Comment - skip it

               first = mark;  last = 0;  lentrim = 0  ! Force a new line to be read

               call next_token ()

            case default ! Must be data proper for current zone (or an unknown keyword such as STRANDID):

               if (.not. number (buffer(first:mark))) then ! Non-numeric token must be an unknown keyword

                  call unknown_keyword ()  ! Ignore it, and skip its value
                  call next_token ()

               else ! A numeric token - we're through with this zone header

                  nzoneheaderlines(nb) = nzoneheaderlines(nb) - 1

                  do ! Until EOF or another zone is encountered, skip the data

                     read (lun, '(a)', iostat=ios) buffer

                     if (ios < 0) exit  ! Or go to 999 directly

                     line = line + 1

                     i = index (buffer(1:4), 'Z')

                     if (i == 0) i = index (buffer(1:4), 'z') ! Assume no leading blanks, or just a few

                     if (i  > 0) exit

                  end do

                  if (ios < 0) exit  ! EOF

!                 Prepare for the next zone:

                  first = i;  mark = i + 3;  lentrim = len_trim (buffer);  last = lentrim

                  call upcase (buffer(first:mark))

                  nb = nb + 1

                  if (nb > max_blocks) then
                     nb  = max_blocks
                     write (*, '(3a, i6)') routine, 'Too many zones in file ', trim (header%filename), '.  Limit:', max_blocks
                     exit
                  end if

                  nzoneheaderlines(nb) = 1  ! Ready for incrementing

               end if

            end select

         end do  ! Process next keyword

  999    continue

         nblocks        = nb
         header%nblocks = nb
         if (ios /= -len_buffer) ios = 0

         if (verbose) write (*, '(a, i5)') ' # zones found:', nb, ' # variables:  ', nvar

         end subroutine count_blocks

!        ---------------------
         subroutine next_token ()  ! Intended for processing zone headers only.  It spans input lines if necessary.
!        ---------------------

!        Use this if converting to upper case is appropriate (not for zone titles or auxiliary data values).
!        Global variables first, last, lentrim, and mark are assumed on input to apply to the previous token.
!        Variable "nb" points to the current zone.  Variable "line" refers to the entire file (for error messages only).
!        Upon return, the next token is in buffer(first:mark), in upper case.
!        EOF is NOT considered a possibility here.

         implicit none

         integer :: first_local, last_local, mark_local

         first = mark + 2

         if (first > last) then

            lentrim = 0
            do while (lentrim == 0)  ! Skip any blank lines
               read (lun, '(a)') buffer
               line = line + 1;  nzoneheaderlines(nb) = nzoneheaderlines(nb) + 1
               lentrim = len_trim (buffer)
            end do

            first = 1;  last = lentrim

            if (verbose) then  ! Avoid echoing the first data line of each zone
               first_local = 1;  last_local = lentrim
               call scan2 (buffer(1:lentrim), blank, first_local, last_local, mark_local)
               if (.not. number (buffer(first_local:mark_local))) write (*, '(a)') buffer(1:lentrim)
            end if

         end if

         call scan2  (buffer(1:lentrim), ' =,', first, last, mark)
         call upcase (buffer(first:mark))

         end subroutine next_token

!        --------------------------
         subroutine unknown_keyword ()
!        --------------------------

         implicit none

         write (*, '(3a)') ' *** Unknown keyword: ', buffer(first:mark), '.  Proceeding.'

         call next_token () ! Get the value presumably associated with the unknown keyword, and ignore it

         end subroutine unknown_keyword

!        --------------------------
         subroutine char_to_integer (text, i1, i2, number)  ! Internal read from text(i1:i2) to integer number
!        --------------------------

         implicit none

         character, intent (in)  :: text * (*)
         integer,   intent (in)  :: i1, i2
         integer,   intent (out) :: number

         character :: format_string * 4

         format_string = '(i?)'
         write (format_string(3:3), '(i1)') i2 - i1 + 1
         read  (text(i1:i2), format_string) number

         end subroutine char_to_integer

!        ------------------------
         subroutine reread_header ()  ! Store any auxiliary dataset info. in header; leave the file ready to read the first zone.
!        ------------------------

         implicit none

         integer :: line, naux

         header%ndatasetaux = ndatasetaux

         if (ndatasetaux > 0) then
            allocate (header%datasetauxname(ndatasetaux), header%datasetauxvalue(ndatasetaux))
         end if

         rewind (lun);  naux = 0

         if (verbose) write (*, '(/, a, /)') ' Rereading input file header:'

         do line = 1, nfileheaderlines

            read (lun, '(a)') buffer
            if (verbose) write (*, '(a)') trim (buffer)

            if (naux < ndatasetaux) then

               first = 1;  lentrim = len_trim (buffer);  last = lentrim

               call scan2  (buffer(1:lentrim), blank, first, last, mark)
               call upcase (buffer(first:mark))

               if (buffer(first:mark) == 'DATASETAUXDATA') then

                  naux  = naux + 1
                  first = mark + 2

                  call scan2 (buffer(1:last), ' =', first, last, mark)

                  header%datasetauxname(naux) = buffer(first:mark)
                  if (verbose) write (*, '(a, a)') ' Dataset auxiliary name:  ', header%datasetauxname(naux)

                  first = mark + 2

                  call scan4 (buffer(1:last), quotes, first, last, mark)

                  header%datasetauxvalue(naux) = buffer(first:mark)
                  if (verbose) write (*, '(a, a)') ' Dataset auxiliary value: ', header%datasetauxvalue(naux)

               end if

            end if

         end do

         end subroutine reread_header

      end subroutine Tec_header_read

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine Tec_header_write (lun, header, ios)
!
!     Open a Tecplot data file and write the header records for 2- or 3-space structured Tecplot data.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Arguments:

      integer,            intent (in)  :: lun     ! Logical unit number, used if the file is ASCII
      type (grid_header), intent (in)  :: header  ! Data structure containing dataset header information
      integer,            intent (out) :: ios     ! 0 = no error

!     Local constants:

      character, parameter :: colon * 1 = ':', comma * 1 = ',', routine * 19 = ' Tec_header_write: '

!     Local variables:

      integer :: k, l, m, n, naux, nvars


!     Execution:

      verbose = ios == 1
      ndim    = header%ndim
      numq    = header%numq
      two_d   = ndim == 2
      nvars   = ndim + numq

      if (header%formatted) then

         open (lun, file=trim (header%filename), status='unknown', iostat=ios)

         if (ios /= 0) then
            write (*, '(3a)') routine, 'Trouble opening file ', trim (header%filename)
            go to 999
         end if

         write (lun, '(4a)') 'TITLE = ', quotes, trim (header%title), quotes, 'VARIABLES ='
         write (lun, '(3a)') (quotes, trim (header%varname(n)), quotes, n = 1, nvars)

         do l = 1, header%ndatasetaux
            write (lun, '(5a)') 'DATASETAUXDATA ', trim (header%datasetauxname(l)), '="',  trim (header%datasetauxvalue(l)), quotes
         end do

      else ! Use the Tecplot I/O function for initializing writing to a binary file

         l = 0
         do n = 1, nvars
            l = len_trim (header%varname(n)) + l
         end do
         l = l + nvars + 1 ! For the commas and the final null

         if (l > max_length) then
            write (*, '(2a, i5)') routine, 'Packed variable names are too long.  Length required:', l
            ios = 999
            go to 999
         end if

         m = 0
         do n = 1, nvars
            l = len_trim (header%varname(n))
            packed_names(m+1:m+l) = header%varname(n)(1:l);  l = m + l + 1
            do k = m+1, m+l
               if (packed_names(k:k) == comma) packed_names(k:k) = colon  ! Double quotes still don't allow embedded commas
            end do
            packed_names(l:l) = comma;              m = l
         end do
         l = m + 1
         packed_names(l:l) = null

!!!      write (*, '(a, i6)') ' Packed names(1:l) passed to TecIni110, with l =', l
!!!      write (*, '(a)') packed_names(1:l)

         eps = epsilon (eps)
         if (eps < 1.e-10) then
            IsDouble = 1
         else
            IsDouble = 0
         end if

         ios = TecIni110 (trim (header%title) // null, packed_names(1:l),       &
                          trim (header%filename) // null, '.' // null, 0, IsDouble)    ! 0 = No debug

         if (ios /= 0) then
            write (*, '(3a)') routine, 'Trouble initializing Tecplot file ', trim (header%filename)
         end if

         do l = 1, header%ndatasetaux
            ios = TecAuxStr110 (trim (header%datasetauxname(l)) // null,  trim (header%datasetauxvalue(l)) // null)
            if (ios /= 0) then
               write (*, '(2a)') routine, 'Trouble writing dataset auxiliary data:'
               write (*, '(5a)') 'DATASETAUXDATA ', trim (header%datasetauxname(l)), '="',  trim (header%datasetauxvalue(l)), quotes
            end if
         end do

      end if

  999 return

      end subroutine Tec_header_write

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine Tec_block_allocate (block, ndim, numq, ios)
!
!     Allocate the x, y[, z] arrays, optional q array, and any auxiliary data for one block or zone.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Arguments:

      type (grid_type), intent (inout) :: block
      integer,          intent (in)    :: ndim
      integer,          intent (in)    :: numq
      integer,          intent (out)   :: ios

!     Local constants:

      character, parameter :: routine_etc * 33 = ' Tec_block_allocate trouble with '

!     Local variables:

      integer :: naux, ni, nj, nk

!     Execution:

      ni = block%ni;  nj = block%nj

      if (ndim == 2) then
         nk = 1
         allocate (block%x(ni,nj,1), block%y(ni,nj,1), stat=ios)
      else
         nk = block%nk
         allocate (block%x(ni,nj,nk), block%y(ni,nj,nk), block%z(ni,nj,nk), stat=ios)
      end if

      if (ios /= 0) then
         write (*, '(2a)') routine_etc, 'x, y[, z].'
         go to 999
      end if

      if (numq > 0) then

         allocate (block%q(numq,ni,nj,nk), stat=ios)

      else ! Avoid undefined arguments elsewhere

         allocate (block%q(1,1,1,1), stat=ios)

      end if

      if (ios /= 0) then
         write (*, '(2a)') routine_etc, 'q array.'
         go to 999
      end if

      block%solutiontime = undefined

      naux = block%nzoneaux

      if (naux > 0) then

         allocate (block%zoneauxname(naux), block%zoneauxvalue(naux), stat=ios)

         if (ios /= 0) write (*, '(2a)') routine_etc, 'zone auxiliary info.'

      end if

  999 if (ios /= 0) write (*, '(a, 5i5)') '  ni, nj, nk, numq, nzoneaux:', ni, nj, nk, numq, naux

      end subroutine Tec_block_allocate

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine Tec_block_read (lun, header, block, ios)
!
!     Read one 2- or 3-space block (zone) of structured data from a Tecplot file (ASCII, POINT or BLOCK packing order).
!     Binary reads have yet to be implemented.
!     The block dimensions have already been determined as part of reading the file header, and array fields for this block
!     should have been allocated via Tec_block_allocate.
!
!     A Tecplot 360 ASCII zone looks like something like the following.
!
!        ZONE T="Zone 1 title"                   [optional zone 1 title]
!        ZONETYPE=Ordered                        [understood]
!        I=idim, J=jdim, K=kdim                  [commas optional, all lines]
!        DATAPACKING=BLOCK                       [or POINT; also: F=POINT]
!        DT=(SINGLE SINGLE .... )                [or (DOUBLE DOUBLE ... )]
!        SOLUTIONTIME=time                       [optional real number]
!        AUXDATA name="value"                    [0+ per zone, 1 per line]
!        xxx.xx  xxx.xx                          [Data proper]
!        ::::::::::::::                          ::::::  ::::::
!        ZONE T="Zone 2 title"                   [Repeat for further zones]
!
!     Alternatively, output from DPLR's Postflow has single-line zone headers that look something like this:
!
!      zone t="flow3d" F=point, i=  97 j=  33 k= 129
!
!     The keywords can be in any order, possibly with a few leading blanks or blank lines.
!     The file is presumed to be positioned ready to read a 'zone' keyword as the first token.
!     Since Tec_header_read has determined the number of lines per zone header, the logic is slightly different from
!     that for the count_blocks procedure.  The loop over keywords here starts with finding a keyword as opposed to
!     already having encountered one, because we don't want to buffer in the first line of numeric data when reparsing
!     the zone header.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Arguments:

      integer,            intent (in)    :: lun     ! Logical unit number
      type (grid_header), intent (inout) :: header  ! Dataset header with ndim, etc.
      type (grid_type),   intent (inout) :: block   ! Block to be read, with %ni/nj/nk/naux set & other fields allocated
      integer,            intent (out)   :: ios     ! 0 = no error

!     Local variables:

      integer   :: first, last, lentrim, line, mark, naux  ! Shared by internal procedures
      integer   :: i, j, k, n, ni, nj, nk

!     Execution:

      two_d = header%ndim == 2
      numq  = header%numq

      if (header%formatted) then

         call read_zone_header ()

      else

!        Not implemented

      end if

!     Read the numeric data for this zone:

      ni = block%ni;  nj = block%nj
      if (two_d) then
         nk = 1  ! For the block-order function data read
      else
         nk = block%nk
      end if

      select case (header%datapacking)

         case (0) ! POINT

            if (header%formatted) then

               if (numq == 0) then
                  if (two_d) then
                     read (lun, *, iostat=ios)  ((block%x(i,j,1), block%y(i,j,1), i = 1, ni), j = 1, nj)
                  else
                     read (lun, *, iostat=ios) (((block%x(i,j,k), block%y(i,j,k), block%z(i,j,k), i = 1, ni), j = 1, nj), k = 1, nk)
                  end if
               else
                  if (two_d) then
                     read (lun, *, iostat=ios) &
                         ((block%x(i,j,1), block%y(i,j,1), block%q(:,i,j,1), i = 1, ni), j = 1, nj)
                  else
                     read (lun, *, iostat=ios) &
                        (((block%x(i,j,k), block%y(i,j,k), block%z(i,j,k), block%q(:,i,j,k), i = 1, ni), j = 1, nj), k = 1, nk)
                  end if
               end if

               if (ios /= 0) then
                  write (*, '(/, a, i4)') ' Trouble reading data in POINT order.  Logical unit:', lun
                  go to 99
               end if

            else ! Not implemented yet
!!!            if (numq == 0) then
!!!               read (lun, iostat=ios) (((x(i,j,k), y(i,j,k), z(i,j,k), i = 1, ni), j = 1, nj), k = 1, nk)
!!!            else
!!!               read (lun, iostat=ios) (((x(i,j,k), y(i,j,k), z(i,j,k), q(:,i,j,k), i = 1, ni), j = 1, nj), k = 1, nk)
!!!            end if
            end if

         case (1) ! BLOCK

            if (header%formatted) then

               if (two_d) then
                  if (numq == 0) then
                     read (lun, *, iostat=ios) block%x, block%y
                  else
                     read (lun, *, iostat=ios) block%x, block%y, (block%q(n,1:ni,1:nj,1), n = 1, numq)
                  end if
               else
                  if (numq == 0) then
                     read (lun, *, iostat=ios) block%x, block%y, block%z
                  else
                     read (lun, *, iostat=ios) block%x, block%y, block%z, (block%q(n,1:ni,1:nj,1:nk), n = 1, numq)
                  end if
               end if

               if (ios /= 0) then
                  write (*, '(/, a, 3i4)') ' Trouble reading data in BLOCK order.  numq, var # and lun:', numq, n, lun
                  go to 99
               end if

            else ! Not implemented yet
!!!            read (lun, iostat=ios) x, y, z
!!!            if (numq > 0) read (lun, iostat=ios) q
            end if

         case default

            ios = 1

      end select

   99 return

!     Local procedures for subroutine Tec_block_read:

      contains

!        ------------------------------
         subroutine read_zone_header ()  ! Read the current zone header with known number of header lines and auxiliaries
!        ------------------------------

         implicit none

!        Execution:

         line = 0;  mark = 0;  last = 0;  naux = 0

!        Always (re)start the loop over possible keywords by locating the next keyword if any (in contrast to count_blocks).

         do ! Until all header lines for this zone have been fully processed - see next_zone_token.

            call next_zone_token ()

            select case (buffer(first:first))

            case ('Z')  ! ZONE or ZONETYPE

               if (mark - first == 3) then ! ZONE (start of a new zone);  zone # is incremented after reading data proper

                  if (buffer(first:mark) == 'ZONE') then
!                    Nothing to do
                  else  ! Unknown keyword: ignore it, and skip its value as well
                     call next_zone_token ()
                  end if

               else  ! ZONETYPE is assumed to be ordered

                  call next_zone_token ()  ! Skip the "ordered" value;  same action if it's an unknown keyword

               end if

            case ('T')  ! Title keyword for this zone; its value may have embedded blanks; assume it's on the same line as T

               if (first == mark) then

                  first = mark + 2;  call scan4 (buffer(1:lentrim), quotes, first, last, mark)

                  if (mark < 0) then ! Blank title signaled with ""
                     mark = -mark
                     block%zone_title = blank
                  else
                     block%zone_title = buffer(first:mark)
                  end if

               else  ! Unknown keyword - ignore it, and skip its value as well
                  call next_zone_token ()
               end if

            case ('D', 'F')  ! Data types or DATAPACKING (or F, used by Postflow to mean DATAPACKING)

               if (mark - first == 1) then  ! DT = (....LE ....LE ...LE )

                  if (buffer(first:mark) == 'DT') then

                     call scan4 (buffer(1:lentrim), '(', first, last, mark)

                  else  ! Unknown keyword and value to ignore
                     call next_zone_token ()
                  end if

               else  ! DATAPACKING or F

                  call next_zone_token ()  ! BLOCK or POINT is already known;  same action for an unknown keyword

               end if

            case ('I', 'J', 'K')  ! I = idim, J = jdim, or K = kdim are already known

               call next_zone_token ()  ! Same action for an unknown keyword

            case ('S')  ! SOLUTIONTIME = xxx

               if (mark - first == 11) then
                  if (buffer(first:mark) == 'SOLUTIONTIME') then
                     call next_zone_token ()
                     call char_to_real (buffer, first, mark, block%solutiontime)
                  else
                     call next_zone_token ()  ! Unknown keyword - ignore it, and skip its value
                  end if
               else
                  call next_zone_token ()  ! Unknown keyword - ignore it, and skip its value
               end if

            case ('A')  ! AUXDATA

               if (mark - first == 6) then
                  if (buffer(first:mark) == 'AUXDATA') then
                     first = mark + 2                     ! Avoid upcasing the name
                     call scan2 (buffer(1:lentrim), ' =,',  first, last, mark)
                     naux  = naux + 1;  block%zoneauxname(naux) = buffer(first:mark)
                     first = mark + 2
                     call scan4 (buffer(1:lentrim), quotes, first, last, mark);  block%zoneauxvalue(naux) = buffer(first:mark)
                  else
                     call next_zone_token ()  ! Unknown keyword - ignore it, and skip its value
                  end if
               else
                  call next_zone_token ()  ! Unknown keyword - ignore it, and skip its value
               end if

            case ('E')  ! End of zone header (see next_zone_token)

               if (mark - first == 2) then
                  if (buffer(first:mark) == 'END') exit  ! Done with this zone header
               else
                  call next_zone_token ()  ! Unknown keyword - ignore, and skip its value
               end if

            case default ! Must be an unknown keyword such as STRANDID:

               call next_zone_token ()  ! Get and skip its value and keep going

            end select

         end do  ! Process next keyword

         end subroutine read_zone_header

!        --------------------------
         subroutine next_zone_token ()  ! This is a variation of "next_token" used by count_blocks
!        --------------------------

!        Use this if converting to upper case is appropriate (not for zone titles or auxiliary data values).
!        Global variables first, last, lentrim, and mark are assumed on input to apply to the previous token.
!        Upon return, the next token is in buffer(first:mark), in upper case.
!        It may have been fudged as 'END' if the known number of header lines has been processed.
!        EOF is NOT considered a possibility here.

         implicit none

         first = mark + 2

!!!      write (6, '(a, 6i6)') ' next_zone_token: nzh, line, first, mark, last, lentrim:', &
!!!                            block%nzoneheaderlines, line, first, mark, last, lentrim

         if (first > last) then

            if (line < block%nzoneheaderlines) then  ! Get a new zone header line

               lentrim = 0

               do while (lentrim == 0)  ! Skip any blank lines
                  read (lun, '(a)') buffer
                  line = line + 1
                  lentrim = len_trim (buffer)
                  if (line == block%nzoneheaderlines) exit
               end do

               first = 1;  last = lentrim  ! Handling a blank last zone header line is awkward

            end if

            if (line == block%nzoneheaderlines .and. first > last) then  ! Done with this zone header

               buffer(1:3) = 'END'
               first = 1;  mark = 3;  last = 0

            end if

         end if

         if (first <= last) then

            call scan2  (buffer(1:lentrim), ' =,', first, last, mark)
            call upcase (buffer(first:mark))

         end if

         end subroutine next_zone_token

!        -----------------------
         subroutine char_to_real (text, i1, i2, real_number)  ! Internal read from text(i1:i2) to real number
!        -----------------------

         implicit none

         character, intent (in)  :: text * (*)
         integer,   intent (in)  :: i1, i2
         real,      intent (out) :: real_number

         integer   :: n
         character :: format_string * 7

         format_string = '(???.0)'
         n = i2 - i1 + 1
         if (n < 10) then
            format_string(2:3) = ' f';  write (format_string(4:4), '(i1)') n
         else ! n < 100
            format_string(2:2) = 'f';   write (format_string(3:4), '(i2)') n
         end if
         read  (text(i1:i2), format_string) real_number

         end subroutine char_to_real

      end subroutine Tec_block_read

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine Tec_block_write (lun, header, block, ios)
!
!     Write one 2- or 3-space block of structured data from a Tecplot file (binary or ASCII, POINT or BLOCK packing order).
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Arguments:

      integer,            intent (in)  :: lun     ! Logical unit number used if the file is ASCII
      type (grid_header), intent (in)  :: header  ! Dataset header with ndim, etc.
      type (grid_type),   intent (in)  :: block   ! Block to be written
      integer,            intent (out) :: ios     ! 0 = no error

!     Local constants:

      integer,   parameter :: ijk_limit = 5  ! Limit on zone dimensions
      character, parameter :: routine * 18 = ' Tec_block_write: '

!     Local variables:

      integer :: i, j, k, li, lj, lk, n, naux, ni, nj, nk, npts, nvars
      integer :: ZoneType, StrandID, ParentZone, IsBlock, NumFaceConnections, FaceNeighborMode
      real    :: SolutionTime

      character  :: aformat*6, eformat*14
      character (ijk_limit) :: istring, jstring, kstring

      integer,   allocatable, dimension (:)     :: PassiveVarList, ValueLocation, ShareVarFromZone, ShareConnectivityFromZone
      real,      allocatable, dimension (:)     :: one_point
      real,      allocatable, dimension (:,:,:) :: one_var
      character, allocatable, dimension (:)     :: chars * 7

!     Execution:

      ndim  = header%ndim
      numq  = header%numq
      two_d = ndim == 2
      nvars = ndim + numq
      ni = block%ni;  nj = block%nj;  naux = block%nzoneaux

      if (two_d) then
         nk = 1
      else
         nk = block%nk
      end if

      npts = ni * nj * nk

      if (header%formatted) then

         write (lun, '(4a)') 'ZONE T=', quotes, trim (block%zone_title), quotes
         call int_to_char (ni, istring, li)  ! Convert ni to string(1:li)
         call int_to_char (nj, jstring, lj)
         if (two_d) then
            write (lun, '(5a)') ' I=', istring(1:li), ', J=', jstring(1:lj),                        ', ZONETYPE=Ordered'
         else
            call int_to_char (nk, kstring, lk)
            write (lun, '(7a)') ' I=', istring(1:li), ', J=', jstring(1:lj), ', K=', kstring(1:lk), ', ZONETYPE=Ordered'
         end if

         if (block%solutiontime /= undefined) write (lun, '(a, f12.6)') 'SOLUTIONTIME=', block%solutiontime

         do n = 1, naux
            write (lun, '(5a)') 'AUXDATA ', trim (block%zoneauxname(n)), '="', trim (block%zoneauxvalue(n)), quotes
         end do

         aformat = '(nnna)'
         write (aformat(2:4), '(i3)') nvars + 2

         allocate (chars(nvars))

         eps = epsilon (eps)
         if (eps < 1.e-10) then
            IsDouble = 1
         else
            IsDouble = 0
         end if

         if (IsDouble == 1) then
            chars(:) = 'DOUBLE '
            eformat = '(1p, nnne17.9)'
         else ! 0
            chars(:) = 'SINGLE '
            eformat = '(1p, nnne15.7)'
         end if

         write (eformat(6:8), '(i3)') nvars

         select case (header%datapacking)

         case (0) ! POINT

            write (lun, '(a)') ' DATAPACKING=POINT'
            write (lun, aformat) ' DT=(', chars, ')'

            allocate (one_point(nvars))

            if (two_d) then

                  do j = 1, nj
                     do i = 1, ni
                        one_point(1) = block%x(i,j,1)
                        one_point(2) = block%y(i,j,1)
                        if (nvars > 2) one_point(3:nvars) = block%q(:,i,j,1)
                        write (lun, eformat) one_point
                     end do
                  end do

            else

               do k = 1, nk
                  do j = 1, nj
                     do i = 1, ni
                        one_point(1) = block%x(i,j,k)
                        one_point(2) = block%y(i,j,k)
                        one_point(3) = block%z(i,j,k)
                        if (nvars > 3) one_point(4:nvars) = block%q(:,i,j,k)
                        write (lun, eformat) one_point
                     end do
                  end do
               end do

            end if

            deallocate (one_point)

         case (1) ! BLOCK

            write (lun, '(a)') ' DATAPACKING=BLOCK'
            write (lun, aformat) ' DT=(', chars, ')'
            write (lun, '(1p, 6e17.9)') block%x
            write (lun, '(1p, 6e17.9)') block%y
            if (.not. two_d) write (lun, '(1p, 6e17.9)') block%z

            if (nvars > ndim) then
               do n = 1, numq
                  write (lun, '(1p, 6e17.9)') block%q(n,1:ni,1:nj,1:nk)
               end do
            end if

         case default

            ios = 999

         end select

         deallocate (chars)

      else ! Binary output

         ZoneType           = 0
         SolutionTime       = block%solutiontime
         StrandID           = 0
         ParentZone         = 0
         IsBlock            = header%datapacking
         NumFaceConnections = 0
         FaceNeighborMode   = 0

         allocate (PassiveVarList(nvars), ValueLocation(nvars), ShareVarFromZone(nvars), ShareConnectivityFromZone(nvars))

         do n = 1, nvars
            PassiveVarList(n)            = 0
            ValueLocation(n)             = 1  ! Vertex-centered
            ShareVarFromZone(n)          = 0
            ShareConnectivityFromZone(n) = 0
         end do

         ios = TecZne110 (trim (block%zone_title) // null, &
                          ZoneType,                        &
                          ni, nj, nk,                      &
                          ni, nj, nk,                      &             ! Not used
                          SolutionTime,                    &
                          StrandID,                        &
                          ParentZone,                      &
                          IsBlock,                         &
                          NumFaceConnections,              &
                          FaceNeighborMode,                &
                          PassiveVarList,                  &
                          ValueLocation,                   &
                          ShareVarFromZone,                &
                          ShareConnectivityFromZone)

         deallocate (PassiveVarList, ValueLocation, ShareVarFromZone, ShareConnectivityFromZone)

         if (ios /= 0) then
            write (*, '(/, 2a, i4)') routine, 'Trouble writing zone info.  Logical unit:', lun
            go to 999
         end if

         do n = 1, naux
            ios = TecZAuxStr110 (trim (block%zoneauxname(n)) // null,  trim (block%zoneauxvalue(n)) // null)
            if (ios /= 0) then
               write (*, '(2a, i5)') routine, 'Trouble writing zone auxiliary data.  Zone #:', n
               write (*, '(5a)') 'AUXDATA ', trim (block%zoneauxname(n)), '="',  trim (block%zoneauxvalue(n)), quotes
            end if
         end do

         select case (header%datapacking)

         case (0) ! POINT

            allocate (one_point(nvars))

            if (two_d) then

                  do j = 1, nj
                     do i = 1, ni
                        one_point(1) = block%x(i,j,1)
                        one_point(2) = block%y(i,j,1)
                        if (nvars > 2) one_point(3:nvars) = block%q(:,i,j,1)
                        ios = TecDat110 (nvars, one_point, IsDouble)
                     end do
                  end do

            else

               do k = 1, nk
                  do j = 1, nj
                     do i = 1, ni
                        one_point(1) = block%x(i,j,k)
                        one_point(2) = block%y(i,j,k)
                        one_point(3) = block%z(i,j,k)
                        if (nvars > 3) one_point(4:nvars) = block%q(:,i,j,k)
                        ios = TecDat110 (nvars, one_point, IsDouble)
                     end do
                  end do
               end do

            end if

            deallocate (one_point)

         case (1) ! BLOCK

            ios = TecDat110 (npts, block%x, IsDouble)  ! Assume x(:,:,:) has been dynamically allocated to fit the data
            ios = TecDat110 (npts, block%y, IsDouble)
            if (.not. two_d) &
            ios = TecDat110 (npts, block%z, IsDouble)

            if (nvars > ndim) then
               allocate (one_var(ni,nj,nk))
               do n = 1, numq
                  one_var = block%q(n,1:ni,1:nj,1:nk)
                  ios = TecDat110 (npts, one_var, IsDouble)
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

      subroutine clone_header (header1, numq, datapacking, header2)

!     Derive one dataset header from another, possibly with a different number of functions.
!     The dataset file name and formatting are not touched, though, as they may be read directly.
!     The dataset title (copied here) may need adjustment upon return.
!     The first ndim + numq variable names are transcribed but may also need adjustment by the application.
!     Any dataset auxiliaries are allocated and transcribed here.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Arguments:

      type (grid_header), intent (in)  :: header1      ! Dataset header to be cloned
      integer,            intent (in)  :: numq         ! Desired number of functions in addition to x, y[, z]
      integer,            intent (in)  :: datapacking  ! Desired ordering: 0 => POINT; 1 => BLOCK
      type (grid_header), intent (out) :: header2      ! Cloned dataset header

!     Local variables:

      integer :: naux, nvars, nvcopy

!     Execution:

      header2%ndim        = header1%ndim
      nvars               = header1%ndim + numq
      nvcopy              = min (nvars, header1%ndim + header1%numq)
      header2%numq        = numq
      header2%nblocks     = header1%nblocks
      header2%datapacking = datapacking
      header2%title       = header1%title   ! May need adjusting by the application

      allocate (header2%varname(nvars))

      header2%varname(1:nvcopy) = header1%varname(1:nvcopy)

      naux = header1%ndatasetaux;  header2%ndatasetaux = naux

      if (naux > 0) then
         allocate (header2%datasetauxname(naux), header2%datasetauxvalue(naux))
         header2%datasetauxname (:) = header1%datasetauxname (:)
         header2%datasetauxvalue(:) = header1%datasetauxvalue(:)
      end if

      end subroutine clone_header

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine clone_zone (block1, ndim, numq, block2, ios)

!     Derive one zone (block) header from another.
!     Any zone auxiliaries are allocated and transcribed here.
!     The variable arrays are also allocated but not copied.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Arguments:

      type (grid_type), intent (in)  :: block1  ! Dataset header to be cloned
      integer,          intent (in)  :: ndim    ! Needed for Tec_block_allocate
      integer,          intent (in)  :: numq    ! Desired number of functions in addition to x, y[, z]
      type (grid_type), intent (out) :: block2  ! Cloned dataset header
      integer,          intent (out) :: ios     ! 0 => no error detected

!     Local variables:

      integer :: naux

!     Execution:

      block2%ni = block1%ni
      block2%nj = block1%nj
      block2%nk = block1%nk
      block2%mi = block1%mi
      block2%mj = block1%mj
      block2%mk = block1%mk

      block2%zone_title = block1%zone_title

      naux = block1%nzoneaux;  block2%nzoneaux = naux

!     Allocate x, y[, z], q if any, and zone auxiliaries if any:

      call Tec_block_allocate (block2, ndim, numq, ios)

      if (ios == 0) then
         block2%solutiontime = block1%solutiontime  ! Not earlier, because Tec_block_allocate initializes it as undefined
         if (naux > 0) then
            block2%zoneauxname (:) = block1%zoneauxname (:)
            block2%zoneauxvalue(:) = block1%zoneauxvalue(:)
         end if
      end if

      end subroutine clone_zone

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine deallocate_header (header)

!     Deallocate any dataset header auxiliaries.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Argument:

      type (grid_header), intent (inout) :: header

!     Execution:

      if (header%ndatasetaux > 0) deallocate (header%datasetauxname, header%datasetauxvalue)

      end subroutine deallocate_header

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine deallocate_blocks (ib1, ib2, ndim, numq, grid, ios)

!     Deallocate the indicated blocks.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Arguments:

      integer, intent (in)  :: ib1, ib2              ! Block range to deallocate
      integer, intent (in)  :: ndim                  ! > 0 means deallocate %z too
      integer, intent (in)  :: numq                  ! > 0 means deallocate %q too
      type (grid_type), intent (inout) :: grid(ib2)  ! At least this many
      integer, intent (out) :: ios                   ! 0 means no problem

!     Local variables:

      integer :: ib

!     Execution:

      do ib = ib1, ib2
         if (ndim == 2) then
            deallocate (grid(ib)%x, grid(ib)%y,             stat=ios)
         else
            deallocate (grid(ib)%x, grid(ib)%y, grid(ib)%z, stat=ios)
         end if
         if (ios /= 0) go to 90
         if (numq > 0) then
            deallocate (grid(ib)%q, stat=ios)
            if (ios /= 0) go to 90
         end if
         if (grid(ib)%nzoneaux > 0) deallocate (grid(ib)%zoneauxname, grid(ib)%zoneauxvalue)
      end do

      go to 99

   90 write (*, '(/, a, 3i5)') ' Trouble deallocating grid_type variable. ib1, ib2, nf: ', ib1, ib2, numq

   99 return

      end subroutine deallocate_blocks

   end module Tecplot_io_module
