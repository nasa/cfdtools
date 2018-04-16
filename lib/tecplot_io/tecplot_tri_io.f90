!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module Tecplot_tri_io

!  This package modularizes certain I/O for triangulated surface data in Tecplot
!  format (3-space only).  See also the earlier Tecplot_io.f90 module, which is
!  restricted to structured (ordered) data.  The two should be kept compatible
!  in their use of the variable name length limit and other constants.
!
!  Restrictions:
!
!    o  The Tecplot "zone" concept is employed here to contain one surface.
!       dataset (such as one surface solution for one set of flight conditions).
!    o  The file may contain more than one such zone, but ...
!    o  Just a single zone or surface dataset is read or written at one time,
!       with all zones containing the same variables.
!    o  Variable names should all appear on the same line.
!    o  Binary reads are not an option; binary or ASCII writes are supported.
!    o  Single and double precision data (-r8 compile or not) are both handled.
!    o  Tecplot ASCII files of triangulated data look like this:
!
!     TITLE = "This line may be omitted."
!     VARIABLES = "X", "Y", "Z", "TEMP", "PRESS", "CP", "MACHE", "ASOUNDE", ...
!     ZONE [T=optional zone title, ]N=96000, E=32000, F=FEPOINT, ET=TRIANGLE
!     0.000000  0.000000 0.000000 2838.386719 51330.925781 1.552663 0.609412 ...
!     0.000883 -0.007150 0.113643 2838.386719 51330.925781 1.552663 0.609412 ...
!     0.000883  0.000000 0.113868 2838.386719 51330.925781 1.552663 0.609412 ...
!     ::::::::::::::::
!     ::::::::::::::::
!     4.882953  0.000000 0.011285 950.867676 16.506409 -0.001166 5.062649 ...
!     1 2 3
!     4 5 6
!     7 8 9
!     10 11 12
!     ::::::::
!     ::::::::
!     95992 95993 95994
!     95995 95996 95997
!     95998 95999 96000
!    [ZONE ...] [Optional further dataset, including connectivity.]
!
!  History:
!
!     03/01/05  D. Saunders  Initial implementation of Tecplot_read_tri_data.
!     04/20/05    "     "    Initial attempt to adapt as an I/O package that
!                            reads or writes one surface dataset at a time.
!     08/25/05    "     "    Variables last & l were not reset when a title line
!                            was read before the variables line.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   private

!  Constants used by this package (consistent with Tecplot_io.f90):

   integer,   parameter :: max_zones  = 200       ! Limit for a file being read
   integer,   parameter :: max_length = 1000      ! For the string into which
                                                  ! variable names are packed
   integer,   parameter :: name_limit = 32        ! Variable name length limit

   character, parameter :: null * 1   = char (0), &
                           quotes * 1 = '"'

!  Awkward variables used by the package:

   integer   :: isdouble                      ! 0 = single precision data (-r4);
                                              ! 1 = double precision data (-r8)
   real      :: eps                           ! epsilon (eps) allows isdouble to
                                              ! be set for binary reads & writes
   character :: packed_names * (max_length)   ! F90 doesn't support dynamically
                                              ! variable string lengths

!  Parsing utility needed from another package:

   external     scan2

!  Utilities provided for applications:

   public :: Tecplot_tri_read   ! Reads  1 surface dataset (+ header on call 1)
   public :: Tecplot_tri_write  ! Writes "     "     "     "     "     "     "

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Tecplot_tri_read (num_call, lun, filename, formatted,         &
                                   geom_title, variable_names, zone_title,     &
                                   nnode, ntri, nf, tri_xyz, tri_f, conn,      &
                                   iostatus)

!     Read the next zone from a Tecplot surface triangulation with associated
!     function values (0 or more beyond the initial "x", "y", "z" inputs).
!     On the first call, open the file and read an optional geometry title and
!     the names of the variables.  There is presently no option to count the
!     number of zones; arguments are intended to handle one zone per call.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer,   intent (in)  :: num_call ! 1 = first call for this file; read
                                          ! an optional title + variable names;
                                          ! 2+ = read the indicated next zone
      integer,   intent (in)  :: lun      ! Logical unit for the input data file
                                          ! which is opened here if num_call = 1
      character, intent (in)  :: filename * (*) ! The binary option has to do
                                          ! the open here, so do it regardless
      logical,   intent (in)  :: formatted! T = ASCII; F = binary (incomplete)
      character, intent (out) :: geom_title * (*) ! Optional geom. title in file
      character, pointer ::   &
        variable_names (:) * (name_limit) ! Variable names, left-justified;
                                          ! array size will be = 3 + nf
      character, intent (out) :: zone_title * (*) ! Optional zone title in file
      integer,   intent (out) :: nnode    ! # surface points found
      integer,   intent (out) :: ntri     ! # triangles found
      integer,   intent (out) :: nf       ! # functions found beyond x,y,z
      real,      pointer :: tri_xyz(:,:)  ! tri_xyz(1:3,*) contains the (x,y,z)
                                          ! coordinates of the triangle vertices
      real,      pointer :: tri_f(:,:)    ! Any additional functions at vertices
      integer,   pointer :: conn(:,:)     ! conn(1:3,n) contains pointers into
                                          ! tri_xyz(:,*) for the vertices of the
                                          ! nth triangle
      integer,   intent (out) :: iostatus !   0 means no error detected;
                                          ! < 0 means no more zones found;
                                          ! > 0 means an I/O error
!     Local constants:

      integer,   parameter :: lunerr         = 6
      character, parameter :: delimiters * 4 = '",= '
      character, parameter :: quotes * 1     = '"'
      character, parameter :: routine * 27   = ' Tecplot_tri_read: Trouble '

!     Local variables:

      integer   :: first, i, l, last, mark
      logical   :: counting
      character :: buffer * 256, format_string * 4

!     Execution:
!     ----------

      if (num_call == 1 ) then ! First call for this file

         open (lun, file=filename, status='old', iostat=iostatus)

         if (iostatus /= 0) then
            write (lunerr, '(/, 3a)') routine, 'opening ', trim (filename)
            go to 99
         end if

         read (lun, '(a)') buffer
         last = len_trim (buffer);  l = last

         if (buffer(1:5) == 'TITLE') then
            first = index (buffer, quotes)
            geom_title = buffer(first+1:last-1)
            read (lun, '(a)') buffer
            last = len_trim (buffer);  l = last
         else
            geom_title = 'Unnamed geometry'
         end if

!        We should now have a line like this:
!        VARIABLES = "X", "Y", "Z", "Var1", ...

!        Count the function names:

         counting = .true.

   10    continue ! Repeat to transfer variable names after allocating the array

            first = index (buffer(1:last), quotes) + 1 ! Avoid VARIABLES = "
            nf    = 0

            do while (first > 2) ! mark = 0 means scan2 reached end of string

               call scan2 (buffer(1:last), delimiters, first, last, mark)

               nf = nf + 1
               if (.not. counting) then
                  if (mark > 0) variable_names(nf) = buffer(first:mark)
               end if
               first = mark + 2
            end do

            nf = nf - 4 ! Don't count x,y,z or the final increment

            if (counting) then
               counting = .false.
               allocate (variable_names(3 + nf))
               write (lunerr, '(a, i3)') ' # variables found beyond x,y,z:', nf
               last = l  ! Lost by last call to scan2
               go to 10
            end if

      end if ! End first call handling of title and variable names

!     Read one "zone" = one triangulated surface dataset.
!     ---------------------------------------------------

!     The next line (if any) should look like this:

!     ZONE [T="Zone title",] N=96000, E=32000, F=FEPOINT, ET=TRIANGLE

      read (lun, '(a)', iostat=iostatus) buffer
      if (iostatus < 0) go to 99  ! EOF

      last = len_trim (buffer)
      first = index (buffer(1:last), quotes) + 1 ! Look for optional zone title

      if (first > 1) then
         call scan2 (buffer(1:last), quotes, first, last, mark)
         zone_title = buffer(first:mark)
      else
         write (zone_title(1:4), '(a1, bz, i3)') 'Z', num_call
         zone_title(5:) = ' '
      end if

      first = index (buffer(1:last), 'N=') + 1

      call scan2 (buffer(1:last), delimiters, first, last, mark)

      format_string = '(i?)'
      write (format_string(3:3), '(i1)') mark - first + 1
      read (buffer(first:mark), format_string, iostat=iostatus) nnode
      if (iostatus /= 0) then
         write (lunerr, '(/, 2a)') routine, 'decoding nnode.'
         go to 99
      end if

      first = index (buffer(1:last), 'E=') + 1

      call scan2 (buffer(1:last), delimiters, first, last, mark)

      write (format_string(3:3), '(i1)') mark - first + 1
      read (buffer(first:mark), format_string, iostat=iostatus) ntri
      if (iostatus /= 0) then
         write (lunerr, '(/, 2a)') routine, 'decoding ntri.'
         go to 99
      end if

      allocate (tri_xyz(3,nnode), tri_f(nf,nnode), conn(3,ntri), stat=iostatus)

      if (iostatus /= 0) then
         write (lunerr, '(/, 2a)') routine, 'allocating triangle space.'
         go to 99
      end if

!!!   read (lun, *, iostat=iostatus) (tri_xyz(:,i), tri_f(:,i), i = 1,nnode)

      do i = 1, nnode
         read (lun, *, iostat=iostatus) tri_xyz(:,i), tri_f(:,i)
      end do

      if (iostatus /= 0) then
         write (lunerr, '(/, 2a)') routine, 'reading triangle data.'
         go to 99
      end if

      read (lun, *, iostat=iostatus) conn

      if (iostatus /= 0) then
         write (lunerr, '(/, 2a)') routine, 'reading connectivity data.'
      end if

   99 continue

      end subroutine Tecplot_tri_read

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Tecplot_tri_write (num_call, lun, filename, formatted,        &
                                    geom_title, variable_names, zone_title,    &
                                    nnode, ntri, nf, tri_xyz, tri_f, conn,     &
                                    iostatus)

!     Write the next zone of a Tecplot surface triangulation with associated
!     function values (0 or more beyond the initial "x", "y", "z" variables).
!     On the first call only, write a geometry title unless it is blank, along
!     with the names of the variables.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer,   intent (in)  :: num_call ! 1 = first call for this file; write
                                          ! an optional title + variable names;
                                          ! 2+ = write the indicated next zone
      integer,   intent (in)  :: lun      ! Logical unit for the output file
                                          ! which is opened here if num_call = 1
      character, intent (in)  :: filename * (*) ! The binary option has to do
                                          ! the open here, so do it regardless
      logical,   intent (in)  :: formatted! T = ASCII; F = binary (incomplete)
      character, intent (in)  :: geom_title * (*) ! File title, unless blank
      character, intent (in)  :: &
        variable_names (*) * (name_limit) ! Variable names, left-justified;
                                          ! array size is 3 + nf
      character, intent (in)  :: zone_title * (*) ! Zone title
      integer,   intent (in)  :: nnode    ! # surface points found
      integer,   intent (in)  :: ntri     ! # triangles found
      integer,   intent (in)  :: nf       ! # functions found beyond x,y,z
      real,      intent (in)  :: &
        tri_xyz(3,nnode)                  ! tri_xyz(1:3,*) contains the (x,y,z)
                                          ! coordinates of the triangle vertices
      real,      intent (in)  :: &
        tri_f(nf,nnode)                   ! Any additional functions at vertices
      integer,   intent (in)  :: &
        conn(3,ntri)                      ! conn(1:3,n) contains pointers into
                                          ! tri_xyz(:,*) for the vertices of the
                                          ! nth triangle
      integer,   intent (out) :: iostatus !   0 means no error detected;
                                          ! < 0 means no more zones found;
                                          ! > 0 means an I/O error
!     Local constants:

      integer,   parameter :: lunerr            = 6
      character, parameter :: comma * 1         = ','
      character, parameter :: delimiters * 4    = '",= '
      character, parameter :: element_type * 24 = ', F=FEPOINT, ET=TRIANGLE'
      character, parameter :: i1_format * 4     = '(i1)'
      character, parameter :: quotes * 1        = '"'
      character, parameter :: routine * 28      = ' Tecplot_tri_write: Trouble '

!     Local variables:

      integer   :: l, m, n, ndigits, nvars
      character :: eformat * 14
      character :: zone_format_no_title   * 21
      character :: zone_format_with_title * 28

!     Data:

      data zone_format_no_title   / '(a7, in, a4, ie, a24)' /
!                                    1234567890123456789012
      data zone_format_with_title / '(a8, a, a5, in, a4, ie, a24)' /

!     Execution:
!     ----------

!     THE BINARY OPTION IS NOT IMPLEMENTED YET.
!     WE ALSO DEFAULT THE PRECISION TO SINGLE FOR NOW.

!     However, adopt the following from Tecplot_io.f90 for formatted output use:

      eps = epsilon (eps)
      if (eps < 1.e-10) then
         isdouble = 1
         eformat  = '(1p, nnne17.9)'
      else
         isdouble = 0
         eformat  = '(1p, nnne15.7)'
      end if

      write (eformat(6:8), '(i3)') 3 + nf

      if (num_call == 1 ) then ! First call for this file

         open (lun, file=filename, status='unknown', iostat=iostatus)

         if (iostatus /= 0) then
            write (lunerr, '(/, 3a)') routine, 'opening ', trim (filename)
            go to 99
         end if

         l = len_trim (geom_title)
         if (l > 0) write (lun, '(3a)') 'TITLE="', geom_title(1:l), quotes

!        The variable names need to be packed for the binary option.
!        Employ the buffer similarly, but the double quotes are needed too:

         nvars = 3 + nf
         packed_names(1:1) = quotes
         m = 1
         do n = 1, nvars
            l = len_trim (variable_names(n))
            packed_names(m+1:m+l) = variable_names(n)(1:l);  l = m + l + 1
            packed_names(l:l+3)   = '", "';         m = l + 3
         end do

         write (lun, '(2a)') 'VARIABLES = ', packed_names(1:l)

      end if ! End first call handling of title and variable names

!     Write one "zone" = one triangulated surface dataset.
!     ----------------------------------------------------

!     The next line should look like this:

!     ZONE [T="Zone title",] N=96000, E=32000, F=FEPOINT, ET=TRIANGLE

      l = len_trim (zone_title)

      if (l == 0) then ! Suppress the zone title

         call count_digits (nnode, ndigits)

         write (zone_format_no_title(7:7),   i1_format) ndigits

         call count_digits (ntri,  ndigits)

         write (zone_format_no_title(15:15), i1_format) ndigits

         write (lun, zone_format_no_title) &
            'ZONE N=', nnode, ', E=', ntri, element_type 
      else

         call count_digits (nnode, ndigits)

         write (zone_format_with_title(14:14), i1_format) ndigits

         call count_digits (ntri,  ndigits)

         write (zone_format_with_title(22:22), i1_format) ndigits

         write (lun, zone_format_with_title) &
            'ZONE T="', zone_title(1:l), '", N=', nnode, ', E=', ntri, &
            element_type
      end if

      if (nf > 0) then
         write (lun, eformat, iostat=iostatus) &
            (tri_xyz(:,n), tri_f(:,n), n = 1, nnode)
      else
         write (lun, eformat, iostat=iostatus) tri_xyz
      end if

      if (iostatus /= 0) then
         write (lunerr, '(/, 2a)') routine, 'writing triangle data.'
         go to 99
      end if

      write (lun, '(i5, 2i6)', iostat=iostatus) conn

      if (iostatus /= 0) then
         write (lunerr, '(/, 2a)') routine, 'writing connectivity data.'
      end if

   99 continue

!     Internal procedure for Tecplot_tri_write:

      contains

         subroutine count_digits (n, ndigits)  ! For 0 <= n < 10 million

         integer, intent (in)  :: n
         integer, intent (out) :: ndigits

         if (n < 10) then
            ndigits = 1
         else if (n < 100) then
            ndigits = 2
         else if (n < 1000) then
            ndigits = 3
         else if (n < 10000) then
            ndigits = 4
         else if (n < 100000) then
            ndigits = 5
         else if (n < 1000000) then
            ndigits = 6
         else
            ndigits = 7
         end if
 
         end subroutine count_digits

      end subroutine Tecplot_tri_write

   end module Tecplot_tri_io
