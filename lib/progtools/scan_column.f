!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine scan_column (lun, icolumn, maxlen, nlines, ier)

!     One-liner:  Scan a text column to EOF; return line count & max. token size
!
!     Description:
!
!        This utility facilitates dynamic allocation of a character array.
!     Unfortunately, there seems to be no way to treat variable-length strings.
!     At best, we can verify that the given limit on string length suffices to
!     ensure that unintended duplicate short strings are not stored when the
!     array is actually read following the scan.  This is important when setting
!     up dictionaries for keyword look-ups.
!
!     Usage:
!
!        >  integer, parameter :: max_length = 20 ! Say
!        >  integer :: icolumn, maxlen, nlines, ier
!        >  character, allocatable, dimension (:) :: dictionary * (max_length)
!        >  Position the file at the first line to be counted.
!        >  call scan_column (lun, icolumn, maxlen, nlines, ier)
!        >  if (ier /= 0) then ! Handle I/O error
!        >  if (maxlen > max_length) then ! Handle non-unique strings [fatal?]
!        >  allocate (dictionary (nlines))
!        >  reposition the file and read the known number of lines appropriately
!
!     09/09/03  D.Saunders  Initial implementation.
!
!     David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Arguments:

      integer, intent (in) ::
     >   lun,                ! Logical unit assumed positioned at first line
     >   icolumn             ! Number of the column to scan; icolumn >= 1

      integer, intent (out) ::
     >   maxlen,             ! Length of longest string found in that column
     >   nlines,             ! Number of lines found
     >   ier                 !   0 = no problems;
                             ! > 0 = read error on line ier
                             !       (probably too few columns)

!     Local constants:

      integer, parameter ::
     >   max_length = 64     ! Must be more than enough for any application;
                             ! length of local string used to read column

!     Local variables:

      integer
     >   ios

      character
     >   token * (max_length)

      character
     >   skip_text (icolumn) * 1  ! 1 suffices because ios=0 for excess chars.

!     Execution:

      maxlen = 0
      nlines = 0

      do ! Until EOF

         read (lun, *, iostat=ios) skip_text (1 : icolumn - 1), token

         if (ios < 0) then ! EOF
            ier = 0
            exit
         end if

         if (ios > 0) then
            ier = nlines + 1
            exit
         end if

         nlines = nlines + 1
         maxlen = max (maxlen, len_trim (token))

      end do

      end subroutine scan_column
