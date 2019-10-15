!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine count_records (lun, formatted, n, ios)
!
!  Count the number of lines in a formatted file, or the number of records in an
!  unformatted file, as an aid to work-space allocation.   Intended applications
!  are to files with contents suited to reading as an array. The calling program
!  should handle any header records appropriately.
!
!  04/11/2018  D.A.Saunders  Initial implementation, with processing of large
!                            files output by NEQAIR in mind.  Long overdue!
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: lun        ! Logical unit; the file is assumed to be
                                       ! open, and it is rewound here upon exit
   logical, intent (in)  :: formatted  ! T/F appropriately
   integer, intent (out) :: n          ! The requested count of lines/records
   integer, intent (out) :: ios        ! ios = 0 means no read error;
                                       ! ios > 0 means a read error
                                       ! ios < 0 is not an option

!  Execution:

   n = 0

   if (formatted) then
      do  ! Until EOF
         read (lun, *, iostat=ios)
         if (ios /= 0) exit
         n = n + 1
      end do
   else
      do  ! Until EOF
         read (lun, iostat=ios)
         if (ios /= 0) exit
         n = n + 1
      end do
   end if

   rewind (lun)

   if (ios < 0) ios = 0  ! ios < 0 means EOF as expected

   end subroutine count_records
