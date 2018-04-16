!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine check_existence (filename, lunerr, ios)

!  Combine Fortran 90's INQUIRE statement with a diagnostic if the indicated
!  file is not found.
!
!  07/07/2014  D.A.Saunders  Initial implementation, for PREPARE_NEQAIR_DATA.
!
!  Author:  David Saunders  ERC, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   character, intent (in)  :: filename*(*)  ! File name of interest
   integer,   intent (in)  :: lunerr        ! Used if the file is not found
   integer,   intent (out) :: ios           ! 0 => file was found
                                            ! 1 => file was not found
!  Local variables:

   logical :: exist

!  Execution:

   inquire (file=filename, exist=exist, iostat=ios)

   if (exist) then
      ios = 0  ! Make sure of it
   else
      write (lunerr, '(2a)') 'File not found: ', trim (filename)
      ios = 1
   end if

   end subroutine check_existence
