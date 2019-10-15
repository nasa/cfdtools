!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine determine_grid_dim (grid_file_name, lun, formatted, ndim, ios)
!
!  Determine whether a PLOT3D-type multiblock grid file is 2- or 3-dimensional.
!  Its form (formatted or unformatted) is assumed to have been determined
!  already (possibly via determine_grid_form), and the existence of the file
!  is also assumed to have been verified.
!
!  If the input file name is 'none', the file is assumed to be open, ready
!  for reading the first record.  It is then rewound upon return.
!
!  Otherwise, the indicated file name is opened here, and the file is closed
!  upon return.
!
!  The main difficulty, at least with an Intel compiler, is that reading the
!  header records with a 3-D utility does not necessarily fail if the file is
!  actually 2-D and formatted.
!
!  Therefore, the (formatted) strategy is to buffer line 2, assuming separate
!  lines for the dimensions of each block, and use an internal read on the
!  buffer that may or may not encounter end-of-data.  The unformatted case
!  can safely try to read the whole header as 3-D and rely on iostat (but see
!  07/01/2018 history entry).
!
!  History:
!
!     12/12/2013  D.A.Saunders  Adaptation of determine_grid_form.
!     04/09/2014    "     "     Opening and closing the file was inconsistent
!                               with the file_prompt utility in xyzq_io.f90.
!                               Workaround: enter grid_file_name = 'none' to
!                               indicate that the file is open and stays so.
!                               Also: Mike Olsen pointed out a better way than
!                               counting tokens on line 2: just do an internal
!                               read of three integers, and trap any I/O error.
!     07/09/2014    "     "     Glitch: ndim wasn't being set for 2-D/unformat.
!     07/01/2018    "     "     An unformatted 2-D grid was wrongly determined
!                               to be 3-D. Check what is found for the third
!                               dimension of all blocks, and switch if an
!                               unlikely value is found.  (Intel compiler.) :(
!
!  Author:  David Saunders, ERC Inc./NASA Ames Research Center, CA
!           Now with AMA, Inc. at NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   character, intent (in)  :: grid_file_name*(*)  ! See above for use of 'none'
   integer,   intent (in)  :: lun
   logical,   intent (in)  :: formatted  ! T|F
   integer,   intent (out) :: ndim       ! 2 or 3
   integer,   intent (out) :: ios        ! ios /= 0 suggests a bad file name

!  Local constants:

   integer,        parameter :: is_unlikely = 100000000
   character (11), parameter :: format = 'unformatted'
   character (1),  parameter :: blank  = ' '

!  Local variables:

   integer        :: i1, ib, l, l2, nb, ni, nj, nk
   logical        :: already_open
   character (64) :: buffer
   integer, allocatable :: nijk(:,:)

!  Execution:

   l =  len_trim (grid_file_name)
   already_open = grid_file_name(1:l) == 'none'

   if (.not. already_open) then
      i1 = 1
      if (formatted) i1 = 3

      open (lun, file=grid_file_name(1:l), form=format(i1:11), &
            status='old', iostat=ios)

      if (ios /= 0) then
         write (*, '(2a)') ' Unable to open ', grid_file_name(1:l)
         go to 99
      end if
   end if

   if (formatted) then  ! Look for 2 or 3 dimensions on line 2

      read (lun, *, iostat=ios) nb
      if (ios /= 0) go to 90

      read (lun, '(a)', iostat=ios) buffer
      if (ios /= 0) then
         write (*, '(2a)') 'Trouble reading line 2 of ', grid_file_name(1:l)
         go to 95
      end if

      l2 = len_trim (buffer)
      read (buffer(1:l2), *, iostat=ios) ni, nj, nk
      if (ios /= 0) then
         ndim = 2
         ios  = 0
      else
         ndim = 3
      end if

   else  ! Unformatted file:  try to read a full 3-D header

      read (lun, iostat=ios) nb
      if (ios /= 0) go to 90

      allocate (nijk(3,nb))

      read (lun, iostat=ios) nijk
      if (ios /= 0) then
         if (.not. already_open) then
            write (*, '(1x, 2a)') grid_file_name(1:l), ' appears to be 2-D.'
         else
            write (*, '(a)') ' The file appears to be 2-D.'
         end if
         ndim = 2
         ios  = 0
      else  ! Make sure that all the third dimensions look plausible
         ndim = 0
         do ib = 1, nb
            if (nijk(3,ib) < 1 .or. nijk(3,ib) > is_unlikely) ndim = 2
         end do
         if (ndim /= 2) ndim = 3  ! Don't diagnose the more likely case
      end if

      deallocate (nijk)

   end if

   go to 95

90 write (*, '(2a)') 'Trouble reading # blocks in ', grid_file_name(1:l)

95 continue

   if (already_open) then
      rewind (lun)
   else
      close (lun)
   end if

99 return

   end subroutine determine_grid_dim
