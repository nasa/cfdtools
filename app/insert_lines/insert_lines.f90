!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program insert_lines
!
!  Starting at a specified line, insert the lines from a specified file at a
!  specified interval.
!
!  02/02/09  David Saunders  Initial implementation, prompted by needing to
!                            convert a list of surface data points back to
!                            multiple structured patches as needed for
!                            Kriging interpolation at CFD grid points.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Constants:

   integer, parameter :: lunin1 = 1  ! Dataset to be modified
   integer, parameter :: lunin2 = 2  ! Dataset containing inserted lines
   integer, parameter :: lunout = 3  ! Output result
   integer, parameter :: lunkbd = 5  ! Keyboard
   integer, parameter :: luncrt = 6  ! screen

!  Variables:

   integer :: i, interval, ios, iskip, ninsert, nskip
   character :: buffer * 132, filename1 * 64, filename2 * 64, filename3 * 64
   character, allocatable, dimension (:) :: insertions * 132

!  Execution:

   ios = 1
   do while (ios /= 0)
      write (luncrt, '(/, a)', advance='no') 'File to be modified: '
      read  (lunkbd, *) filename1
      open  (lunin1, file=filename1, status='old', iostat=ios)
   end do

   ios = 1
   do while (ios /= 0)
      write (luncrt, '(a)', advance='no') 'File to be inserted: '
      read  (lunkbd, *) filename2
      open  (lunin2, file=filename2, status='old', iostat=ios)
   end do

   ios = 1
   do while (ios /= 0)
      write (luncrt, '(a)', advance='no') 'Output file name:    '
      read  (lunkbd, *) filename3
      open  (lunout, file=filename3, status='unknown', iostat=ios)
   end do

   write (luncrt, '(a)', advance='no') '# lines above first insertion:  '
   read  (lunkbd, *) nskip

   write (luncrt, '(a)', advance='no') '# lines between the insertions: '
   read  (lunkbd, *) interval

   ninsert = -1
   do while (ios == 0)
      ninsert = ninsert + 1
      read (lunin2, '(a)', iostat=ios)
   end do

   write (luncrt, '(a, i6)') '# lines found to insert:', ninsert

   allocate (insertions(ninsert))

   rewind (lunin2)

   read (lunin2, '(a)') insertions(:);  close (lunin2)

   iskip = nskip
   ios = 0
   do while (ios == 0)
      do i = 1, iskip
         read (lunin1, '(a)', iostat=ios) buffer
         if (ios /= 0) exit
         write (lunout, '(a)') trim (buffer)
      end do
      if (ios /= 0) exit

      write (lunout, '(a)') (trim (insertions(i)), i = 1, ninsert)
      iskip = interval
   end do

   close (lunin1);  close (lunout);  deallocate (insertions)

   end program insert_lines
