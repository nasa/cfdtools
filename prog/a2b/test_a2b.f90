!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program test_a2b

!  Simple-minded way of checking output from A2B: prompt for a stream file
!  position number and data type (known to the tester), read the data unit,
!  and print it to the screen in a loop over prompts.
!
!  History:
!
!  12/28/2017  David Saunders  Initial way of testing A2B on some simple file.
!  12/29/2017  DAS             Had to include "nchar" for alphanumeric tokens.
!  01/03/2018  DAS             Handling of logicals had been overlooked.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   integer, parameter :: &
      lunin   = 1,       &
      lunkbd  = 5,       &
      luncrt  = 6

   integer :: &
      inumber, ios, itype, n, nchar

   real :: &
      rnumber

   logical :: &
      cr, eof, TorF

   character (132) :: &
      cstring, filename

   call reads (luncrt, 'Stream file name: ', &
               lunkbd, filename, cr, eof)
   open (lunin, file=filename, access='stream', form='unformatted')

!  Prompt for a known position and data type:

   do  ! Until ^D

      write (luncrt, '(a)', advance='no') &
         'Enter position n and data type 1=int, 2=real, 3=string, 4=logical: '
      read  (lunkbd, *, iostat=ios) n, itype
      if (ios < 0) exit

      select case (itype)
         case (1)  ! Integer scalar
            read  (lunin, pos=n, iostat=ios) inumber
            if (ios /= 0) then
               write (luncrt, 101)
               cycle
            end if
            write (luncrt, '(a, i12)') 'Integer found:', inumber
         case (2)  ! Real scalar
            read  (lunin, pos=n) rnumber
            if (ios /= 0) then
               write (luncrt, 101)
               cycle
            end if
            write (luncrt, '(a, es25.15)', iostat=ios) 'Real found:', rnumber
         case (3)  ! Character string
            read  (lunin, pos=n, iostat=ios) nchar, cstring(1:nchar)
            if (ios /= 0) then
               write (luncrt, 101)
               cycle
            end if
            write (luncrt, '(2a)') 'String found: ', cstring(1:nchar)
         case (4)  ! Logical value
            read  (lunin, pos=n, iostat=ios) TorF
            if (ios /= 0) then
               write (luncrt, 101)
               cycle
            end if
            write (luncrt, '(a, l2)') 'Logical found:', TorF
      end select

   end do

101 format ('Bad read; likely a wrong entry n for pos=n keyword. Try again.')

   end program test_a2b
