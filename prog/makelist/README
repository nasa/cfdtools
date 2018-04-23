      program makelist

!     Make a list of integers, probably for entering to another program.
!
!     07/23/04  David Saunders, ELORET/NASA Ames Research Center.

      implicit none

      integer :: i, i1, i2, inc, l

      write (*, '(/, a)') ' Enter i1, i2, inc: '
      read  (*, *) i1, i2, inc

      l = max (abs (i1), abs (i2))
      if (i1 < 0 .or. i2 < 0) l = 10 * l  ! Allow for the -ve sign

      select case (l)

      case (1:9)
         write (*, '(i1)') (i, i = i1, i2, inc)
      case (10:99)
         write (*, '(i2)') (i, i = i1, i2, inc)
      case (100:999)
         write (*, '(i3)') (i, i = i1, i2, inc)
      case (1000:9999)
         write (*, '(i4)') (i, i = i1, i2, inc)
      case (10000:99999)
         write (*, '(i5)') (i, i = i1, i2, inc)
      case (100000:999999)
         write (*, '(i6)') (i, i = i1, i2, inc)
      case (1000000:9999999)
         write (*, '(i7)') (i, i = i1, i2, inc)
      case (10000000:99999999)
         write (*, '(i8)') (i, i = i1, i2, inc)
      case (100000000:999999999)
         write (*, '(i9)') (i, i = i1, i2, inc)
      case default
         write (*, *) (i, i = i1, i2, inc)
      end select

      end program makelist
