!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program makelist
!
!     Make a list of integers, probably for entering to another program, or
!     inserting ordinals as a column in a dataset.  (See utility COLUMNEDIT
!     for doing that.)
!
!  History:
!
!     07/23/04  David Saunders  Original simple utility (output to screen).
!     01/06/21    "      "      Send the list to a file unless it's fairly
!                               short, when cutting and pasting is easy.
!
!  Author:
!     David Saunders, ELORET/NASA Ames Research Center.
!                     Now with AMA, Inc. at NASA Ames.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   integer :: i, i1, i2, inc, l, lun, ndigits, nlist
   character (4) cformat

   cformat = '(in)'
   write (*, '(/, a)') ' Enter i1, i2, inc: '
   read  (*, *) i1, i2, inc

   l = max (abs (i1), abs (i2))
   ndigits = min (max (1, int (alog10 (real (10*l)))), 9)
   if (i1 < 0 .or. i2 < 0) ndigits = ndigits + 1  ! Allow for the -ve sign
   write (cformat(3:3), '(i1)') ndigits
   lun = 6
   nlist = nint (real (i2 - i1) / real (inc)) + 1

   if (nlist > 100) then
      lun = 1
      open (lun, file='ordinals.dat', status='unknown')
   end if

   write (6, '(a)') cformat
   write (lun, cformat) (i, i = i1, i2, inc)

   end program makelist
