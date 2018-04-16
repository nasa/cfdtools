!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine get_coordinates (case_name, ndim, coords, ios)
!
!  Purpose:
!
!     Extract one or more numerical values from a character string.  The initial
!  application is to a file name of the type <ID>.xxx.v10.0h71.dat where the
!  desired output pair of coordinates (ndim = 2) is (10.0, 71.0).  Given that
!  the <ID> (identifier) indicated here may also contain digits, we choose to
!  stipulate that the extracted coordinates are the LAST ndim numerical values
!  in the string, not the first.
!
!  Sample string:          antenna1_17.attenuation_v10.0_h71.dat
!
!  Further constraints:  o Each "coordinate" is returned as real.
!
!                        o Negative values are not handled (yet).
!                        o Embedded decimal points are permitted, but not
!                          trailing decimal points.  (Integers are OK.)
!                        o Each coordinate value is preceded by a nonnumeric
!                          character or characters, but the actual coordinate
!                          names (v, h in this example) are not important.
!
!  Strategy:               Isolate the last two numeric substrings using a
!                          backward scan of the characters.
!
!     This fairly general functionality may have other applications too.  Making
!  the utility an argument to the application routine also allows for different
!  conventions to apply, but the argument list should not need changing.
!
!  History:  D.A.Saunders  02/14/16  Initial implementation for Radio_Blackout.
!              "     "     02/10/17  The initial history date was missing.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   character (*), intent (in)  :: case_name    ! String containing numerical
                                               ! coordinates to be extracted
   integer,       intent (in)  :: ndim         ! Number of coordinate dimensions
                                               ! such as 2 for V and h
   real,          intent (out) :: coords(ndim) ! Desired ntuple of coordinate
                                               ! values
   integer,       intent (out) :: ios          ! 0 => success;
                                               ! 1 => less than ndim valid real
                                               !      values found
!  Local constants:

   integer, parameter :: notset = -9999999
   real, parameter    :: unset  = -9999999.
   character (1)      :: period = '.'

!  Local variables:

   integer :: i, i1, i2, inumber, j, l, m, n
   real    :: rnumber

!  Procedures:

   logical  :: number
   external :: decode_number, number

!  Execution:

   coords(:) = unset

   l = len_trim (case_name);  m = l

   do n = ndim, 1, -1
      i1 = notset
      do i = m, 1, -1 
         if (number (case_name(i:i))) then
            i2 = i
            do j = i - 1, 1, -1
               if (.not. number (case_name(j:j))) then
                  if (case_name(j:j) == period) then  ! Keep going
                  else
                     i1 = j + 1
                     exit
                  end if
               end if
            end do
            m = max (1, j - 1)
            exit
         end if
      end do
      if (i1 == notset) then
         ios = 1
      else
         call decode_number (case_name(i1:i2), inumber, rnumber, ios)
         if (ios == 0) then  ! Neither integer nor real (shouldn't be possible)
             ios = 1
         else
             ios = 0;  coords(n) = rnumber
         end if
      end if
      if (ios /= 0) exit
   end do

   if (ios == 0) then  ! Make sure ndim coordinates were found:
      do n = ndim, 1, -1
         if (coords(n) == unset) then
            ios = 1
            exit
         end if
      end do
   end if

   end subroutine get_coordinates
