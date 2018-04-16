!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  These forms of the cosd, sind, and tand trigonometric intrinsics are an
!  unfortunate requirement for some compilers, even in the year 2007.
!
!  "Degrees" forms of the intrinsics are often preferable, but if supported
!  at all are most likely to be done so via the "radians" forms that are now
!  normally implemented in hardware.
!
!  Link to the present functions only if your compiler does not provide them.
!  Declarations such as "real :: cosd" in calling routines should suffice not
!  to clash with the intrinsics when they are available.
!
!  Note that the equivalent of the "-r8" switch for compiling "real"s with 64-
!  bit precision is assumed.  Failing that, an even more distasteful workaround
!  is required by Fortran 90 for all real variables and constants.
!
!  10/24/2007  David Saunders  Response to Todd White's finding that g95 does
!              ELORET Corp.    not support the "degrees" intrinsics.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module degrees_to_radians

!     Hard-coding the conversion constant may be preferable, but this approach
!     costs little and always produces the relevant precision:

      implicit none

      real,    public,  save :: radians_per_degree
      logical, private, save :: first = .true.
      public                 :: calculate_conversion_factor
      
      contains

         subroutine calculate_conversion_factor ()

         if (first) then
            first = .false.
            radians_per_degree = atan (1.0) / 45.0
         end if

         end subroutine calculate_conversion_factor

   end module degrees_to_radians

!  ----------------------------
   real function cosd (degrees)
!  ----------------------------

   use degrees_to_radians

   real, intent (in) :: degrees

   call calculate_conversion_factor ()

   cosd = cos (degrees * radians_per_degree)

   end function cosd

!  ----------------------------
   real function sind (degrees)
!  ----------------------------

   use degrees_to_radians

   real, intent (in) :: degrees

   call calculate_conversion_factor ()

   sind = sin (degrees * radians_per_degree)

   end function sind

!  ----------------------------
   real function tand (degrees)
!  ----------------------------

   use degrees_to_radians

   real, intent (in) :: degrees

   call calculate_conversion_factor ()

   tand = tan (degrees * radians_per_degree)

   end function tand
