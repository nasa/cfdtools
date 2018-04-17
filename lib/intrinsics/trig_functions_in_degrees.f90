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

module trigd

    implicit none

    real, parameter :: radians_per_degree = acos(-1.0d0)/180.0d0
    real, parameter :: degrees_per_radian = 180.0d0/acos(-1.0d0)

contains

real function cosd (degrees)
    real, intent (in) :: degrees
    cosd = cos (degrees * radians_per_degree)
end function cosd

real function sind (degrees)
    real, intent (in) :: degrees
    sind = sin (degrees * radians_per_degree)
end function sind

real function tand (degrees)
    real, intent (in) :: degrees
    tand = tan (degrees * radians_per_degree)
end function tand

real function acosd (value)
    real, intent (in) :: value
    acosd = acos (value) * degrees_per_radian
end function acosd

real function asind (value)
    real, intent (in) :: value
    asind = asin (value) * degrees_per_radian
end function asind

real function atand (value)
    real, intent (in) :: value
    atand = atan(value) * degrees_per_radian
end function atand

real function atan2d (x,y)
    real, intent (in) :: x, y
    atan2d = atan2(x,y) * degrees_per_radian
end function atan2d

end module trigd
