!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine getargs (maxargs, maxchars, progname, nargs, args)

!  Description:
!
!     Use Fortran 90 intrinsics "iargc" and "getarg" to return any command line
!  arguments as an array of character strings.  Each string is terminated here
!  with a null character for possible passing to a C function.  The program name
!  entered as an argument is inserted as element 0 of the arguments array, since
!  this is what C also expects.
!
!     Note that the null characters inserted here count as an extra character
!  when the argument lengths are recovered with the Fortran intrinsic len_trim.
!
!     Intel compilers require the -Vaxlib switch to access the intrinsics.
!
!  11/03/05  D.A.Saunders  Initial implementation.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer,   intent (in)  :: maxargs             ! 0:maxargs are allowed for
   integer,   intent (in)  :: maxchars            ! Max. length of any arg. + 1
   character, intent (in)  :: progname*(*)        ! Program name --> args(0)
   integer,   intent (out) :: nargs               ! # args. found; 0 => none
   character, intent (out) :: args(0:maxargs)*(*) ! Array of null-terminated
                                                  ! command line arguments, with
                                                  ! args(0) = progname;
                                                  ! args(0:nargs) are returned
!  Local variables:

   integer   :: i, l
   character :: null*1

!  Fortran intrinsic function:

   integer :: iargc                               ! Counts command line args.

!  Execution:

   null    = char (0)
   nargs   = min (iargc (), maxargs)              ! Omits program name, unlike C
   args(0) = progname
   l       = min (len_trim (progname) + 1, maxchars)
   args(0)(l:l) = null

   do i = 1, nargs                                ! May be zero

      call getarg (i, args(i))                    ! Unix intrinsic

      l = min (len_trim (args(i)) + 1, maxchars)
      args(i)(l:l) = null

   end do 

   end subroutine getargs
