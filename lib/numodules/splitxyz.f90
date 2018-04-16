!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine splitxyz (coord, test, cutoff, n, x, y, z, nsplit)

!  Split an (x,y,z) dataset in-place according to the indicated coordinate and
!  splitting criterion.  The input arrays are repacked with elements 1:nsplit
!  satisfying the comparison test against the cutoff value for x, y, or z.
!
!  08/13/2008  D.A.Saunders  Initial implementation for splitting a cloud of
!                            laser scan data so that a forebody and aft body
!                            of a recovered entry vehicle can be processed
!                            separately.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   character, intent (in)    :: coord*1       ! 'X', 'Y', or 'Z' (either case)
   character, intent (in)    :: test*(2)      ! '< ', '<=', '==', '> ', or '>='
   real,      intent (in)    :: cutoff        ! Cut-off value to compare with
   integer,   intent (in)    :: n             ! # input data points 
   real,      intent (inout) :: x(n),  &      ! Data coordinates, repacked upon
                                y(n),  &      ! return such that the elements
                                z(n)          ! 1:nsplit satisfy the cutoff test
   integer,   intent (out)   :: nsplit        ! # points satisfying the test

!  Local variables:

   integer :: i, icoord, ikeep

!  Execution:

   select case (coord)
      case ('x', 'X')
         icoord = 1
      case ('y', 'Y')
         icoord = 2
      case ('z', 'Z')
         icoord = 3
   end select

   ikeep = 0

   select case (test)

      case ('< ')

         select case (icoord)
            case (1)  ! X
               do i = 1, n
                  if (x(i) <  cutoff) then
                     ikeep = ikeep + 1
                     x(ikeep) = x(i)
                     y(ikeep) = y(i)
                     z(ikeep) = z(i)
                  end if
               end do
            case (2)  ! Y
               do i = 1, n
                  if (y(i) <  cutoff) then
                     ikeep = ikeep + 1
                     x(ikeep) = x(i)
                     y(ikeep) = y(i)
                     z(ikeep) = z(i)
                  end if
               end do
            case (3)  ! Z
               do i = 1, n
                  if (z(i) <  cutoff) then
                     ikeep = ikeep + 1
                     x(ikeep) = x(i)
                     y(ikeep) = y(i)
                     z(ikeep) = z(i)
                  end if
               end do
         end select

      case ('<=')

         select case (icoord)
            case (1)  ! X
               do i = 1, n
                  if (x(i) <= cutoff) then
                     ikeep = ikeep + 1
                     x(ikeep) = x(i)
                     y(ikeep) = y(i)
                     z(ikeep) = z(i)
                  end if
               end do
            case (2)  ! Y
               do i = 1, n
                  if (y(i) <= cutoff) then
                     ikeep = ikeep + 1
                     x(ikeep) = x(i)
                     y(ikeep) = y(i)
                     z(ikeep) = z(i)
                  end if
               end do
            case (3)  ! Z
               do i = 1, n
                  if (z(i) <= cutoff) then
                     ikeep = ikeep + 1
                     x(ikeep) = x(i)
                     y(ikeep) = y(i)
                     z(ikeep) = z(i)
                  end if
               end do
         end select

      case ('==')

         select case (icoord)
            case (1)  ! X
               do i = 1, n
                  if (x(i) == cutoff) then
                     ikeep = ikeep + 1
                     x(ikeep) = x(i)
                     y(ikeep) = y(i)
                     z(ikeep) = z(i)
                  end if
               end do
            case (2)  ! Y
               do i = 1, n
                  if (y(i) == cutoff) then
                     ikeep = ikeep + 1
                     x(ikeep) = x(i)
                     y(ikeep) = y(i)
                     z(ikeep) = z(i)
                  end if
               end do
            case (3)  ! Z
               do i = 1, n
                  if (z(i) == cutoff) then
                     ikeep = ikeep + 1
                     x(ikeep) = x(i)
                     y(ikeep) = y(i)
                     z(ikeep) = z(i)
                  end if
               end do
         end select

      case ('> ')

         select case (icoord)
            case (1)  ! X
               do i = 1, n
                  if (x(i) >  cutoff) then
                     ikeep = ikeep + 1
                     x(ikeep) = x(i)
                     y(ikeep) = y(i)
                     z(ikeep) = z(i)
                  end if
               end do
            case (2)  ! Y
               do i = 1, n
                  if (y(i) >  cutoff) then
                     ikeep = ikeep + 1
                     x(ikeep) = x(i)
                     y(ikeep) = y(i)
                     z(ikeep) = z(i)
                  end if
               end do
            case (3)  ! Z
               do i = 1, n
                  if (z(i) >  cutoff) then
                     ikeep = ikeep + 1
                     x(ikeep) = x(i)
                     y(ikeep) = y(i)
                     z(ikeep) = z(i)
                  end if
               end do
         end select

      case ('>=')

         select case (icoord)
            case (1)  ! X
               do i = 1, n
                  if (x(i) >= cutoff) then
                     ikeep = ikeep + 1
                     x(ikeep) = x(i)
                     y(ikeep) = y(i)
                     z(ikeep) = z(i)
                  end if
               end do
            case (2)  ! Y
               do i = 1, n
                  if (y(i) >= cutoff) then
                     ikeep = ikeep + 1
                     x(ikeep) = x(i)
                     y(ikeep) = y(i)
                     z(ikeep) = z(i)
                  end if
               end do
            case (3)  ! Z
               do i = 1, n
                  if (z(i) >= cutoff) then
                     ikeep = ikeep + 1
                     x(ikeep) = x(i)
                     y(ikeep) = y(i)
                     z(ikeep) = z(i)
                  end if
               end do
         end select

   end select

   nsplit = ikeep

   end subroutine splitxyz
