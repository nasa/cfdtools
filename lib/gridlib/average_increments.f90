!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine average_increments (ni, nj, nk, x, y, z, average_ds)

!  Calculate the average off-face grid spacings for faces 1:6 of a grid block.
!  To identify the most likely solid wall face, use iwall = minloc (average_ds).
!
!  08/23/05  David Saunders  Initial implementation, for identifying wall faces.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)                       :: ni, nj, nk    ! Block dims.
   real,    intent (in), dimension (ni,nj,nk) :: x, y, z       ! Block coords.
   real,    intent (out)                      :: average_ds(6) ! (1) => i = 1,
                                                               ! (2) => i = ni,
                                                               ! and so on
!  Local constants:

   real, parameter :: zero = 0.

!  Local variables:

   integer :: i, j, k

!  Execution:

   average_ds(:) = zero

!  The i direction:

   do k = 1, nk
      do j = 1, nj
         average_ds(1) = average_ds(1) + sqrt ( &
            (x(2,j,k) - x(1,j,k))**2 + &
            (y(2,j,k) - y(1,j,k))**2 + &
            (z(2,j,k) - z(1,j,k))**2   )
         average_ds(2) = average_ds(2) + sqrt ( &
            (x(ni,j,k) - x(ni-1,j,k))**2 + &
            (y(ni,j,k) - y(ni-1,j,k))**2 + &
            (z(ni,j,k) - z(ni-1,j,k))**2   )
      end do
   end do

   average_ds(1:2) = average_ds(1:2) / real (nj * nk)

!  The j direction:

   do k = 1, nk
      do i = 1, ni
         average_ds(3) = average_ds(3) + sqrt ( &
            (x(i,2,k) - x(i,1,k))**2 + &
            (y(i,2,k) - y(i,1,k))**2 + &
            (z(i,2,k) - z(i,1,k))**2   )
         average_ds(4) = average_ds(4) + sqrt ( &
            (x(i,nj,k) - x(i,nj-1,k))**2 + &
            (y(i,nj,k) - y(i,nj-1,k))**2 + &
            (z(i,nj,k) - z(i,nj-1,k))**2   )
      end do
   end do

   average_ds(3:4) = average_ds(3:4) / real (nk * ni)

!  The k direction:

   do j = 1, nj
      do i = 1, ni
         average_ds(5) = average_ds(5) + sqrt ( &
            (x(i,j,2) - x(i,j,1))**2 + &
            (y(i,j,2) - y(i,j,1))**2 + &
            (z(i,j,2) - z(i,j,1))**2   )
         average_ds(6) = average_ds(6) + sqrt ( &
            (x(i,j,nk) - x(i,j,nk-1))**2 + &
            (y(i,j,nk) - y(i,j,nk-1))**2 + &
            (z(i,j,nk) - z(i,j,nk-1))**2   )
      end do
   end do

   average_ds(5:6) = average_ds(5:6) / real (ni * nj)

   end subroutine average_increments
