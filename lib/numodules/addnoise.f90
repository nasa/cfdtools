!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine addnoise (nrow, ncol, table, percent)

!     Add Gaussian noise to all columns of the given data table.  The mean of
!  the noise is zero and the standard deviation is the given percentage of the
!  data range for each column.
!
!  06/24/2017  D.A.Saunders  Initial implementation for thermocouple and iso-
!                            therm sensor data.
!
!  Author:  David Saunders, AMA, Inc at NASA Ames Research Center, CA.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)    :: nrow, ncol        ! The table is rectangular
   real,    intent (inout) :: table(nrow,ncol)  ! The noise is added in-place
   real,    intent (in)    :: percent           ! % of each column data range
                                                ! defining noise std. deviation
!  Local variables:

   integer :: i, iseed, j
   real    :: cnoise, colmax, colmin, fraction, sdev

!  Procedures:

   real, external :: gauss  ! Standard normal distribn. utility (mu=0, sigma=1)

!  Execution:

   fraction = percent * 0.01
   do j = 1, ncol
      colmax = maxval (table(:,j))
      colmin = minval (table(:,j))
      sdev   = fraction*(colmax - colmin)
      iseed  = -123457           ! < 0 initializes a sequence
      do i = 1, nrow
         cnoise = gauss (iseed)
         table(i,j) = table(i,j) + sdev*cnoise
      end do
   end do

   end subroutine addnoise
