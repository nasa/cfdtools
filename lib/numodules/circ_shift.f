!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE CIRC_SHIFT (N, KEY_IN, KEY_OUT, X)
!
!     Circular shift elements of a real 1-D array so that input element KEY_IN
!     becomes output element KEY_OUT (which is probably 1).  Note that shifting
!     left or right cannot reverse the inherent circular order.  The shift is
!     performed in-place.
!
!     12/27/01  D. Saunders  Initial implementation, for reordering points
!                            formed by slicing a surface triangulation.
!
!     Author:  David Saunders, ELORET/NASA Ames Research Center, Mtn. View, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

!     Arguments:

      INTEGER, INTENT (IN) ::
     >   N,                ! Number of elements
     >   KEY_IN,           ! Input element KEY_IN becomes output element KEY_OUT
     >   KEY_OUT           ! but the inherent circular order remains the same

      REAL, INTENT (INOUT) ::
     >   X(N)              ! Elements circular-shifted left or right, in-place

!     Local variables:

      INTEGER
     >   KDIFF

      REAL
     >   T(N)              ! Always enough (easy way out)

!     Execution:

      KDIFF = KEY_OUT - KEY_IN

      IF (KDIFF > 0) THEN      ! Right shift

         T(1:KDIFF)     = X(N-KDIFF+1:N)
         X(KDIFF+1:N)   = X(1:N-KDIFF)
         X(1:KDIFF)     = T(1:KDIFF)

      ELSE IF (KDIFF < 0) THEN ! Left shift

         KDIFF = -KDIFF
         T(1:KDIFF)     = X(1:KDIFF)
         X(1:N-KDIFF)   = X(1+KDIFF:N)
         X(N-KDIFF+1:N) = T(1:KDIFF)

      END IF

      END SUBROUTINE CIRC_SHIFT
