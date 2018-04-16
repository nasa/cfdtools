C+----------------------------------------------------------------------
C
      SUBROUTINE DISTPTOL (N, P0, P1, P2, T, DIST)
C
C     One-liner:  Distance from a point to a line in N-space.
C
C     DISTPTOL computes the perpendicular distance from point P0 to the
C     line defined by points P1 and P2, using the method which is valid
C     for any number of dimensions N >= 2.
C
C     Method:
C
C     If P is the foot of the perpendicular, we have
C
C        P = P1 + t (P2 - P1)  or  P1 + t C  where  C = P2 - P1
C
C     Then, from the fact that (P0 - P)'(P2 - P1) = 0, we can show
C
C            C'D
C        t = ---       where  D = P0 - P1  and  ' = transpose
C            C'C
C
C     and the squared distance is given by
C                                                        (C'D) ** 2
C        || P0 - P || ** 2  =  (P0 - P)'(P0 - P) = D'D - ----------
C                                                           C'C
C     History:
C
C     02/15/98  DAS  Initial implementation.
C
C     Author:  D.A.Saunders  Sterling Software/NASA Ames, Mtn. View, CA
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER N        ! I Number of dimensions; N >= 2
      REAL    P0(N)    ! I Target point coordinates
      REAL    P1(N),   ! I Two points defining target line;
     >        P2(N)    !   P1 .NE. P2 is assumed
      REAL    T        ! O Distance from P1 of the foot of the perpendicular
                       !   along the line through P1 and P2
      REAL    DIST     ! O The shortest distance from P0 to this line

C     Local  constants:

      REAL, PARAMETER :: ZERO = 0.

C     Local variables:

      INTEGER I
      REAL    CI, DI, CTC, CTD, DTD

C     Execution:

      CTC = ZERO
      CTD = ZERO
      DTD = ZERO

      DO I = 1, N
         CI  = P2(I) - P1(I)
         DI  = P0(I) - P1(I)
         CTC = CTC + CI * CI
         CTD = CTD + CI * DI
         DTD = DTD + DI * DI
      END DO

      T = CTD / CTC
      DIST = SQRT (MAX (DTD - T * CTD, ZERO))

      END SUBROUTINE DISTPTOL
