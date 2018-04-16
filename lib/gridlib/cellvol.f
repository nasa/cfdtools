C**********************************************************************
C
      SUBROUTINE CELLVOL (IL, JL, KL, I1, I2, J1, J2, K1, K2,
     >                    X, Y, Z, VOL, HAND)
C
C     Argument-driven portion of FLO87's METRIC routine for calculating
C     cell volumes, as needed for a grid quality check. HAND = +/-1 for
C     RH/LH xyz resp. VP1, VP3, & VP5 negative signs added empirically.
C
C     c. 1992   David Saunders  Adaptation from FLO87 (A. Jameson).
C     c. 1999     "      "      Fortran 90 internal procedure is not
C                               efficient (on some systems).  Retain
C                               the statement function for now.
C     04/22/04    "      "      Adapted for use on one or more cells.
C
C**********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   IL, JL, KL              ! Grid array dimensions

      INTEGER, INTENT (IN) ::
     >   I1, I2, J1, J2, K1, K2  ! Sub-ranges to calculate volumes for

      REAL, DIMENSION (IL,JL,KL), INTENT (IN) ::
     >   X, Y, Z

      REAL, INTENT (IN) ::
     >   HAND                    ! Enter +1 if you're looking for bad
                                 ! cells in a RH system, or trying to
                                 ! tell if the system is RH or LH
      REAL, INTENT (OUT) ::
     >   VOL(I1+1:I2,J1+1:J2,K1+1:K2)

C     Local constants:

      REAL, PARAMETER ::
     >   FOURTH = 1./ 4., SIXTH = 1./ 6., R8TH = 1./ 8.

C     Local variables:

      INTEGER
     >   I, J, K, L, M, N
      REAL
     >   FACTOR, XP, YP, ZP, VP1, VP2, VP3, VP4, VP5, VP6,
     >   VOLPYM, XA, YA, ZA, XB, YB, ZB, XC, YC, ZC, XD, YD, ZD

C     Statement function (more efficient than internal procedure):

      VOLPYM (XA, YA, ZA, XB, YB, ZB, XC, YC, ZC, XD, YD, ZD) =
     >      (XP - FOURTH * (XA + XB + XC + XD))
     >   * ((YA - YC) * (ZB - ZD) - (ZA - ZC) * (YB - YD))
     >   +  (YP - FOURTH * (YA + YB + YC + YD))
     >   * ((ZA - ZC) * (XB - XD) - (XA - XC) * (ZB - ZD))
     >   +  (ZP - FOURTH * (ZA + ZB + ZC + ZD))
     >   * ((XA - XC) * (YB - YD) - (YA - YC) * (XB - XD))

C     VOLPYM = Volume of a pyramid with four-sided base, times 6
C     (adapted from Jameson code in FLO87).

C     Execution:

      FACTOR = HAND * SIXTH

      DO K = K1 + 1, K2
         N = K - 1
         DO J = J1 + 1, J2
            M = J - 1
            DO I = I1 + 1, I2
               L = I - 1
               XP =(X(I,J,K) +X(I,M,K) +X(I,M,N) +X(I,J,N)
     >             +X(L,J,K) +X(L,M,K) +X(L,M,N) +X(L,J,N))*R8TH
               YP =(Y(I,J,K) +Y(I,M,K) +Y(I,M,N) +Y(I,J,N)
     >             +Y(L,J,K) +Y(L,M,K) +Y(L,M,N) +Y(L,J,N))*R8TH
               ZP =(Z(I,J,K) +Z(I,M,K) +Z(I,M,N) +Z(I,J,N)
     >             +Z(L,J,K) +Z(L,M,K) +Z(L,M,N) +Z(L,J,N))*R8TH
               VP1 = -VOLPYM (X(I,J,K), Y(I,J,K), Z(I,J,K),
     >                        X(I,M,K), Y(I,M,K), Z(I,M,K),
     >                        X(I,M,N), Y(I,M,N), Z(I,M,N),
     >                        X(I,J,N), Y(I,J,N), Z(I,J,N))
               VP2 =  VOLPYM (X(L,J,K), Y(L,J,K), Z(L,J,K),
     >                        X(L,M,K), Y(L,M,K), Z(L,M,K),
     >                        X(L,M,N), Y(L,M,N), Z(L,M,N),
     >                        X(L,J,N), Y(L,J,N), Z(L,J,N))
               VP3 = -VOLPYM (X(I,J,K), Y(I,J,K), Z(I,J,K),
     >                        X(I,J,N), Y(I,J,N), Z(I,J,N),
     >                        X(L,J,N), Y(L,J,N), Z(L,J,N),
     >                        X(L,J,K), Y(L,J,K), Z(L,J,K))
               VP4 =  VOLPYM (X(I,M,K), Y(I,M,K), Z(I,M,K),
     >                        X(I,M,N), Y(I,M,N), Z(I,M,N),
     >                        X(L,M,N), Y(L,M,N), Z(L,M,N),
     >                        X(L,M,K), Y(L,M,K), Z(L,M,K))
               VP5 = -VOLPYM (X(I,J,K), Y(I,J,K), Z(I,J,K),
     >                        X(L,J,K), Y(L,J,K), Z(L,J,K),
     >                        X(L,M,K), Y(L,M,K), Z(L,M,K),
     >                        X(I,M,K), Y(I,M,K), Z(I,M,K))
               VP6 =  VOLPYM (X(I,J,N), Y(I,J,N), Z(I,J,N),
     >                        X(L,J,N), Y(L,J,N), Z(L,J,N),
     >                        X(L,M,N), Y(L,M,N), Z(L,M,N),
     >                        X(I,M,N), Y(I,M,N), Z(I,M,N))
               VOL(I,J,K) = (VP1 + VP2 + VP3 + VP4 + VP5 + VP6) * FACTOR
            END DO
         END DO
      END DO

      END SUBROUTINE CELLVOL
