      SUBROUTINE DTILC1(XKNOTS,NKNOTS,KORD,IFAIL)
C***********************************************************************
C
C  PURPOSE  DTILC1 CHECKS IF A GIVEN KNOT VECTOR IS IN ASCENDING ORDER
C           WITH NO KNOT HAVING MULTIPLICITY GREATER THAN THE
C           ORDER.
C
C  USAGE    REAL XKNOTS(NKNOTS)
C           CALL DTILC1(XKNOTS,NKNOTS,KORD,IFAIL)
C
C  INPUT    XKNOTS  ARRAY OF NKNOTS VALUES CONTIAINING THE KNOTS.
C
C           NKNOTS  NUMBER OF VALUES IN THE NKNOTS ARRAY.
C
C           KORD    ORDER OF THE B-SPLINES BASED ON THE KNOTS.
C
C  OUTPUT   IFAIL   PARAMETER WHICH INDICATES IF THE KNOT SEQUENCE
C                   IS VALID.
C
C                   IFAIL .EQ. 0  KNOT SEQUENCE VALID.
C
C                   IFAIL .EQ. I  XKNOTS(I) .GT. XKNOTS(I+1) OR THE
C                                 THE MULTIPLICITY OF XKNOTS(I)
C                                 EXCEEDS KORD.
C
C***********************************************************************
C
C
C=======================================================================
C     PARAMETERS
C=======================================================================
C
      DOUBLE PRECISION XKNOTS(*)
      INTEGER IFAIL, KORD, NKNOTS
C
C=======================================================================
C     INTERNAL VARIABLES
C=======================================================================
C
      INTEGER I, IL, IL1, IL2, IS
C
C=======================================================================
C     CHECK ASCENDING ORDER AND MULTIPLICITY
C=======================================================================
C
      IFAIL = 0
      IL = NKNOTS - KORD
      DO 10 I = 1, IL
          IL1 = I + 1
          IL2 = I + KORD
          IF( XKNOTS(I) .LE. XKNOTS(IL1) .AND.
     1        XKNOTS(I) .LT. XKNOTS(IL2) ) GO TO 10
          IFAIL = I
          GO TO 30
10    CONTINUE
      IF( KORD .EQ. 1 ) GO TO 30
      IS = IL + 1
      IL = IL + KORD - 1
      DO 20 I = IS, IL
          IF( XKNOTS(I) .LE. XKNOTS(I+1) ) GO TO 20
          IFAIL = I
          GO TO 30
20    CONTINUE
30    RETURN
      END
