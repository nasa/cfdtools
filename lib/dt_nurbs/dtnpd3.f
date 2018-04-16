C<*>
      SUBROUTINE DTNPD3 ( X,   KNOTS,  NBS,    KORD,   IPOS,   IER    )
C     ==================================================================
C
C     --------------
C     ... PARAMETERS
C     --------------
C
      INTEGER      NBS,      KORD,     IPOS,     IER
C
      DOUBLE PRECISION       X,        KNOTS(*)
C
C     ==================================================================
C
C     ---------------------
C     ... SET IPOS IN RANGE
C     ---------------------
C
      IER    = 0
      IPOS   = MAX0 ( KORD, IPOS )
      IPOS   = MIN0 ( NBS,  IPOS )
C
C     --------------------
C     ... CHECK X IN RANGE
C     --------------------
C
      IF ( X .LT. KNOTS(KORD) .OR. X .GT. KNOTS(NBS+1) ) THEN
C
          IER = -50
          RETURN
C
      END IF
C
C     ----------------------------------
C     ... CHECK IF X EQUALS KNOTS(NBS+1)
C     ----------------------------------
C
      IF ( X .EQ. KNOTS(NBS+1) ) THEN
C
          IPOS = NBS
  5       CONTINUE
          IF ( KNOTS(IPOS) .LT. KNOTS(IPOS+1) ) RETURN
          IPOS = IPOS - 1
          GOTO 5
C
      END IF
C
C     ------------------------------------------------
C     ... KNOTS(KORD) .LE. X .AND. X .LT. KNOTS(NBS+1)
C     ------------------------------------------------
C
      IF ( X .LT. KNOTS(IPOS) ) THEN
C
C     ----------------------
C     ... X .LT. KNOTS(IPOS)
C     ----------------------
C
 10       IPOS = IPOS - 1
C
          IF ( X .GE. KNOTS(IPOS) ) THEN
C
              RETURN
C
          END IF
          GO TO 10
C
      ELSE
C
C     ----------------------
C     ... X .GE. KNOTS(IPOS)
C     ----------------------
C
 20       IPOS = IPOS + 1
C
          IF ( X .LT. KNOTS(IPOS) ) THEN
C
              IPOS   = IPOS - 1
              RETURN
C
          END IF
          GO TO 20
C
      END IF
C
C     ==================================================================
C
C     ----------
C     ... RETURN
C     ----------
C
      END
