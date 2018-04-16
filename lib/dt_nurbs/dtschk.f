      SUBROUTINE DTSCHK ( C, IER )
C
C     ===================================================================
C
C     WRITTEN BY RICH MASTRO
C
C     MODIFIED TO ALLOW RATIONAL SPLINES
C         BY DEBORAH PARSONS
C            5/26/89
C
C     --------------
C     ... PARAMETERS
C     --------------
C
      INTEGER           IER
C
      DOUBLE PRECISION  C(*)
      CHARACTER*8       SUBNAM
C
C     ----------------------
C     ... INTERNAL VARIABLES
C     ----------------------
C
      INTEGER           I            , IDOM         , ILK          ,
     1                  KORD         , NBS          , NDOM
C
C     -------------
C     ... EXTERNALS
C     -------------
C
      DOUBLE PRECISION  DTMCON
C
C     ===================================================================
C
      SUBNAM = 'DTSCHK  '
      IER = 0
C
      IF ( C(1) .LT. 1.0D0 .OR. C(1) .GT. DTMCON(11) ) THEN
C
          IER = -1
          GOTO 9900
C
      END IF
C
      IF ( (C(2) .GT. DTMCON(11)) .OR. (C(2) .LT. -DTMCON(11))
     1     .OR. ((C(2) .GE. -1.0) .AND. (C(2) .LE. 0.0))) THEN
C
          IER = -2
          GOTO 9900
C
      END IF
C
C     --------------------------------------
C     ... LOOP THROUGH INDEPENDENT VARIABLES
C     --------------------------------------
C
      NDOM = INT ( C(1) )
      ILK  = 2 + 3 * NDOM
C
      DO 40 IDOM = 1, NDOM
C
C         -------------------------
C         ... CHECK FOR VALID ORDER
C         -------------------------
C
          IF ( C(2+IDOM) .LT. 1.0D0 .OR.
     1         C(2+IDOM) .GT. DTMCON(11) ) THEN
C
              IER = -3
              GOTO 9900
C
          END IF
C
C         ---------------------------------------
C         ... CHECK FOR VALID NUMBER OF B-SPLINES
C         ---------------------------------------
C
          IF ( C(2+NDOM+IDOM) .LT. C(2+IDOM) .OR.
     1         C(2+NDOM+IDOM) .GT. DTMCON(11) ) THEN
C
              IER = -4
              GOTO 9900
C
          END IF
C
          KORD = INT ( C(2+IDOM) )
          NBS  = INT ( C(2+NDOM+IDOM) )
C
C         ------------------------------------------
C         ... FIRST KNOT MUST HAVE MULTIPLICITY KORD
C         ------------------------------------------
C
          DO 10 I = 1, KORD - 1
C
              IF ( C(ILK+I) .EQ. C(ILK+I+1) ) GO TO 10
C
                  IER = -5
                  GOTO 9900
C
 10       CONTINUE
C
          IF ( C(ILK+KORD) .GE. C(ILK+KORD+1) ) THEN
C
              IER = -5
              GOTO 9900
C
          END IF
C
C         ------------------------------------
C         ... INTERIOR KNOTS MUST BE ASCENDING
C         ------------------------------------
C
          DO 20 I = KORD + 1, NBS
C
              IF ( C(ILK+I) .LE. C(ILK+I+1) .AND.
     1             C(ILK+I) .LT. C(ILK+I+KORD) ) GO TO 20
C
                  IER = -5
                  GOTO 9900
C
 20       CONTINUE
C
          IF ( C(ILK+NBS) .EQ. C(ILK+NBS+1) ) THEN
C
              IER = -5
              GOTO 9900
C
          END IF
C
C         --------------------------------------------
C         ... LAST KNOT MUST HAVE MULTIPLICITY OF KORD
C         --------------------------------------------
C
          DO 30 I = NBS + 1, NBS + KORD - 1
C
              IF ( C(ILK+I) .NE. C(ILK+I+1) ) THEN
C
                  IER = -5
                  GOTO 9900
C
              END IF
C
 30       CONTINUE
C
          ILK = ILK + KORD + NBS
C
 40   CONTINUE
C
C     ==================================================================
C
 9900 CONTINUE
      IF (IER .LT. 0) THEN
          CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF
C
      RETURN
      END
