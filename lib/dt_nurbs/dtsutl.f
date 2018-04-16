      SUBROUTINE DTSUTL (C, IOPT, OUTPUT, IER)

      DOUBLE PRECISION C(*), OUTPUT(*), XMIN, XMAX
      INTEGER IOPT, IER, NDOM, IDOM

      IF (IOPT .EQ. 1) THEN
        CALL DTSUT1 (C, OUTPUT)
        RETURN
        ENDIF
      IF (IOPT .EQ. 2) THEN
        NDOM = C(1)
        DO 10 IDOM = 1,NDOM
          CALL DTSUT2 (C, IDOM, XMIN, XMAX)
          OUTPUT(2*IDOM-1) = XMIN
          OUTPUT(2*IDOM) = XMAX
10        CONTINUE
        RETURN
        ENDIF
      IF (IOPT .EQ. 3) THEN
        CALL DTSUT3 (C, OUTPUT)
        RETURN
        ENDIF
      END

