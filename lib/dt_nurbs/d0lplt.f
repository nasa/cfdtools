      SUBROUTINE D0LPLT (CMEM, IMEM, DMEM, LP)

C     Locate pointers in Library data types.  Used by D2WRIT.  Not permitted
C     to call any other Library subroutines because D2WRIT temporarily 
C     "corrupts" the header blocks.
C
C   DYNAMIC MEMORY INPUT-OUTPUT
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT-OUTPUT:
C     LP      Array of context data.
C             LP(1) = Header index for entity
C             LP(2) = MEM data type number of entity.
C             LP(3) = 0  Initial call for this entity.
C                   > 0  Beginning absolute index of previous pointer block.
C             LP(4) = Ending absolute index of previous pointer block.
C             LP(5) = Absolute index of end of entity's integer data space.
C             Upon return, LP(3) and LP(4) have been updated to locate the
C             beginning and end of the next pointer block in the entity's
C             integer data space, or to LP(5)+1 if no more pointers.
C
C   CALLS:    none
C
C   HISTORY:
C     20Sep93 P. Kraushar   Created
C*****
      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER LP(14)

C     LIBTYP is the lowest type number used by the Library, minus one.
      INTEGER LIBTYP
      PARAMETER (LIBTYP=233)

      INTEGER I
C****
      IF (LP(2) .LE. LIBTYP) THEN
C        Unknown type (should never happen)
         LP(3) = -1
         LP(4) = LP(5)+1
         RETURN
      ENDIF
C     Begin Case:
C     type:  DM, IM, CM, GI, GE, PJ, JS, TS, LP, EG, CS, BQ, BF,
C           234,235,236,237,238,239,240,241,242,243,244,245,246,
      GOTO (100,100,100,237,238,239,240,241,242,243,244,100,100,
C            DQ, IQ, CQ, DB, IB, CB, DA, IA, CA
C           247,248,249,250,251,252,253,254,255
     +      100,100,100,100,100,100,100,100,100), LP(2)-LIBTYP

C     Case of no internal pointers to other entities
  100 CONTINUE
      LP(3) = LP(5) + 1
      LP(4) = LP(5) + 1
      GOTO 999

C     IGES Index
  237 CONTINUE
      IF (LP(3) .EQ. 0) THEN
C        Use LP(6) to record the end of the second block of pointers
         I = IMEM(LP(4)+2)
         IF (I .GT. 0) THEN
            LP(6) = MIN (LP(4) + 2 + I, LP(5))
         ELSE IF (I .EQ. 0) THEN
            LP(6) = 0
         ELSE
            LP(6) = LP(5)
         ENDIF
         LP(3) = LP(4) + 1
         LP(4) = LP(4) + 1
      ELSE IF (LP(4) .LT. LP(6)) THEN
         LP(3) = LP(4) + 2
         LP(4) = LP(6)
      ELSE
         LP(3) = LP(5) + 1
         LP(4) = LP(5) + 1
      ENDIF
      GOTO 999

C     IGES Entity
  238 CONTINUE
      IF (LP(3) .EQ. 0) THEN
C        LP(9) is index to 0th position of integer data space
         LP(9) = LP(4)
         LP(3) = LP(9) + 1
         LP(4) = LP(9) + 2
C        LP(6) is last scanned parameter
         LP(6) = 0
C        LP(7) is number of last parameter
         LP(7) = IMEM(LP(9)+3)
C        LP(8) is index to 0th parameter TYPCOD in character data space
         LP(8) = IMEM(LP(1)+1) + 23
C        LP(9) changes to index of 0th parameter value or ref in integer space
         LP(9) = LP(9) + 19
      ELSE IF (LP(6) .LT. LP(7)) THEN
C        Scan parameter TYPCODs for the next T, if any.
C        START LOOP
 2380       CONTINUE
            LP(6) = LP(6) + 1
            I = LP(8) + LP(6)
            IF (LP(6) .LT. LP(7) .AND. CMEM(I:I) .NE. 'T') GOTO 2380
C        END LOOP
         IF (LP(6) .LT. LP(7)) THEN
            LP(3) = LP(9) + LP(6)
C           START LOOP
 2381          CONTINUE
               LP(6) = LP(6) + 1
               I = LP(8) + LP(6)
               IF (LP(6) .LT. LP(7) .AND. CMEM(I:I) .EQ. 'T') GOTO 2381
C           END LOOP
            IF (LP(6) .LT. LP(7)) THEN
               LP(4) = LP(9) + LP(6) - 1
            ELSE IF (CMEM(I:I) .EQ. 'T' .OR. CMEM(I:I) .EQ. 'R') THEN
C              Also test if last parameter is TYPCOD R.
               LP(4) = LP(9) + LP(6)
            ELSE
               LP(4) = LP(9) + LP(6) - 1
            ENDIF
         ELSE IF (CMEM(I:I) .EQ. 'T' .OR. CMEM(I:I) .EQ. 'R') THEN
C           Also test if last parameter is TYPCOD R.
            LP(3) = LP(9) + LP(6)
            LP(4) = LP(3)
         ELSE
            LP(3) = LP(5) + 1
            LP(4) = LP(5) + 1
         ENDIF
      ELSE
         LP(3) = LP(5) + 1
         LP(4) = LP(5) + 1
      ENDIF
      GOTO 999

C     Point on Joined Surface
  239 CONTINUE
      IF (LP(3) .EQ. 0) THEN
         LP(3) = LP(4) + 1
         LP(4) = LP(4) + 2
      ELSE
         LP(3) = LP(5) + 1
         LP(4) = LP(5) + 1
      ENDIF
      GOTO 999

C     Joined Surface
  240 CONTINUE
      IF (LP(3) .EQ. 0) THEN
         I = IMEM(LP(4)+2)
         IF (I .GT. 0) THEN
            LP(3) = LP(4) + 3
            LP(4) = MIN (LP(4) + 2 + I, LP(5))
         ELSE IF (I .EQ. 0) THEN
            LP(3) = LP(5) + 1
            LP(4) = LP(5) + 1
         ELSE
            LP(3) = LP(4) + 3
            LP(4) = LP(5)
         ENDIF
      ELSE
         LP(3) = LP(5) + 1
         LP(4) = LP(5) + 1
      ENDIF
      GOTO 999

C     Trimmed Surface of a Joined Surface
  241 CONTINUE
      IF (LP(3) .EQ. 0) THEN
C        Use LP(6) to record the end of the second block of pointers
         I = IMEM(LP(4)+5)
         IF (I .GT. 0) THEN
            LP(6) = MIN (LP(4) + 5 + I, LP(5))
         ELSE IF (I .EQ. 0) THEN
            LP(6) = 0
         ELSE
            LP(6) = LP(5)
         ENDIF
         LP(3) = LP(4) + 1
         LP(4) = LP(4) + 2
      ELSE IF (LP(4) .LT. LP(6)) THEN
         LP(3) = LP(4) + 4
         LP(4) = LP(6)
      ELSE
         LP(3) = LP(5) + 1
         LP(4) = LP(5) + 1
      ENDIF
      GOTO 999

C     Loop of a Trimmed Surface
  242 CONTINUE
      IF (LP(3) .EQ. 0) THEN
C        Use LP(6) to record the end of the second block of pointers
         I = IMEM(LP(4)+3)
         IF (I .GT. 0) THEN
            LP(6) = MIN (LP(4) + 3 + I, LP(5))
         ELSE IF (I .EQ. 0) THEN
            LP(6) = 0
         ELSE
            LP(6) = LP(5)
         ENDIF
         LP(3) = LP(4) + 1
         LP(4) = LP(4) + 1
      ELSE IF (LP(4) .LT. LP(6)) THEN
         LP(3) = LP(4) + 3
         LP(4) = LP(6)
      ELSE
         LP(3) = LP(5) + 1
         LP(4) = LP(5) + 1
      ENDIF
      GOTO 999

C     Edge of a Trimmed Surface
  243 CONTINUE
      IF (LP(3) .EQ. 0) THEN
         LP(3) = LP(4) + 1
         LP(4) = LP(4) + 2
      ELSE IF (LP(4) .LT. LP(5)) THEN
         LP(3) = LP(5)
         LP(4) = LP(5)
      ELSE
         LP(3) = LP(5) + 1
         LP(4) = LP(5) + 1
      ENDIF
      GOTO 999

C     Curve on a Surface
  244 CONTINUE
C     This type not yet fully defined
      IF (LP(3) .EQ. 0) THEN
         LP(3) = LP(4) + 1
         LP(4) = LP(5)
      ELSE
         LP(3) = LP(5) + 1
         LP(4) = LP(5) + 1
      ENDIF
      GOTO 999

C     END CASE
  999 CONTINUE
      RETURN
      END
