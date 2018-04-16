      SUBROUTINE D0CQNF (CMEM, IMEM, DMEM, IHDR, ICQNDX, IER)
C
C     Find the next free (deleted or unitialized) position in a
C     Character Sequence Entity.  If none found, ICQNDX is set to
C     point to 1 greater than the number of Sequences currently
C     allocated.
C
C     INPUT
C        IHDR     Header index for a Character Sequence entity
C
C     OUTPUT
C        ICQNDX   The next free Sequence (may be a new Sequence
C                 if all are currently in use.)
C
C        IER      Error code
C                  0 no problems
C                 -1 Invalid or corrupted Character Sequence entity
C
C     Calls
C        D1FEEI   Fetch sequence pointers
C
C     04Jul93  D. Parsons  Created
C ***

      CHARACTER*(*)     CMEM
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER  IHDR, ICQNDX, IER, IERX
      INTEGER  MXSEQ, I, ENDPOS

      IER    = 0
      IERX   = 0
      ICQNDX = 0

      CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, 1, MXSEQ, IERX)
      IF (IERX .NE. 0) THEN
         IER = -1
         GOTO 9999
      ENDIF

      DO 10, I = 1, MXSEQ
         CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, MXSEQ+1, ENDPOS, IERX)
         IF (IERX .NE. 0) THEN
            IER = -1
            GOTO 9999
         ENDIF
         IF (ENDPOS .LT. 0) THEN
            ICQNDX = I
            GOTO 9999
         ENDIF
   10 CONTINUE

      ICQNDX = MXSEQ+1

 9999 CONTINUE
      RETURN
      END

