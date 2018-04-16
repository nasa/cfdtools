      SUBROUTINE D0IGND (CMEM, IMEM, DMEM, IHDR, DMEMX, IER)
C
C     Get next free DMEM position in IGES Entity parameter space
C
C     Subroutine to pop the next free DMEM space from the IGES
C     DMEM free stack.
C
C     INPUT:
C        IHDR     Header index for an IGES entity.
C
C     OUTPUT:
C        DMEMX    Next free DMEM index
C        IER      Error return
C              0     no errors
C              -1    Invalid or corrupted IGES entity
C              -2    no DMEM index positions are available
C
C     Calls:
C        D1FEEI   Fetch the pointer (next free) from the free stack
C        D1STEI   Store the new first free position to the stack
C        D1FEED   Fetch the next pointer from the stack
C
C     04Jul93  D. Parsons     Created.
C ***
      CHARACTER*(*)     CMEM
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IHDR, IER, ISTART, INEXT, DMEMX, IERX
      DOUBLE PRECISION  DNEXT

      IER   = 0
      IERX  = 0
      DMEMX = 0

      CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, 4, ISTART, IERX)
      IF (IERX .NE. 0) THEN
         IER = -1
      ELSE
         IF (ISTART .EQ. 0) THEN
C           No free positions
            IER = -2
         ELSE
            DMEMX = ISTART
            CALL D1FEED (CMEM, IMEM, DMEM, IHDR, ISTART, DNEXT, IERX)
            IF (IERX .NE. 0) THEN
               IER = -1
            ELSE
               INEXT = INT(DNEXT)
               CALL D1STEI (CMEM, IMEM, DMEM, INEXT, 4, IHDR, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -1
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      RETURN
      END
