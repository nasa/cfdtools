      SUBROUTINE D0IGFD (CMEM, IMEM, DMEM, IHDR, IFREE, IER)
C
C     Free a double precision parameter storage location from DMEM.
C
C     Subroutine to push a newly freed DMEM space into the IGES
C     DMEM free stack.
C
C     INPUT:
C        IHDR     Header index for an IGES entity.
C        IFREE    The position in DMEM being freed
C
C     OUTPUT:
C        IER      Error return
C              0     no errors
C              -1    Invalid or corrupted IGES entity
C              -2    INEXT < 1 or > LEND
C
C     Calls:
C        D1FEEI   Fetch the pointer (next free) from the free stack
C        D1STEI   Store the newly freed position to the free stack
C        D1STED   Store the old pointer (next free) to the stack
C
C     04Jul93  D. Parsons     Created.
C ***
      CHARACTER*(*)     CMEM
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER  IHDR, IFREE, IER, ISTART, IERX

      IER  = 0
      IERX = 0

      CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, 4, ISTART, IERX)
      IF (IERX .NE. 0) THEN
         IER = -1
      ELSE
         CALL D1STEI (CMEM, IMEM, DMEM, IFREE, 4, IHDR, IERX)
         CALL D1STED (CMEM, IMEM, DMEM, ISTART, IFREE, IHDR, IERX)
         IF (IERX .NE. 0) THEN
            IER = -2
         ENDIF
      ENDIF

      RETURN
      END
