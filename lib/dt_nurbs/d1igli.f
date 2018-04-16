      SUBROUTINE D1IGLI (CMEM, IMEM, DMEM, IGIHDR, NSL, NDE, IER)

C     PURPOSE:
C        Fetch IGES Index Length (shadow routine)
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IGIHDR   Header index of IGES INDEX entity
C
C     OUTPUT:
C        NSL      Number of Start Section Lines
C        NDE      Number of Directory Section Lines
C        IER      Error flag.  (none returned)
C
C     CALLS:
C        D1FEEI   Fetch Integer Element
C
C
C     HISTORY:
C        17Jun93  D. Parsons   Created.
C
C     ------

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER         NSL, NDE, IER
      INTEGER         IGIHDR, IERX, LENC
      
      IER  = 0
      IERX = 0
      NSL  = 0
      NDE  = 0

C     Compute the NSL from the length of the allocated CMEM

      LENC = ABS(IMEM(IGIHDR+4))
      NSL = (LENC - 3)/72
      
C     Fetch the NDE from the IMEM space.
      
      CALL D1FEEI (CMEM, IMEM, DMEM, IGIHDR, 2, NDE, IERX)
      IF (IERX .NE. 0) NDE = 0

 9000 CONTINUE

      IF (IER .NE. 0) THEN
         NSL = 0
         NDE = 0
      ENDIF

      RETURN
      END
