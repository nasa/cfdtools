      LOGICAL FUNCTION D2TSTP (CMEM, IMEM, DMEM, IDE)

C     PURPOSE:
C        Test pointer validity.
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IDE      Pointer to be checked
C
C     OUTPUT:
C        D2TSTP   True -> pointer seems valid
C                 False-> pointer is invalid (no reason returned)
C
C     CALLS:
C        D0PTR
C
C     HISTORY:
C        12Feb92  D. Parsons   Created.
C
C     ------

C     Long name alias:
C        ENTRY D2_TEST_POINTER (CMEM, IMEM, DMEM, IDE)

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*), IDE
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IHDR, ITYP, IER

C     ------

C     Call the general pointer-check utility routine

      CALL D0PTR (CMEM, IMEM, DMEM, IDE, -1, 0, IHDR, ITYP, IER)

      D2TSTP = (IER .EQ. 0)

      RETURN
      END
