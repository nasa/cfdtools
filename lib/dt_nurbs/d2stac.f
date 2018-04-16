      SUBROUTINE D2STAC (CMEM, IMEM, DMEM, CA, JLO, JHI, IDE, IER)

C     PURPOSE:
C        Store subarray of character data into the data space of
C        entity IDE
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        CA       Character array to put into CMEM
C        JLO      Relative index into the character space
C        JHI      Relative index into the character space
C        IDE      Pointer to entity to store
C
C     OUTPUT:
C        IER      Error flag.  If negative, the space remains unchanged
C                 -1  = Dynamic memory is corrupt or uninitialized
C                 -2  = Null pointer (IDE = 0)
C                 -3  = Garbage value in pointer (points outside valid
C                       range in IMEM)
C                 -4  = Either a deleted entity or a garbage value in
C                       pointer
C                 -5  = Pointer refers to a deleted entity
C                 -8  = JLO <= 0
C                 -9  = JHI > LENC
C                 -10 = JLO > JHI
C
C     CALLS:
C        D0PTR    Pointer check
C        D1STAC   Store character array (shadow)
C        DTERR    Error handler
C
C     HISTORY:
C        11Mar92  D. Parsons   Created.
C        28May93  D. Parsons   Extracted guts to D1STAC
C
C     ------

C     Long name alias:
C        ENTRY D2_STORE_CMEM_SUBARRAY (CMEM, IMEM, DMEM, CA, JLO,
C    +         JHI, IDE, IER)

      CHARACTER         CMEM*(*), CA*(*)
      INTEGER           IMEM(*), IDE, JLO, JHI
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IHDR, ITYP, IER

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2STAC'/

C     ------

      IER = 0

C     Call the general pointer-check utility routine

      CALL D0PTR (CMEM, IMEM, DMEM, IDE, 0, 0, IHDR, ITYP, IER)
      IF (IER .NE. 0) GOTO 9000

      CALL D1STAC (CMEM, IMEM, DMEM, CA, JLO, JHI, IHDR, IER)
      IF (IER .NE. 0) THEN
         IER = IER - 7
         GOTO 9000
      ENDIF

      GOTO 9999

C     Error return

 9000 CONTINUE

      CALL DTERR (1, SUBNAM, IER, 0)

 9999 CONTINUE

      RETURN
      END
