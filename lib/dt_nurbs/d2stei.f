      SUBROUTINE D2STEI (CMEM, IMEM, DMEM, IELM, JSUB, IDE, IER)

C     PURPOSE:
C        Store element of integer data into the data space of
C        entity IDE
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IELM     Integer to put into IMEM
C        JSUB     Relative index into the integer space
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
C                 -8 =  JSUB <= 0
C                 -9 =  JSUB > LENI
C
C     CALLS:
C        D0PTR    Pointer check
C        D1STEI   Store integer element (shadow)
C        DTERR    Error handler
C
C     HISTORY:
C        12Mar92  D. Parsons   Created.
C        01Jun93  D. Parsons   Extracted guts to D1STEI
C
C     ------

C     Long name alias:
C        ENTRY D2_STORE_IMEM_ELEMENT (CMEM, IMEM, DMEM, IELM,
C    +         JSUB, IDE, IER)

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*), IDE, JSUB, IELM
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IHDR, ITYP, IER

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2STEI'/

C     ------

      IER = 0

C     Call the general pointer-check utility routine

      CALL D0PTR (CMEM, IMEM, DMEM, IDE, 0, 0, IHDR, ITYP, IER)
      IF (IER .NE. 0) GOTO 9000

C     Store the element to IMEM
      CALL D1STEI (CMEM, IMEM, DMEM, IELM, JSUB, IHDR, IER)
      IF (IER .NE. 0) THEN
          IER = IER -7
          GOTO 9000
      ENDIF

      RETURN

C     Error return

 9000 CONTINUE

      CALL DTERR (1, SUBNAM, IER, 0)

 9999 CONTINUE

      RETURN
      END
