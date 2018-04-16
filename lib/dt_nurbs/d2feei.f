      SUBROUTINE D2FEEI (CMEM, IMEM, DMEM, IDE, JSUB, IELM, IER)

C     PURPOSE:
C        Fetch element of integer data from the data space of
C        entity IDE
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IDE      Pointer to entity to fetch
C        JSUB     Relative index into the integer space
C
C     OUTPUT:
C        IELM     Integer to put data into
C        IER      Error flag.  If negative, IELM is unchanged
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
C        D1FEEI   Fetch integer element
C        DTERR    Error handler
C
C     HISTORY:
C        11Mar92  D. Parsons   Created.
C        28May93  D. Parsons   Extracted guts to D1FEEI
C
C     ------

C     Long name alias:
C        ENTRY D2_FETCH_IMEM_ELEMENT (CMEM, IMEM, DMEM, IDE, JSUB,
C    +         IELM, IER)

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*), IDE, JSUB, IELM
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IHDR, ITYP, IER

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2FEEI'/

C     ------

      IER = 0

C     Call the general pointer-check utility routine

      CALL D0PTR (CMEM, IMEM, DMEM, IDE, 0, 0, IHDR, ITYP, IER)
      IF (IER .NE. 0) GOTO 9000

C     Get the element from IMEM
      CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, JSUB, IELM, IER)
      IF (IER .NE. 0) THEN
         IER = IER-7
         GOTO 9000
      ENDIF

      RETURN

C     Error return

 9000 CONTINUE

      CALL DTERR (1, SUBNAM, IER, 0)

 9999 CONTINUE

      RETURN
      END
