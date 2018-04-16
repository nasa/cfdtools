      SUBROUTINE D2FEED (CMEM, IMEM, DMEM, IDE, JSUB, DELM, IER)

C     PURPOSE:
C        Fetch element of double precision data from the data
C        space of entity IDE
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IDE      Pointer to entity to fetch
C        JSUB     Relative index into the double precision space
C
C     OUTPUT:
C        DELM     Double precision to put data into
C        IER      Error flag.  If negative, DELM is unchanged
C                 -1  = Dynamic memory is corrupt or uninitialized
C                 -2  = Null pointer (IDE = 0)
C                 -3  = Garbage value in pointer (points outside valid
C                       range in IMEM)
C                 -4  = Either a deleted entity or a garbage value in
C                       pointer
C                 -5  = Pointer refers to a deleted entity
C                 -8 =  JSUB <= 0
C                 -9 =  JSUB > LEND
C
C     CALLS:
C        D0PTR, DTERR, D1FEED
C
C     HISTORY:
C        11Mar92  D. Parsons   Created.
C        29Jun93  D. Parsons   Extracted guts to D1FEED.
C
C     ------

C     Long name alias:
C        ENTRY D2_FETCH_DMEM_ELEMENT (CMEM, IMEM, DMEM, IDE,
C    +         JSUB, DELM, IER)

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*), IDE, JSUB
      DOUBLE PRECISION  DMEM(*), DELM

      INTEGER           IHDR, ITYP, IER

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2FEED'/

C     ------

      IER = 0

C     Call the general pointer-check utility routine

      CALL D0PTR (CMEM, IMEM, DMEM, IDE, 0, 0, IHDR, ITYP, IER)
      IF (IER .EQ. 0) THEN
         CALL D1FEED (CMEM, IMEM, DMEM, IHDR, JSUB, DELM, IER)
      ENDIF

      IF (IER .NE. 0) CALL DTERR (1, SUBNAM, IER, 0)

 9999 CONTINUE

      RETURN
      END
