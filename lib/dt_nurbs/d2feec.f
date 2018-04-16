      SUBROUTINE D2FEEC (CMEM, IMEM, DMEM, IDE, JSUB, CELM, IER)

C     PURPOSE:
C        Fetch element of character data from the data space of
C        entity IDE
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IDE      Pointer to entity to fetch
C        JSUB     Relative index into the character space
C
C     OUTPUT:
C        CELM     Character to put data into
C        IER      Error flag.  If negative, CELM is unchanged
C                 -1  = Dynamic memory is corrupt or uninitialized
C                 -2  = Null pointer (IDE = 0)
C                 -3  = Garbage value in pointer (points outside valid
C                       range in IMEM)
C                 -4  = Either a deleted entity or a garbage value in
C                       pointer
C                 -5  = Pointer refers to a deleted entity
C                 -8 =  JSUB <= 0
C                 -9 =  JSUB > LENC
C
C     CALLS:
C        D0PTR, DTERR, D1FEEC
C
C     HISTORY:
C        11Mar92  D. Parsons   Created.
C
C     ------

C     Long name alias:
C        ENTRY D2_FETCH_CMEM_ELEMENT (CMEM, IMEM, DMEM, IDE, JSUB,
C    +         CELM, IER)

      CHARACTER         CMEM*(*), CELM
      INTEGER           IMEM(*), IDE, JSUB
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IHDR, ITYP, IER

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2FEEC'/

C     ------

      IER = 0

C     Call the general pointer-check utility routine

      CALL D0PTR (CMEM, IMEM, DMEM, IDE, 0, 0, IHDR, ITYP, IER)
      IF (IER .EQ. 0) THEN
         CALL D1FEEC (CMEM, IMEM, DMEM, IHDR, JSUB, CELM, IER)
      ENDIF

      IF (IER .NE. 0) CALL DTERR (1, SUBNAM, IER, 0)

      RETURN
      END
