      SUBROUTINE D2STED (CMEM, IMEM, DMEM, DELM, JSUB, IDE, IER)

C     PURPOSE:
C        Store element of double precision data into the data space of
C        entity IDE
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        DELM     Double precision to put into DMEM
C        JSUB     Relative index into the double precision space
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
C                 -9 =  JSUB > LEND
C
C     CALLS:
C        D0PTR, DTERR
C
C     HISTORY:
C        12Mar92  D. Parsons   Created.
C
C     ------

C     Long name alias:
C        ENTRY D2_STORE_DMEM_ELEMENT (CMEM, IMEM, DMEM, DELM, JSUB,
C    +         IDE, IER)

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*), IDE, JSUB
      DOUBLE PRECISION  DMEM(*), DELM

      INTEGER           IHDR, ITYP, IER

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2STED'/

C     ------

      IER = 0

C     Call the general pointer-check utility routine

      CALL D0PTR (CMEM, IMEM, DMEM, IDE, 0, 0, IHDR, ITYP, IER)
      IF (IER .EQ. 0) THEN
         CALL D1STED (CMEM, IMEM, DMEM, DELM, JSUB, IHDR, IER)
      ENDIF

C     Error return

 9000 CONTINUE

      IF (IER .NE. 0) CALL DTERR (1, SUBNAM, IER, 0)

      RETURN
      END
