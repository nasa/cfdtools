      SUBROUTINE D2IGLI (CMEM, IMEM, DMEM, IGI, NSL, NDE, IER)

C     PURPOSE:
C        Fetch IGES Index Length
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IGI      Pointer to Index entity
C
C     OUTPUT:
C        NSL      Number of Start Section Lines
C        NDE      Number of Directory Section Lines
C        IER      Error flag.
C                 -1  = Dynamic memory is corrupt or uninitialized.
C                 -2  = IGI does not point to a valid IGES Index Entity
C
C     CALLS:
C        D0PTR    Pointer check
C        D1FEEI   Fetch Integer Element
C        DTERR    Error handler
C
C
C     HISTORY:
C        16Jun93  D. Parsons   Created.
C
C     ------

C     Long name alias:
C        ENTRY D2_IGES_INDEX_LENGTH (CMEM, IMEM, DMEM, IGI, NSL, NDE, IER)

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      EXTERNAL D1MBAD
      LOGICAL  D1MBAD

C     Get the Reserved Entity Type ID Numbers

C =====================================================================
C
C     INCLUDE FILE DITYID.INC
C
C     Reserved Entity Type ID Numbers
C
C     HISTORY
C        12May92     D. Parsons     created
C
C                          IGES Index
      INTEGER     ENTGI
      PARAMETER  (ENTGI = 237)

C =====================================================================

      INTEGER         IGI, NSL, NDE, IER
      INTEGER         IGIHDR, ITYP, IERX
      CHARACTER*6     SUBNAM
      DATA            SUBNAM /'D2IGLI'/ 
      
      IER  = 0
      IERX = 0
      NSL  = 0
      NDE  = 0

C     Error checking

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
            IER = -1
            GOTO 9000
         ENDIF
      ENDIF
      
      CALL D0PTR (CMEM, IMEM, DMEM, IGI, ENTGI, 0, IGIHDR,
     +    ITYP, IERX)
      IF (IERX .NE. 0) THEN
          IER = -2
          GOTO 9000
      ENDIF 
      
      CALL D1IGLI (CMEM, IMEM, DMEM, IGIHDR, NSL, NDE, IER)

 9000 CONTINUE

      IF (IER .NE. 0) THEN
         CALL DTERR (1, SUBNAM, IER, 0)
         NSL = 0
         NDE = 0
      ENDIF

      RETURN
      END
