      SUBROUTINE D2IGXI (CMEM, IMEM, DMEM, IGI, MSL, MDE, IER)

C     PURPOSE:
C        IGES Extend Index.  Extend an IGES Index entity.
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IGI      Index entity ID
C        MSL      Amount by which to increase the number of Start Section 
C                 lines allocated
C        MDE      Amount by which to increase the number of Directory Section 
C                 lines allocated, which is twice the number of IGES Entities 
C                 that can be placed in the index.
C
C     OUTPUT:
C        IER      Error flag.  If negative, the entity is not defined.
C                 -1  = Dynamic memory is corrupt or uninitialized
C                 -2  = IGI is a null pointer
C                 -3  = Garbage value found in IGI
C                 -4  = Ambiguous--garbage pointer or deleted entity
C                 -5  = IGI points to a deleted entity
C                 -6  = IGI does not point to a IGES Index
C                       entity
C                 -7  = Character space is locked.
C                 -8  = Integer space is locked.
C                 -9  = Double Precision space is locked (shouldn't 
C                       happen).
C                 -10 = Insufficient space available for character data.
C                 -11 = Insufficient space available for integer data.
C                 -12 = Insufficient space available for double 
C                       precision data.  (shouldn't happen)
C                 -13 = MSL < 0.
C                 -14 = MDE < 0.
C                 -15 = MDE not even.
C                 -999= Unexpected error in a called subroutine.
C
C     CALLS:
C        D0PTR    pointer check
C        D1DEFX   entity definition extend
C        D1FEEI   fetch entity element
C        D1STEI   store entity element
C        DTERR    error handler
C
C
C     HISTORY:
C        21Jun93  D. Parsons   Created.
C
C     ------

C     Long name alias:
C        ENTRY D2_IGES_EXTEND_INDEX (CMEM, IMEM, DMEM, IGI, MSL, MDE, IER)

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

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

      INTEGER     IGI, MSL, MDE, IER
      INTEGER     IERX, IHDR, ITYP, MOREC, MOREI, MORED, NEED, NDE
      LOGICAL     INITIZ
      CHARACTER*6 SUBNAM

      DATA SUBNAM /'D2IGXI'/

C     ------

      IER  = 0
      IERX = 0
      NEED = 0

C     Check pointer and type, then locate relevant data spaces
      CALL D0PTR (CMEM, IMEM, DMEM, IGI, ENTGI, 0, IHDR, ITYP, IER)
      IF (IER .LT. 0) GOTO 9000

      
C     Check the validity of input arguments

      IF (MSL .LT. 0) THEN
         IER = -13
         GOTO 9000
      ENDIF

      IF (MDE .LT. 0) THEN
         IER = -14
         GOTO 9000
      ENDIF

      IF (MOD(MDE,2) .NE. 0) THEN
          IER = -15
          GOTO 9000
      ENDIF                           
      
C     Extend the Character and integer data

      MOREC = MSL*72
      MOREI = MDE/2
      MORED = 0
      INITIZ = .TRUE.

      IF ((MDE .GT. 0) .OR. (MSL .GT. 0)) THEN
         CALL D1DEFX (CMEM, IMEM, DMEM, IHDR, MOREC, MOREI, MORED,
     +      INITIZ, NEED, IERX)
         IF (IERX .NE. 0) THEN
            IF ((IERX .LE. -8) .AND. (IERX .GT. -13)) THEN
               IER = IERX + 1
            ELSE
               IER = -999
            ENDIF
            GOTO 9000
         ENDIF
      ENDIF
      
C     Store the updated NDE into IGI
      IF (MDE .GT. 0) THEN
         CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, 2, NDE, IERX)
         CALL D1STEI (CMEM, IMEM, DMEM, NDE+MDE, 2, IHDR, IERX)
      ENDIF

C     Error Handling

 9000 CONTINUE

      IF (IER .EQ. -999) THEN
         CALL DTERR (5, SUBNAM, IER, 0)
      ELSE IF ((IER .EQ. -10) .OR. 
     +         (IER .EQ. -11) .OR. 
     +         (IER .EQ. -12)) THEN
         CALL DTERR (2, SUBNAM, IER, NEED)
      ELSE IF (IER .NE. 0) THEN
         CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

 9999 CONTINUE

      RETURN
      END
