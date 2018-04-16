      SUBROUTINE D2IGXE (CMEM, IMEM, DMEM, IGE, MPAR, MSTR, MSTRCH, 
     +                   MREAL, IER)

C     PURPOSE:
C        IGES Extend Entity.  Extend an IGES entity.
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IGE      IGES entity ID
C        MPAR     Amount by which to increase number of Parameters
C        MSTR     Amount by which to increase number of string 
C                 parameters to allocate
C        MSTRCH   Amount by which to increase total number of 
C                 characters in strings
C        MREAL    Amount by which to increase number of double 
C                 precision parameters
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
C                 -9  = Double Precision space is locked.
C                 -10 = Insufficient space available for character data.
C                 -11 = Insufficient space available for integer data.
C                 -12 = Insufficient space available for double
C                       precision data.
C                 -13 = MPAR < 0.
C                 -14 = MSTR < 0.
C                 -15 = MSTRCH < 0.
C                 -16 = MREAL < 0.
C                 -17 = Entity does not point to a valid Character Sequence
C                       Entity
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
C        06Aug93  D. Parsons   Add support of DMEM "Free" stack.
C
C     ------

C     Long name alias:
C        ENTRY D2_IGES_EXTEND_ENTITY (CMEM, IMEM, DMEM, IGE, MPAR, MSTR, 
C     +                   MSTRCH, MREAL, IER)

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
C                          IGES Entity Data
      INTEGER     ENTGE
      PARAMETER  (ENTGE = 238)

C =====================================================================

      INTEGER     IGE, MPAR, MSTR, MSTRCH, MREAL, IER
      INTEGER     IERX, IGEHDR, ITYP, NEED
      CHARACTER*6 SUBNAM

      DATA SUBNAM /'D2IGXE'/

C     ------

      IER  = 0
      IERX = 0
      NEED = 0

C     Check pointer and type, then locate relevant data spaces
      CALL D0PTR (CMEM, IMEM, DMEM, IGE, ENTGE, 0, IGEHDR, ITYP, IER)
      IF (IER .LT. 0) GOTO 9000

C     Check the validity of input arguments

      IF (MPAR .LT. 0) THEN
         IER = -13
         GOTO 9000
      ENDIF

      IF (MSTR .LT. 0) THEN
         IER = -14
         GOTO 9000
      ENDIF

      IF (MSTRCH .LT. 0) THEN
         IER = -15
         GOTO 9000
      ENDIF

      IF (MREAL .LT. 0) THEN
         IER = -16
         GOTO 9000
      ENDIF

      CALL D1IGXE (CMEM, IMEM, DMEM, IGEHDR, MPAR, MSTR, MSTRCH, 
     +                   MREAL, NEED, IER)
     
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
