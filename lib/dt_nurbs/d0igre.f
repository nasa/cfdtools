      SUBROUTINE D0IGRE (CMEM, IMEM, DMEM, LUNIT, IENT, IGI, IETYPE, 
     +                   IDVAL, CDVAL, IPARST, IER)

C     Read an IGES entity into the DTRC Dynamic Memory, without
C     regard to its Entity Type (for now).
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        LUNIT    The I/O unit to use when working on the IGES file. 
C        IENT     The Entity Number (always odd).
C        IGI      The IGES Index with which this entity belongs.
C        IETYPE   The entity type number.
C        IDVAL    Integer parameters from the Directory Section
C        CDVAL    Character parameters from the Directory Section
C        IPARST   The starting line number of all Parameter section.
C
C     OUTPUT:
C        IER      Returned error value.  Zero implies no errors.
C                 Otherwise,
C                 -n An error occurred while reading or processing line 
C                    'n' of the IGES File.
C
C   CALLS:
C
C   HISTORY:
C     8/5/93      D. Parsons   Created
C*****

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)
      
      INTEGER        LUNIT, IENT, IGI, IETYPE, IDVAL(15), IPARST, IER
      CHARACTER*8    CDVAL(3)
      
      CHARACTER*(80) LINE, STRING
      INTEGER        IERX, IOCHK, I
      INTEGER        NPAR, NSTR, NSTRCH, NREAL, IGE
      INTEGER        IPLAST, IPCURP, ICURLN, IVAL, ILEN, IPAR
      DOUBLE PRECISION  RVAL
      CHARACTER      EOP, EOR, TYPE
      LOGICAL        EORFLG 

C****
      IER  = 0
      IERX = 0

      ICURLN = IPARST+IDVAL(2)-1
      IPLAST = ICURLN+IDVAL(13)-1
      IPCURP = 1
            
C     Allocate an IGES Entity
C        Initial guess at parameter counts:
C           NPAR = 5 * number of parameter lines.
C           NREAL = NPAR
C           NSTR = 0
C           NSTRCH = 0 

C        The STORE routines will expand these values, if necessary.

      IETYPE = IDVAL(1)
      
      NPAR = IDVAL(13)
      NSTR = 0
      NSTRCH = 0
      NREAL  = NPAR
      
      CALL D2IGDE (CMEM, IMEM, DMEM, NPAR, NSTR, NSTRCH, NREAL,
     +                   IGI, IENT, IGE, IERX) 
      IF (IERX .NE. 0) THEN
         IER = -ICURLN
         GOTO 9000
      ENDIF
      
C     Store the Directory information to the new entity

      DO 10 I = 1, 19
         IF (( I .EQ. 1) .AND. (IDVAL(1) .NE. IETYPE)) THEN
            IER = -ICURLN
            GOTO 9000
         ENDIF
         IF (I .LE. 10) THEN
            CALL D2IGSI (CMEM, IMEM, DMEM, 'I', IDVAL(I), IGE, 'D',
     +         I, IERX)
         ELSE IF (I .EQ. 11) THEN
C           (Skipped parameter)
         ELSE IF (I .LE. 15) THEN
            CALL D2IGSI (CMEM, IMEM, DMEM, 'I', IDVAL(I-1), IGE, 'D',
     +         I, IERX)
         ELSE IF (I .LE. 18) THEN
            CALL D2IGSC (CMEM, IMEM, DMEM, 'C', CDVAL(I-15), IGE, 'D',
     +         I, IERX)
         ELSE
            CALL D2IGSI (CMEM, IMEM, DMEM, 'I', IDVAL(I-4), IGE, 'D',
     +         I, IERX)
         ENDIF
         IF (IERX .NE. 0) THEN
            IER = -ICURLN
            GOTO 9000
         ENDIF
   10 CONTINUE
         
C     Read and store the Parameter Lines

      READ (LUNIT, '(A80)', REC=ICURLN, IOSTAT=IOCHK) LINE
      IF ((LINE(73:73) .NE. 'P') .OR. (IOCHK .NE. 0)) THEN
         IER = -ICURLN
         GOTO 9000
      ENDIF
      
      CALL D2IGFC (CMEM, IMEM, DMEM, IGI, 'G', 1, EOP, 1, IERX)
      CALL D2IGFC (CMEM, IMEM, DMEM, IGI, 'G', 2, EOR, 1, IERX)
      
      CALL D0IGRP (LINE, LUNIT, ICURLN, IPCURP, EOP, EOR, 
     +             EORFLG, TYPE, IVAL, RVAL, STRING, ILEN)
     
C     Verify the entity type number one last time.

      IF ((TYPE .NE. 'I') .OR. (IVAL .NE. IETYPE)) THEN
         IER = -ICURLN
         GOTO 9000
      ENDIF

C     The entity type number is not stored again; 
C        start with parameter 1

      IPAR = 1

C     (WHILE .NOT. EOR...)

      IF (EORFLG) GOTO 9000
      
  100 CONTINUE
         IF (ICURLN .GT. IPLAST) THEN
            IER = -ICURLN
            GOTO 9000
         ENDIF
         
         CALL D0IGRP (LINE, LUNIT, ICURLN, IPCURP, EOP, EOR, 
     +                EORFLG, TYPE, IVAL, RVAL, STRING, ILEN)
     
         IF (TYPE .EQ. 'C') THEN
            CALL D2IGSC (CMEM, IMEM, DMEM, TYPE, STRING(1:ILEN), IGE, 
     +                   'P', IPAR, IERX)
         ELSE IF (TYPE .EQ. 'R') THEN
            CALL D2IGSD (CMEM, IMEM, DMEM, 'D', RVAL, IGE, 'P',
     +                   IPAR, IERX)
         ELSE
            CALL D2IGSI (CMEM, IMEM, DMEM, TYPE, IVAL, IGE, 'P',
     +                   IPAR, IERX)
         ENDIF
         IF (IERX .NE. 0) THEN
            IER = -ICURLN
            GOTO 9000
         ENDIF 
         
         IPAR = IPAR+1
         
      IF (.NOT. EORFLG) GOTO 100
      
 9000 CONTINUE
      RETURN
      END
