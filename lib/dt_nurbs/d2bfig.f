      SUBROUTINE D2BFIG (CMEM, IMEM, DMEM, IBF, IGI, JDE, IGE, IER)
      
C   Convert a B-Spline Function Entity to an IGES Entity 
C   type 126 (B-Spline Curve) or 128 (B-Spline Surface).
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM     Dynamically managed character data
C     IMEM     Dynamically managed integer data
C     DMEM     Dynamically managed double precision data
C
C   INPUT:
C     IBF      Entity Pointer to B-Spline Function Entity.
C     IGI      Pointer to IGES Index Entity.
C     JDE      DE number for this new IGES Entity.
C
C   OUTPUT:
C     IGE      Pointer to the IGES Spline Entity (type 126 or 128)
C     IER      Returned error value.  Zero implies no errors.  
C              -1  = Dynamic memory is corrupt or uninitialized
C              -2  = IBF is a null pointer
C              -3  = Garbage value found in IBF
C              -4  = Ambiguous--IBF is garbage pointer or deleted entity
C              -5  = IBF points to a deleted entity
C              -6  = IBF does not point to a B-spline Function entity.
C              -7  = IGI is a null pointer
C              -8  = Garbage value found in IGI
C              -9  = Ambiguous--IGI is garbage pointer or deleted entity
C              -10 = IGI points to a deleted entity
C              -11 = IGI does not point to an IGES Index entity.
C              -12 = The JDEth directory entry is not empty.
C              -13 = JDE < 1 or JDE > NDE.
C              -14 = Number of independent variables < 1 or > 2
C              -15 = Number of dependent variables < 1 or > 3
C              -16 = Rational Spline; Encountered a Weight <= 0.0
C              -999 = Error in a lower-level subroutine
C
C   CALLS:
C     D0PTR
C     DTERR
C     D2IGFI
C     D2IGDE 
C     D2IGSI
C     D2IPSI 
C     D2IPSD
C     D2IGSD
C
C   HISTORY:
C     9/20/93  D. Parsons  Created
C
C*****
C
C   Long Name Alias:
C     ENTRY D2_BSPLINE_TO_IGES (CMEM, IMEM, DMEM, IBF, IGI, JDE, IGE, IER)

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      
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

C                          B-spline Function
      INTEGER     ENTBF
      PARAMETER  (ENTBF = 246)

C =====================================================================

      INTEGER  IGI, JDE, IGE, IBF, IER
      INTEGER  IERX
      INTEGER  IHDR, IGEHDR, ITYP, JBLO, JBHI
       
      CHARACTER*6    SUBNAM
      DATA SUBNAM   /'D2BFIG'/
      
C****

      IER = 0
      IERX = 0
      
C     Check pointer and type, then locate relevant data spaces
      CALL D0PTR (CMEM, IMEM, DMEM, IBF, ENTBF, 0, IHDR, ITYP, IER)
      IF (IER .LT. 0) GOTO 9900

C     Check pointer and type, then locate relevant data spaces
      CALL D0PTR (CMEM, IMEM, DMEM, IGI, ENTGI, -5, IHDR, ITYP, IER)
      IF (IER .LT. 0) GOTO 9900
      
C     Fetch the bounds of the "C-Array"

      CALL D2FEBD (CMEM, IMEM, DMEM, IBF, JBLO, JBHI, IERX)
      
C     Do the conversion

      CALL D2CYIG (CMEM, IMEM, DMEM, DMEM(JBLO), IGI, JDE, IGE, IERX)
      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9900
      ENDIF
      
C     Unlock the B-Spline Function Entity

      CALL D2UNLD (CMEM, IMEM, DMEM, IBF, IERX)
      
C     Generate the label for the Entity

      CALL D0LBIG (CMEM, IMEM, DMEM, IBF, IGE)
      
C     Error reporting section

 9900 CONTINUE

      IF (IER .NE. 0) THEN
         IF (IER .EQ. -999) THEN
            CALL DTERR (5, SUBNAM, IER, 0) 
         ELSE
            CALL DTERR (1, SUBNAM, IER, 0) 
         ENDIF
          
C        Try to erase the entity; ignore any error

         IF (IGE .NE. 0) THEN
            IGEHDR = IGE/256
            CALL D0GEER (CMEM, IMEM, DMEM, IGEHDR, 0, IERX)
         ENDIF
         IGE = 0
      ENDIF

      RETURN
      END
    
