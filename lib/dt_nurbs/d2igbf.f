      SUBROUTINE D2IGBF (CMEM, IMEM, DMEM, IGE, IBF, IER)
      
C   Convert an IGES Entity type 126 (B-Spline Curve) or 128 (B-Spline
C   Surface) to a B-Spline Function Entity.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM     Dynamically managed character data
C     IMEM     Dynamically managed integer data
C     DMEM     Dynamically managed double precision data
C
C   INPUT:
C     IGE      Pointer to the IGES Spline Entity (type 126 or 128)
C
C   OUTPUT:
C     IBF      Entity Pointer to B-Spline Function Entity.
C     IER      Returned error value.  Zero implies no errors.  
C              -1  = Dynamic memory is corrupt or uninitialized
C              -2  = IGE is a null pointer
C              -3  = Garbage value found in IGE
C              -4  = Ambiguous--garbage pointer or deleted entity
C              -5  = IGE points to a deleted entity
C              -6  = IGE does not point to an IGES Entity
C              -7  = IGE is not a B-Spline (type 126 or 128)
C              -8  = Data in Entity is inconsistent with IGES Standard
C              -9  = Insufficient storage for B-Spline Function Entity
C
C
C   CALLS:
C     D0PTR
C     DTERPT
C     DTERR
C
C   HISTORY:
C     9/20/93  D. Parsons  Created from HSRIGS.
C
C*****
C
C   Long Name Alias:
C     ENTRY D2_IGES_TO_BSPLINE (CMEM, IMEM, DMEM, IGE, IBF, IER)

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
C                          IGES Entity
      INTEGER     ENTGE
      PARAMETER  (ENTGE = 238)

C =====================================================================

      INTEGER  IGE, IBF, IER
      INTEGER  ITYPE, IHDR, ITYP, IERX
      INTEGER  IDATA(9)
       
      INTEGER  LN, M, JBLO, JBHI
      INTEGER  INDEX1, INDEX2, IDEG1, IDEG2, K1, K2, N1, N2
      INTEGER  IPROP3, NEED, LENLBL
      CHARACTER*16 LABEL 
      
      CHARACTER*6    SUBNAM
      DATA SUBNAM   /'D2IGBF'/
      
C****

      IER = 0
      IERX = 0
      
C     Check pointer and type

      CALL D0PTR (CMEM, IMEM, DMEM, IGE, ENTGE, 0, IHDR, ITYP, IER)
      IF (IER .NE. 0) GOTO 9900

C     Generate a label for the B-Spline Function Entity 

      CALL D0LABL (CMEM, IMEM, DMEM, IHDR, LABEL, LENLBL)
      
C     Calculate the needed length of the "C-Vector"
      
      CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'D', 1, ITYPE, IERX)
      
      IF (ITYPE .EQ. 126) THEN
         CALL DTERPT(0)
         CALL D2IPFI (CMEM, IMEM, DMEM, 6, IGE, 1, 1, 1, IDATA, IERX)
         CALL DTERPT(1)
         IF (IERX .NE. 0) THEN
            IER = -8
            GOTO 9900
         ENDIF
         
         INDEX1 = IDATA(1)
         IDEG1  = IDATA(2)
         K1     = IDEG1 + 1
         N1     = INDEX1 + 1
         LN     = 1
         IPROP3 = IDATA(5)
         
         IF (IPROP3 .EQ. 0) THEN
            M = 4
         ELSE
            M = 3
         ENDIF

         NEED = 5 + K1 + N1 + M * N1

      ELSE IF (ITYPE .EQ. 128) THEN
         CALL DTERPT(0)
         CALL D2IPFI (CMEM, IMEM, DMEM, 9, IGE, 1, 1, 1, IDATA, IERX)
         CALL DTERPT(1)
         IF (IERX .NE. 0) THEN
            IER = -8
            GOTO 9900
         ENDIF

         INDEX1 = IDATA(1)
         INDEX2 = IDATA(2)
         IDEG1  = IDATA(3)
         IDEG2  = IDATA(4)
         K1     = IDEG1 + 1
         K2     = IDEG2 + 1
         N1     = INDEX1 + 1
         N2     = INDEX2 + 1
         LN     = 2
         IPROP3 = IDATA(7)

         IF (IPROP3 .EQ. 0) THEN
            M = 4
         ELSE
            M = 3
         ENDIF
            
         NEED = 8 + K1 + N1 + K2 + N2 + M * N1 * N2
            
      ELSE
         IER = -7
         GOTO 9900
      ENDIF
      
C     Define the B-Spline Function Entity

      CALL D2BFDF (CMEM, IMEM, DMEM, NEED, LABEL, IBF, IERX)
      IF (IERX .NE. 0) THEN
         IER = -9
         GOTO 9900
      ENDIF
      
C     Fetch the bounds of the "C-Array"

      CALL D2FEBD (CMEM, IMEM, DMEM, IBF, JBLO, JBHI, IERX)
      
C     Convert the IGES Entity to the B-Spline Function Entity

      CALL DTERPT(0)
      CALL D2IGCY (CMEM, IMEM, DMEM, IGE, NEED, DMEM(JBLO), IER)
      CALL DTERPT(1)
      
      IF (IER .NE. 0) GOTO 9900
      
C     Unlock the B-Spline Function Entity

      CALL D2UNLD (CMEM, IMEM, DMEM, IBF, IERX)
      
C     Error reporting section

 9900 CONTINUE

      IF (IER .NE. 0) THEN
          CALL DTERR (1, SUBNAM, IER, 0) 
          CALL DTERPT(0)
          CALL D2ERAS (CMEM, IMEM, DMEM, IBF, IERX)
          CALL DTERPT(1)
      ENDIF

      RETURN
      END
           
