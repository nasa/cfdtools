      SUBROUTINE D2CYIG (CMEM, IMEM, DMEM, C, IGI, JDE, IGE, IER)
      
C   Convert a standard Library Spline Array (C) to an IGES Entity 
C   type 126 (B-Spline Curve) or 128 (B-Spline Surface).
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM     Dynamically managed character data
C     IMEM     Dynamically managed integer data
C     DMEM     Dynamically managed double precision data
C
C   INPUT:
C     C        Spline array.
C     IGI      Pointer to IGES Index Entity.
C     JDE      DE number for this new IGES Entity.
C
C   OUTPUT:
C     IGE      Pointer to the IGES Spline Entity (type 126 or 128)
C     IER      Returned error value.  Zero implies no errors.  
C              -1  = Dynamic memory is corrupt or uninitialized
C              -2  = IGI is a null pointer
C              -3  = Garbage value found in IGI
C              -4  = Ambiguous--garbage pointer or deleted entity
C              -5  = IGI points to a deleted entity
C              -6  = IGI does not point to an IGES Index entity.
C              -7  = The JDEth directory entry is not empty.
C              -8  = JDE < 1 or JDE > NDE.
C              -9  = Number of independent variables < 1 or > 2
C              -10 = Number of dependent variables < 1 or > 3
C              -11 = Rational Spline; Encountered a Weight <= 0.0
C              -999= Error in a lower-level subroutine
C
C   CALLS:
C     D0PTR
C     DTERPT
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
C     ENTRY D2_CARRAY_TO_IGES (CMEM, IMEM, DMEM, C, IGE, IER)

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

C =====================================================================

      INTEGER  IGI, JDE, IGE, IER
      DOUBLE PRECISION C(*)
      INTEGER  IERX, IDATA(9), ITYPE, IGEHDR
      DOUBLE PRECISION NORM(3)
      LOGICAL  RATNL
      INTEGER  IA, IB, IC, N1, N2, NPAR, NREAL, NSTR, NSTRCH, NUM
      INTEGER  IHDR, ITYP, IPOLD, I, IP
      INTEGER  INDEP, DEP, IPWT, K1, K2, IPKT1, IPKT2, NKTS1, NKTS2
      INTEGER  DEG1, DEG2, IPPT, NPTS
      DOUBLE PRECISION DVAL, WEIGHT
       
      CHARACTER*6    SUBNAM
      DATA SUBNAM   /'D2CYIG'/
      
C****

      IER = 0
      IERX = 0
      
C     Check pointer and type, then locate relevant data spaces
      CALL D0PTR (CMEM, IMEM, DMEM, IGI, 0, 0, IHDR, ITYP, IER)
      IF (IER .LT. 0) GOTO 9900

C     Check that the JDEth pointer is empty 
      CALL DTERPT(0)
      CALL D2IGFI (CMEM, IMEM, DMEM, IGI, 'I', JDE, IPOLD, IERX)
      CALL DTERPT(1)
      IF (IERX .NE. 0) THEN
         IF (IERX .EQ. -8) THEN 
            IER = -8
         ELSE
            CALL DTERR (1, 'D2IGFI', IERX, 0)
            IER = -999
         ENDIF
         GOTO 9900
      ENDIF
      
      IF (IPOLD .NE. 0) THEN
         IER = -7
         GOTO 9900
      ENDIF 
         
      INDEP = C(1)
      
      IF ((INDEP.LT.1) .OR. (INDEP.GT.2)) THEN
         IER = -9
         GOTO 9900
      ENDIF
      
      IF(C(2) .GT. 0.0) THEN
         DEP=C(2)
         RATNL = .FALSE.
      ELSE
         DEP=ABS(C(2))-1
         RATNL = .TRUE.
      ENDIF
      
      IF ((DEP.LT.1) .OR. (DEP.GT.3)) THEN
         IER = -10
         GOTO 9900
      ENDIF

      IF (INDEP.EQ.1) THEN
         ITYPE = 126
         DEG1 = C(3) - 1
         K1 = C(4) - 1
C
         IDATA(1) = K1
         IDATA(2) = DEG1               
C        Non-planar versus planar
         IF (DEP .EQ. 3) THEN
            IDATA(3) = 0
         ELSE
            IDATA(3) = 1
         ENDIF
C        Closure?
         IDATA(4) = 0
C        Rational versus Polynomial
         IF (RATNL) THEN
            IDATA(5) = 0
         ELSE
            IDATA(5) = 1
         ENDIF
C        Periodicity?
         IDATA(6) = 0
         
C        Starting position for Knots in C vector:
         IPKT1=6

C        Number of Knots:
         NKTS1=(K1+1)+(DEG1+1)

C        Starting position for Control Points in C vector:
         IPPT = IPKT1 + NKTS1

C        Number of Control points for each dependent variable
         NPTS = K1+1

C        Starting position of weights in C vector
         IPWT = IPPT+DEP*NPTS
         
C        Determine the number of IGES parameters needed
C        
         N1 = 1 + K1-DEG1
         IA = N1 + 2*DEG1
         
         NPAR  = 16 + IA + 4*K1
         NREAL = NPAR-6
         NSTRCH = 0

C        Create the IGES Entity

         CALL D2IGDE (CMEM, IMEM, DMEM, NPAR, NSTR, NSTRCH, NREAL,
     +                   IGI, JDE, IGE, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9900
         ENDIF
         
C        Store the Entity type

         CALL D2IGSI (CMEM, IMEM, DMEM, 'I', ITYPE, IGE, 'D', 1, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9900
         ENDIF

C        Store the Integer data

         CALL D2IPSI (CMEM, IMEM, DMEM, 6, IDATA, 1, IGE, 'I',
     +      1, 1, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9900
         ENDIF

C        Copy the knots

         NUM = NKTS1
         IP = 7
         CALL D2IPSD (CMEM, IMEM, DMEM, NUM, C(IPKT1), 1, IGE, 'D',
     +      IP, 1, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9900
         ENDIF
         
         IP = IP+NUM
     
C        If rational, copy the weights. Otherwise, store 1.0D0's

         NUM = NPTS
         IF (RATNL) THEN
            CALL D2IPSD (CMEM, IMEM, DMEM, NUM, C(IPWT), 1, 
     +         IGE, 'D', IP, 1, IERX)
         ELSE
            CALL D2IPSD (CMEM, IMEM, DMEM, NUM, 1.0D0, 0, 
     +         IGE, 'D', IP, 1, IERX)
         ENDIF
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9900
         ENDIF
         
         IP = IP+NUM

C        Calculate and store the Control Points

         NUM = NPTS
         
C        X's
         IF (RATNL) THEN         
            DO 100 I = 0, NPTS-1
               WEIGHT = C(IPWT+I)
               IF (WEIGHT .LE. 0.0D0) THEN
                  IER = -11
                  GOTO 9900
               ENDIF
               DVAL = C(IPPT+I) / WEIGHT
               CALL D2IGSD (CMEM, IMEM, DMEM, 'D', DVAL, IGE, 'P',
     +            IP+3*I, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -999
                  GOTO 9900
               ENDIF
  100       CONTINUE
         ELSE
            CALL D2IPSD (CMEM, IMEM, DMEM, NUM, C(IPPT), 1, IGE, 'D',
     +         IP, 3, IERX)
            IF (IERX .NE. 0) THEN
               IER = -999
               GOTO 9900
            ENDIF
         ENDIF
         
C        Y's         
         IF (DEP .GT. 1) THEN
            IF (RATNL) THEN         
               DO 200 I = 0, NPTS-1
                  WEIGHT = C(IPWT+I)
                  IF (WEIGHT .LE. 0.0D0) THEN
                     IER = -11
                     GOTO 9900
                  ENDIF
                  DVAL = C(IPPT+NPTS+I) / WEIGHT
                  CALL D2IGSD (CMEM, IMEM, DMEM, 'D', DVAL, IGE, 'P',
     +               IP+1+3*I, IERX)
                  IF (IERX .NE. 0) THEN
                     IER = -999
                     GOTO 9900
                  ENDIF
  200          CONTINUE
            ELSE
               CALL D2IPSD (CMEM, IMEM, DMEM, NUM, C(IPPT+NPTS), 1, 
     +            IGE, 'D', IP+1, 3, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -999
                  GOTO 9900
               ENDIF
            ENDIF
         ELSE
            CALL D2IPSD (CMEM, IMEM, DMEM, NUM, 0.0D0, 0, IGE, 'D', 
     +         IP+1, 3, IERX)
            IF (IERX .NE. 0) THEN
               IER = -999
               GOTO 9900
            ENDIF
         ENDIF

C        Z's         
         IF (DEP .GT. 2) THEN
            IF (RATNL) THEN         
               DO 300 I = 0, NPTS-1
                  WEIGHT = C(IPWT+I)
                  IF (WEIGHT .LE. 0.0D0) THEN
                     IER = -11
                     GOTO 9900
                  ENDIF
                  DVAL = C(IPPT+2*NPTS+I) / WEIGHT
                  CALL D2IGSD (CMEM, IMEM, DMEM, 'D', DVAL, IGE, 'P',
     +               IP+2+3*I, IERX)
                  IF (IERX .NE. 0) THEN
                     IER = -999
                     GOTO 9900
                  ENDIF
  300          CONTINUE
            ELSE
               CALL D2IPSD (CMEM, IMEM, DMEM, NUM, C(IPPT+2*NPTS), 1, 
     +            IGE, 'D', IP+2, 3, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -999
                  GOTO 9900
               ENDIF
            ENDIF
         ELSE
            CALL D2IPSD (CMEM, IMEM, DMEM, NUM, 0.0D0, 0, IGE, 'D', 
     +         IP+2, 3, IERX)
            IF (IERX .NE. 0) THEN
               IER = -999
               GOTO 9900
            ENDIF
         ENDIF
         
         IP = IP + 3*NPTS

C        V(0) and V(1) are set to the end knots 

         CALL D2IGSD (CMEM, IMEM, DMEM, 'D', C(6), IGE, 'P',
     +         IP, IERX)               
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9900
         ENDIF

         CALL D2IGSD (CMEM, IMEM, DMEM, 'D', C(5+NKTS1), IGE, 'P',
     +         IP+1, IERX)               
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9900
         ENDIF

         IP = IP+2
         
C        Set the Normal vector to (0,0,1) if planar, (0,0,0) otherwise

         NORM(1) = 0.0D0
         NORM(2) = 0.0D0
         IF (DEP .EQ. 3) THEN
            NORM(3) = 0.0D0
         ELSE
            NORM(3) = 1.0D0
         ENDIF
         
         NUM = 3
         CALL D2IPSD (CMEM, IMEM, DMEM, NUM, NORM, 1, 
     +      IGE, 'D', IP, 1, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9900
         ENDIF
      ELSE
         ITYPE = 128
         DEG1 = C(3) - 1
         DEG2 = C(4) - 1
         K1 = C(5) - 1
         K2 = C(6) - 1
C
         IDATA(1) = K1
         IDATA(2) = K2
         IDATA(3) = DEG1
         IDATA(4) = DEG2               
         IDATA(5) = 0
         IDATA(6) = 0
         
         IF (RATNL) THEN
            IDATA(7) = 0
         ELSE
            IDATA(7) = 1
         ENDIF
         
         IDATA(8) = 0
         IDATA(9) = 0
         
C        Number of Knots:
         NKTS1=(K1+1)+(DEG1+1)
         NKTS2=(K2+1)+(DEG2+1)

C        Starting position for Knots in C vector:
         IPKT1=9
         IPKT2=9+NKTS1 

C        Starting position for Control Points in C vector:
         IPPT = IPKT2 + NKTS2

C        Number of Control points for each dependent variable
         NPTS = (K1+1)*(K2+1)

C        Starting position of weights in C vector
         IPWT = IPPT+DEP*NPTS
         
C        Determine the number of IGES parameters needed
C        
         N1 = 1 + K1 - DEG1
         N2 = 1 + K2 - DEG2
         IA = N1 + 2*DEG1
         IB = N2 + 2*DEG2
         IC = NPTS
         
         NPAR  = 15 + IA + IB + 4*IC
         NREAL = NPAR-9
         NSTRCH = 0

C        Create the IGES Entity

         CALL D2IGDE (CMEM, IMEM, DMEM, NPAR, NSTR, NSTRCH, NREAL,
     +                   IGI, JDE, IGE, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9900
         ENDIF
         
C        Store the Entity type

         CALL D2IGSI (CMEM, IMEM, DMEM, 'I', ITYPE, IGE, 'D', 1, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9900
         ENDIF

C        Store the Integer data

         CALL D2IPSI (CMEM, IMEM, DMEM, 9, IDATA, 1, IGE, 'I',
     +      1, 1, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9900
         ENDIF

C        Copy the knots

         NUM = NKTS1 + NKTS2
         IP = 10
         CALL D2IPSD (CMEM, IMEM, DMEM, NUM, C(IPKT1), 1, IGE, 'D',
     +      IP, 1, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9900
         ENDIF
         
         IP = IP+NUM
     
C        If rational, copy the weights. Otherwise, store 1.0D0's

         NUM = NPTS
         IF (RATNL) THEN
            CALL D2IPSD (CMEM, IMEM, DMEM, NUM, C(IPWT), 1, 
     +         IGE, 'D', IP, 1, IERX)
         ELSE
            CALL D2IPSD (CMEM, IMEM, DMEM, NUM, 1.0D0, 0, 
     +         IGE, 'D', IP, 1, IERX)
         ENDIF
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9900
         ENDIF
         
         IP = IP+NUM

C        Calculate and store the Control Points

         NUM = NPTS
         
C        X's
         IF (RATNL) THEN         
            DO 400 I = 0, NPTS-1
               WEIGHT = C(IPWT+I)
               IF (WEIGHT .LE. 0.0D0) THEN
                  IER = -11
                  GOTO 9900
               ENDIF
               DVAL = C(IPPT+I) / WEIGHT
               CALL D2IGSD (CMEM, IMEM, DMEM, 'D', DVAL, IGE, 'P',
     +            IP+3*I, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -999
                  GOTO 9900
               ENDIF
  400       CONTINUE
         ELSE
            CALL D2IPSD (CMEM, IMEM, DMEM, NUM, C(IPPT), 1, IGE, 'D',
     +         IP, 3, IERX)
            IF (IERX .NE. 0) THEN
               IER = -999
               GOTO 9900
            ENDIF
         ENDIF
         
C        Y's         
         IF (DEP .GT. 1) THEN
            IF (RATNL) THEN         
               DO 500 I = 0, NPTS-1
                  WEIGHT = C(IPWT+I)
                  IF (WEIGHT .LE. 0.0D0) THEN
                     IER = -11
                     GOTO 9900
                  ENDIF
                  DVAL = C(IPPT+NPTS+I) / WEIGHT
                  CALL D2IGSD (CMEM, IMEM, DMEM, 'D', DVAL, IGE, 'P',
     +               IP+1+3*I, IERX)
                  IF (IERX .NE. 0) THEN
                     IER = -999
                     GOTO 9900
                  ENDIF
  500          CONTINUE
            ELSE
               CALL D2IPSD (CMEM, IMEM, DMEM, NUM, C(IPPT+NPTS), 1, 
     +            IGE, 'D', IP+1, 3, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -999
                  GOTO 9900
               ENDIF
            ENDIF
         ELSE
            CALL D2IPSD (CMEM, IMEM, DMEM, NUM, 0.0D0, 0, IGE, 'D', 
     +         IP+1, 3, IERX)
            IF (IERX .NE. 0) THEN
               IER = -999
               GOTO 9900
            ENDIF
         ENDIF

C        Z's         
         IF (DEP .GT. 2) THEN
            IF (RATNL) THEN         
               DO 600 I = 0, NPTS-1
                  WEIGHT = C(IPWT+I)
                  IF (WEIGHT .LE. 0.0D0) THEN
                     IER = -11
                     GOTO 9900
                  ENDIF
                  DVAL = C(IPPT+2*NPTS+I) / WEIGHT
                  CALL D2IGSD (CMEM, IMEM, DMEM, 'D', DVAL, IGE, 'P',
     +               IP+2+3*I, IERX)
                  IF (IERX .NE. 0) THEN
                     IER = -999
                     GOTO 9900
                  ENDIF
  600          CONTINUE
            ELSE
               CALL D2IPSD (CMEM, IMEM, DMEM, NUM, C(IPPT+2*NPTS), 1, 
     +            IGE, 'D', IP+2, 3, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -999
                  GOTO 9900
               ENDIF
            ENDIF
         ELSE
            CALL D2IPSD (CMEM, IMEM, DMEM, NUM, 0.0D0, 0, IGE, 'D', 
     +         IP+2, 3, IERX)
            IF (IERX .NE. 0) THEN
               IER = -999
               GOTO 9900
            ENDIF
         ENDIF
         
         IP = IP + 3*NPTS

C        U(0), U(1), V(0) and V(1) are set to the end knots 

         CALL D2IGSD (CMEM, IMEM, DMEM, 'D', C(9), IGE, 'P',
     +         IP, IERX)               
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9900
         ENDIF

         CALL D2IGSD (CMEM, IMEM, DMEM, 'D', C(8+NKTS1), IGE, 'P',
     +         IP+1, IERX)               
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9900
         ENDIF

         CALL D2IGSD (CMEM, IMEM, DMEM, 'D', C(9+NKTS1), IGE, 'P',
     +         IP+2, IERX)               
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9900
         ENDIF

         CALL D2IGSD (CMEM, IMEM, DMEM, 'D', C(8+NKTS1+NKTS2), IGE, 'P',
     +         IP+3, IERX)               
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9900
         ENDIF
      ENDIF

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
    
