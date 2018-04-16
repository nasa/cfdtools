      SUBROUTINE D2IGEG (CMEM, IMEM, DMEM, IGE, ILP, IXEG, IEG, IER)

C     Translate an IGES type 126 B-Spline Curve into an Edge entity.  This
C     is only valid if the curve is planar with all z-coordinates zero.
C     Optionally, connect the Edge to a Loop.  Note that any joined
C     edges must be connected separately using D2CNEE.
C
C   DYNAMIC MEMORY INPUT-OUTPUT
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IGE     MEM pointer to an IGES Entity entity
C     ILP     MEM pointer to a Loop to which this Edge belongs, if any.
C             Otherwise, ILP must be zero.
C     IXEG    Index of this Edge in the list of Edge entities in ILP.
C             Ignored if ILP is zero.
C   OUTPUT:
C     IEG     MEM pointer to new Edge entity
C     IER     Returned error code.  Zero implies no errors.  If negative,
C             the IEG pointer is set null (0).
C             +1  = Translation successful, but requested connection to
C                   joined surface was unsuccessful
C             -1  = Dynamic memory is corrupt or uninitialized
C             -2  = IGE is a null pointer
C             -3  = Garbage value found in IGE
C             -4  = Ambiguous - garbage pointer or deleted entity
C             -5  = IGE points to a deleted entity
C             -6  = IGE does not point to an IGES Entity entity
C             -7  = IGE is not IGES Entity type 126 (B-spline curve)
C             -8  = IGE is not planar (i.e. it does not map into a surface's
C                   parameter space)
C             -9  = IGE does not have all z-coordinates zero.
C             -10 = Unsuccessful at allocating memory for Edge entity
C             -11 = Curve declared polynomial but weights not all equal
C             -999= Unexpected error from lower-level routine
C
C   CALLS:
C     D0LABL
C     D0PTR
C     D2BFDF
C     D2CNEL
C     D2EGER
C     D2ERAS
C     D2IGBF
C     D2IGFI
C     D2IPFD
C     D2UNLD
C     DTERR
C
C   HISTORY:
C     24Aug93 P. Kraushar   Created
C*****
C
C   Long Name Alias:
C     ENTRY D2_IGES_TO_EDGE (CMEM, IMEM, DMEM, IGE, ILP, IXEG,
C    +                       IEG, IER)
C
      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER IGE, ILP, IXEG, IEG, IER

      INTEGER ENTGE
      PARAMETER (ENTGE=238)

      CHARACTER LABEL*16
      INTEGER IHDR, ITYP, IERX, ICBF, I, KORD, NC, JP, JX, JY, JW, J0
      INTEGER IPOLY, LEN
      DOUBLE PRECISION W
      CHARACTER SUBNAM*6

      DATA SUBNAM /'D2IGEG'/

C****
      IER = 0
      IEG = 0
      ICBF = 0
C     Check pointer and type of IGE
      CALL D0PTR (CMEM, IMEM, DMEM, IGE, ENTGE, 0, IHDR, ITYP, IER)
      IF (IER .LT. 0) GOTO 9000
      IER = -999
      CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'D', 1, ITYP, IERX)
      IF (IERX .NE. 0) GOTO 9000
      IF (ITYP .NE. 126) THEN
         IER = -7
         GOTO 9000
      ENDIF
      CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'P', 3, ITYP, IERX)
      IF (IERX .NE. 0) GOTO 9000
      IF (ITYP .NE. 1) THEN
         IER = -8
         GOTO 9000
      ENDIF

C     Get max index and degree and convert to ncoef and order
      CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'P', 1, NC, IERX)
      IF (IERX .NE. 0) GOTO 9000
      CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'P', 2, KORD, IERX)
      IF (IERX .NE. 0) GOTO 9000
      NC = NC + 1
      KORD = KORD + 1

C     Is it polynomial or rational?
      CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'P', 5, IPOLY, IERX)
      IF (IERX .NE. 0) GOTO 9000

C     Allocate B-spline function
      LEN = 5 + (NC + KORD) + (3-IPOLY)*NC
      CALL D2BFDF (CMEM, IMEM, DMEM, LEN, ' ', ICBF, IERX)
      IF (IERX .NE. 0) THEN
         IERX = -10
         GOTO 9000
      ENDIF

C     Verify all z coordinates are zero
      CALL D2FEBD (CMEM, IMEM, DMEM, ICBF, J0, JW, IERX)
C     Note that the DMEM space of the B-spline is now explicitly locked.
C     The lock permits us to work with DMEM  directly.
      IF (IERX .NE. 0) GOTO 9000
      JP = 9 + NC + KORD + NC
      CALL D2IPFD (CMEM, IMEM, DMEM, NC, IGE, JP, 3, 1, DMEM(J0), IERX)
      IF (IERX .NE. 0) GOTO 9000
      DO 100 I=J0,J0+NC-1
         IF (DMEM(I) .NE. 0.0D0) THEN
            IER = -9
            GOTO 9000
         ENDIF
  100 CONTINUE

C     Start constructing B-spline
      DMEM(J0) = 1.0D0
      IF (IPOLY .EQ. 1) THEN
         DMEM(J0+1) = 2.0D0
      ELSE
         DMEM(J0+1) = -3.0D0
      ENDIF
      DMEM(J0+2) = KORD
      DMEM(J0+3) = NC
      DMEM(J0+4) = KORD
C     Fetch knots
      CALL D2IPFD (CMEM, IMEM, DMEM, NC+KORD, IGE, 7, 1, 1, DMEM(J0+5),
     +             IERX)
      IF (IERX .NE. 0) GOTO 9000
      JX = J0 + 5 + NC + KORD
      JY = JX + NC
      JW = JY + NC
      JP = 7 + NC + KORD
      IF (IPOLY .EQ. 1) THEN
C        If polynomial, verify all weights are equal
         CALL D2IPFD (CMEM, IMEM, DMEM, NC, IGE, JP, 1, 1, DMEM(JX), 
     +                IERX)
         IF (IERX .NE. 0) GOTO 9000
         W = DMEM(JX)
         DO 200 I=JX+1,JX+NC-1
            IF (DMEM(I) .NE. W) THEN
               IER = -11
               GOTO 9000
            ENDIF
  200    CONTINUE
      ELSE
C        Fetch weights for rational spline
         CALL D2IPFD (CMEM, IMEM, DMEM, NC, IGE, JP, 1, 1, DMEM(JW),
     +                IERX)
         IF (IERX .NE. 0) GOTO 9000
      ENDIF
C     Fetch x coordinates    
      JP = 7 + NC + KORD + NC
      CALL D2IPFD (CMEM, IMEM, DMEM, NC, IGE, JP, 3, 1, DMEM(JX), IERX)
      IF (IERX .NE. 0) GOTO 9000
C     Fetch y coordinates
      JP = 8 + NC + KORD + NC
      CALL D2IPFD (CMEM, IMEM, DMEM, NC, IGE, JP, 3, 1, DMEM(JY), IERX)
      IF (IERX .NE. 0) GOTO 9000

C     Correct for Library  version of coeficients
      IF (IPOLY .EQ. 0) THEN
         DO 300 I=0,NC-1
            DMEM(JX+I) = DMEM(JX+I) * DMEM(JW+I)
            DMEM(JY+I) = DMEM(JY+I) * DMEM(JX+I)
  300    CONTINUE
      ELSE IF (W .NE. 1.0D0) THEN
         DO 400 I=0,NC-1
            DMEM(JX+I) = DMEM(JX+I) * W
            DMEM(JY+I) = DMEM(JY+I) * W
  400    CONTINUE
      ENDIF

C     Do not fail to unlock the DMEM space of the B-spline
      CALL D2UNLD (CMEM, IMEM, DMEM, ICBF, IERX)
      IF (IERX .NE. 0) GOTO 9000

C     Fetch the label and allocate the Edge entity
      CALL D0LABL (CMEM, IMEM, DMEM, IHDR, LABEL, LEN)
      CALL D2EGDF (CMEM, IMEM, DMEM, ICBF, 0, 0, 0, LABEL, IEG, IERX)
      IF (IERX .NE. 0) THEN
         IERX = -10
         GOTO 9000
      ENDIF
      IER = 0

C     Finally, connect to Loop, if selected.
      IF (ILP .NE. 0) THEN
         CALL D2CNEL (CMEM, IMEM, DMEM, IXEG, IEG, ILP, IERX)
         IF (IERX .NE. 0) THEN
            IER = 1
            CALL DTERR (0, SUBNAM, IER, 0)
         ENDIF
      ENDIF
      RETURN

C     Error Handling

 9000 CONTINUE
C     Report error
      IF (IER .LE. -999) THEN
         CALL DTERR (4, SUBNAM, IER, 0)
      ELSE
         CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF
C     Attempt to erase any partially generated entity
      IF (ICBF .NE. 0) THEN
         CALL DTERPT (0)
         CALL D2ERAS (CMEM, IMEM, DMEM, ICBF, IERX)
         CALL DTERPT (1)
      ENDIF
      IF (IEG .NE. 0) THEN
         CALL DTERPT (0)
         CALL D2EGER (CMEM, IMEM, DMEM, IEG, IERX)
         CALL DTERPT (1)
         IEG = 0
      ENDIF
C     End of error handling
      RETURN
      END
