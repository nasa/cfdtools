      SUBROUTINE D0DMPC (CMEM, IMEM, DMEM, LUNIT, IER)

C   Dynamic Memory Pointer Check
C     Tests pointer consistency for the following entity types
C
C        BF    B-Spline Function
C        EG    Edge
C        LP    Loop
C        TS    Trimmed Surface
C        JS    Joined Surface
C        PJ    Point on a Joined Surface
C        CQ    Character Sequence
C        GI    IGES Index
C        GE    IGES
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM     Dynamically managed character data
C     IMEM     Dynamically managed integer data
C     DMEM     Dynamically managed double precision data
C
C   INPUT:
C     LUNIT    Unit to write error messages to
C              If LUNIT = -1, no messages are written
C
C   OUTPUT:
C     IER      Returned error value.  Zero implies no errors.
C              If LUNIT >= 0, IER = -999 for any errors.
C
C     (errors -1 through -71 from D1MBAD)
C
C     -1,3   = Internal length doesn't match given length
C     -4,6   = Internal length not in valid range
C     -7     = Check value at end of IMEM is changed
C     -8     = Error check level flag at beginning of CMEM has invalid value
C     -9     = Check value at end of CMEM is changed
C     -10    = Check value at beginning of DMEM is changed
C     -11    = Check value at end of DMEM is changed
C     -12    = Lowest header (IMEM(16)) is out of range
C     -13    = Lowest header (IMEM(16)) is invalid (MOD 8 test)
C     -14,16 = Lowest free location is out of range
C     -17,19 = Header of highest lock is out of range
C     -20,22 = Header of highest lock is invalid (MOD 8 test)
C     -23,25 = Header of highest lock is marked deleted or zero
C     -26,28 = Lowest unlocked location does not match end of highest lock
C     -29,31 = Lowest search location is out of range
C     -32,34 = Lowest unlocked location is not at default when no locks
C     -35,37 = Lowest search location is not at default when no locks
C     -38    = Top of active list is out of range
C     -39    = Top of active list is invalid (MOD 8 test)
C     -40    = Pointer check value for top of active list is wrong
C     -41    = Top of free list is out of range
C     -42    = Top of free list is invalid (MOD 8 test)
C     -43    = Pointer check value for top of free list is wrong
C     -44    = Last new header sequence number is out of range (0-255)
C     -45    = Pointer check value in header block is wrong
C     -46    = List link value in header block is negative
C     -47,49 = Data space beginning location is out of range
C     -50,52 = Locked data space overlaps free space
C     -53,55 = Locked data space below lowest search location
C     -56,58 = Deleted data space below lowest search location
C     -59    = List link in header out of range
C     -60    = List link in header invalid (MOD 8 test)
C     -61    = Free header block found in active list
C     -62    = Active list loops back on itself
C     -63    = Character data spaces out of order or have gaps
C     -64    = Integer data spaces out of order or have gaps
C     -65    = Double precision data spaces out of order or have gaps
C     -66    = Active list ends before reaching all active header blocks
C     -67    = Active list continues too far
C     -68    = Active header block found in free list
C     -69    = Free list loops back on itself or enters active list
C     -70    = Free list ends before reaching all free header blocks
C     -71    = Free list continues too far
C
C     -99    = Unexpected error from D0PTR
C
C     -999   = LUNIT >= 0 and some error has occurred, see messages
C              for details.
C
C        ENTITY SPECIFIC ERRORS:
C
C     -246xx = B-spline Function
C
C        -24601 = LENC .LT. 0
C        -24602 = LENI .NE. 0
C        -24603 = LEND .LT. 0
C        -24604 = Spline contained in DMEM is invalid
C
C     -243xx = Edge of a Trimmed Surface
C
C        -24301 = LENC .LT. 0
C        -24302 = LENI .NE. 4
C        -24303 = LEND .NE. 0
C        -24304 = IBFC .NE. 0 and does not point to a B-Spline function
C        -24305 = ILP .NE. 0 and does not point to a Loop entity
C        -24306 = ILP .NE. 0 and IXEG .EQ. 0
C        -24307 = ILP .EQ. 0 and IXEG .NE. 0
C        -24308 = IXEG .LT. 0 or IXEG .GT. ILP(MAXEG)
C        -24309 = ILP(IXEG) does not point to this entity
C        -24310 = JEG .NE. 0 and does not point to an Edge entity
C        -24311 = IEG(JEG(JEG)) does not point back to this IEG
C
C     -242xx = Loop of a Trimmed Surface
C
C        -24201 = LENC .LT. 0
C        -24202 = LENI .NE. MAXEG+3
C        -24203 = LEND .NE. 0
C        -24204 = ITS .NE. 0 and does not point to a Trimmed Surface
C        -24205 = ITS .NE. 0 and IXLP .EQ. 0
C        -24206 = ITS .EQ. 0 and IXLP .NE. 0
C        -24207 = IXLP .LT. 0 or IXLP .GT. ITS(MAXLP)
C        -24208 = ITS(IXLP) does not point to this entity
C        -24209 = MAXEG .LT. 1
C        -24210 = At least one non-zero Edge pointer does not point to
C                 an Edge
C        -24211 = At least one pointed to Edge pointer does not point
C                 back to this Loop
C        -24212 = At least one pointed to Edge pointer points back to
C                 this Loop, but has an incorrect IXEG
C
C     -241xx = Trimmed Surface of a Joined Surface
C
C        -24101 = LENC .LT. 0
C        -24102 = LENI .NE. MAXLP+4
C        -24103 = LEND .NE. 0
C        -24104 = IBFS does not point to a surface B-Spline
C        -24105 = IJS .NE. 0 and does not point to a Joined Surface
C        -24106 = IJS .NE. 0 and IXTS .EQ. 0
C        -24107 = IJS .EQ. 0 and IXTS .NE. 0
C        -24108 = IXTS .LT. 0 or IXTS .GT. IJS(MAXTS)
C        -24109 = IJS(IXTS) does not point to this entity
C        -24110 = IORIEN .NE. +1 or -1
C        -24111 = MAXLP .LT. 1
C        -24112 = At least one non-zero Loop pointer does not point to
C                 a Loop
C        -24113 = At least one pointed to Loop pointer does not point
C                 back to this Trimmed surface
C        -24114 = At least one pointed to Loop pointer points back to
C                 this Loop, but has an incorrect IXLP
C
C     -240xx = Joined Surface
C
C        -24001 = LENC .LT. 0
C        -24002 = LENI .NE. MAXTS+2
C        -24003 = LEND .NE. 0
C        -24004 = ICLS .LT. -1 .OR. ICLS .GT. 1
C        -24005 = MAXTS .LT. 1
C        -24006 = At least one non-zero Trimmed Surface pointer does
C                 not point to a Trimmed Surface
C        -24007 = At least one pointed to Trimmed Surface pointer
C                 does not point back to this Joined surface
C        -24008 = At least one pointed to Trimmed Surface pointer points
C                 back to this Joined Surface, but has an incorrect IXTS
C
C     -239xx = Point on a Joined Surface
C
C        -23901 = LENC .LT. 0
C        -23902 = LENI .NE. 2
C        -23903 = LEND .NE. 6
C        -23904 = IEG .NE. 0 and does not point to an Edge entity
C        -23905 = ITS .NE. 0 and does not point to a Trimmed Surface
C        -23906 = IEG->IBFC does not point to a B-spline Function
C                 entity.
C        -23907 = IEG->ILP does not point to a Loop entity
C        -23908 = IEG->ILP->ITS does not point to a Trimmed Surface
C                 entity
C        -23909 = IEG->ILP->ITS->IBFS does not point to a B-Spline
C                 function entity
C        -23910 = ITS is not zero and IEG->ILP->ITS are different.
C        -23911 = Insufficient DMEM for working storage.
C        -23912 = The Curve B-spline function entity has C(3) .LE. 0
C        -23913 = The Curve B-spline function entity has C(4) < C(3)
C        -23914 = The Curve B-spline function entity has invalid knot
C                 set
C        -23915 = The Curve B-spline function entity has denominator
C                 = 0
C        -23916 = T is inside too small an interval
C        -23917 = T out of range.
C        -23918 = The Curve B-spline function entity has NDOM .NE. 1
C        -23919 = The Curve B-spline function entity has NRNG .LT. 1
C        -23920 = The Curve B-spline function entity does not have 2
C                 dependent variables
C        -23921 = The Surface B-spline function entity has K(i)
C                 .LE. 0 for some i
C        -23922 = The Surface B-spline function entity has a number
C                 of B-spline coefficients with respect to the ith
C                 independent variable is less than Ki for some i
C        -23923 = The Surface B-spline function entity has invalid knot
C                 set
C        -23924 = The Surface B-spline function entity has denominator
C                 = 0
C        -23925 = U,V is inside too small an interval
C        -23926 = U or V is out of range.
C        -23927 = The Surface B-spline function entity has NDOM .NE. 2
C        -23928 = The Surface B-spline function entity has NRNG .LT. 1
C        -23929 = The Surface B-spline function entity does not have
C                 3 dependent variables
C        -23930 = IBFC evaluated at (T) does not match given (U,V)
C        -23931 = IBFS evaluated at (U,V) does not match given (X,Y,Z)
C
C     -238xx = IGES Entity
C
C        -23801 = Unable to get parameter lengths
C        -23802 = LENC != 23 (Global) or NPAR+23 (Normal)
C        -23803 = (Global IGES entity) NPAR != 24
C        -23804 = (Global IGES entity) NSTR != 10
C        -23805 = (Global IGES entity) NSTRCH != ??
C        -23806 = (Global IGES entity) NREAL != 4
C        -23807 = (Global IGES entity) LENI != 28
C        -23808 = LENI != NPAR+18
C        -23809 = LEND != NREAL
C        -23810 = IGI not a valid IGES Index
C        -23811 = ICQ not a valid Character Sequence
C        -23812 = (Global IGES entity) TYPCODs incorrect
C        -23813 = Invalid TYPCOD found
C        -23814 = TYPCOD = 'C', but index > LENC
C        -23815 = TYPCOD = 'P' or '*', but IPAR not a valid IGES entity
C        -23816 = TYPCOD = 'L', but IVAL not 0 or 1
C        -23817 = TYPCOD = 'D' or 'E', but index > LEND
C        -23818 = TYPCOD = 'T' or 'R', but IPAR not a valid ICQ entity
C        -23819 = DMEM Free-stack pointer > LEND
C        -23820 = More DMEM Free-stack pointers than LEND
C        -23821 = Number of real TYPCODs > NREAL
C        -23822 = Number of character TYPCODs > NSTR
C
C     -237xx = IGES Index Entity
C
C        -23701 = LENC < 3
C        -23702 = LENI < 2
C        -23703 = LEND != 0
C        -23704 = NSL < 0
C        -23705 = NDE < 0
C        -23706 = LENC != NSL*73+3
C        -23707 = LENI != NDE/2+2
C        -23708 = Global IGE invalid
C        -23709 = IGES IGE invalid
C
C     -249xx = Character Sequence
C
C        -24901 = LENC < 0
C        -24902 = LENI < 1
C        -24903 = LEND != 0
C        -24904 = Non-increasing Ending Positions
C        -24905 = LENC < Last ending position
C
C   HISTORY:
C     16Jun92  D. Parsons     Created.
C     22Jul92  D. Parsons     Added detailed Point checking
C     01Jul93  D. Parsons     Added Character Sequence and IGES entities
C
C*****

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      EXTERNAL D1MBAD
      LOGICAL  D1MBAD

C =====================================================================
C
C     INCLUDE FILE DITYID.INC
C
C     Reserved Entity Type ID Numbers
C
C     HISTORY
C        12May92     D. Parsons     created
C
C                          B-spline Function
      INTEGER     ENTBF
      PARAMETER  (ENTBF = 246)

C                          Edge of a Trimmed Surface
      INTEGER     ENTEG
      PARAMETER  (ENTEG = 243)

C                          Loop of a Trimmed Surface
      INTEGER     ENTLP
      PARAMETER  (ENTLP = 242)

C                          Trimmed Surface of a Joined Surface
      INTEGER     ENTTS
      PARAMETER  (ENTTS = 241)

C                          Joined Surface
      INTEGER     ENTJS
      PARAMETER  (ENTJS = 240)

C                          Point on a Joined Surface
      INTEGER     ENTPJ
      PARAMETER  (ENTPJ = 239)

C                          IGES Entity Data
      INTEGER     ENTGE
      PARAMETER  (ENTGE = 238)

C                          IGES Index
      INTEGER     ENTGI
      PARAMETER  (ENTGI = 237)

C                          Character Sequence Entity
      INTEGER     ENTCQ
      PARAMETER  (ENTCQ = 249)
C =====================================================================

      INTEGER  LUNIT, IER, IERX, IDMY
      INTEGER  ILHDR, ITYP, IADDR, IDE, I1, I2

C     *****

      IER   = 0
      IERX  = 0
      I1   = IDE/256
      I2   = MOD(ABS(IDE),256)

      CALL DTERPT (0)

      IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT,100) I1, I2, IERX
            IER = -999
         ELSE
            IER = IERX
         ENDIF
         GOTO 9999
      ENDIF

C     Loop through entities

      ILHDR = IMEM(16)

      DO 1000, IADDR = IMEM(2)-8, ILHDR, -8

         IDE = IMEM(IADDR)

         CALL D0PTR (CMEM, IMEM, DMEM, IDE, 0, 0, IDMY, ITYP, IERX)
         IF (IERX .NE. 0) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT,110) I1, I2, IERX
               IER = -999
            ELSE
               IER = -99
            ENDIF
            GOTO 9999
         ENDIF

C        Is this one of the recognized entities?

         IF      (ITYP .EQ. ENTBF) THEN
C           B-Spline function
            CALL D0DMBF (CMEM, IMEM, DMEM, IDE, LUNIT, IERX)
            IF (IERX .NE. 0) THEN
               IF (IERX .EQ. -999) THEN
                  IER = -999
               ELSE
                  IER = -100*ENTBF-IERX
               ENDIF
            ENDIF

         ELSE IF (ITYP .EQ. ENTEG) THEN
C           Edge entity
            CALL D0DMEG (CMEM, IMEM, DMEM, IDE, LUNIT, IERX)
            IF (IERX .NE. 0) THEN
               IF (IERX .EQ. -999) THEN
                  IER = -999
               ELSE
                  IER = -100*ENTEG-IERX
               ENDIF
            ENDIF


         ELSE IF (ITYP .EQ. ENTLP) THEN
C           Loop entity
            CALL D0DMLP (CMEM, IMEM, DMEM, IDE, LUNIT, IERX)
            IF (IERX .NE. 0) THEN
               IF (IERX .EQ. -999) THEN
                  IER = -999
               ELSE
                  IER = -100*ENTLP-IERX
               ENDIF
            ENDIF


         ELSE IF (ITYP .EQ. ENTTS) THEN
C           Trimmed Surface entity
            CALL D0DMTS (CMEM, IMEM, DMEM, IDE, LUNIT, IERX)
            IF (IERX .NE. 0) THEN
               IF (IERX .EQ. -999) THEN
                  IER = -999
               ELSE
                  IER = -100*ENTTS-IERX
               ENDIF
            ENDIF


         ELSE IF (ITYP .EQ. ENTJS) THEN
C           Joined Surface entity
            CALL D0DMJS (CMEM, IMEM, DMEM, IDE, LUNIT, IERX)
            IF (IERX .NE. 0) THEN
               IF (IERX .EQ. -999) THEN
                  IER = -999
               ELSE
                  IER = -100*ENTJS-IERX
               ENDIF
            ENDIF


         ELSE IF (ITYP .EQ. ENTPJ) THEN
C           Point on a Joined Surface
            CALL D0DMPJ (CMEM, IMEM, DMEM, IDE, LUNIT, IERX)
            IF (IERX .NE. 0) THEN
               IF (IERX .EQ. -999) THEN
                  IER = -999
               ELSE
                  IER = -100*ENTPJ-IERX
               ENDIF
            ENDIF


         ELSE IF (ITYP .EQ. ENTCQ) THEN
C           Character Sequence
            CALL D0DMCQ (CMEM, IMEM, DMEM, IDE, LUNIT, IERX)
            IF (IERX .NE. 0) THEN
               IF (IERX .EQ. -999) THEN
                  IER = -999
               ELSE
                  IER = -100*ENTCQ-IERX
               ENDIF
            ENDIF


         ELSE IF (ITYP .EQ. ENTGI) THEN
C           IGES Index
            CALL D0DMGI (CMEM, IMEM, DMEM, IDE, LUNIT, IERX)
            IF (IERX .NE. 0) THEN
               IF (IERX .EQ. -999) THEN
                  IER = -999
               ELSE
                  IER = -100*ENTGI-IERX
               ENDIF
            ENDIF


         ELSE IF (ITYP .EQ. ENTGE) THEN
C           IGES
            CALL D0DMGE (CMEM, IMEM, DMEM, IDE, LUNIT, IERX)
            IF (IERX .NE. 0) THEN
               IF (IERX .EQ. -999) THEN
                  IER = -999
               ELSE
                  IER = -100*ENTGE-IERX
               ENDIF
            ENDIF


         ENDIF

C        An error was found, either a message was printed, go on;
C        or return IER if LUNIT < 0

         IF (IER .NE. 0) THEN
            IF (LUNIT .GE. 0) THEN
               IER = -999
            ELSE
               GOTO 9999
            ENDIF
         ENDIF

 1000 CONTINUE

 9999 CONTINUE

      CALL DTERPT(1)

      RETURN

  100 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' ERROR NUMBER',I5,' RETURNED FROM D1MBAD -- ABORTING')
  110 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' ERROR NUMBER',I5,' RETURNED FROM D0PTR  -- ABORTING')

      END

