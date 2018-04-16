      SUBROUTINE D2FETN (CMEM, IMEM, DMEM, ITYP, TYPNAM, NAMLEN, IER)

C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     ITYP    Type id number
C   OUTPUT:
C     TYPNAM  Type name string
C     NAMLEN  Length of used part of TYPNAM
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C             >0  Name was truncated to fit in TYPNAM.  IER value is 
C                 true length.
C             -1  Dynamic memory is corrupt or uninitialized
C             -2  ITYP out of the valid range of 1 to 255.
C             -3  No name defined for this user type number.
C             -4  No name defined for this Library type number.
C
C   CALLS:    D2CQFE, D0FITS, DTERR, DTERPT
C   HISTORY:
C     28Jan92 P Kraushar    Created.
C*****
C
C   Long Name Alias:
C     ENTRY D2_FETCH_TYPE_NAME (CMEM, IMEM, DMEM, ITYP, TYPNAM,
C    +    NAMLEN, IER)

      EXTERNAL D1MBAD
      LOGICAL  D1MBAD

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)

      INTEGER ITYP, NAMLEN, IER
      CHARACTER TYPNAM*(*)

C     LIBTYP is the lowest type number currently used by the Library, minus one.
      INTEGER LIBTYP
      PARAMETER (LIBTYP=236)
C****

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IER)) THEN
            IER = -1
            GOTO 9900
         ENDIF
      ENDIF

      IF (ITYP .LE. 0 .OR. ITYP .GT. 255) THEN
C       Not a legitimate type number
        CALL D0FITS ('** NOT A TYPE', TYPNAM, NAMLEN, IER)
        IER = -2
      ELSE IF (ITYP .GT. LIBTYP) THEN
C       Implement a large case statement on ITYP values for Library types.
C       Begin Case:
        GO TO (237,238,239,240,241,242,243,244,245,246,247,248,249,
     +      250,251,252,253,254,255), ITYP-LIBTYP
  237   CALL D0FITS ('GI IGES Index',
     +      TYPNAM, NAMLEN, IER)
        GO TO 999
  238   CALL D0FITS ('GE IGES Entity Data',
     +      TYPNAM, NAMLEN, IER)
        GO TO 999
  239   CALL D0FITS ('PJ Point on a Joined Surface',
     +      TYPNAM, NAMLEN, IER)
        GO TO 999
  240   CALL D0FITS ('JS Joined Surface',
     +      TYPNAM, NAMLEN, IER)
        GO TO 999
  241   CALL D0FITS ('TS Trimmed Surface of a Joined Surface',
     +      TYPNAM, NAMLEN, IER)
        GO TO 999
  242   CALL D0FITS ('LP Loop of a Trimmed Surface',
     +      TYPNAM, NAMLEN, IER)
        GO TO 999
  243   CALL D0FITS ('EG Edge of a Loop',
     +      TYPNAM, NAMLEN, IER)
        GO TO 999
  244   CALL D0FITS ('CS Curve on a Surface',
     +      TYPNAM, NAMLEN, IER)
        GO TO 999
  245   CALL D0FITS ('BQ B-spline Function Sequence',
     +      TYPNAM, NAMLEN, IER)
        GO TO 999
  246   CALL D0FITS ('BF B-spline Function',
     +      TYPNAM, NAMLEN, IER)
        GO TO 999
  247   CALL D0FITS ('DQ DoubleP Array Sequence',
     +      TYPNAM, NAMLEN, IER)
        GO TO 999
  248   CALL D0FITS ('IQ Integer Array Sequence',
     +      TYPNAM, NAMLEN, IER)
        GO TO 999
  249   CALL D0FITS ('CQ Character String Sequence',
     +      TYPNAM, NAMLEN, IER)
        GO TO 999
  250   CALL D0FITS ('DB DoubleP Array w Bounds',
     +      TYPNAM, NAMLEN, IER)
        GO TO 999
  251   CALL D0FITS ('IB Integer Array w Bounds',
     +      TYPNAM, NAMLEN, IER)
        GO TO 999
  252   CALL D0FITS ('CB Character Array w Bounds',
     +      TYPNAM, NAMLEN, IER)
        GO TO 999
  253   CALL D0FITS ('DA DoubleP Array',
     +      TYPNAM, NAMLEN, IER)
        GO TO 999
  254   CALL D0FITS ('IA Integer Array',
     +      TYPNAM, NAMLEN, IER)
        GO TO 999
  255   CALL D0FITS ('CA Character String',
     +      TYPNAM, NAMLEN, IER)
        GO TO 999
C       end case
  999   CONTINUE

      ELSE IF (ITYP .LE. 192) THEN
C       Should be a user defined type named in the call to D2INIT
        CALL DTERPT (0)
        CALL D2CQFE (CMEM, IMEM, DMEM, IMEM(IMEM(2)-8), ITYP, TYPNAM,
     +      NAMLEN, IER)
        CALL DTERPT (1)
        IF (IER .LT. 0) THEN
C         Either no user types were named or this one was out of range
          CALL D0FITS ('** UNNAMED USER TYPE', TYPNAM, NAMLEN, IER)
          IER = -3
        END IF
      ELSE
C       Currently undefined Library type
        CALL D0FITS ('** UNDEFINED LIB TYPE',TYPNAM, NAMLEN, IER)
        IER = -4
      END IF

 9900 CONTINUE

      IF (IER .LT. 0) THEN
        CALL DTERR (1, 'D2FETN  ', IER, 0)
        NAMLEN = -1
      END IF
      RETURN
      END
