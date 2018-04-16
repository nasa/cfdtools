C+----------------------------------------------------------------------
C
      FUNCTION NUMBER( STRING )
C
C     Acronym:
C
C        if ( NUMBER( string ) ) then
C           string very likely represents a number - integer or real
C
C     Description and usage:
C
C        A simple(-minded) test for numeric data is implemented by
C        searching the input string for legitimate characters:
C                digits 0 to 9, D, E, -, + and .
C        Insurance is provided by requiring that a numeric string
C        have at least one digit, at most one D, E or .
C        and at most two -s or +s.  Note that a few ambiguities remain:
C
C           (a)  A string might have the form of numeric data but be
C                intended as text.  No general test can hope to detect
C                such cases.
C
C           (b)  There is no check for correctness of the data format.
C                For example a meaningless string such as 'E1.+2-'
C                will be accepted as numeric.
C
C        Despite these weaknesses, the method should work in the
C        majority of cases.
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        NUMBER      -       L      O    Set .TRUE. if STRING appears
C                                        to be numerical data, else
C                                        set .FALSE.
C        STRING      *       C    I      Input data to be tested,
C                                        assumed to be in upper case.
C
C     Notes:
C
C        (1)  It is assumed that STRING has been extracted by
C             a "token" utility - hence the upper case assumption.
C
C        (2)  The scan of STRING stops at the first blank.
C
C        (3)  COMPLEX data with parentheses will not look numeric.
C
C     Environment:  ANSI FORTRAN 77.
C
C     Michael Saunders, Systems Optimization Lab., Stanford University.
C     12 Nov  1985    Initial design and coding, starting from the
C                     routine ALPHA from Informatics General, Inc.
C     23 May  1986    OPNUMB name changed to NUMBER, analogous to ALPHA,
C                     for use at NASA Ames - D. Saunders, Informatics.
C
C-----------------------------------------------------------------------

C  Arguments:

      LOGICAL          NUMBER
      CHARACTER*(*)    STRING

C  Local variables:

      LOGICAL         NUM
      INTEGER         J, LENGTH, NDIGIT, NEXP, NMINUS, NPLUS, NPOINT
      CHARACTER*1     ATOM

C  Executable statements:

      NDIGIT = 0
      NEXP   = 0
      NMINUS = 0
      NPLUS  = 0
      NPOINT = 0
      NUM    = .TRUE.
      LENGTH = LEN (STRING)
      J      = 0

   10    J    = J + 1
         ATOM = STRING (J:J)
         IF      (ATOM .GE. '0'  .AND.  ATOM .LE. '9') THEN
            NDIGIT = NDIGIT + 1
         ELSE IF (ATOM .EQ. 'D'  .OR.   ATOM .EQ. 'E') THEN
            NEXP   = NEXP   + 1
         ELSE IF (ATOM .EQ. '-') THEN
            NMINUS = NMINUS + 1
         ELSE IF (ATOM .EQ. '+') THEN
            NPLUS  = NPLUS  + 1
         ELSE IF (ATOM .EQ. '.') THEN
            NPOINT = NPOINT + 1
         ELSE IF (ATOM .EQ. ' ') THEN
            J      = LENGTH
         ELSE
            NUM    = .FALSE.
         END IF

         IF (NUM  .AND.  J .LT. LENGTH) GO TO 10

      NUMBER = NUM
     $         .AND.  NDIGIT .GE. 1
     $         .AND.  NEXP   .LE. 1
     $         .AND.  NMINUS .LE. 2
     $         .AND.  NPLUS  .LE. 2
     $         .AND.  NPOINT .LE. 1

      RETURN
      END
