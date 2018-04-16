C+------------------------------------------------------------------------------
C
      SUBROUTINE GRIDSPEC (CASE, LUNDAT, LUNERR, IFIRST, ILAST, NPTS,
     >                     METHOD, CPARAM, IPARAM, RPARAM)
C
C  PURPOSE:
C        GRIDSPEC is one solution to the problem of specifying multiple
C     1-dimensional point distributions noninteractively.  It reads the
C     parameters for one distribution at a time from a keyword-style control
C     file which is assumed to be open.
C
C        GRIDSPEC was introduced with the ARCDIS/DISTRIB combination in mind,
C     and prompted by a double-delta wing surface grid generation application,
C     where the wing edges have indefinite numbers of discontinuities requiring
C     distributions within various subintervals.
C
C  METHOD:
C        A sample control file will help:
C
C     Planform edge distribution controls for demo. wing       ! Title (ignored)
C     Edge=Leading                                             ! Case
C     FIRST=1  LAST=2  NPTS=30  METHOD=12  C1=L R1=-.1  R2=-1. ! One spec./line
C     FIRST=2  LAST=3  NPTS=20  METHOD=12  C1=L R1=-1.  R2=-.1
C     FIRST=3  LAST=83 NPTS=20  METHOD=12  C1=B R1=-.1  R2=-.1 ! B=Bessel(loose)
C     Edge=Root                                                ! Another case
C     FIRST=1  LAST=21 NPTS=30  METHOD=12  C1=B R1=-.1  R2=-1.
C     FIRST=21 LAST=36 NPTS=20  METHOD=12  C1=B R1=-1.  R2=-.1
C     FIRST=36 LAST=51 NPTS=20  METHOD=12  C1=B R1=-.1  R2=-.1   
C     Edge=Trailing
C     FIRST=1  LAST=2  NPTS=50  METHOD=3  C1=L  I1=2  R1=.5  ! Modif. sine dist.
C
C     Explanation:
C
C        The first line of parameters specifies a 30-point distribution using
C     method 12 (HTDIS2) on the arc defined by points 1 and 2 (of some external
C     set of discrete points), with the initial and final increments (-ve=%)
C     .1% and 1% of the interval length.  C1=L controls the type of parametric
C     fit (piecewise-linear) to be used by ARCDIS's calls to PLSFIT3D.  Any
C     C2, C3, ... would apply to the distribution defined by METHOD, as do
C     R1, R2, ... and I1, I2, ...
C
C        Keyword/value strings are assumed to be no longer than 12 characters.
C     Delimiters can be any of those supported by the TOKENS utility: blank,
C     tab, comma, ':', or '='.  All specs. for a distribution must be on the
C     same line (<= 80 characters).
C
C  USAGE:
C
C        Trying to return all of the sub-distributions for one case was
C     discarded in favor of returning one specification at a time (because
C     of the indefinite numbers of integer, real, and string parameters
C     associated with any one distribution).  Therefore, the application
C     must call GRIDSPEC in a loop over an indefinite number of distributions
C     for each case.
C
C        If NPTS = 0 on input (as it should be on the first call for a given
C     case), the control file is rewound and searched for CASE from the second
C     line on (title ignored).  Then (or if NPTS > 0 as it normally will be
C     after the first call for a case), one line is parsed for distribution
C     parameters.  NPTS = 0 on return is used to indicate EOF or an unidentified
C     distribution-spec. keyword (presumably the start of the next case).
C     NPTS = 0 on return also if an input case is not found.  This may be
C     OK in the application.
C
C  ARGUMENTS:
C     ARG     DIM    TYPE   I/O/S   DESCRIPTION
C     CASE     -       C    I       Keyword/value string defining requested case
C                                   in uppercase.  E.g. 'EDGE=LEADING'.  Both
C                                   keyword and value are checked for a match.
C                                   No other keywords should appear on the line.
C                                   The keyword should not clash with the well-
C                                   defined possibilities common to all specs.
C     LUNDAT,  -       I    I       Logical units for the grid spec. file and
C     LUNERR                        for diagnostics (assumed open already).
C     IFIRST,  -       I      O     First and last indices in (external) X,Y,Z
C     ILAST                         arrays defining subarc to be redistributed.
C     NPTS     -       I    I/O     Number of pts. to distribute on the subarc.
C                                   See USAGE above for meaning of NPTS=0 on
C                                   input and output.
C     METHOD   -       I      O     Method number - see summary in DISTRIB.
C     CPARAM   *      C*(*)   O     CHARACTER parameter(s) packed in order.
C                                   All methods require at least 1 for ARCDIS's
C                                   PMETHOD argument.  No others are foreseen
C                                   but if DISTRIB ever needs any, they will be
C                                   appended to this string.
C     IPARAM   *       I      O     INTEGER parameter(s) used by some methods,
C                                   packed in the order described in DISTRIB.
C     RPARAM   *       R      O     REAL parameters used by certain methods.
C                                   Usually 1 to 3 values - see summary in
C                                   DISTRIB.  Any NEGATIVE value intended as an
C                                   increment indicates a PERCENTAGE of the
C                                   interval rather than an absolute increment.
C                                   DISTRIB does any conversion from relative
C                                   to absolute in-place.
C
C  EXTERNAL REFERENCES:
C     GETLINE  Utility for reading text files
C     LOOKUP   Keyword dictionary utility
C     PAIRS    Keyword/value parsing utility
C
C  FILES USED:
C     See LUNDAT and LUNERR arguments.
C
C  ENVIRONMENT:
C     VAX/VMS; FORTRAN 77 with extensions:
C     >  IMPLICIT NONE
C     >  A few names longer than 6 characters
C     >  ! comments
C
C  HISTORY:
C     04/12/91  D.Saunders  Initial implementation of a scheme for dealing
C                           with an indefinite number of sub-distributions
C                           without a lot of application-specific prompting.
C                           Jeff Trosin suggested the file-driven approach.
C     06/24/91    "    "    Made a missing case non-fatal: return NPTS=0
C                           with a "Proceeding" diagnostic.
C     08/16/2013  "    "    Fixed a couple of typos in the header (only).
C
C  AUTHOR:   David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   LUNDAT, LUNERR, IFIRST,  ILAST, NPTS, METHOD, IPARAM (*)
      REAL
     >   RPARAM (*)
      CHARACTER
     >   CASE * (*), CPARAM * (*)

C     Local constants:

      INTEGER
     >   MAXKEYS, MAXLEN
      CHARACTER
     >   COMMENT * 1, SUBNAME * 11

      PARAMETER
     >  (COMMENT = '!',
     >   MAXLEN  = 12,   ! Max. length of any token. Beware of FORMATS below.
     >   MAXKEYS = 13,   ! Max. number of keyword/value pairs on one line
     >   SUBNAME = ' GRIDSPEC: ')

C     Local variables:

      INTEGER
     >   ENTRY, I, IOS, N, NCHARS, NPAIRS
      CHARACTER
     >   CASEKEY * (MAXLEN), CASEVAL * (MAXLEN),
     >   DICTRY (MAXKEYS) * (MAXLEN), KEYS (MAXKEYS) * (MAXLEN),
     >   LINE * 80, VALUES (MAXKEYS) * (MAXLEN)

C     Procedures:

      EXTERNAL
     >   GETLINE, LOOKUP, PAIRS

C     Storage:

      DATA      ! Keep the dictionary alphabetic uppercase if adding I4, etc.
     >   DICTRY /
     >      'C1', 'C2', 'C3', 'FIRST', 'I1', 'I2', 'I3', 'LAST',
     >      'METHOD', 'NPTS', 'R1', 'R2', 'R3'/

C     Execution:

C     NPTS = 0 on input means we have to locate CASE within the control file:

      IF (NPTS .EQ. 0) THEN

C        Tokenize the CASE string:

         NPAIRS = 1
         CALL PAIRS (CASE, NPAIRS, CASEKEY, CASEVAL)

C        Look for this case in the control file:

         REWIND LUNDAT
         READ (LUNDAT, 1003)  ! Skip title

  100    CONTINUE
            CALL GETLINE (LUNDAT, COMMENT, LINE, NCHARS, IOS)
            IF (IOS .LT. 0) GO TO 900    ! EOF: case is absent
            IF (IOS .GT. 0) GO TO 910    ! Read error
            IF (NCHARS .EQ. 0) GO TO 100 ! Blank line

            NPAIRS = 1                   ! Tokenize the first keyword/value
            CALL PAIRS (LINE (1 : NCHARS), NPAIRS, KEYS, VALUES)
            IF (KEYS (1) .NE. CASEKEY .OR. VALUES (1) .NE. CASEVAL)
     >   GO TO 100
      END IF

C     Control file is positioned within the target case.  Read the next line:

  200 CALL GETLINE (LUNDAT, COMMENT, LINE, NCHARS, IOS)
      IF (IOS .LT. 0) GO TO 800          ! Normal EOF
      IF (IOS .GT. 0) GO TO 910          ! Read error
      IF (NCHARS .EQ. 0) GO TO 200       ! Blank line

      NPAIRS = MAXKEYS                   ! Tokenize the whole line
      CALL PAIRS (LINE (1 : NCHARS), NPAIRS, KEYS, VALUES)
      IF (NPAIRS .EQ. 1) GO TO 800       ! Must be bumping into the next case.

C     Identify all the distribution parameters:

      DO 300, I = 1, NPAIRS

         CALL LOOKUP (MAXKEYS, DICTRY, .TRUE., KEYS (I), ENTRY)
         IF (ENTRY .LE. 0) GO TO 920     ! Invalid keyword

         IF (KEYS (I) (1 : 5) .EQ. 'FIRST') THEN
            READ (VALUES (I), 1001, ERR=930) IFIRST

         ELSE IF (KEYS (I) (1 : 4) .EQ. 'LAST') THEN
            READ (VALUES (I), 1001, ERR=930) ILAST

         ELSE IF (KEYS (I) (1 : 4) .EQ. 'NPTS') THEN
            READ (VALUES (I), 1001, ERR=930) NPTS

         ELSE IF (KEYS (I) (1 : 6) .EQ. 'METHOD') THEN
            READ (VALUES (I), 1001, ERR=930) METHOD

         ELSE
C           Find output element n from the In, Rn, or Cn keyword:

            READ (KEYS (I) (2 : 2), 1000) N

            IF (KEYS (I) (1 : 1) .EQ. 'I') THEN
               READ (VALUES (I), 1001, ERR=930) IPARAM (N)

            ELSE IF (KEYS (I) (1 : 1) .EQ. 'R') THEN
               READ (VALUES (I), 1002, ERR=930) RPARAM (N)

            ELSE IF (KEYS (I) (1 : 1) .EQ. 'C') THEN
               CPARAM (N : N) = VALUES (I)

            END IF
         END IF

  300 CONTINUE

      GO TO 999


  800 NPTS = 0     ! No more distributions for the present case
      GO TO 999


C     Error handling:

  900 WRITE (LUNERR, 1003)
     >   SUBNAME, 'Absent grid spec. case: ', CASE, ' Proceeding ...'
      GO TO 999

  910 WRITE (LUNERR, 1003) SUBNAME, 'Read error.'
      GO TO 990

  920 WRITE (LUNERR, 1003) SUBNAME, 'Invalid keyword = ', KEYS (I)
      GO TO 990

  930 WRITE (LUNERR, 1003) SUBNAME, 'Invalid value = ', VALUES (I),
     >                     ' Corresponding keyword: ', KEYS (I)
      GO TO 990


  990 STOP 'Aborting.'

  999 RETURN

C     Formats:

 1000 FORMAT (I1)
 1001 FORMAT (BN, I12)
 1002 FORMAT (BN, F12.0)
 1003 FORMAT (A, A, A)

      END     ! End of GRIDSPEC
